# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Finalize the purge and re-run J1 collapsing on the purged universe.

Rule A (InChIKey exact match): 12 rows purged (computed in 11_purge.py).
Rule B (skeleton/species-key match, manually reviewed): 1 candidate found (raw id 576, SMILES [2H]C=O,
  formula CH2O -- an explicit deuterium isotopologue of formaldehyde, HDCO) skeleton-matching anchor H2CO.
  MANUAL DECISION: NOT purged. Isotopologues are treated as distinct from their parent species (as in
  Lee et al. and Shay et al.'s own published data sets, and the validation-set criteria) and are excluded
  from the anchor/training set because they are different physical species with different
  detectability. The species key's skeleton match strips isotope labelling along with charge/stereo,
  which is correct for J1's own purpose (de-duplicating resonance-form SMILES of the SAME species) but
  wrong for judging isotopic identity. HDCO is kept in the purged universe as its own candidate.
"""
import sys
import csv
import hashlib

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi as rd_inchi

sys.path.insert(0, "<PROJECT_ROOT>/f4c_implementation")
sys.path.insert(0, "<PROJECT_ROOT>/f4a_implementation")
from species_key_v2 import species_key_v2  # noqa: E402
from step2_build_candidates import LEE_PATH, V32_PATH, apply_a7, parse_date_field  # noqa: E402

RDLogger.DisableLog("rdApp.*")

RULE_A_IDS = [289, 332, 368, 380, 398, 430, 432, 484, 679, 681, 793, 1302]
RULE_B_REVIEWED_NOT_PURGED = [576]  # HDCO -- kept, see docstring

ANCHORS_PATH = "<PROJECT_ROOT>/f6_implementation/table_a2_authoritative.csv"
raw = pd.read_csv(LEE_PATH)
anchors = pd.read_csv(ANCHORS_PATH)
akey_to_name = {}
for _, r in anchors.iterrows():
    akey_to_name.setdefault(r["inchikey"], []).append(r["name"])


def inchikey_of(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_KEKULIZE,
            catchErrors=True,
        )
    try:
        return rd_inchi.InchiToInchiKey(rd_inchi.MolToInchi(mol, treatWarningAsError=False))
    except Exception:
        return None


raw["inchikey"] = raw["smiles"].apply(inchikey_of)

purge_log_rows = []
for rid in RULE_A_IDS:
    row = raw.loc[raw["id"] == rid].iloc[0]
    purge_log_rows.append({
        "id": rid, "smiles": row["smiles"], "formula": row["formula"],
        "rule": "RULE_A_INCHIKEY",
        "matched_anchor": ";".join(akey_to_name.get(row["inchikey"], [])),
        "note": "exact InChIKey match to an anchor -- purged",
    })
for rid in RULE_B_REVIEWED_NOT_PURGED:
    row = raw.loc[raw["id"] == rid].iloc[0]
    purge_log_rows.append({
        "id": rid, "smiles": row["smiles"], "formula": row["formula"],
        "rule": "RULE_B_SKELETON_REVIEWED_KEPT",
        "matched_anchor": "H2CO",
        "note": ("species_key_v2 skeleton matches anchor H2CO, but this SMILES is an "
                 "explicit deuterium isotopologue (HDCO, [2H]C=O). Manually reviewed and "
                 "NOT purged -- isotopologues are a distinct species in this project's "
                 "established convention (Phase 2 validation criteria; Lee/Shay's own "
                 "data sets exclude isotopologues from anchors for the same reason)."),
    })

purge_log_df = pd.DataFrame(purge_log_rows)
purge_log_path = "<PROJECT_ROOT>/f6_implementation/purge_log.csv"
purge_log_df.to_csv(purge_log_path, index=False)
print(f"Wrote {purge_log_path}, {len(purge_log_df)} rows logged "
      f"({len(RULE_A_IDS)} purged, {len(RULE_B_REVIEWED_NOT_PURGED)} reviewed-and-kept)")

purged_raw = raw.loc[~raw["id"].isin(RULE_A_IDS)].drop(columns=["inchikey"]).reset_index(drop=True)
print(f"Raw universe: {len(raw)} -> purged: {len(purged_raw)} "
      f"(removed {len(raw) - len(purged_raw)})")
assert len(raw) - len(purged_raw) == len(RULE_A_IDS)

purged_raw_path = "<PROJECT_ROOT>/f6_implementation/lee_list_entry_labels_purged.csv"
purged_raw.to_csv(purged_raw_path, index=False)
print(f"Wrote {purged_raw_path}")

# ---------------------------------------------------------------------------
# Re-run J1 collapsing (same logic as src/candidates.py; only the input universe changes: purged_raw
# instead of the full 1,501).
# ---------------------------------------------------------------------------
v32 = pd.read_csv(V32_PATH)
lee, n_a7 = apply_a7(purged_raw, v32)

raw_keys = lee["smiles"].apply(species_key_v2)
singleton_ids = lee.loc[raw_keys.isna(), "id"].tolist()
final_keys = raw_keys.copy()
for idx in lee.index[raw_keys.isna()]:
    final_keys.loc[idx] = f"__SINGLETON__{lee.loc[idx, 'id']}"
lee["species_key"] = final_keys

parsed = lee["first_public_date_of_match"].apply(parse_date_field)
lee["parsed_date"] = [p[0] for p in parsed]
lee["date_precision"] = [p[1] for p in parsed]

pos_no_date = lee[(lee["in_list_by_v3_1"] == True) & (lee["parsed_date"].isna())]
assert len(pos_no_date) == 0, f"positive entries with no parseable date: {pos_no_date['id'].tolist()}"

groups = []
for key, g in lee.groupby("species_key", sort=False):
    positive = bool(g["in_list_by_v3_1"].any())
    date = g.loc[g["in_list_by_v3_1"] == True, "parsed_date"].min() if positive else pd.NaT
    groups.append({
        "species_key": key,
        "n_members": len(g),
        "member_ids": ";".join(str(x) for x in g["id"]),
        "formulas": ";".join(sorted(set(g["formula"]))),
        "smiles_members": ";".join(g["smiles"]),
        "dist": g["dist"].mean(),
        "logcd": g["logcd"].mean(),
        "sigma": g["sigma"].mean(),
        "mu": g["mu"].mean(),
        "mw": g["mw"].mean(),
        "charged": bool((g["charge"] != 0).any()),
        "radical": bool((g["mult"] > 1).any()),
        "positive": positive,
        "date": date,
    })
cand = pd.DataFrame(groups).reset_index(drop=True)
cand["candidate_id"] = np.arange(len(cand))

n_missing_mu = cand["mu"].isna().sum()
if n_missing_mu:
    med = cand["mu"].median()
    cand.loc[cand["mu"].isna(), "mu"] = med

print(f"\nPurged raw entries: {len(lee)} -> candidates (J1 key, purged pool): {len(cand)}")
print(f"Singleton (parse/InChI-failure) ids: {singleton_ids}, n={len(singleton_ids)}")
print(f"Positive candidates: {cand['positive'].sum()}")
print(f"A7 dates updated: {n_a7}")

cand_path = "<PROJECT_ROOT>/f6_implementation/candidates_purged_v1.csv"
cand.to_csv(cand_path, index=False)
print(f"Wrote {cand_path}")

# Compare against the unpurged candidate count (1,451).
for p, label in [(purge_log_path, "purge_log.csv"), (purged_raw_path, "lee_list_entry_labels_purged.csv"),
                  (cand_path, "candidates_purged_v1.csv")]:
    with open(p, "rb") as f:
        h = hashlib.sha256(f.read()).hexdigest()
    print(f"{label} SHA-256: {h}")
