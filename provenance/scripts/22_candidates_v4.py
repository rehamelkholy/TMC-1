# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Purge (RULE_A/RULE_B, unchanged from 11_purge.py / 12_collapse_candidates.py -- the same 12 RULE_A ids + 1
RULE_B reviewed-and-kept id, since these depend only on the raw 1,501-entry Lee universe vs. the anchor
table, not on which detections file is used) + unchanged J1 species-key collapsing logic (species_key,
apply_a7, parse_date_field) on lee_list_entry_labels_v4.csv in place of lee_list_entry_labels_v3_1.csv.
"""
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "<PROJECT_ROOT>/f4a_implementation")
sys.path.insert(0, "<PROJECT_ROOT>/f4c_implementation")
from step2_build_candidates import apply_a7, parse_date_field  # noqa: E402
from species_key_v2 import species_key_v2  # noqa: E402

RULE_A_IDS = [289, 332, 368, 380, 398, 430, 432, 484, 679, 681, 793, 1302]
RULE_B_REVIEWED_NOT_PURGED = [576]  # HDCO

LEE_V4_PATH = "<PROJECT_ROOT>/f9_implementation/lee_list_entry_labels_v4.csv"
V4_PATH = "<PROJECT_ROOT>/f9_implementation/post2021_detections_v4.csv"

raw = pd.read_csv(LEE_V4_PATH)
assert len(raw) == 1501

purged_raw = raw.loc[~raw["id"].isin(RULE_A_IDS)].reset_index(drop=True)
assert len(raw) - len(purged_raw) == len(RULE_A_IDS)
purged_path = "<PROJECT_ROOT>/f9_implementation/lee_list_entry_labels_purged_v4.csv"
purged_raw.to_csv(purged_path, index=False)
print(f"Raw universe: {len(raw)} -> purged: {len(purged_raw)} (removed {len(RULE_A_IDS)}, "
      f"RULE_B id {RULE_B_REVIEWED_NOT_PURGED} reviewed-and-kept, unchanged from F6)")

v4 = pd.read_csv(V4_PATH)
lee, n_a7 = apply_a7(purged_raw, v4)  # v4 already carries A7-corrected dates; idempotent

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
        "species_key": key, "n_members": len(g),
        "member_ids": ";".join(str(x) for x in g["id"]),
        "formulas": ";".join(sorted(set(g["formula"]))),
        "smiles_members": ";".join(g["smiles"]),
        "dist": g["dist"].mean(), "logcd": g["logcd"].mean(), "sigma": g["sigma"].mean(),
        "mu": g["mu"].mean(), "mw": g["mw"].mean(),
        "charged": bool((g["charge"] != 0).any()), "radical": bool((g["mult"] > 1).any()),
        "positive": positive, "date": date,
    })
cand = pd.DataFrame(groups).reset_index(drop=True)
cand["candidate_id"] = np.arange(len(cand))

n_missing_mu = cand["mu"].isna().sum()
if n_missing_mu:
    med = cand["mu"].median()
    cand.loc[cand["mu"].isna(), "mu"] = med

print(f"\nPurged raw entries: {len(lee)} -> candidates (J1 key): {len(cand)}")
print(f"Singleton ids: {singleton_ids}, n={len(singleton_ids)}")
print(f"Positive candidates: {cand['positive'].sum()}")
print(f"A7 dates updated (idempotent, already correct in v4): {n_a7}")
print(f"Missing mu imputed: {n_missing_mu}")

cand_path = "<PROJECT_ROOT>/f9_implementation/candidates_purged_v4.csv"
cand.to_csv(cand_path, index=False)
print(f"Wrote {cand_path}")

# Compare against the purged candidate table of the v3.2-based universe
cand_f6 = pd.read_csv("<PROJECT_ROOT>/f6_implementation/candidates_purged_v1.csv")
print(f"\nF6 (v3.2-based) candidates: {len(cand_f6)}, positive: {cand_f6['positive'].sum()}")
print(f"F9 (v4-based)   candidates: {len(cand)}, positive: {cand['positive'].sum()}")

import hashlib
for p, label in [(purged_path, "lee_list_entry_labels_purged_v4.csv"), (cand_path, "candidates_purged_v4.csv")]:
    with open(p, "rb") as f:
        h = hashlib.sha256(f.read()).hexdigest()
    print(f"{label} SHA-256: {h}")
