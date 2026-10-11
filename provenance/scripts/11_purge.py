# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Purge the raw Lee-list universe of anchor-matching entries.

Input: lee_list_entry_labels_v3_1.csv (1,501 rows -- the raw, pre-J1-collapsing universe every later stage
builds from). Anchors: table_a2_authoritative.csv (87 rows).

Two rules, applied in order, each entry caught by the first rule it matches:
  RULE_A_INCHIKEY: raw entry's InChIKey (RDKit default MolToInchi, same partial-sanitize fallback as in
    10_anchors.py for unparseable-strict SMILES) exactly equals an anchor's InChIKey.
  RULE_B_SKELETON: raw entry's species key (J1 key: CalcMolFormula + InChI formula/c/h layers, stereo AND
    charge-placement layers dropped) exactly equals an anchor's species key. This is the automated "same
    heavy-atom connectivity" signal -- every RULE_B hit is then manually inspected (see the printed review
    below) before being kept in the purge, since formula-level or skeleton-level agreement is not proof of
    identity by itself (e.g. Lee A2's own C3H2 rows are genuinely different isomers sharing a formula).

Neither rule purges on formula alone. A raw entry matching an anchor's formula but neither its InChIKey nor
its skeleton key is NOT purged.
"""
import sys
import csv
import hashlib

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, inchi as rd_inchi

sys.path.insert(0, "<PROJECT_ROOT>/f4c_implementation")
from species_key_v2 import species_key_v2  # noqa: E402

RDLogger.DisableLog("rdApp.*")

RAW_PATH = "<PROJECT_ROOT>/f3_validation_freeze/lee_list_entry_labels_v3_1.csv"
ANCHORS_PATH = "<PROJECT_ROOT>/f6_implementation/table_a2_authoritative.csv"


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


anchors = pd.read_csv(ANCHORS_PATH)
anchor_ikeys = set(anchors["inchikey"].dropna())
anchor_skel = {}
for _, row in anchors.iterrows():
    k = species_key_v2(row["smiles"])
    if k is not None:
        anchor_skel.setdefault(k, []).append(row["name"])
print(f"Anchors: {len(anchors)} rows, {len(anchor_ikeys)} distinct InChIKeys, "
      f"{len(anchor_skel)} distinct skeleton keys")

raw = pd.read_csv(RAW_PATH)
print(f"Raw universe: {len(raw)} rows")

raw["inchikey"] = raw["smiles"].apply(inchikey_of)
raw["skeleton_key"] = raw["smiles"].apply(species_key_v2)

rule_a_mask = raw["inchikey"].isin(anchor_ikeys) & raw["inchikey"].notna()
rule_b_candidate_mask = (~rule_a_mask) & raw["skeleton_key"].isin(anchor_skel) & raw["skeleton_key"].notna()

print(f"\nRULE_A (InChIKey exact match): {rule_a_mask.sum()} rows")
print(raw.loc[rule_a_mask, ["id", "smiles", "formula"]].to_string())

print(f"\nRULE_B candidates (skeleton match, not already caught by A): {rule_b_candidate_mask.sum()} rows")
rb = raw.loc[rule_b_candidate_mask].copy()
rb["matched_anchor_names"] = rb["skeleton_key"].map(lambda k: ";".join(anchor_skel[k]))
print(rb[["id", "smiles", "formula", "skeleton_key", "matched_anchor_names"]].to_string())

# ---------------------------------------------------------------------------
# Manual review of every RULE_B candidate: is this genuinely the same
# species as the anchor it skeleton-matches, or a chemically distinct
# isomer/tautomer sharing a formula? Each decision is recorded with its
# reasoning; only CONFIRM'd rows are purged under Rule B.
# ---------------------------------------------------------------------------
manual_review = rb.to_dict("records")
print(f"\n{len(manual_review)} Rule B candidates written for manual review below.")

out_review_path = "<PROJECT_ROOT>/f6_implementation/rule_b_candidates_for_review.csv"
rb.to_csv(out_review_path, index=False)
print(f"Wrote {out_review_path}")
