# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Canonicalize the reference lists with the same full-sanitization
convention used for the pool, so the cross-list join in the next step is
apples-to-apples.

detections_labeled_vs_lee.csv (76 detections) was not available when this
script was run -- this script canonicalizes what IS available: Lee's
1,510-entry list and the 40 pos_primary in-list detections from
lee_list_entry_labels.csv. The missing 76-detections file is flagged, not
substituted or guessed.
"""
import os
import sys
import warnings

warnings.filterwarnings("ignore")

ROOT = "<PROJECT_ROOT>"
out_dir = os.path.join(ROOT, "lee_broader_pool_feasibility")

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def canon_full(smi):
    if pd.isna(smi):
        return None
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        return Chem.MolToSmiles(m)
    except Exception:
        return None


def classify(smi):
    """closed_shell / radical / ion / unknown_sanitize_failed"""
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return "unknown_sanitize_failed"
        charged = any(a.GetFormalCharge() != 0 for a in m.GetAtoms())
        radical = sum(a.GetNumRadicalElectrons() for a in m.GetAtoms()) > 0
        if charged:
            return "ion"
        if radical:
            return "radical"
        return "closed_shell"
    except Exception:
        return "unknown_sanitize_failed"


# ---------------------------------------------------------------------
# Lee's 1,510-entry recommendation list
# ---------------------------------------------------------------------
lee1510 = pd.read_csv(os.path.join(ROOT, "lee_recommendations_parsed.csv"))
lee1510["canonical_smiles_full"] = lee1510["smiles"].apply(canon_full)
lee1510["class"] = lee1510["smiles"].apply(classify)
n_fail = lee1510["canonical_smiles_full"].isna().sum()
print(f"Lee's 1,510-entry list: {len(lee1510)} rows, {n_fail} fail full-sanitize "
      f"canonicalization ({n_fail} 'unknown_sanitize_failed').")
lee1510.to_csv(os.path.join(out_dir, "lee1510_canonicalized.csv"), index=False)

# ---------------------------------------------------------------------
# The 40 pos_primary in-list detections (from lee_list_entry_labels.csv)
# ---------------------------------------------------------------------
labels = pd.read_csv(os.path.join(ROOT, "lee_list_entry_labels.csv"))
pos40 = labels[labels["pos_primary"] == 1].copy()
assert len(pos40) == 40
pos40["canonical_smiles_full"] = pos40["smiles"].apply(canon_full)
pos40["class"] = pos40["smiles"].apply(classify)
n_fail40 = pos40["canonical_smiles_full"].isna().sum()
print(f"40 pos_primary in-list detections: {n_fail40} fail full-sanitize canonicalization.")
pos40.to_csv(os.path.join(out_dir, "pos40_canonicalized.csv"), index=False)

print("\nClass breakdown, Lee's 1,510:")
print(lee1510["class"].value_counts())
print("\nClass breakdown, 40 pos_primary detections:")
print(pos40["class"].value_counts())

print("\nNOTE: detections_labeled_vs_lee.csv (76 detections) was not found on "
      "Reham's connected PC -- item 2's coverage-of-76-detections sub-check "
      "could not be run. Flagged for Reham in the final report.")
