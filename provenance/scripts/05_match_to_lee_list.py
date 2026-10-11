# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Match the 76 v3 detections to Lee et al. (2021)'s 1510-entry released list, using the same A/B/C rule as
protocol_v1.md:
  A: identical canonical SMILES
  B: identical after removing stereo (only if the Lee entry itself carries no stereo marks, per
     protocol_v1.md's exact wording)
  C: same formula (incl. charge) + same connectivity/H-skeleton via InChI (only if the Lee entry carries no
     stereo marks)
This determines in-list vs out-of-list for the half-year table.
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

lee = pd.read_csv("<PROJECT_ROOT>/lee_broader_pool_feasibility/lee1510_canonicalized.csv")
v3 = pd.read_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3.csv")

def has_stereo_marks(smi):
    return ("@" in str(smi)) or ("/" in str(smi)) or ("\\" in str(smi))

def inchi_skeleton_key(mol):
    s = inchi.MolToInchi(mol, logLevel=None)
    if not s:
        return None
    parts = s.split("/")
    layers = {}
    for p in parts[2:]:
        if p and p[0].isalpha():
            layers[p[0]] = p
    key_layers = {k: v for k, v in layers.items() if k not in ("b", "t", "m", "s")}
    formula_token = parts[1] if len(parts) > 1 else ""
    return formula_token + "|" + "|".join(f"{k}:{v}" for k, v in sorted(key_layers.items()))

# Build Lee-side lookup structures
lee["has_stereo"] = lee["smiles"].apply(has_stereo_marks)
lee_A = set(lee["canonical_smiles_full"].dropna())

lee_noStereo = lee[~lee["has_stereo"]].copy()
def stereo_free(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)
lee_noStereo["stereo_free_smiles"] = lee_noStereo["smiles"].apply(stereo_free)
lee_B = set(lee_noStereo["stereo_free_smiles"].dropna())

def lee_inchi_key(row):
    smi = row["canonical_smiles_full"]
    if not isinstance(smi, str) or not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return inchi_skeleton_key(m)
lee_noStereo["inchi_key"] = lee_noStereo.apply(lee_inchi_key, axis=1)
lee_C = set(lee_noStereo["inchi_key"].dropna())

results = []
for _, r in v3.iterrows():
    smi = r["smiles"]
    m = Chem.MolFromSmiles(smi)
    if m is None:
        results.append({"formula": r["formula"], "smiles": smi, "parse_ok": False,
                          "level_A": False, "level_B": False, "level_C": False,
                          "in_lee_list": False})
        continue
    can_full = Chem.MolToSmiles(m)
    m2 = Chem.MolFromSmiles(smi)
    Chem.RemoveStereochemistry(m2)
    sf = Chem.MolToSmiles(m2)
    ik = inchi_skeleton_key(m)

    A = can_full in lee_A
    B = sf in lee_B
    C = (ik in lee_C) if ik else False
    results.append({"formula": r["formula"], "smiles": smi, "parse_ok": True,
                      "level_A": A, "level_B": B, "level_C": C,
                      "in_lee_list": A or B or C})

res = pd.DataFrame(results)
res.to_csv("<PROJECT_ROOT>/f3_validation_freeze/v3_vs_lee1510_match.csv", index=False)
print(res["in_lee_list"].value_counts())
print(f"\nA:{res['level_A'].sum()} B:{res['level_B'].sum()} C:{res['level_C'].sum()} "
      f"union(in_lee_list):{res['in_lee_list'].sum()} / {len(res)}")
