# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Re-run A/B/C matching of detections v3.1 (77 rows) against Lee's 1,501-entry universe
(lee_list_entry_labels.csv: lee1510_canonicalized.csv's 1510 rows minus 9 already-known-anchor duplicates,
same ids, same SMILES for the 1501 retained). Build lee_list_entry_labels_v3_1.csv (all 1501 + first_public_date
of the matched detection(s)) and reconcile against the original 40/76 labeling.

Ties are handled explicitly: in_lee_list is computed PER DETECTION directly against the full Lee-side key
sets (independent of any other detection), and ALL matching detections are listed on the Lee-entry side.
Reversing a Lee-entry -> single-detection map would silently drop ties (e.g. cis- and trans-CH3CHCHCN both
reduce to the same stereo-free skeleton and match the SAME Lee entry).
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

lee1510 = pd.read_csv("<PROJECT_ROOT>/lee_broader_pool_feasibility/lee1510_canonicalized.csv")
lee1501 = pd.read_csv("<PROJECT_ROOT>/lee_list_entry_labels.csv")
v3_1 = pd.read_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3_1.csv")

lee = lee1510[lee1510["id"].isin(lee1501["id"])].copy()
assert len(lee) == 1501

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

def stereo_free(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)

lee["has_stereo"] = lee["smiles"].apply(has_stereo_marks)

# ---- Lee-side key sets (for A/B/C membership testing) ----
lee_A_set = set(lee["canonical_smiles_full"].dropna())

lee_noStereo = lee[~lee["has_stereo"]].copy()
lee_noStereo["stereo_free_smiles"] = lee_noStereo["smiles"].apply(stereo_free)
lee_B_set = set(lee_noStereo["stereo_free_smiles"].dropna())

def lee_inchi_key(smi):
    if not isinstance(smi, str) or not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return inchi_skeleton_key(m)

lee_noStereo["inchi_key"] = lee_noStereo["canonical_smiles_full"].apply(lee_inchi_key)
lee_C_set = set(lee_noStereo["inchi_key"].dropna())

# ---- v3.1-side keys, and DIRECT per-detection in_lee_list (independent test) ----
v3rows = []
for _, r in v3_1.iterrows():
    smi = r["smiles"]
    m = Chem.MolFromSmiles(smi)
    if m is None:
        v3rows.append(dict(formula=r["formula"], parse_ok=False, level_A=False,
                            level_B=False, level_C=False, in_lee_list=False))
        continue
    can_full = Chem.MolToSmiles(m)
    m2 = Chem.MolFromSmiles(smi)
    Chem.RemoveStereochemistry(m2)
    sf = Chem.MolToSmiles(m2)
    ik = inchi_skeleton_key(m)
    A = can_full in lee_A_set
    B = sf in lee_B_set
    C = (ik in lee_C_set) if ik else False
    v3rows.append(dict(formula=r["formula"], parse_ok=True, can_full=can_full, sf=sf,
                        ik=ik, level_A=A, level_B=B, level_C=C, in_lee_list=A or B or C))
v3keys = pd.DataFrame(v3rows)

v3_1 = v3_1.merge(v3keys[["formula", "level_A", "level_B", "level_C", "in_lee_list"]],
                   on="formula", how="left")
v3_1.to_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3_1.csv", index=False)
print(f"v3.1 detections in Lee's list (direct per-detection test): "
      f"{v3_1['in_lee_list'].sum()} / {len(v3_1)}")

# ---- Lee-entry side: list ALL matching detections (not just one) ----
det_fpd = dict(zip(v3_1["formula"], v3_1["first_public_date"]))

lee_A_map, lee_B_map, lee_C_map = {}, {}, {}
for _, d in v3keys[v3keys["parse_ok"]].iterrows():
    lee_A_map.setdefault(d["can_full"], []).append(d["formula"])
    lee_B_map.setdefault(d["sf"], []).append(d["formula"])
    if d["ik"]:
        lee_C_map.setdefault(d["ik"], []).append(d["formula"])

lee = lee.merge(lee_noStereo[["id", "stereo_free_smiles", "inchi_key"]], on="id", how="left")

def match_all(row):
    matches = []
    can = row["canonical_smiles_full"]
    if isinstance(can, str):
        matches += [(f, "A") for f in lee_A_map.get(can, [])]
    sf = row.get("stereo_free_smiles")
    if isinstance(sf, str):
        matches += [(f, "B") for f in lee_B_map.get(sf, []) if f not in [m[0] for m in matches]]
    ik = row.get("inchi_key")
    if isinstance(ik, str):
        matches += [(f, "C") for f in lee_C_map.get(ik, []) if f not in [m[0] for m in matches]]
    if not matches:
        return "", "", ""
    formulas = [m[0] for m in matches]
    levels = [m[1] for m in matches]
    fpds = [str(det_fpd.get(f, "")) for f in formulas]
    return ";".join(formulas), ";".join(levels), ";".join(fpds)

matched = lee.apply(match_all, axis=1, result_type="expand")
matched.columns = ["matched_detection_formula", "match_level", "first_public_date_of_match"]
lee = pd.concat([lee, matched], axis=1)

out = lee1501.merge(
    lee[["id", "matched_detection_formula", "match_level", "first_public_date_of_match"]],
    on="id", how="left")
out["in_list_by_v3_1"] = out["matched_detection_formula"].str.len() > 0
out.to_csv("<PROJECT_ROOT>/f3_validation_freeze/lee_list_entry_labels_v3_1.csv", index=False)

print(f"Lee entries matched to >=1 v3.1 detection: {out['in_list_by_v3_1'].sum()} / {len(out)}")

# check for Lee entries matched by >1 detection (ties, e.g. cis/trans case)
ties = out[out["matched_detection_formula"].str.contains(";", na=False)]
print(f"\nLee entries matched by >1 v3.1 detection (documented, not dropped): {len(ties)}")
if len(ties):
    print(ties[["id", "formula", "matched_detection_formula", "match_level"]].to_string())

# ---- reconciliation: original 40/76 vs new v3.1-based count ----
orig = pd.read_csv("<PROJECT_ROOT>/detections_labeled_vs_lee.csv")
orig["orig_in_list"] = orig["level"].isin(
    ["A_exact", "B_exact_stereo_agnostic", "C_same_formula_and_H_skeleton"])
comp = orig[["formula", "level", "orig_in_list"]].merge(
    v3_1[["formula", "in_lee_list"]], on="formula", how="outer", indicator=True)
comp = comp[comp["formula"] != "C4H8"]
diff = comp[comp["orig_in_list"] != comp["in_lee_list"]]
print(f"\nReconciliation vs original 40/76 labeling (excluding dropped C4H8, "
      f"{len(comp)} comparable formulas):")
print(diff.to_string())
diff.to_csv("<PROJECT_ROOT>/f3_validation_freeze/reconcile_v3_1_diffs.csv", index=False)
print(f"\nFinal v3.1 in-list count: {v3_1['in_lee_list'].sum()} / {len(v3_1)}")
