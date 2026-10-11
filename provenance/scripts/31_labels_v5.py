# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
lee_list_entry_labels_v4.csv -> lee_list_entry_labels_v5.csv. Matching logic identical to
21_labels_v4.py (rules A/B/C: exact canonical SMILES, stereo-free canonical SMILES, InChI-skeleton); the
ONLY change is the detection-side input: post2021_detections_v5.csv (85 rows) in place of
post2021_detections_v4.csv (81 rows). Does not touch v4. Also refreshes post2021_detections_v5.csv's own
level_A/level_B/level_C/in_lee_list columns against the TRUE 1,501-row raw universe.
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

lee1510 = pd.read_csv("<PROJECT_ROOT>/lee_broader_pool_feasibility/lee1510_canonicalized.csv")
lee1501 = pd.read_csv("<PROJECT_ROOT>/lee_list_entry_labels.csv")
v5 = pd.read_csv("<PROJECT_ROOT>/f12_implementation/post2021_detections_v5.csv")

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

# ---- v5-side keys, direct per-detection in_lee_list ----
v5rows = []
for _, r in v5.iterrows():
    smi = r["smiles"]
    m = Chem.MolFromSmiles(smi)
    if m is None:
        v5rows.append(dict(formula=r["formula"], parse_ok=False, level_A=False,
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
    v5rows.append(dict(formula=r["formula"], parse_ok=True, can_full=can_full, sf=sf,
                        ik=ik, level_A=A, level_B=B, level_C=C, in_lee_list=A or B or C))
v5keys = pd.DataFrame(v5rows)

v5_out = v5.drop(columns=["level_A", "level_B", "level_C", "in_lee_list"]).merge(
    v5keys[["formula", "level_A", "level_B", "level_C", "in_lee_list"]], on="formula", how="left")
v5_out.to_csv("<PROJECT_ROOT>/f12_implementation/post2021_detections_v5.csv", index=False)
print(f"v5 detections in Lee's list (direct per-detection test, TRUE 1501-row universe): "
      f"{v5_out['in_lee_list'].sum()} / {len(v5_out)}")
print(v5keys[["formula", "level_A", "level_B", "level_C", "in_lee_list"]].to_string(index=False))

# ---- Lee-entry side: list ALL matching detections (not just one) ----
det_fpd = dict(zip(v5["formula"], v5["first_public_date"]))

lee_A_map, lee_B_map, lee_C_map = {}, {}, {}
for _, d in v5keys[v5keys["parse_ok"]].iterrows():
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
out["in_list_by_v3_1"] = out["matched_detection_formula"].str.len() > 0  # column name kept
# for downstream code compatibility (32_candidates_v5.py etc. hard-reference it)
out.to_csv("<PROJECT_ROOT>/f12_implementation/lee_list_entry_labels_v5.csv", index=False)

print(f"\nLee entries matched to >=1 v5 detection: {out['in_list_by_v3_1'].sum()} / {len(out)}")

ties = out[out["matched_detection_formula"].str.contains(";", na=False)]
print(f"Lee entries matched by >1 v5 detection (documented, not dropped): {len(ties)}")
if len(ties):
    print(ties[["id", "formula", "matched_detection_formula", "match_level"]].to_string())

# ---- reconciliation vs v4-based lee labels: which raw ids changed status ----
v4_labels = pd.read_csv("<PROJECT_ROOT>/f9_implementation/lee_list_entry_labels_v4.csv")
comp = v4_labels[["id", "in_list_by_v3_1"]].rename(columns={"in_list_by_v3_1": "in_list_v4"}).merge(
    out[["id", "in_list_by_v3_1"]].rename(columns={"in_list_by_v3_1": "in_list_v5"}), on="id")
newly_positive = comp[(~comp["in_list_v4"]) & (comp["in_list_v5"])]
print(f"\nRaw Lee entries newly positive under v5 (were negative under v4): {len(newly_positive)}")
if len(newly_positive):
    detail = newly_positive.merge(out[["id", "formula", "matched_detection_formula"]], on="id")
    print(detail.to_string(index=False))
newly_positive.merge(out[["id", "formula", "matched_detection_formula"]], on="id").to_csv(
    "<PROJECT_ROOT>/f12_implementation/newly_positive_raw_ids_v5.csv", index=False)
