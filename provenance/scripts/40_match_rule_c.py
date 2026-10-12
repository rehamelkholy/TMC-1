# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Match every member SMILES of all 1,439 candidates in
candidates_purged_v5.csv to (a) the 87 Lee anchors (table_a2_authoritative_v2.csv)
and (b) the TMC-1 inventory (item 1), using rule C (this project's rule_c.py,
validated against all 5 provided anchor examples and the general
isotopologue/isomer differentiation tests -- see rule_c.py's __main__ block).

Each candidate's member SMILES (smiles_members, ';'-separated for n_members>1)
is matched independently -- NOT the pre-collapsed species_key column -- per
the requirement to match EVERY member SMILES, so this does not
just trust the existing (buggy, per the v5 audit) collapsing/purging logic.
"""
import pandas as pd
from rule_c import rule_c_key

CAND_PATH = "<PROJECT_ROOT>/f12_implementation/candidates_purged_v5.csv"
ANCHORS_PATH = "<PROJECT_ROOT>/f7_implementation/table_a2_authoritative_v2.csv"

# ---- Inventory (item 1): TMC-1 species compiled from Kaifu et al. 2004 /
# Gratier et al. 2016 / Agundez & Wakelam 2013 / astromol census / ADS,
# restricted (after a formula-presence pre-filter against candidates_purged_v5.csv)
# to species whose formula actually occurs among the 1,439 candidates' members.
# name -> (SMILES, first_public_date, date_basis, reference)
INVENTORY = {
    "HCO+": ("C#[O+]", None, "unresolved", "bounded pre-1978 by Langer et al. 1978, ApJ 225, L139 (H13CO+ toward TMC-1); exact TMC-1 discovery paper not pinned -- see note"),
    "H2S": ("S", "1989-10", "online", "Minh, Irvine & Ziurys 1989, ApJ 345, L63"),
    "NS": ("[N]=[S]", "1994-01", "online", "McGonagle, Irvine & Ohishi 1994, ApJ 422, 621"),
    "HNCS": ("N=C=S", "2010-12", "online", "Adande, Halfen, Ziurys, Quan & Herbst 2010, ApJ 725, 561"),
    "HC(O)CN": ("O=CC#N", "2021-05-27", "arXiv v1 (arXiv:2105.12996)", "Cernicharo et al. 2021, A&A 650, L14"),
    "CH3O": ("[O]C", None, "unresolved", "methoxy radical -- TMC-1 reference not yet located; only checked because its formula occurs among the candidates"),
    "HCNH+": ("C#[NH+]", None, "unresolved", "protonated HCN -- TMC-1 reference not yet located; only checked because its formula occurs among the candidates"),
    "H2O": ("O", None, "unresolved", "ubiquitous ISM species -- TMC-1 reference not yet located; checked only because its formula occurs among the candidates"),
    "O2": ("O=O", None, "unresolved", "TMC-1 (tentative in general ISM; checked only because its formula occurs among the candidates)"),
    "CO": ("[C-]#[O+]", None, "unresolved", "ubiquitous ISM species -- TMC-1 reference not yet located; checked only because its formula occurs among the candidates"),
    "NO": ("[N]=O", None, "unresolved", "TMC-1 reference not yet located; checked only because its formula occurs among the candidates"),
}

# ---- Load anchors, compute rule C keys ----
anchors = pd.read_csv(ANCHORS_PATH)
anchor_keys = {}
for _, row in anchors.iterrows():
    k, charge, formula = rule_c_key(row["smiles"])
    if k is not None:
        anchor_keys.setdefault(k, []).append(row["name"])
print(f"Anchors: {len(anchors)} rows -> {len(anchor_keys)} distinct rule-C keys")

# ---- Inventory keys ----
inv_keys = {}
for name, (smi, date, basis, ref) in INVENTORY.items():
    k, charge, formula = rule_c_key(smi)
    if k is None:
        print(f"WARNING: inventory species {name!r} SMILES {smi!r} failed to parse under rule C")
        continue
    inv_keys.setdefault(k, []).append(name)
print(f"Inventory: {len(INVENTORY)} species -> {len(inv_keys)} distinct rule-C keys")

# ---- Explode candidates into member SMILES ----
cand = pd.read_csv(CAND_PATH)
rows = []
for _, row in cand.iterrows():
    ids = str(row["member_ids"]).split(";")
    smis = str(row["smiles_members"]).split(";")
    assert len(ids) == len(smis), (row["candidate_id"], ids, smis)
    for mid, smi in zip(ids, smis):
        rows.append({"candidate_id": row["candidate_id"], "member_id": mid.strip(),
                      "member_smiles": smi.strip(), "candidate_positive": row["positive"],
                      "candidate_formula": row["formulas"]})
members = pd.DataFrame(rows)
print(f"Candidates: {len(cand)} rows -> {len(members)} member SMILES exploded")

members["rule_c_key"] = members["member_smiles"].apply(lambda s: rule_c_key(s)[0])
n_unparsed = members["rule_c_key"].isna().sum()
print(f"Unparsed members (rule C key = None): {n_unparsed}")
if n_unparsed:
    print(members.loc[members["rule_c_key"].isna(), ["candidate_id", "member_id", "member_smiles"]].to_string())

# ---- Match against anchors ----
anchor_matches = []
for _, r in members.iterrows():
    if r["rule_c_key"] in anchor_keys:
        anchor_matches.append({
            "candidate_id": r["candidate_id"], "member_id": r["member_id"],
            "member_smiles": r["member_smiles"], "candidate_formula": r["candidate_formula"],
            "candidate_already_positive": bool(r["candidate_positive"]),
            "matched_anchor_names": ";".join(anchor_keys[r["rule_c_key"]]),
            "rule_c_key": r["rule_c_key"],
        })
anchor_matches_df = pd.DataFrame(anchor_matches)
print(f"\nAnchor matches (member level): {len(anchor_matches_df)}")
if len(anchor_matches_df):
    print(anchor_matches_df.to_string(index=False))

# ---- Match against inventory (excluding any member that already matched an anchor) ----
anchor_matched_member_keys = set(anchor_matches_df["rule_c_key"]) if len(anchor_matches_df) else set()
inv_matches = []
for _, r in members.iterrows():
    if r["rule_c_key"] in inv_keys and r["rule_c_key"] not in anchor_matched_member_keys:
        inv_matches.append({
            "candidate_id": r["candidate_id"], "member_id": r["member_id"],
            "member_smiles": r["member_smiles"], "candidate_formula": r["candidate_formula"],
            "candidate_already_positive": bool(r["candidate_positive"]),
            "matched_inventory_names": ";".join(inv_keys[r["rule_c_key"]]),
            "rule_c_key": r["rule_c_key"],
        })
inv_matches_df = pd.DataFrame(inv_matches)
print(f"\nInventory matches (member level): {len(inv_matches_df)}")
if len(inv_matches_df):
    print(inv_matches_df.to_string(index=False))

out_dir = "<PROJECT_ROOT>/f13_implementation"
anchor_matches_df.to_csv(f"{out_dir}/item2_anchor_matches.csv", index=False)
inv_matches_df.to_csv(f"{out_dir}/item2_inventory_matches.csv", index=False)
print(f"\nWrote {out_dir}/item2_anchor_matches.csv, {out_dir}/item2_inventory_matches.csv")

# ---- Confirm the 5+5 expected matches are present ----
EXPECTED_ANCHOR_CIDS = [293, 376, 396, 490, 523]
EXPECTED_INV_NAMES = ["HCO+", "H2S", "NS", "HNCS", "HC(O)CN"]
got_anchor_cids = set(anchor_matches_df["candidate_id"]) if len(anchor_matches_df) else set()
print("\nExpected anchor candidate_ids found:", {c: (c in got_anchor_cids) for c in EXPECTED_ANCHOR_CIDS})
got_inv_names = set()
if len(inv_matches_df):
    for s in inv_matches_df["matched_inventory_names"]:
        got_inv_names.update(s.split(";"))
print("Expected inventory names found:", {n: (n in got_inv_names) for n in EXPECTED_INV_NAMES})
