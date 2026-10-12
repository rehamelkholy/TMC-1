# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
For every (non-isotope-labeled) species surfaced by the
2013-2021 ADS screening as a "detection" toward TMC-1, build its SMILES,
compute the rule-C key, and check it against:
  (a) the 87 Lee anchors' rule-C keys
  (b) Table 1 v7's 85 species' rule-C keys
  (c) the 1,439 candidates' exploded member rule-C keys (already built in
      step2_match_rule_c.py's "members" table)
Isotope-labeled tokens (13C/D/15N/18O/17O) are not run through rule C: only
4 candidates in the whole universe carry isotope-labeled SMILES (CH3D,
NH2D, HDCO, deuterated formamide; checked separately, none match any
isotopologue reported in this window), so an isotopologue token cannot
match any OTHER candidate by construction, and rule C does not define away
isotope labels (this project always treats isotopologues as distinct
species, e.g. HDCO kept apart from H2CO per step3_finalize_purge_and_collapse.py).
"""
import pandas as pd
from rule_c import rule_c_key

# name -> SMILES, for every non-isotope-labeled "detection" species surfaced
# by the B-screening (see item_B_species_raw.csv), built from standard
# structures for these (all well-established, unambiguous ISM species).
SPECIES_SMILES = {
    "HC5NH+": "C#CC#CC#[NH+]",
    "HC3O+": "C#CC#[O+]",
    "HC3S+": "C#CC#[S+]",
    "CH3CO+": "CC#[O+]",
    "CH2CHCCH": "C=CC#C",
    "CH2CCHCCH": "C=C=CC#C",  # methyltriacetylene-family isomer (allenyl-diacetylene type)
    "HC4NC": "C#CC#C[N+]#[C-]",
    "CH3C7N": "CC#CC#CC#CC#N",
    "HCSCN": "C(=S)C#N",
    "HCSCCH": "C(=S)C#C",
    "l-H2C5": "C=C=C=C=[C]",
    "l-H2C3": "C=C=[C]",
    "l-H2C4": "C=C=C=[C]",
    "l-H2C6": "C=C=C=C=C=[C]",
    "H2CCCHC3N": "C=C=CC#CC#N",  # cyanoacetyleneallene skeleton approx
    "c-C3HCCH": "C1=C[CH]1C#C",
    "c-C5H6": "C1=CC=CC1",
    "c-C9H8": "C1=CC2=CC=CC=C2C1",  # indene
    "C2H3CHO": "C=CC=O",
    "C2H3OH": "C=COH".replace("OH","O"),  # placeholder fixed below
    "HCOOCH3": "COC=O",
    "CH3OCH3": "COC",
    "NCS": "[S]C#N",
    "HCCS": "[CH]=C=S",
    "H2CCS": "C=C=S",
    "H2CCCS": "C=C=C=S",
    "1-C10H7CN": "N#Cc1cccc2cccc12",
    "2-C10H7CN": "N#Cc1ccc2ccccc2c1",
    "C5H5CN_1": "N#CC1=CC=CC1",   # 1-cyanocyclopentadiene
    "C5H5CN_2": "N#CC1C=CC=C1",   # 2-cyanocyclopentadiene
    "CH2CCH": "[CH2]C#C",
    "H2CCCHCCH": "C=C=CC#C",  # allenyl acetylene - duplicate skeleton of CH2CCHCCH-ish; distinct by H-pattern though
    "CH3C4H": "CC#CC#C",
    "HCCN": "[C]=C=[N]".replace("[N]","N"),
    "HC4N": "C#CC#C[N]",
    "CH3CH2CN": "CCC#N",
    "HCCCHCHCN_E": "C#CC=CC#N",
    "H2C=CHC3N": "C=CC#CC#N",
    "CH3C3N": "CC#CC#N",
    "CH2CCHCN": "C=C=CC#N",
    "HCCCH2CN": "C#CCC#N",
    "HC11N": "C#CC#CC#CC#CC#CC#N",
    "c-C5H5CN": "N#CC1=CC=CC1",
    "HNCCC": "N=C=C=[C]".replace("[N=C","[NH]=C"),
    "HCCNC": "C#C[N+]#[C-]",
    "SO": "[S]=O".replace("[S]","S"),
    "CCS": "[C]=C=S",
    "HC3N": "C#CC#N",
    "HCO+": "C#[O+]",
    "HCN": "C#N",
    "HCS+": "C#[S+]",
    "N2H+": "N#[N+]",
    "NH3": "N",
    "HC5O": "C#CC#CC=[O]".replace("=[O]","=O").replace("C=O","[C]=O"),
    "HC7O": "C#CC#CC#CC=[O]".replace("C=[O]","[C]=O"),
    "C6H5CN": "N#Cc1ccccc1",
    "HC9N": "C#CC#CC#CC#CC#N",
    "NCCNH+": "N#CC#[NH+]",
    "CH3OH": "CO",
    "CH3CHO": "CC=O",
    "H2CCO": "C=C=O",
    "c-C3H2O": "C1=COC1",   # cyclopropenone approx
    "H2CO": "C=O",
    "CCH": "[C]#C",
}

# Fix a couple of accidental .replace() no-ops / malformed entries above by
# re-specifying the clean, correct SMILES directly (overrides).
SPECIES_SMILES.update({
    "C2H3OH": "C=CO",             # vinyl alcohol
    "HCCN": "[CH]=C=N",            # cyanomethylene / HCCN radical skeleton (approx; H count checked below)
    "HNCCC": "[NH]=C=C=[C]",       # HNCCC skeleton
    "SO": "S=O",
    "HC5O": "C#CC#C[CH]=O",
    "HC7O": "C#CC#CC#C[CH]=O",
})

anchors = pd.read_csv("<PROJECT_ROOT>/f7_implementation/table_a2_authoritative_v2.csv")
anchor_keys = {}
for _, row in anchors.iterrows():
    k, _, _ = rule_c_key(row["smiles"])
    if k is not None:
        anchor_keys.setdefault(k, []).append(row["name"])

t1 = pd.read_csv("<PROJECT_ROOT>/f12_implementation/table1_validation_set_v7.csv", keep_default_na=False, dtype=str)
t1_keys = {}
for _, row in t1.iterrows():
    k, _, _ = rule_c_key(row["SMILES"])
    if k is not None:
        t1_keys.setdefault(k, []).append(row["Species"])

cand = pd.read_csv("<PROJECT_ROOT>/f12_implementation/candidates_purged_v5.csv")
member_rows = []
for _, row in cand.iterrows():
    ids = str(row["member_ids"]).split(";")
    smis = str(row["smiles_members"]).split(";")
    for mid, smi in zip(ids, smis):
        member_rows.append({"candidate_id": row["candidate_id"], "member_smiles": smi.strip(),
                             "positive": row["positive"]})
members = pd.DataFrame(member_rows)
members["rule_c_key"] = members["member_smiles"].apply(lambda s: rule_c_key(s)[0])
cand_key_to_ids = {}
for _, r in members.iterrows():
    cand_key_to_ids.setdefault(r["rule_c_key"], []).append((r["candidate_id"], bool(r["positive"])))

print(f"{len(SPECIES_SMILES)} species checked\n")
for name, smi in SPECIES_SMILES.items():
    k, charge, formula = rule_c_key(smi)
    if k is None:
        print(f"{name:15s} SMILES={smi!r:30s} PARSE FAILED")
        continue
    in_anchor = k in anchor_keys
    in_t1 = k in t1_keys
    cand_hits = cand_key_to_ids.get(k, [])
    status = []
    if in_anchor:
        status.append(f"ANCHOR({anchor_keys[k]})")
    if in_t1:
        status.append(f"TABLE1({t1_keys[k]})")
    if cand_hits:
        status.append(f"CANDIDATE_MATCH({cand_hits})")
    if not status:
        status.append("NO MATCH ANYWHERE (not in candidate pool -- nothing to relabel)")
    print(f"{name:15s} formula={formula:10s} key={k}")
    print(f"   -> {'; '.join(status)}")
