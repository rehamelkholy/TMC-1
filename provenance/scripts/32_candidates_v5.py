# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
candidates_purged_v4.csv (1440 rows) -> candidates_purged_v5.csv (1439 rows).

NOT a from-scratch rebuild via 22_candidates_v4.py's groupby/np.arange(len(cand)) machinery: candidate_id
there is assigned POSITIONALLY after grouping, so purging raw id 378 (candidate_id 366's sole member) from
scratch would silently renumber every candidate_id after position 366, invalidating every fixed-ID
reference used throughout (171, 651, 263, 113, 891, the rankings, Table 1's LeeID column, etc.).

Safe to patch in place instead: candidate_id 366 (species_key C2HN|C2HN|c1-2-3|h1H) is a SINGLETON group,
n_members=1, member_ids="378" only -- so removing raw id 378 does not touch any other candidate's
membership, and removing one row does not reorder the remaining rows a from-scratch groupby(sort=False)
would iterate over, so every other candidate_id is provably unchanged either way. The four relabelled
candidates keep every column except positive/date, which are set to exactly what 22_candidates_v4.py's own
logic would produce (parse_one_date's month-only rule: day set to 1) from post2021_detections_v5.csv /
lee_list_entry_labels_v5.csv.

Actions:
  1. DROP candidate_id 366 (C2HN, [CH+]=C=[N-]) -- anchor charge-audit item 1 (HCCN): training-inventory
     duplicate of the charge-corrected HCCN anchor, removed exactly like the 21 previously-removed anchor
     duplicates (RULE_A/RULE_B), not relabelled.
  2. RELABEL positive=True, date=<matched detection's date> for:
       candidate_id 171 (C5N)    -> date 1998-05-01
       candidate_id 651 (C4HN)   -> date 2021-02-01
       candidate_id 263 (NS+)    -> date 2018-01-01
       candidate_id 113 (C3HS+)  -> date 2021-01-01
     All four dated well before every test origin (earliest 2022-01-01) -> known at every origin per A3,
     never a test positive.
"""
import pandas as pd

v4 = pd.read_csv("<PROJECT_ROOT>/f9_implementation/candidates_purged_v4.csv")
assert len(v4) == 1440

RELABEL = {
    171: "1998-05-01",  # C5N
    651: "2021-02-01",  # C4HN (l-HC4N)
    263: "2018-01-01",  # NS+
    113: "2021-01-01",  # C3HS+ (HC3S+)
}
DROP_ID = 366  # C2HN, [CH+]=C=[N-] -- HCCN charge-encoding duplicate

# sanity: confirm the pre-patch state matches the expected values
for cid, expected_smiles, expected_charged in [
    (171, "[C]#CC#CC#N", False),
    (651, "C#C[C+]=C=[N-]", False),
    (263, "N#[S+]", True),
    (113, "C#CC#[S+]", True),
    (366, "[CH+]=C=[N-]", False),
]:
    row = v4[v4["candidate_id"] == cid].iloc[0]
    assert row["smiles_members"] == expected_smiles, (cid, row["smiles_members"])
    assert bool(row["charged"]) == expected_charged, (cid, row["charged"])
    assert row["positive"] == False, (cid, "expected positive=False pre-patch")
    assert row["n_members"] == 1, (cid, "expected singleton")

v5 = v4[v4["candidate_id"] != DROP_ID].copy()
assert len(v5) == 1439

for cid, date_str in RELABEL.items():
    mask = v5["candidate_id"] == cid
    assert mask.sum() == 1
    v5.loc[mask, "positive"] = True
    v5.loc[mask, "date"] = date_str

n_positive = int(v5["positive"].sum())
print(f"candidates_purged_v5.csv: {len(v5)} rows (was {len(v4)}, dropped candidate_id {DROP_ID})")
print(f"Positive candidates: {n_positive} (was {int(v4['positive'].sum())} in v4)")

out_path = "<PROJECT_ROOT>/f12_implementation/candidates_purged_v5.csv"
v5.to_csv(out_path, index=False)
print(f"Wrote {out_path}")

print(v5[v5["candidate_id"].isin(RELABEL.keys())][
    ["candidate_id", "species_key", "member_ids", "smiles_members", "positive", "date"]
].to_string(index=False))

# ---- also build lee_list_entry_labels_purged_v5.csv (the raw purged universe, for the reseal scripts'
# data_summary block: n_raw_lee_entries_purged etc.) ----
RULE_A_IDS = [289, 332, 368, 380, 398, 430, 432, 484, 679, 681, 793, 1302]
HCCN_DUP_RAW_ID = 378  # candidate 366's sole member (anchor charge-audit)

raw_v5 = pd.read_csv("<PROJECT_ROOT>/f12_implementation/lee_list_entry_labels_v5.csv")
assert len(raw_v5) == 1501
purged_raw_v5 = raw_v5.loc[~raw_v5["id"].isin(RULE_A_IDS + [HCCN_DUP_RAW_ID])].reset_index(drop=True)
assert len(raw_v5) - len(purged_raw_v5) == len(RULE_A_IDS) + 1
purged_out = "<PROJECT_ROOT>/f12_implementation/lee_list_entry_labels_purged_v5.csv"
purged_raw_v5.to_csv(purged_out, index=False)
print(f"\nWrote {purged_out}: {len(purged_raw_v5)} rows "
      f"(1501 - {len(RULE_A_IDS)} RULE_A - 1 HCCN-charge-dup [raw id {HCCN_DUP_RAW_ID}])")
