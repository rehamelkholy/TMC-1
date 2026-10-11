# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
lee_list_entry_labels_purged_v5.csv (1488 rows) ->
lee_list_entry_labels_purged_v6.csv (1483 rows).

Removes the 5 raw-universe entries that are the sole members of the 5
anchor-duplicate candidates purged in item 2/3 (rule C: connectivity/
charge-placement variants of 5 of the 87 Lee anchors), exactly as the v5 update's
own HCCN-charge-duplicate removal (raw id 378) was handled in
step2_build_candidates_purged_v5.py:
  raw id 302 -- candidate_id 293 (H2COH+/CH3O+ anchor dup)
  raw id 389 -- candidate_id 376 (HNC3/C3HN anchor dup)
  raw id 411 -- candidate_id 396 (HOCO+/CHO2+ anchor dup)
  raw id 509 -- candidate_id 490 (HC3NH+/C3H2N+ anchor dup)
  raw id 544 -- candidate_id 523 (NCCNH+/C2HN2+ anchor dup)

All 5 are singleton candidates (n_members=1), confirmed via
candidates_purged_v6.csv showing 0 rows for all 5 candidate_ids -- so this
removal does not touch any other candidate's membership, matching the
same "verified safe to patch in place" logic as the v4->v5 step.

Does not touch lee_list_entry_labels_purged_v5.csv. This file is the
load_purged_raw() input for the v6 leave-self-out feature rebuild,
in place of the v5 version.
"""
import pandas as pd

v5 = pd.read_csv("<PROJECT_ROOT>/f12_implementation/lee_list_entry_labels_purged_v5.csv")
assert len(v5) == 1488

REMOVE_RAW_IDS = [302, 389, 411, 509, 544]

for rid, expected_smiles in [
    (302, "[CH2+]O"), (389, "[C-]#C[C+]=N"), (411, "O=[C+]O"),
    (509, "C#C[C+]=N"), (544, "N#C[C+]=N"),
]:
    row = v5[v5["id"] == rid].iloc[0]
    assert row["smiles"] == expected_smiles, (rid, row["smiles"])

v6 = v5[~v5["id"].isin(REMOVE_RAW_IDS)].copy()
assert len(v6) == 1488 - len(REMOVE_RAW_IDS) == 1483

out_path = "<PROJECT_ROOT>/f13_implementation/lee_list_entry_labels_purged_v6.csv"
v6.to_csv(out_path, index=False)
print(f"lee_list_entry_labels_purged_v6.csv: {len(v6)} rows (was {len(v5)}, removed "
      f"{len(REMOVE_RAW_IDS)} anchor-dup raw ids: {REMOVE_RAW_IDS})")
print(f"Wrote {out_path}")
