# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Detections v3.1 (77 rows) -> v3.2; v3.1 is not overwritten.
Apply protocol amendment 2's A7 date corrections, which supersede v3.1's previous "year only" / unresolved
dates for three detections (all had DOI-only citations with no arXiv id located in the original
research pass):
  C3H+, C5H+ (Cernicharo et al. 2022a) -> arXiv:2201.02434, first-public-date 2022-01 (month-level, per
    A7's explicit allowance: "Month-level arXiv-ID dates suffice for annual origins")
  H2CCCH+ (Silva et al. 2023) -> arXiv:2307.01733, first-public-date 2023-07 (month-level)
No other rows touched.
"""
import pandas as pd

v3_1 = pd.read_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3_1.csv")
df = v3_1.copy()
for c in ["arxiv_id", "arxiv_v1_date", "first_public_date", "first_public_date_basis",
          "date_confidence", "flags"]:
    df[c] = df[c].astype("object")

A7_CORRECTIONS = {
    "C3H+": dict(arxiv_id="2201.02434", first_public_date="2022-01",
                 first_public_date_basis="arxiv_id_month_only (amendments_2 A7)"),
    "C5H+": dict(arxiv_id="2201.02434", first_public_date="2022-01",
                 first_public_date_basis="arxiv_id_month_only (amendments_2 A7)"),
    "H2CCCH+": dict(arxiv_id="2307.01733", first_public_date="2023-07",
                     first_public_date_basis="arxiv_id_month_only (amendments_2 A7)"),
}

change_log_rows = []
for formula, corr in A7_CORRECTIONS.items():
    idx = df.index[df["formula"] == formula]
    assert len(idx) == 1, f"{formula} not found or duplicated in v3.1"
    i = idx[0]
    old_fpd = df.loc[i, "first_public_date"]
    old_arxiv = df.loc[i, "arxiv_id"]
    old_conf = df.loc[i, "date_confidence"]
    df.loc[i, "arxiv_id"] = corr["arxiv_id"]
    df.loc[i, "arxiv_v1_date"] = corr["first_public_date"]
    df.loc[i, "first_public_date"] = corr["first_public_date"]
    df.loc[i, "first_public_date_basis"] = corr["first_public_date_basis"]
    df.loc[i, "date_confidence"] = "medium"  # month-level arXiv ID, not a day-level read
    df.loc[i, "flags"] = (str(df.loc[i, "flags"]) if pd.notna(df.loc[i, "flags"]) else "") + \
        (f" [F3c A7: date/arXiv id superseded per protocol_v1_amendments_2.md -- "
         f"was '{old_fpd}' (unresolved, no arXiv id), now arXiv:{corr['arxiv_id']} "
         f"({corr['first_public_date']}).]")
    change_log_rows.append({
        "formula": formula, "action": "CORRECT_DATE_AND_ARXIV_ID",
        "old_first_public_date": old_fpd, "old_arxiv_id": old_arxiv,
        "old_date_confidence": old_conf,
        "new_first_public_date": corr["first_public_date"],
        "new_arxiv_id": corr["arxiv_id"], "new_date_confidence": "medium",
        "evidence": f"protocol_v1_amendments_2.md A7 (Reham-supplied): arXiv:{corr['arxiv_id']}",
    })

df.to_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3_2.csv", index=False)
pd.DataFrame(change_log_rows).to_csv(
    "<PROJECT_ROOT>/f3_validation_freeze/v3_2_change_log.csv", index=False)

print(f"v3.2 written: {len(df)} rows (unchanged count from v3.1, only 3 rows' "
      f"date/arxiv_id/date_confidence fields updated)")
print(pd.DataFrame(change_log_rows)[["formula", "old_first_public_date", "new_first_public_date",
                                       "new_arxiv_id"]].to_string(index=False))
