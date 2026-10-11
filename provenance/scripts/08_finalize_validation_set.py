# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
import pandas as pd
import re

v3 = pd.read_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3.csv")
match = pd.read_csv("<PROJECT_ROOT>/f3_validation_freeze/v3_vs_lee1510_match.csv")

v3 = v3.merge(match[["formula", "level_A", "level_B", "level_C", "in_lee_list"]],
              on="formula", how="left")

def extract_year_half(fpd):
    """Best-effort half-year bucket from the messy first_public_date strings."""
    if not isinstance(fpd, str) or not fpd:
        return "UNKNOWN"
    m = re.match(r"(\d{4})-(\d{2})", fpd)
    if m:
        year, month = int(m.group(1)), int(m.group(2))
        half = "H1" if month <= 6 else "H2"
        return f"{year}-{half}"
    m = re.match(r"(\d{4})\s*\(year only", fpd)
    if m:
        return f"{m.group(1)}-UNKNOWN_HALF"
    return "UNKNOWN"

v3["half_year"] = v3["first_public_date"].apply(extract_year_half)

v3.to_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3.csv", index=False)

# dates_evidence.csv -- the focused date-provenance deliverable
dates_ev = v3[["formula", "first_public_date", "first_public_date_basis",
               "arxiv_v1_date", "journal_received", "journal_accepted",
               "journal_online", "arxiv_id", "doi", "date_confidence",
               "half_year", "flags"]].copy()
dates_ev.to_csv("<PROJECT_ROOT>/f3_validation_freeze/dates_evidence.csv", index=False)

# Step 5: counts per half-year, in-list vs out-of-list
tbl = v3.groupby(["half_year", "in_lee_list"]).size().unstack(fill_value=0)
tbl.columns = ["out_of_lee_list", "in_lee_list"] if False else tbl.columns
tbl = tbl.rename(columns={True: "in_lee_list", False: "out_of_lee_list"})
tbl["total"] = tbl.sum(axis=1)
tbl = tbl.sort_index()
tbl.to_csv("<PROJECT_ROOT>/f3_validation_freeze/step5_halfyear_counts.csv")

print(tbl.to_string())
print(f"\nTotal: in_lee_list={v3['in_lee_list'].sum()}, "
      f"out_of_lee_list={(~v3['in_lee_list'].astype(bool)).sum()}")
print(f"\nUNKNOWN half_year bucket rows (need attention):")
print(v3.loc[v3['half_year'].str.contains('UNKNOWN'), ['formula', 'first_public_date', 'date_confidence']].to_string())
