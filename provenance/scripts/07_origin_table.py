# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Origin table: for each origin 2022-01-01, 2023-01-01, 2024-01-01, 2025-01-01 (protocol_v1's origins), report:
  - in-list detections KNOWN before the origin (amendment A3: includes detections dated before Aug-2021,
    e.g. C3N-/C5N-, not just Lee's 87 anchors)
  - test positives inside that origin's window [origin_k, origin_k+1), the last window ending 2026-09-21
  - of those test positives, how many have low/unresolved date_confidence, and how many lie within 60 days
    of a window boundary (amendment A4 sensitivity flags)
No modeling -- this is a data-preparation summary table.
"""
import pandas as pd
import re
from datetime import datetime, timedelta

v3_1 = pd.read_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3_1.csv")

def parse_date(fpd):
    """Best-effort parse. Returns (datetime_or_None, precision) where
    precision in {'day','month','year',None}."""
    if not isinstance(fpd, str) or not fpd:
        return None, None
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})", fpd)
    if m:
        try:
            return datetime(int(m.group(1)), int(m.group(2)), int(m.group(3))), "day"
        except ValueError:
            pass
    m = re.match(r"(\d{4})-(\d{2})\b", fpd)
    if m:
        return datetime(int(m.group(1)), int(m.group(2)), 15), "month"
    m = re.match(r"(\d{4})\s*\(year only", fpd)
    if m:
        return datetime(int(m.group(1)), 7, 1), "year"
    return None, None

v3_1["parsed_date"], v3_1["date_precision"] = zip(*v3_1["first_public_date"].apply(parse_date))

ORIGINS = [datetime(2022, 1, 1), datetime(2023, 1, 1), datetime(2024, 1, 1), datetime(2025, 1, 1)]
FREEZE_DATE = datetime(2026, 9, 21)

rows = []
for i, origin in enumerate(ORIGINS):
    window_end = ORIGINS[i + 1] if i + 1 < len(ORIGINS) else FREEZE_DATE
    known_before = v3_1[(v3_1["in_lee_list"] == True) & (v3_1["parsed_date"].notna()) &
                        (v3_1["parsed_date"] < origin)]
    in_window = v3_1[(v3_1["in_lee_list"] == True) & (v3_1["parsed_date"].notna()) &
                      (v3_1["parsed_date"] >= origin) & (v3_1["parsed_date"] < window_end)]
    low_or_unresolved = in_window[in_window["date_confidence"].isin(["low", "unresolved"])]
    near_start = in_window[abs(in_window["parsed_date"] - origin) <= timedelta(days=60)]
    near_end = in_window[abs(in_window["parsed_date"] - window_end) <= timedelta(days=60)]
    near_boundary = pd.concat([near_start, near_end]).drop_duplicates(subset="formula")
    rows.append({
        "origin": origin.strftime("%Y-%m-%d"),
        "window_end": window_end.strftime("%Y-%m-%d"),
        "in_list_known_before_origin": len(known_before),
        "test_positives_in_window": len(in_window),
        "test_positives_low_or_unresolved_date_confidence": len(low_or_unresolved),
        "test_positives_within_60d_of_boundary": len(near_boundary),
        "test_positive_formulas": ";".join(sorted(in_window["formula"])),
        "low_or_unresolved_formulas": ";".join(sorted(low_or_unresolved["formula"])),
        "near_boundary_formulas": ";".join(sorted(near_boundary["formula"])),
    })

origin_df = pd.DataFrame(rows)
origin_df.to_csv("<PROJECT_ROOT>/f3_validation_freeze/origin_table.csv", index=False)

# no-date-at-all check (shouldn't happen given the earlier fallbacks, but verify)
no_date = v3_1[v3_1["parsed_date"].isna() & (v3_1["in_lee_list"] == True)]
print("in-list detections with NO parseable date at all:", len(no_date))
if len(no_date):
    print(no_date[["formula", "first_public_date"]].to_string())

# first-origin drop check (protocol: drop first origin if <12 in-list known then)
first_known = rows[0]["in_list_known_before_origin"]
print(f"\nIn-list detections known before first origin (2022-01-01): {first_known} "
      f"({'>= 12, KEEP first origin' if first_known >= 12 else '< 12, DROP first origin per protocol'})")

print(origin_df.drop(columns=["test_positive_formulas", "low_or_unresolved_formulas",
                                "near_boundary_formulas"]).to_string(index=False))

v3_1.to_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3_1.csv", index=False)
