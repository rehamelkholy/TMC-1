"""
Sensitivity A4(b) as an explicit date shift (protocol amendment 4, J2). The shift is applied once
to produce a shifted-date candidate table, which is then run through the unchanged
evaluation_engine.run_origin_evaluation(..., sensitivity=None), so the A4(b) re-run reuses the
primary evaluation code path with different dates.
"""
import pandas as pd
from evaluation_engine import ORIGINS  # unchanged: [2022-01-01, 2023-01-01, 2024-01-01, 2025-01-01]

BASELINE_START = pd.Timestamp(2022, 1, 1)


def apply_j2_shift(cand):
    """Returns a COPY of cand with the 'date' column replaced per J2. Only
    positive candidates with a real date >= 2022-01-01 are eligible to move;
    a candidate within 60 days of exactly one origin boundary is moved to
    the other side of it (origin - 1 day if it was on/after; origin + 1 day
    if it was before)."""
    cand2 = cand.copy()
    new_dates = []
    shift_log = []
    for _, row in cand2.iterrows():
        d = row["date"]
        if row["positive"] and pd.notna(d) and d >= BASELINE_START:
            moved = False
            for o in ORIGINS:
                delta_days = (d - o).days
                if abs(delta_days) <= 60:
                    if d >= o:
                        new_d = o - pd.Timedelta(days=1)
                    else:
                        new_d = o + pd.Timedelta(days=1)
                    new_dates.append(new_d)
                    shift_log.append({
                        "candidate_id": row["candidate_id"], "formulas": row["formulas"],
                        "orig_date": d.strftime("%Y-%m-%d"), "origin": o.strftime("%Y-%m-%d"),
                        "new_date": new_d.strftime("%Y-%m-%d"),
                    })
                    moved = True
                    break
            if not moved:
                new_dates.append(d)
        else:
            new_dates.append(d)
    cand2["date"] = new_dates
    return cand2, shift_log


if __name__ == "__main__":
    from candidates import build_candidates_v2
    cand, _ = build_candidates_v2()
    cand2, log = apply_j2_shift(cand)
    print("candidates moved under J2:", len(log))
    for l in log:
        print(" ", l)
