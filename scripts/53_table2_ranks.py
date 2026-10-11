"""
Table 2 ("would we have caught it?" rank table), released regime (data/frozen/candidates.csv):
one row per test positive (29, summed over the 4 origins), rank under Combined ranking, Detectability
and Similarity (low_dist) among that origin's test pool (ties: mean rank). Species names are resolved to
the exact strings of tables/table1_validation_set.csv (keyed by LeeID/member_ids).
"""
import os
import sys
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from paths import TABLES
from evaluation_engine import known_mask_at, ORIGINS, FREEZE_DATE
from regimes import build_released_v7

OUT_DIR = TABLES

cand, Xs, per_origin = build_released_v7()
dates = cand["date"].to_numpy()
positive = cand["positive"].to_numpy()

# species-name resolution: tables/table1_validation_set.csv, keyed by LeeID (candidate's own
# member_ids) -- a candidate with n_members > 1 resolves via its member ids jointly.
t1 = pd.read_csv(os.path.join(TABLES, "table1_validation_set.csv"),
                  dtype=str, keep_default_na=False)
leeid_to_species = {}
for _, r in t1.iterrows():
    leeid = r["LeeID"].strip()
    if not leeid:
        continue
    for part in leeid.split(";"):
        leeid_to_species.setdefault(part.strip(), []).append(r["Species"])


def resolve_species_name(cand_row):
    member_ids = [str(int(x)) for x in str(cand_row["member_ids"]).split(";")]
    names = []
    for mid in member_ids:
        names.extend(leeid_to_species.get(mid, []))
    names = sorted(set(names))
    if names:
        return "/".join(names)
    return cand_row["formulas"]  # fallback: raw candidate formula, if no Table 1 match found


rows = []
unresolved = []
for oi, o in enumerate(per_origin):
    origin = ORIGINS[oi]
    next_o = ORIGINS[oi + 1] if oi + 1 < len(ORIGINS) else FREEZE_DATE

    known = known_mask_at(dates, origin) & positive
    test_pool = ~known
    known_next = known_mask_at(dates, next_o) & positive
    test_positive = test_pool & known_next
    idx_test = np.where(test_pool)[0]
    y_test = test_positive[idx_test]

    assert len(y_test) == o["test_pool_size"]
    assert int(y_test.sum()) == o["test_positives"]
    assert np.array_equal(y_test.astype(bool), o["_test_positive_bool"])

    combined_score = o["_combined_score"]
    detect_score = o["_raw_scores_test"]["detectability"]
    lowdist_score = o["_raw_scores_test"]["low_dist"]

    rank_combined = rankdata(-combined_score, method="average")
    rank_detect = rankdata(-detect_score, method="average")
    rank_lowdist = rankdata(-lowdist_score, method="average")

    test_pos_local = np.where(y_test)[0]
    for local_i in test_pos_local:
        global_i = idx_test[local_i]
        cand_row = cand.iloc[global_i]
        species = resolve_species_name(cand_row)
        if species == cand_row["formulas"] and str(cand_row["member_ids"]).split(";")[0] not in leeid_to_species:
            unresolved.append((cand_row["candidate_id"], cand_row["formulas"]))
        rows.append({
            "species": species,
            "species_key": cand_row["species_key"],
            "candidate_id": int(cand_row["candidate_id"]),
            "origin": o["origin"],
            "test_pool_size": o["test_pool_size"],
            "rank_combined": float(rank_combined[local_i]),
            "rank_detectability": float(rank_detect[local_i]),
            "rank_similarity": float(rank_lowdist[local_i]),
        })

table2_df = pd.DataFrame(rows, columns=[
    "species", "species_key", "candidate_id", "origin", "test_pool_size",
    "rank_combined", "rank_detectability", "rank_similarity",
])
assert len(table2_df) == 29, f"expected 29 test positives, got {len(table2_df)}"
if unresolved:
    print(f"WARNING: {len(unresolved)} rows fell back to raw candidate formula "
          f"(no Table 1 v9 species-name match found): {unresolved}")
out_path = os.path.join(OUT_DIR, "table2_ranks.csv")
table2_df.to_csv(out_path, index=False)
print(f"Wrote {out_path} ({len(table2_df)} rows)")
print(table2_df.to_string(index=False))
