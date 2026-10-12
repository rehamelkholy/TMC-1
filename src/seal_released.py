"""
Primary evaluation (released regime) under protocol amendments 9 (R1-R3: exact random
expectations, exact tie expectation, sha256 seeds).

Output path is argv[1]; the non-deterministic "generated"/"elapsed_seconds" fields are
removed so repeated runs are byte-identical.

Evaluation logic (evaluation_engine, bootstrap, feedback_learner, boundary_shift) is imported
unchanged. The candidate table is data/frozen/candidates.csv (1,434 rows); date bases are read
from data/frozen/post2021_detections.csv. The "primary" features (dist/logcd/sigma) are
Lee's released values as of 2021-07, baked into the candidate table itself and invariant to
which detections/labels are attached (protocol_v1's definition of "primary"); only the
label/date attached to each candidate_id, and which candidate_ids exist, depend on the
label version.
"""
import hashlib
import json
import os
import sys
import time
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from paths import DATA_FROZEN
from evaluation_engine import (  # unchanged
    standardize_features, run_origin_evaluation, pool_metrics,
    best_gating_baseline, FEATURE_COLS, TOPK_LIST, ORIGINS, BASE_SEED,
    N_TIEBREAK_TRIALS, GATING_BASELINES, NONGATING, raw_scores,
)
from bootstrap import run_bootstrap, N_BOOT, BOOT_SEED  # unchanged
from feedback_learner import run_v3_feedback, obs_to_target_avg, N_SEEDS  # unchanged
from boundary_shift import apply_j2_shift  # unchanged


def load_candidates_purged_v7():
    cand = pd.read_csv(os.path.join(DATA_FROZEN, "candidates.csv"))
    cand["date"] = pd.to_datetime(cand["date"])
    lee = pd.read_csv(os.path.join(DATA_FROZEN, "lee_list_entry_labels_purged.csv"))
    cand.attrs["n_raw_entries"] = len(lee)
    return cand, lee


# Date bases for the six positives added at the final label freeze (all outside Lee's v3.1 list,
# so the list-based lookup below cannot find them): taken from their own rows in
# post2021_detections.csv. The ten list-external positives added at v6 keep basis None.
F14_BASIS_FALLBACK = {"1305": "journal_online_day", "522": "journal_month_only",
                       "498": "arxiv_v1", "318": "arxiv_v1", "512": "arxiv_v1", "1312": "arxiv_v1"}


def attach_date_basis_v7(cand, lee):
    v6 = pd.read_csv(os.path.join(DATA_FROZEN, "post2021_detections.csv"))
    basis_by_formula = dict(zip(v6["formula"], v6["first_public_date_basis"]))
    date_basis = []
    for _, row in cand.iterrows():
        if not row["positive"]:
            date_basis.append(None)
            continue
        member_ids = [int(x) for x in row["member_ids"].split(";")]
        sub = lee[lee["id"].isin(member_ids) & (lee["in_list_by_v3_1"] == True)]
        candidates_detections = []
        for _, srow in sub.iterrows():
            formulas = str(srow["matched_detection_formula"]).split(";")
            for f in formulas:
                if f in basis_by_formula:
                    candidates_detections.append(f)
        bases = sorted(set(basis_by_formula[f] for f in candidates_detections))
        if not bases and str(row["candidate_id"]) in F14_BASIS_FALLBACK:
            bases = [F14_BASIS_FALLBACK[str(row["candidate_id"])]]
        if len(bases) > 1:
            date_basis.append("AMBIGUOUS:" + ";".join(bases))
        else:
            date_basis.append(bases[0] if bases else None)
    cand = cand.copy()
    cand["date_basis"] = date_basis
    return cand


def clean_origin_for_json(o):
    return {k: v for k, v in o.items() if not k.startswith("_")}


def main():
    t_start = time.time()

    cand, lee = load_candidates_purged_v7()
    cand = attach_date_basis_v7(cand, lee)
    Xs, mean, std = standardize_features(cand)
    y_pos = cand["positive"].to_numpy()

    data_summary = {
        "n_raw_lee_entries_purged": int(cand.attrs.get("n_raw_entries", 1483)),
        "n_raw_lee_entries_before_purge": 1501,
        "n_anchor_matched_rows_purged": 18,
        "n_candidates_after_species_key_collapse": len(cand),
        "n_positive_candidates": int(cand["positive"].sum()),
        "feature_cols": FEATURE_COLS,
        "feature_means": {c: float(m) for c, m in zip(FEATURE_COLS, mean)},
        "feature_stds": {c: float(s) for c, s in zip(FEATURE_COLS, std)},
        "v7_vs_v6_positive_delta": int(cand["positive"].sum()) - 57,
        "v7_vs_v6_candidate_count_delta": len(cand) - 1434,
    }

    per_origin_primary = run_origin_evaluation(cand, Xs, sensitivity=None)
    pooled_primary = pool_metrics(per_origin_primary)
    best_baseline = best_gating_baseline(pooled_primary)
    for o in per_origin_primary:
        o["_baseline_score"] = o["_raw_scores_test"][best_baseline]
    bootstrap_res = run_bootstrap(per_origin_primary, best_baseline, n_boot=N_BOOT, seed=BOOT_SEED)
    primary_out = {
        "per_origin": [clean_origin_for_json(o) for o in per_origin_primary],
        "pooled": pooled_primary, "best_gating_baseline": best_baseline,
        "bootstrap_combined_vs_best_baseline": bootstrap_res,
    }

    per_origin_a4a = run_origin_evaluation(cand, Xs, sensitivity="A4a")
    pooled_a4a = pool_metrics(per_origin_a4a)
    best_a4a = best_gating_baseline(pooled_a4a)
    sens_a4a_out = {"per_origin": [clean_origin_for_json(o) for o in per_origin_a4a],
                     "pooled": pooled_a4a, "best_gating_baseline": best_a4a}

    cand_j2, j2_shift_log = apply_j2_shift(cand)
    per_origin_a4b = run_origin_evaluation(cand_j2, Xs, sensitivity=None)
    pooled_a4b = pool_metrics(per_origin_a4b)
    best_a4b = best_gating_baseline(pooled_a4b)
    sens_a4b_out = {"per_origin": [clean_origin_for_json(o) for o in per_origin_a4b],
                     "pooled": pooled_a4b, "best_gating_baseline": best_a4b,
                     "j2_shift_log": j2_shift_log}

    v3_res = run_v3_feedback(Xs, y_pos, n_seeds=N_SEEDS)
    raw = raw_scores(cand)
    fixed_orderings = {}
    for m in GATING_BASELINES:
        avg, _ = obs_to_target_avg(raw[m], y_pos, v3_res["target_positives"],
                                    N_TIEBREAK_TRIALS, seed=BASE_SEED)
        fixed_orderings[m] = avg
    best_fixed = min(fixed_orderings, key=lambda m: fixed_orderings[m])
    v3_out = {
        **v3_res, "fixed_ordering_expected_obs_to_50pct": fixed_orderings,
        "best_fixed_ordering": best_fixed, "best_fixed_ordering_value": fixed_orderings[best_fixed],
        "share_of_seeds_learner_beats_best_fixed": float(
            np.mean(np.array(v3_res["obs_to_50pct_per_seed"]) < fixed_orderings[best_fixed])),
    }

    elapsed = time.time() - t_start
    results = {
        "task": "F15 (Amendments 9 re-seal) of F14 (F6/F9/F12c primary rerun on v7)",
        "run_config": {
            "base_seed": BASE_SEED, "n_tiebreak_trials": N_TIEBREAK_TRIALS,
            "n_bootstrap": N_BOOT, "bootstrap_seed": BOOT_SEED, "n_v3_seeds": N_SEEDS,
            "topk_list": TOPK_LIST, "origins": [o.strftime("%Y-%m-%d") for o in ORIGINS],
        },
        "implementation_notes": [
            "Identical to f12_implementation/seal_results_f6_v5.py's evaluation logic "
            "(step4_engine.py, step5_bootstrap.py, step6_v3_feedback.py, j2_shift.py, all "
            "unchanged). Only the candidate table and date-basis source changed: "
            "post2021_detections_v7.csv / candidates_purged_v7.csv (this task's item 4 "
            "output) in place of the v5 files. candidates_purged_v7.csv = "
            "candidates_purged_v5.csv with 5 anchor-duplicate candidates (293/376/396/490/"
            "523) dropped and 9 candidate_ids (268/272/282/298/383/484/524/535/765) "
            "relabelled positive; every other candidate_id, and every dist/logcd/sigma "
            "value, is byte-identical to v5.",
        ],
        "data_summary": data_summary,
        "primary": primary_out,
        "sensitivity_A4a_non_arxiv_dates_dropped": sens_a4a_out,
        "sensitivity_A4b_60day_boundary_flip": sens_a4b_out,
        "v3_feedback_learner": v3_out,
    }

    out_path = os.path.abspath(sys.argv[1])
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating,)) else
                   (int(x) if isinstance(x, (np.integer,)) else str(x)))
    with open(out_path, "rb") as f:
        h = hashlib.sha256(f.read()).hexdigest()
    with open(out_path + ".sha256", "w") as f:
        f.write(h + "\n")
    print(os.path.basename(out_path), "written:", out_path)
    print("SHA-256:", h)
    print("elapsed:", round(elapsed, 1), "s")
    print("known_before_origin/test_pool_size/test_positives per origin:")
    for o in per_origin_primary:
        print(" ", o["origin"], o["known_before_origin"], o["test_pool_size"], o["test_positives"])


if __name__ == "__main__":
    main()
