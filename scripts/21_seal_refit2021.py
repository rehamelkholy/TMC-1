"""
Refit-2021 regime (protocol amendment 8, N1) under protocol amendments 9 (R1-R3). Output path is
argv[1]; the non-deterministic "generated"/"elapsed_seconds" fields are removed so repeated runs
are byte-identical.

data/frozen/features_refit2021.csv is REUSED VERBATIM, not recomputed: the anchors-only fit (87
anchors + 643 shipped training points, no post-2021 detections) does not read positive/date, and
the candidate skeleton (data/frozen/candidates.csv, 1,434 rows) is cell-identical to earlier label
versions on every column the fit depends on (species_key, member_ids, dist, logcd, sigma, mu, mw,
charged, radical). Only the evaluation is re-run.

Additionally writes data/derived/features_refit2021_relabelled.csv: the per-origin, per-candidate
(dist, logcd, sigma, is_positive, date) table used for this evaluation, i.e. the unchanged
features_refit2021.csv feature values joined onto the candidate skeleton/labels.
"""
import os
import sys
import time
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import hashlib
import json
import numpy as np
import pandas as pd

from paths import DATA_FROZEN, DATA_DERIVED
from evaluation_engine import (
    standardize_features, run_origin_evaluation, pool_metrics,
    best_gating_baseline, FEATURE_COLS, TOPK_LIST, ORIGINS, BASE_SEED,
    N_TIEBREAK_TRIALS, GATING_BASELINES, NONGATING, raw_scores,
)
from bootstrap import run_bootstrap, N_BOOT, BOOT_SEED
from boundary_shift import apply_j2_shift
from seal_released import load_candidates_purged_v7, attach_date_basis_v7, clean_origin_for_json

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


features = pd.read_csv(os.path.join(DATA_FROZEN, "features_refit2021.csv"))  # REUSED, unchanged
cand_skel, lee = load_candidates_purged_v7()
cand_skel = attach_date_basis_v7(cand_skel, lee)
log(f"Loaded features_refit_2021.csv ({len(features)} rows, REUSED unchanged) and the v6 candidate "
    f"skeleton ({len(cand_skel)} candidates, {int(cand_skel['positive'].sum())} positive).")

cand_skel["date"] = pd.to_datetime(cand_skel["date"])
cand_j2_skel, j2_shift_log = apply_j2_shift(cand_skel)
log(f"Applied J2 date shift to the v6 skeleton once ({len(j2_shift_log)} candidates moved).")


def build_cand_for_origin(origin_str, base_skel):
    sub = features[features["origin"] == origin_str].set_index("candidate_id")
    cand = base_skel.copy()
    cand["dist"] = cand["candidate_id"].map(sub["dist"])
    cand["logcd"] = cand["candidate_id"].map(sub["logcd"])
    cand["sigma"] = cand["candidate_id"].map(sub["sigma"])
    n_nan = cand[["dist", "logcd", "sigma"]].isna().any(axis=1).sum()
    if n_nan:
        log(f"  WARNING: {n_nan} candidates have NaN features at origin {origin_str}")
    return cand


ORIGIN_STRS = [o.strftime("%Y-%m-%d") for o in ORIGINS]

per_origin_primary, per_origin_a4a, per_origin_a4b = [], [], []
relabelled_rows = []
for oi, ostr in enumerate(ORIGIN_STRS):
    log(f"--- Evaluating origin {oi} ({ostr}) ---")
    cand_i = build_cand_for_origin(ostr, cand_skel)
    Xs_i, mean_i, std_i = standardize_features(cand_i)

    res_primary = run_origin_evaluation(cand_i, Xs_i, sensitivity=None)[oi]
    assert res_primary["origin"] == ostr
    per_origin_primary.append(res_primary)

    res_a4a = run_origin_evaluation(cand_i, Xs_i, sensitivity="A4a")[oi]
    per_origin_a4a.append(res_a4a)

    cand_i_j2 = build_cand_for_origin(ostr, cand_j2_skel)
    Xs_i_j2, _, _ = standardize_features(cand_i_j2)
    res_a4b = run_origin_evaluation(cand_i_j2, Xs_i_j2, sensitivity=None)[oi]
    per_origin_a4b.append(res_a4b)

    log(f"  known_before_origin={res_primary['known_before_origin']}, "
        f"test_pool_size={res_primary['test_pool_size']}, test_positives={res_primary['test_positives']}")

    for _, r in cand_i.iterrows():
        relabelled_rows.append({
            "origin": ostr, "candidate_id": int(r["candidate_id"]), "species_key": r["species_key"],
            "is_positive": bool(r["positive"]),
            "date": r["date"].strftime("%Y-%m-%d") if pd.notna(r["date"]) else "",
            "dist": float(r["dist"]), "logcd": float(r["logcd"]), "sigma": float(r["sigma"]),
            "mu": float(r["mu"]), "mw": float(r["mw"]),
            "charged": bool(r["charged"]), "radical": bool(r["radical"]),
        })

relabelled_path = os.path.join(DATA_DERIVED, "features_refit2021_relabelled.csv")
pd.DataFrame(relabelled_rows).to_csv(relabelled_path, index=False)
log(f"Wrote {relabelled_path}: {len(relabelled_rows)} rows ({len(cand_skel)} candidates x {len(ORIGIN_STRS)} origins)")

pooled_primary = pool_metrics(per_origin_primary)
best_baseline = best_gating_baseline(pooled_primary)
for o in per_origin_primary:
    o["_baseline_score"] = o["_raw_scores_test"][best_baseline]
bootstrap_res = run_bootstrap(per_origin_primary, best_baseline, n_boot=N_BOOT, seed=BOOT_SEED)

pooled_a4a = pool_metrics(per_origin_a4a)
best_a4a = best_gating_baseline(pooled_a4a)
pooled_a4b = pool_metrics(per_origin_a4b)
best_a4b = best_gating_baseline(pooled_a4b)

primary_out = {"per_origin": [clean_origin_for_json(o) for o in per_origin_primary],
                "pooled": pooled_primary, "best_gating_baseline": best_baseline,
                "bootstrap_combined_vs_best_baseline": bootstrap_res}
sens_a4a_out = {"per_origin": [clean_origin_for_json(o) for o in per_origin_a4a],
                "pooled": pooled_a4a, "best_gating_baseline": best_a4a}
sens_a4b_out = {"per_origin": [clean_origin_for_json(o) for o in per_origin_a4b],
                "pooled": pooled_a4b, "best_gating_baseline": best_a4b, "j2_shift_log": j2_shift_log}

elapsed = time.time() - t0
results = {
    "task": "F15 (Amendments 9 re-seal) of F14 (F7e refit-2021 rerun on v7)",
    "run_config": {"base_seed": BASE_SEED, "n_tiebreak_trials": N_TIEBREAK_TRIALS,
                    "n_bootstrap": N_BOOT, "bootstrap_seed": BOOT_SEED, "topk_list": TOPK_LIST,
                    "origins": ORIGIN_STRS},
    "implementation_notes": [
        "features_refit_2021.csv REUSED VERBATIM from f7e_implementation/ (not recomputed): the G1 "
        "anchors-only fit does not read positive/date, and candidates_purged_v7.csv is cell-identical "
        "to candidates_purged_v5.csv on every column the fit depends on (species_key, member_ids, "
        "dist, logcd, sigma, mu, mw, charged, radical) for every surviving candidate_id -- v6 only "
        "drops 5 anchor-duplicate candidate_ids and flips positive/date for 9 ids. Only the reseal "
        "(this script) changed: v6 candidate skeleton (candidates_purged_v7.csv, 1434 rows) joined "
        "onto the unchanged features, same step4_engine.py/step5_bootstrap.py/j2_shift.py evaluation "
        "code as F12c's own seal_results_f7_refit2021_v5.py. features_refit2021_relabelled_v7.csv "
        "(delivered alongside this seal) is the per-origin join actually used here.",
    ],
    "primary": primary_out,
    "sensitivity_A4a_non_arxiv_dates_dropped": sens_a4a_out,
    "sensitivity_A4b_60day_boundary_flip": sens_a4b_out,
}

out_path = os.path.abspath(sys.argv[1])
with open(out_path, "w") as f:
    json.dump(results, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating,)) else
               (int(x) if isinstance(x, (np.integer,)) else (bool(x) if isinstance(x, (np.bool_,)) else str(x))))
with open(out_path, "rb") as f:
    h = hashlib.sha256(f.read()).hexdigest()
with open(out_path + ".sha256", "w") as f:
    f.write(f"{h}  {os.path.basename(out_path)}\n")
log(f"{os.path.basename(out_path)} written: {out_path}")
log(f"SHA-256: {h}")

with open(relabelled_path, "rb") as f:
    h2 = hashlib.sha256(f.read()).hexdigest()
log(f"{os.path.basename(relabelled_path)} SHA-256: {h2}")
log(f"elapsed: {elapsed:.1f}s")
