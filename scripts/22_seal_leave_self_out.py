"""
Leave-self-out regime under protocol amendments 9 (R1-R3). Output path is argv[1]; the
non-deterministic "generated"/"elapsed_seconds" fields are removed so repeated runs are
byte-identical.

Input: data/frozen/features_leave_self_out.csv (not recomputed here). The evaluation must be run
under the locked analysis environment (env/requirements-analysis.txt, scikit-learn 1.8.0): the lbfgs
fit is sensitive to the scikit-learn version, and the leave-self-out features were generated in a
separate embedding environment (scikit-learn 1.0.2).
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

from paths import DATA_FROZEN
from evaluation_engine import (
    standardize_features, run_origin_evaluation, pool_metrics,
    best_gating_baseline, FEATURE_COLS, TOPK_LIST, ORIGINS, BASE_SEED,
    N_TIEBREAK_TRIALS, GATING_BASELINES, NONGATING, raw_scores,
)
from bootstrap import run_bootstrap, N_BOOT, BOOT_SEED
from boundary_shift import apply_j2_shift
from seal_released import load_candidates_purged_v7 as load_candidates_purged_v6, attach_date_basis_v7 as attach_date_basis_v6, clean_origin_for_json

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


features = pd.read_csv(os.path.join(DATA_FROZEN, "features_leave_self_out.csv"))
cand_skel, lee = load_candidates_purged_v6()
cand_skel = attach_date_basis_v6(cand_skel, lee)
log(f"Loaded features_recomputed_lso_v10.csv ({len(features)} rows) and v6 candidate skeleton "
    f"({len(cand_skel)} candidates, {int(cand_skel['positive'].sum())} positive).")

cand_skel["date"] = pd.to_datetime(cand_skel["date"])
cand_j2_skel, j2_shift_log = apply_j2_shift(cand_skel)


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
for oi, ostr in enumerate(ORIGIN_STRS):
    log(f"--- Evaluating origin {oi} ({ostr}) ---")
    cand_i = build_cand_for_origin(ostr, cand_skel)
    Xs_i, mean_i, std_i = standardize_features(cand_i)

    res_primary = run_origin_evaluation(cand_i, Xs_i, sensitivity=None)[oi]
    per_origin_primary.append(res_primary)
    res_a4a = run_origin_evaluation(cand_i, Xs_i, sensitivity="A4a")[oi]
    per_origin_a4a.append(res_a4a)
    cand_i_j2 = build_cand_for_origin(ostr, cand_j2_skel)
    Xs_i_j2, _, _ = standardize_features(cand_i_j2)
    res_a4b = run_origin_evaluation(cand_i_j2, Xs_i_j2, sensitivity=None)[oi]
    per_origin_a4b.append(res_a4b)
    log(f"  known_before_origin={res_primary['known_before_origin']}, "
        f"test_pool_size={res_primary['test_pool_size']}, test_positives={res_primary['test_positives']}")

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

import sklearn, scipy, platform

env_fingerprint = {
    "python": platform.python_version(),
    "python_executable": sys.executable,
    "scikit-learn": sklearn.__version__,
    "numpy": np.__version__,
    "scipy": scipy.__version__,
    "pandas": pd.__version__,
}

elapsed = time.time() - t0
results = {
    "task": "F15 (Amendments 9 re-seal) of F14, step 3: F7 leave-self-out, v7 (F14-frozen) labels, sealed under the locked environment; (F13f text follows) F7 leave-self-out, v6 "
            "labels, features_recomputed_lso_v10.csv -- IDENTICAL input and evaluation code "
            "to results_repo_f7_lso_v8.json, run under the locked environment (system "
            "python3, scikit-learn 1.8.0) instead of the venv (scikit-learn 1.0.2) that "
            "produced results_repo_f7_lso_v8.json. Supersedes results_repo_f7_lso_v8.json "
            "FOR REPORTING; that file is kept on disk unmodified.",
    "environment": env_fingerprint,
    "run_config": {"base_seed": BASE_SEED, "n_tiebreak_trials": N_TIEBREAK_TRIALS,
                    "n_bootstrap": N_BOOT, "bootstrap_seed": BOOT_SEED, "topk_list": TOPK_LIST,
                    "origins": ORIGIN_STRS},
    "implementation_notes": [
        "features_recomputed_lso_v10.csv is read UNCHANGED from disk (SHA-256 "
        "33db98c315867986afe5c11dc44d10fb665da3853d3c5f668ab924771b664850, same file used by "
        "results_repo_f7_lso_v8.json) -- no feature recomputation in this task. The evaluation "
        "code (step4_engine.py, step5_bootstrap.py, j2_shift.py, seal_results_f6_v6.py) is "
        "identical/unmodified; the only variable changed from seal_results_f7_lso_v8.py is the "
        "interpreter this script is run under (system python3 / scikit-learn 1.8.0, matching "
        "the released_2021/refit_2021 regimes' own sealed environment, in place of the venv's "
        "scikit-learn 1.0.2). See 'environment' above for the exact versions and the F11d log "
        "entry for the root-cause finding (lbfgs fits differ slightly by sklearn version on "
        "identical data) that this re-seal exists to control for under a third regime.",
        "results_repo_f7_lso_v8.json is superseded by this file FOR REPORTING purposes (same "
        "features input, single-environment-consistent evaluation) but is kept on disk "
        "unmodified, per explicit instruction.",
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
log(f"elapsed: {elapsed:.1f}s")
