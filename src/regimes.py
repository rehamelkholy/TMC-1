"""
The three per-origin regimes used by Figure 1's bootstrap, Figure 2 and Table 2:
  - released:       data/frozen/candidates.csv via seal_released.load_candidates_purged_v7 /
                    attach_date_basis_v7 (the loaders seal_released.py itself uses).
  - refit_2021:     data/derived/features_refit2021_relabelled.csv joined per-origin onto the candidate skeleton.
  - leave_self_out: data/frozen/features_leave_self_out.csv joined per-origin onto the candidate skeleton.
Evaluation code (evaluation_engine) unchanged. Run under the locked analysis environment.
"""
import os
import sys
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

from paths import DATA_DERIVED, DATA_FROZEN  # noqa: E402
from evaluation_engine import standardize_features, run_origin_evaluation, ORIGINS  # noqa: E402
from seal_released import load_candidates_purged_v7, attach_date_basis_v7  # noqa: E402

ORIGIN_STRS = [o.strftime("%Y-%m-%d") for o in ORIGINS]


def _build_cand_for_origin(features, origin_str, base_skel, log=None):
    sub = features[features["origin"] == origin_str].set_index("candidate_id")
    cand = base_skel.copy()
    cand["dist"] = cand["candidate_id"].map(sub["dist"])
    cand["logcd"] = cand["candidate_id"].map(sub["logcd"])
    cand["sigma"] = cand["candidate_id"].map(sub["sigma"])
    n_nan = cand[["dist", "logcd", "sigma"]].isna().any(axis=1).sum()
    if n_nan and log is not None:
        log(f"  WARNING: {n_nan} candidates have NaN features at origin {origin_str}")
    return cand


def build_released_v7():
    """Regime 1: candidates.csv's own dist/logcd/sigma."""
    cand, lee = load_candidates_purged_v7()
    cand = attach_date_basis_v7(cand, lee)
    Xs, mean, std = standardize_features(cand)
    per_origin = run_origin_evaluation(cand, Xs, sensitivity=None)
    return cand, Xs, per_origin


def build_refit_2021_v7(log=None):
    """Regime 2: features_refit2021_relabelled.csv, joined per-origin onto the candidate
    skeleton exactly like scripts/21_seal_refit2021.py."""
    features = pd.read_csv(os.path.join(DATA_DERIVED, "features_refit2021_relabelled.csv"))
    cand_skel, lee = load_candidates_purged_v7()
    cand_skel = attach_date_basis_v7(cand_skel, lee)
    cand_skel["date"] = pd.to_datetime(cand_skel["date"])

    per_origin = []
    for oi, ostr in enumerate(ORIGIN_STRS):
        cand_i = _build_cand_for_origin(features, ostr, cand_skel, log=log)
        Xs_i, _, _ = standardize_features(cand_i)
        res = run_origin_evaluation(cand_i, Xs_i, sensitivity=None)[oi]
        per_origin.append(res)
    return cand_skel, per_origin


def build_lso_v10(log=None):
    """Regime 3: features_leave_self_out.csv, joined per-origin onto the candidate
    skeleton exactly like scripts/22_seal_leave_self_out.py. Reads only the unsealed features CSV;
    never opens a sealed results file."""
    features = pd.read_csv(os.path.join(DATA_FROZEN, "features_leave_self_out.csv"))
    cand_skel, lee = load_candidates_purged_v7()
    cand_skel = attach_date_basis_v7(cand_skel, lee)
    cand_skel["date"] = pd.to_datetime(cand_skel["date"])

    per_origin = []
    for oi, ostr in enumerate(ORIGIN_STRS):
        cand_i = _build_cand_for_origin(features, ostr, cand_skel, log=log)
        Xs_i, _, _ = standardize_features(cand_i)
        res = run_origin_evaluation(cand_i, Xs_i, sensitivity=None)[oi]
        per_origin.append(res)
    return cand_skel, per_origin
