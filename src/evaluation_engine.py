"""
Core evaluation engine (origin-based evaluation, baselines and bootstrap), implemented from
protocol_v1.md and protocol amendments 1-3 (items I1-I9).

Protocol amendment 9 (protocol_v1_amendments_9.md): R1 random AUC/top-k are exact expectations; R2 top-k
with ties is the exact expectation under uniform random tie-breaking; R3 seed_for() uses sha256, not hash().
"""
import hashlib

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

BASE_SEED = 20260921
N_TIEBREAK_TRIALS = 200
TOPK_LIST = [25, 50, 100, 200]
ORIGINS = [pd.Timestamp(2022, 1, 1), pd.Timestamp(2023, 1, 1),
           pd.Timestamp(2024, 1, 1), pd.Timestamp(2025, 1, 1)]
FREEZE_DATE = pd.Timestamp(2026, 9, 21)
FEATURE_COLS = ["dist", "logcd", "sigma", "mu", "mw", "charged", "radical"]

GATING_BASELINES = ["high_logcd", "detectability", "low_dist"]
NONGATING = ["high_sigma", "ucb"]
ALL_METHODS = ["combined"] + GATING_BASELINES + ["random"] + NONGATING


def _stable_int(x):
    """Stable 64-bit integer from str(x) (Amendments 9, R3): replaces Python's salted hash()."""
    return int.from_bytes(hashlib.sha256(str(x).encode("utf-8")).digest()[:8], "big")


def seed_for(*parts):
    """Deterministic per-(origin,method,...) seed derived from BASE_SEED, stable across processes
    (Amendments 9, R3: sha256-based, not Python's hash())."""
    s = BASE_SEED
    for p in parts:
        s = (s * 1000003 + _stable_int(p)) % (2**31 - 1)
    return s


def standardize_features(cand):
    X = cand[FEATURE_COLS].astype(float).to_numpy()
    mean = X.mean(axis=0)
    std = X.std(axis=0, ddof=0)
    std_safe = np.where(std == 0, 1.0, std)
    Xs = (X - mean) / std_safe
    return Xs, mean, std_safe


def raw_scores(cand):
    """I4 baseline + non-gating comparator scores, from raw (unstandardized)
    features."""
    logcd = cand["logcd"].to_numpy(float)
    dist = cand["dist"].to_numpy(float)
    sigma = cand["sigma"].to_numpy(float)
    mu = cand["mu"].to_numpy(float)
    detect = np.where(mu > 0, logcd + 2 * np.log10(np.where(mu > 0, mu, 1.0)), -99.0)
    return {
        "high_logcd": logcd.copy(),
        "detectability": detect,
        "low_dist": -dist,
        "high_sigma": sigma.copy(),
        "ucb": logcd + sigma,
    }


def topk_exact_counts(score, y, k_list):
    """Amendments 9, R2: exact expectation of the top-k positive count under uniform random
    tie-breaking. Candidates strictly above the cutoff score count fully; the tied block at the
    cutoff contributes (tied positives) x (remaining slots) / (tied candidates). `score`=None means
    the random ordering (R1): k x positives / n."""
    y = np.asarray(y, dtype=bool)
    n = len(y)
    if score is None:
        npos = int(y.sum())
        return {k: (min(k, n) * npos / n if n else 0.0) for k in k_list}
    s = np.asarray(score, dtype=float)
    s_sorted = np.sort(s)[::-1]
    out = {}
    for k in k_list:
        kk = min(k, n)
        if kk <= 0:
            out[k] = 0.0
            continue
        cutoff = s_sorted[kk - 1]
        above = s > cutoff
        tied = s == cutoff
        n_above = int(above.sum())
        pos_above = int((y & above).sum())
        t_n = int(tied.sum())
        t_pos = int((y & tied).sum())
        out[k] = float(pos_above + t_pos * (kk - n_above) / t_n)
    return out


def topk_avg_counts(score, y, k_list, n_trials=None, seed=None):
    """Amendments 9, R2 (replaces I4's 200-draw tie-break average): exact expectation under uniform
    random tie-breaking. `n_trials` and `seed` are accepted for backward compatibility and ignored."""
    return topk_exact_counts(score, y, k_list)


def auc_score(score, y):
    y = np.asarray(y, dtype=int)
    if y.sum() == 0 or y.sum() == len(y):
        return float("nan")
    return float(roc_auc_score(y, score))


def random_auc_and_topk(y, k_list, n_trials=None, seed=None):
    """Amendments 9, R1: the random ordering's AUC is its exact expectation, 0.5 (NaN if the origin has no
    positives or no negatives), and its top-k count is the exact expectation k x (positives / pool size).
    `n_trials` and `seed` are accepted for backward compatibility and ignored."""
    y = np.asarray(y, dtype=bool)
    npos = int(y.sum())
    mean_auc = float("nan") if (npos == 0 or npos == len(y)) else 0.5
    return mean_auc, topk_exact_counts(None, y, k_list)


def known_mask_at(dates, origin, flip_mask=None):
    """dates: array of pd.Timestamp/NaT. known = positive & date < origin.
    flip_mask (A4b): boolean array, True for candidates whose comparison
    against THIS origin boundary is inverted (within 60 days of it)."""
    base = np.array([(pd.notna(d)) and (d < origin) for d in dates])
    if flip_mask is not None:
        base = np.where(flip_mask, ~base, base)
    return base


def build_a4b_flip_masks(cand):
    """For each of the 4 origins, which positive candidates lie within 60 days
    of that boundary (their known-status comparison against that boundary is
    inverted under sensitivity A4(b))."""
    flips = {}
    for o in ORIGINS:
        d = cand["date"]
        near = pd.Series([
            (pd.notna(x) and abs((x - o).days) <= 60) for x in d
        ], index=cand.index)
        flips[o] = near.to_numpy()
    return flips


def run_origin_evaluation(cand, Xs, sensitivity=None):
    """sensitivity in {None, 'A4a', 'A4b'}. Returns per-origin dict + pooled."""
    dates = cand["date"].to_numpy()
    positive = cand["positive"].to_numpy()
    raw = raw_scores(cand)
    flip_masks = build_a4b_flip_masks(cand) if sensitivity == "A4b" else None

    if sensitivity == "A4a":
        is_arxiv = cand["date_basis"].astype(str).str.contains("arxiv", case=False).to_numpy()
        # only meaningful for positives; drop these positives from the test-positive
        # pool entirely (I8: "drop from the test pool test positives whose date
        # is not an arXiv v1 date") -- treat them as if not-yet-detected (no date)
        drop = positive & (~is_arxiv)
    else:
        drop = np.zeros(len(cand), dtype=bool)

    def known_at(boundary):
        """known-at-boundary mask, honoring the A4(b) 60-day flip for real
        origin boundaries only (the trailing FREEZE_DATE is not an origin and
        is never flipped)."""
        flip = flip_masks.get(boundary) if (sensitivity == "A4b" and boundary in flip_masks) else None
        return known_mask_at(dates, boundary, flip) & positive & (~drop)

    per_origin = []
    for oi, origin in enumerate(ORIGINS):
        next_o = ORIGINS[oi + 1] if oi + 1 < len(ORIGINS) else FREEZE_DATE

        known = known_at(origin)
        test_pool = ~known  # I3: test pool = candidates not known at o
        known_next = known_at(next_o)
        test_positive = test_pool & known_next  # known at next origin but not at this one = in this window
        # (candidates dated after the window stay in the pool as negatives: automatically
        # true since they're not "known_next" either, so test_positive=False for them)

        # training labels: y=1 for known, 0 for all others (I3)
        y_train = known.astype(int)
        if y_train.sum() == 0 or y_train.sum() == len(y_train):
            model = None
        else:
            model = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=1000)
            model.fit(Xs, y_train)

        idx_test = np.where(test_pool)[0]
        y_test = test_positive[idx_test]
        Xs_test = Xs[idx_test]

        if model is not None:
            combined_score = model.predict_proba(Xs_test)[:, 1]
        else:
            combined_score = np.zeros(len(idx_test))

        origin_result = {
            "origin": origin.strftime("%Y-%m-%d"),
            "window_end": next_o.strftime("%Y-%m-%d"),
            "known_before_origin": int(known.sum()),
            "test_pool_size": int(len(idx_test)),
            "test_positives": int(y_test.sum()),
            "methods": {},
            # raw arrays kept for bootstrap; stripped before JSON serialization
            "_test_positive_bool": y_test.astype(bool),
            "_combined_score": combined_score,
            "_raw_scores_test": {m: raw[m][idx_test] for m in GATING_BASELINES + NONGATING},
        }

        seed_combined = seed_for(str(sensitivity), origin, "combined")
        origin_result["methods"]["combined"] = {
            "auc": auc_score(combined_score, y_test),
            "topk": topk_avg_counts(combined_score, y_test, TOPK_LIST, N_TIEBREAK_TRIALS, seed_combined),
        }

        for m in GATING_BASELINES + NONGATING:
            s_test = raw[m][idx_test]
            seed_m = seed_for(str(sensitivity), origin, m)
            origin_result["methods"][m] = {
                "auc": auc_score(s_test, y_test),
                "topk": topk_avg_counts(s_test, y_test, TOPK_LIST, N_TIEBREAK_TRIALS, seed_m),
            }

        seed_rand = seed_for(str(sensitivity), origin, "random")
        rand_auc, rand_topk = random_auc_and_topk(y_test, TOPK_LIST, N_TIEBREAK_TRIALS, seed_rand)
        origin_result["methods"]["random"] = {"auc": rand_auc, "topk": rand_topk}

        per_origin.append(origin_result)

    return per_origin


def pool_metrics(per_origin, methods=ALL_METHODS):
    """Pooled counts = sum over origins. Pooled AUC = positive-weighted mean
    of per-origin AUCs (I5)."""
    pooled = {}
    for m in methods:
        aucs = []
        weights = []
        topk_sums = {k: 0.0 for k in TOPK_LIST}
        for o in per_origin:
            npos = o["test_positives"]
            a = o["methods"][m]["auc"]
            if npos > 0 and not np.isnan(a):
                aucs.append(a)
                weights.append(npos)
            for k in TOPK_LIST:
                topk_sums[k] += o["methods"][m]["topk"][k]
        pooled_auc = float(np.average(aucs, weights=weights)) if aucs else float("nan")
        pooled[m] = {"auc": pooled_auc, "topk": topk_sums}
    return pooled


def best_gating_baseline(pooled):
    """I5: best baseline = the gating baseline (not random) with the highest
    pooled top-100 count."""
    return max(GATING_BASELINES, key=lambda m: pooled[m]["topk"][100])
