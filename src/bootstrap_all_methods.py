"""
Figure 1 needs an ABSOLUTE 95% CI on pooled AUC for each of 6 methods (combined, the 3 gating
baselines, high_sigma, ucb) in each of 3 regimes. bootstrap.run_bootstrap() only computes a single
PAIRED DIFFERENCE (combined - best_baseline), not each method's own absolute pooled AUC
distribution, so it cannot be reused as-is.

This module is a wrapper that reuses bootstrap.py's resampling scheme: 2000 resamples, seed
20260921, positives and negatives resampled separately WITH REPLACEMENT within each origin, using
the same origin-level pos_idx/neg_idx split off of `_test_positive_bool`. The only change is what is
computed inside each resample: each method's own AUC on the resampled origin-level test set (no
top-k/tie-break needed for an AUC-only CI), pooled origin-weighted (weight = original len(pos_idx),
matching pool_metrics' positive-weighted pooling), keeping the full 2000-length distribution per
method to take 2.5/97.5 percentiles from.
"""
import numpy as np
from sklearn.metrics import roc_auc_score

N_BOOT = 2000
BOOT_SEED = 20260921


def method_score_arrays(per_origin, method):
    """Per-origin score array for `method` on that origin's test pool, exactly what
    evaluation_engine.run_origin_evaluation already computed and retained (_combined_score for
    'combined', _raw_scores_test[m] for every other method)."""
    out = []
    for o in per_origin:
        if method == "combined":
            out.append(o["_combined_score"])
        else:
            out.append(o["_raw_scores_test"][method])
    return out


def run_bootstrap_all_methods(per_origin, methods, n_boot=N_BOOT, seed=BOOT_SEED):
    """Returns {method: {"auc_ci95": (lo,hi), "auc_boot_mean": m}}."""
    rng = np.random.default_rng(seed)

    origin_data = []
    for o in per_origin:
        y = o["_test_positive_bool"]
        pos_idx = np.where(y)[0]
        neg_idx = np.where(~y)[0]
        scores = {m: method_score_arrays([o], m)[0] for m in methods}
        origin_data.append({"pos_idx": pos_idx, "neg_idx": neg_idx, "scores": scores})

    boot_aucs = {m: np.full(n_boot, np.nan) for m in methods}

    for b in range(n_boot):
        num = {m: 0.0 for m in methods}
        den = {m: 0.0 for m in methods}
        for od in origin_data:
            pos_idx, neg_idx = od["pos_idx"], od["neg_idx"]
            if len(pos_idx) == 0:
                continue
            pos_boot = rng.choice(pos_idx, size=len(pos_idx), replace=True)
            neg_boot = (rng.choice(neg_idx, size=len(neg_idx), replace=True)
                        if len(neg_idx) else np.array([], dtype=int))
            idx_boot = np.concatenate([pos_boot, neg_boot])
            y_boot = np.concatenate([np.ones(len(pos_boot)), np.zeros(len(neg_boot))])
            weight = len(pos_idx)  # original (not resampled) positive count, matches pool_metrics
            if not (0 < y_boot.sum() < len(y_boot)):
                continue
            for m in methods:
                sc = od["scores"][m][idx_boot]
                auc = roc_auc_score(y_boot, sc)
                num[m] += auc * weight
                den[m] += weight
        for m in methods:
            boot_aucs[m][b] = num[m] / den[m] if den[m] > 0 else np.nan

    results = {}
    for m in methods:
        vals = boot_aucs[m]
        vals = vals[~np.isnan(vals)]
        lo = float(np.percentile(vals, 2.5)) if len(vals) else float("nan")
        hi = float(np.percentile(vals, 97.5)) if len(vals) else float("nan")
        mean = float(np.mean(vals)) if len(vals) else float("nan")
        results[m] = {"auc_ci95": (lo, hi), "auc_boot_mean": mean, "n_boot_valid": int(len(vals))}
    return results
