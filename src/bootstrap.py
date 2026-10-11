"""
Paired bootstrap (protocol item I6): 2000 resamples, seed 20260921. Within each
origin, resample positives and negatives separately with replacement,
identical resampled indices used for both methods being compared (combined
vs best gating baseline). Recompute pooled AUC and pooled top-100 count each
resample; report 95% percentile intervals of (combined - best_baseline).

Within each resample, ties are broken by ONE random draw shared between the two
paired methods, so the pairing stays exact. I6 does not specify nesting a
multi-trial tie-break inside every resample, which would require 2000*200
re-rankings per origin per method; a single shared draw is the tractable reading
of "identical indices for both methods".
"""
import numpy as np
from sklearn.metrics import roc_auc_score

N_BOOT = 2000
BOOT_SEED = 20260921


def run_bootstrap(per_origin_primary, best_baseline, n_boot=N_BOOT, seed=BOOT_SEED):
    """per_origin_primary: list of dicts from run_origin_evaluation(..., sensitivity=None),
    but we need the RAW per-candidate score/label arrays, not just the aggregated
    metrics -- caller must attach `_test_positive_bool`, `_combined_score`,
    `_baseline_score` to each origin dict before calling this (see driver)."""
    rng = np.random.default_rng(seed)
    diffs_auc = np.zeros(n_boot)
    diffs_top100 = np.zeros(n_boot)

    # pre-split positive/negative index arrays per origin (fixed across all resamples)
    origin_data = []
    for o in per_origin_primary:
        y = o["_test_positive_bool"]
        pos_idx = np.where(y)[0]
        neg_idx = np.where(~y)[0]
        origin_data.append({
            "pos_idx": pos_idx, "neg_idx": neg_idx,
            "combined": o["_combined_score"], "baseline": o["_baseline_score"],
        })

    for b in range(n_boot):
        pooled_auc_c_num, pooled_auc_c_den = 0.0, 0.0
        pooled_auc_k_num, pooled_auc_k_den = 0.0, 0.0
        top100_c, top100_k = 0.0, 0.0
        for od in origin_data:
            pos_idx, neg_idx = od["pos_idx"], od["neg_idx"]
            if len(pos_idx) == 0:
                continue
            pos_boot = rng.choice(pos_idx, size=len(pos_idx), replace=True)
            neg_boot = rng.choice(neg_idx, size=len(neg_idx), replace=True) if len(neg_idx) else np.array([], dtype=int)
            idx_boot = np.concatenate([pos_boot, neg_boot])
            y_boot = np.concatenate([np.ones(len(pos_boot)), np.zeros(len(neg_boot))])

            sc_c = od["combined"][idx_boot]
            sc_k = od["baseline"][idx_boot]

            if 0 < y_boot.sum() < len(y_boot):
                auc_c = roc_auc_score(y_boot, sc_c)
                auc_k = roc_auc_score(y_boot, sc_k)
                pooled_auc_c_num += auc_c * len(pos_idx)
                pooled_auc_k_num += auc_k * len(pos_idx)
                pooled_auc_c_den += len(pos_idx)
                pooled_auc_k_den += len(pos_idx)

            jitter = rng.random(len(idx_boot))  # shared tie-break for both methods (paired)
            order_c = np.lexsort((jitter, -sc_c))
            order_k = np.lexsort((jitter, -sc_k))
            kk = min(100, len(idx_boot))
            top100_c += y_boot[order_c][:kk].sum()
            top100_k += y_boot[order_k][:kk].sum()

        pooled_auc_c = pooled_auc_c_num / pooled_auc_c_den if pooled_auc_c_den else float("nan")
        pooled_auc_k = pooled_auc_k_num / pooled_auc_k_den if pooled_auc_k_den else float("nan")
        diffs_auc[b] = pooled_auc_c - pooled_auc_k
        diffs_top100[b] = top100_c - top100_k

    ci_auc = (float(np.percentile(diffs_auc, 2.5)), float(np.percentile(diffs_auc, 97.5)))
    ci_top100 = (float(np.percentile(diffs_top100, 2.5)), float(np.percentile(diffs_top100, 97.5)))
    return {
        "n_boot": n_boot, "seed": seed,
        "auc_diff_ci95": ci_auc, "top100_diff_ci95": ci_top100,
        "auc_diff_mean": float(np.mean(diffs_auc)), "top100_diff_mean": float(np.mean(diffs_top100)),
    }
