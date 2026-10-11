"""
Figure 2: budget vs yield, released regime (data/frozen/candidates.csv), computed with the protocol
amendment 9 engine (evaluation_engine: R2 exact expectation under uniform random tie-breaking, R1 exact
random curve) and checked against results/sealed/released.json. Every k = 1..200 value of every curve
is additionally checked against an independent exact-expectation computation (unique-score groups via
pandas), and all 16 sealed counts at k = 25/50/100/200 must match.
"""
import os
import sys
import json
import time
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from paths import SEALED, FIGURES
from evaluation_engine import (
    seed_for, topk_avg_counts, random_auc_and_topk, GATING_BASELINES, N_TIEBREAK_TRIALS, ORIGINS,
)
from regimes import build_released_v7
from plot_style import apply_rcparams, save_fig, METHOD_LABEL, METHOD_COLOR, METHOD_LINESTYLE, TWO_COL_W

OUT_DIR = FIGURES
K_FULL = list(range(1, 201))
N_BOOT = 2000
BOOT_SEED = 20260921
BAND_METHODS = ["combined", "detectability"]
CURVE_METHODS = ["combined"] + GATING_BASELINES + ["random"]

EXPECTED_POOLS = [1400, 1390, 1384, 1375]
EXPECTED_POS = [10, 6, 9, 4]

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


log("Building released regime...")
_, _, per_origin = build_released_v7()
pools = [o["test_pool_size"] for o in per_origin]
positives = [o["test_positives"] for o in per_origin]
log(f"test_pool_size = {pools}")
log(f"test_positives = {positives}")
assert pools == EXPECTED_POOLS, f"pool mismatch: {pools} vs {EXPECTED_POOLS}"
assert positives == EXPECTED_POS, f"positives mismatch: {positives} vs {EXPECTED_POS}"

# ---- point estimates: pooled (summed) count-at-k for every k=1..200, combined + baselines ----
pooled_curve = {m: np.zeros(len(K_FULL)) for m in CURVE_METHODS}
for oi, o in enumerate(per_origin):
    origin_ts = ORIGINS[oi]
    assert origin_ts.strftime("%Y-%m-%d") == o["origin"]
    y = o["_test_positive_bool"]
    for m in ["combined"] + GATING_BASELINES:
        score = o["_combined_score"] if m == "combined" else o["_raw_scores_test"][m]
        seed_m = seed_for(str(None), origin_ts, m)
        counts = topk_avg_counts(score, y, K_FULL, N_TIEBREAK_TRIALS, seed_m)
        pooled_curve[m] += np.array([counts[k] for k in K_FULL])
    log(f"  origin {o['origin']}: point-estimate curves done.")

# ---- random curve: exact analytic expectation ----
rate_sum = sum(p / n for p, n in zip(positives, pools))
random_exact = np.array([k * rate_sum for k in K_FULL])
pooled_curve["random"] = random_exact
k100_val = float(random_exact[K_FULL.index(100)])
log(f"Random curve exact expectation at k=100: {k100_val:.6f}")

# ---- independent exact-expectation check at EVERY k = 1..200 (not using evaluation_engine.topk_exact_counts) ----
def independent_pooled_curve(score_key):
    tot = np.zeros(len(K_FULL))
    for o in per_origin:
        y = o["_test_positive_bool"].astype(int)
        sc = o["_combined_score"] if score_key == "combined" else o["_raw_scores_test"][score_key]
        df = pd.DataFrame({"s": np.asarray(sc, float), "y": y})
        g = df.groupby("s").agg(n=("y", "size"), p=("y", "sum")).sort_index(ascending=False)
        n_cum_before, p_cum_before = 0, 0
        blocks = [(n_cum_before := n_cum_before, ) for _ in []]  # unused; the blocks are computed below
        nb = g["n"].to_numpy(); pb = g["p"].to_numpy()
        cn = np.concatenate([[0], np.cumsum(nb)]); cp = np.concatenate([[0], np.cumsum(pb)])
        for ki, k in enumerate(K_FULL):
            kk = min(k, len(df))
            j = int(np.searchsorted(cn, kk, side="left"))   # block containing rank kk: cn[j-1] < kk <= cn[j]
            if cn[j] == kk:
                val = cp[j]
            else:
                val = cp[j - 1] + pb[j - 1] * (kk - cn[j - 1]) / nb[j - 1]
            tot[ki] += val
    return tot


every_k_rows = []
for m in ["combined"] + GATING_BASELINES:
    ref = independent_pooled_curve(m)
    maxdiff = float(np.max(np.abs(ref - pooled_curve[m])))
    every_k_rows.append({"method": m, "k_range": "1..200", "max_abs_diff_vs_independent_exact": maxdiff,
                         "match": maxdiff < 1e-9})
rand_ref = np.array([sum(min(k, p_n) * p / p_n for p, p_n in zip(positives, pools)) for k in K_FULL])
every_k_rows.append({"method": "random", "k_range": "1..200",
                     "max_abs_diff_vs_independent_exact": float(np.max(np.abs(rand_ref - pooled_curve["random"]))),
                     "match": bool(np.max(np.abs(rand_ref - pooled_curve["random"])) < 1e-9)})
pd.DataFrame(every_k_rows).to_csv(os.path.join(OUT_DIR, "fig2_everyk_check.csv"), index=False)
log("Every-k check (k=1..200) vs independent exact expectation: " + "; ".join(
    f"{r['method']} max|diff|={r['max_abs_diff_vs_independent_exact']:.2e}" for r in every_k_rows))
assert all(r["match"] for r in every_k_rows), "every-k check failed"

# ---- consistency check against sealed pooled topk at k=25/50/100/200 ----
with open(os.path.join(SEALED, "released.json")) as f:
    sealed_f6 = json.load(f)
sealed_pooled = sealed_f6["primary"]["pooled"]
check_ks = [25, 50, 100, 200]
consistency_rows = []
for m in ["combined"] + GATING_BASELINES:
    for k in check_ks:
        computed = float(pooled_curve[m][K_FULL.index(k)])
        sealed_val = float(sealed_pooled[m]["topk"][str(k)])
        match = abs(computed - sealed_val) < 1e-9
        note = ""
        if not match:
            # A mismatch can only arise from a score tie at this k; ties are resolved by the exact
            # expectation (protocol amendment 9, R2), so under amendment 9 every count must match.
            note = ("expected: unpinned PYTHONHASHSEED changes seed_for()'s hash-derived tie-break "
                    "seed across process launches; low_dist has a genuine score tie at this k/origin")
        consistency_rows.append({"method": m, "k": k, "computed": computed, "sealed": sealed_val,
                                  "match": match, "note": note})
n_mismatch = sum(1 for r in consistency_rows if not r["match"])
assert n_mismatch == 0, "sealed-count check at k=25/50/100/200 failed"
log(f"Consistency check: Figure 2 values @ k=25/50/100/200 vs sealed pooled top-k (combined + "
    f"3 baselines) -- {len(consistency_rows)-n_mismatch}/{len(consistency_rows)} match exactly, "
    f"{n_mismatch} mismatched (expected tie-break seed non-reproducibility, see note column).")
pd.DataFrame(consistency_rows).to_csv(os.path.join(OUT_DIR, "fig2_consistency_check.csv"), index=False)

# ---- bootstrap bands (combined, detectability) ----
rng = np.random.default_rng(BOOT_SEED)
origin_data = []
for o in per_origin:
    y = o["_test_positive_bool"]
    pos_idx = np.where(y)[0]
    neg_idx = np.where(~y)[0]
    origin_data.append({
        "pos_idx": pos_idx, "neg_idx": neg_idx,
        "combined": o["_combined_score"], "detectability": o["_raw_scores_test"]["detectability"],
    })

boot_curves = {m: np.zeros((N_BOOT, len(K_FULL))) for m in BAND_METHODS}
for b in range(N_BOOT):
    totals = {m: np.zeros(len(K_FULL)) for m in BAND_METHODS}
    for od in origin_data:
        pos_idx, neg_idx = od["pos_idx"], od["neg_idx"]
        if len(pos_idx) == 0:
            continue
        pos_boot = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        neg_boot = rng.choice(neg_idx, size=len(neg_idx), replace=True) if len(neg_idx) else np.array([], dtype=int)
        idx_boot = np.concatenate([pos_boot, neg_boot])
        y_boot = np.concatenate([np.ones(len(pos_boot)), np.zeros(len(neg_boot))])
        jitter = rng.random(len(idx_boot))
        for m in BAND_METHODS:
            sc = od[m][idx_boot]
            order = np.lexsort((jitter, -sc))
            y_ord = y_boot[order]
            cum = np.cumsum(y_ord)
            n = len(cum)
            for ki, k in enumerate(K_FULL):
                kk = min(k, n)
                totals[m][ki] += cum[kk - 1] if kk > 0 else 0.0
    for m in BAND_METHODS:
        boot_curves[m][b] = totals[m]
    if (b + 1) % 500 == 0:
        log(f"  bootstrap resample {b+1}/{N_BOOT}")

band_ci = {}
for m in BAND_METHODS:
    lo = np.percentile(boot_curves[m], 2.5, axis=0)
    hi = np.percentile(boot_curves[m], 97.5, axis=0)
    band_ci[m] = (lo, hi)
log("Bootstrap bands done.")

# ---- assemble fig2_data.csv ----
rows = []
for m in CURVE_METHODS:
    lo_arr, hi_arr = band_ci.get(m, (None, None))
    for ki, k in enumerate(K_FULL):
        rows.append({
            "method": m, "k": k, "found_mean": float(pooled_curve[m][ki]),
            "ci_lo": float(lo_arr[ki]) if lo_arr is not None else "",
            "ci_hi": float(hi_arr[ki]) if hi_arr is not None else "",
        })
fig2_df = pd.DataFrame(rows, columns=["method", "k", "found_mean", "ci_lo", "ci_hi"])
fig2_path = os.path.join(OUT_DIR, "fig2_data.csv")
fig2_df.to_csv(fig2_path, index=False)
log(f"Wrote {fig2_path} ({len(fig2_df)} rows)")

# ---- plot ----
apply_rcparams()
fig, ax = plt.subplots(figsize=(TWO_COL_W, 3.2))
for m in CURVE_METHODS:
    sub = fig2_df[fig2_df.method == m]
    ax.plot(sub["k"], sub["found_mean"], color=METHOD_COLOR[m], linestyle=METHOD_LINESTYLE[m],
            linewidth=1.3, label=METHOD_LABEL[m])
    if m in BAND_METHODS:
        lo_arr, hi_arr = band_ci[m]
        ax.fill_between(K_FULL, lo_arr, hi_arr, color=METHOD_COLOR[m], alpha=0.18, linewidth=0)

ax.axvline(100, color="gray", linestyle=":", linewidth=0.9, zorder=0)
ax.set_xlabel("Candidates searched per origin, $k$")
ax.set_ylabel("Test positives found\n(summed over 4 origins)")
ax.set_xlim(1, 200)
ax.set_ylim(bottom=0)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.legend(loc="upper left", frameon=False, fontsize=6, ncol=1)

save_fig(fig, os.path.join(OUT_DIR, "fig2"))
plt.close(fig)
log("Wrote fig2.pdf / fig2.png")

log(f"DONE. Total runtime {time.time()-t0:.1f}s")
