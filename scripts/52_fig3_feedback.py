"""
Figure 3: feedback simulation, released regime (data/frozen/candidates.csv), 63 positives,
target = ceil(0.5*63) = 32. Consistency: the learner's median observations-to-target(32) must equal
results/sealed/released.json's v3_feedback_learner.median (200.0).
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

from paths import SEALED, FIGURES
from evaluation_engine import raw_scores, topk_avg_counts, GATING_BASELINES, BASE_SEED, N_TIEBREAK_TRIALS
from feedback_learner import N_SEEDS
from regimes import build_released_v7
from feedback_curve import run_v3_feedback_curve, obs_to_target_from_curve
from plot_style import apply_rcparams, save_fig, METHOD_LABEL, METHOD_COLOR, METHOD_LINESTYLE, TWO_COL_W

OUT_DIR = FIGURES
MAX_OBS = 800
BATCH = 5

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


log("Building released regime (candidates_purged_v7.csv, full 1434-candidate pool)...")
cand, Xs, _ = build_released_v7()
y_pos = cand["positive"].to_numpy()
n_positives_total = int(y_pos.sum())
target = -(-n_positives_total // 2)  # ceil(0.5 * n_positives_total)
log(f"n_positives_total={n_positives_total}, target={target}")
assert n_positives_total == 63 and target == 32, "expected 63 positives / target 32 per F11e spec"

log(f"Running learner feedback curve: {N_SEEDS} seeds x {MAX_OBS//BATCH} checkpoints...")
n_obs_grid, curves = run_v3_feedback_curve(Xs, y_pos, n_seeds=N_SEEDS, max_obs=MAX_OBS, batch_size=BATCH)
learner_median = np.median(curves, axis=0)
learner_q25 = np.percentile(curves, 25, axis=0)
learner_q75 = np.percentile(curves, 75, axis=0)
log("Learner curve done.")

obs_to_target = obs_to_target_from_curve(n_obs_grid, curves, target)
median_obs_to_target = float(np.median(obs_to_target))

with open(os.path.join(SEALED, "released.json")) as f:
    sealed_f6 = json.load(f)
sealed_median = sealed_f6["v3_feedback_learner"]["median"]
match = median_obs_to_target == sealed_median
log(f"CONSISTENCY CHECK: learner median observations-to-target({target}) = {median_obs_to_target} "
    f"vs sealed v3_feedback_learner.median = {sealed_median} -- {'MATCH' if match else 'MISMATCH'}")
assert match, "Figure 3 consistency check failed"

# ---- 3 fixed baseline orderings ----
raw = raw_scores(cand)
k_list_nonzero = [int(k) for k in n_obs_grid if k > 0]
fixed_curves = {}
for m in GATING_BASELINES:
    counts = topk_avg_counts(raw[m], y_pos, k_list_nonzero, N_TIEBREAK_TRIALS, seed=BASE_SEED)
    fixed_curves[m] = np.array([0.0] + [counts[k] for k in k_list_nonzero])
log("Fixed-ordering baseline curves done.")

# ---- Random (expectation): analytic, no seed ----
random_curve = n_obs_grid * (n_positives_total / len(y_pos))

# ---- assemble fig3_data.csv ----
rows = []
for i, n_obs in enumerate(n_obs_grid):
    rows.append({"method": "combined", "n_obs": int(n_obs), "value": float(learner_median[i]),
                 "band_lo": float(learner_q25[i]), "band_hi": float(learner_q75[i])})
for m in GATING_BASELINES:
    for i, n_obs in enumerate(n_obs_grid):
        rows.append({"method": m, "n_obs": int(n_obs), "value": float(fixed_curves[m][i]),
                     "band_lo": "", "band_hi": ""})
for i, n_obs in enumerate(n_obs_grid):
    rows.append({"method": "random", "n_obs": int(n_obs), "value": float(random_curve[i]),
                 "band_lo": "", "band_hi": ""})

fig3_df = pd.DataFrame(rows, columns=["method", "n_obs", "value", "band_lo", "band_hi"])
fig3_path = os.path.join(OUT_DIR, "fig3_data.csv")
fig3_df.to_csv(fig3_path, index=False)
log(f"Wrote {fig3_path} ({len(fig3_df)} rows)")

pd.DataFrame([{
    "check": f"learner median observations-to-target({target})",
    "computed": median_obs_to_target, "expected": sealed_median,
    "match": match,
}]).to_csv(os.path.join(OUT_DIR, "fig3_consistency_check.csv"), index=False)

# ---- plot ----
apply_rcparams()
fig, ax = plt.subplots(figsize=(TWO_COL_W, 3.2))

sub = fig3_df[fig3_df.method == "combined"]
ax.plot(sub["n_obs"], sub["value"], color=METHOD_COLOR["combined"], linewidth=1.5,
        label="Feedback searcher")
ax.fill_between(sub["n_obs"], sub["band_lo"].astype(float), sub["band_hi"].astype(float),
                color=METHOD_COLOR["combined"], alpha=0.18, linewidth=0)

for m in GATING_BASELINES:
    sub = fig3_df[fig3_df.method == m]
    ax.plot(sub["n_obs"], sub["value"], color=METHOD_COLOR[m], linestyle=METHOD_LINESTYLE[m],
            linewidth=1.2, label=METHOD_LABEL[m])

sub = fig3_df[fig3_df.method == "random"]
ax.plot(sub["n_obs"], sub["value"], color=METHOD_COLOR["random"], linestyle=METHOD_LINESTYLE["random"],
        linewidth=1.2, label=METHOD_LABEL["random"])

ax.axhline(target, color="gray", linestyle=":", linewidth=0.9, zorder=0)
ax.set_xlabel("Observations")
ax.set_ylabel("Positives revealed")
ax.set_xlim(0, MAX_OBS)
ax.set_ylim(bottom=0)
ax.legend(loc="upper left", frameon=False, fontsize=6, ncol=1)

save_fig(fig, os.path.join(OUT_DIR, "fig3"))
plt.close(fig)
log("Wrote fig3.pdf / fig3.png")

log(f"DONE. Total runtime {time.time()-t0:.1f}s")
