"""
Figure 1: pooled AUC with 95% bootstrap CI, 3 regimes x 6 methods.
released = data/frozen/candidates.csv; refit-2021 = data/derived/features_refit2021_relabelled.csv;
leave-self-out = data/frozen/features_leave_self_out.csv. Consistency target: the sealed pooled AUCs in
results/sealed/released.json, refit2021.json and leave_self_out.json.
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

from paths import SEALED, FIGURES
from evaluation_engine import pool_metrics, best_gating_baseline, GATING_BASELINES, NONGATING
from regimes import build_released_v7, build_refit_2021_v7, build_lso_v10
from bootstrap_all_methods import run_bootstrap_all_methods, N_BOOT, BOOT_SEED
from plot_style import apply_rcparams, save_fig, METHOD_LABEL, METHOD_COLOR, ONE_COL_W

OUT_DIR = FIGURES

FIG1_METHODS = ["combined", "high_logcd", "detectability", "low_dist", "high_sigma", "ucb"]

REGIME_LABEL_V7 = {
    "released": "Released (pre-registered)",
    "refit_2021": "Refit at 2021 (post hoc)",
    "lso_v10": "Updated leave-self-out (post hoc)",
}
REGIME_POSTHOC_V7 = {"released": False, "refit_2021": True, "lso_v10": True}

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


# ---- already-sealed comparison values (read directly, never written) ----
with open(os.path.join(SEALED, "released.json")) as f:
    sealed_f6 = json.load(f)
with open(os.path.join(SEALED, "refit2021.json")) as f:
    sealed_refit = json.load(f)
with open(os.path.join(SEALED, "leave_self_out.json")) as f:
    sealed_lso = json.load(f)
log("Comparison target for lso_v10: results_repo_f7_lso_v10.json (system python3 / "
    "sklearn 1.8.0 seal) -- NOT the venv-sealed results_repo_f7_lso_v10.json.")

SEALED_POOLED = {
    "released": {m: sealed_f6["primary"]["pooled"][m]["auc"] for m in FIG1_METHODS},
    "refit_2021": {m: sealed_refit["primary"]["pooled"][m]["auc"] for m in FIG1_METHODS},
    "lso_v10": {m: sealed_lso["primary"]["pooled"][m]["auc"] for m in FIG1_METHODS},
}

log("Building regime: released (candidates_purged_v7.csv)")
_, _, po_released = build_released_v7()
log("Building regime: refit_2021 (features_refit2021_relabelled_v7.csv)")
_, po_refit = build_refit_2021_v7(log=log)
log("Building regime: lso_v10 (features_recomputed_lso_v10.csv) -- direct, same process, "
    "no cross-environment bridge needed now")
_, po_lso = build_lso_v10(log=log)

REGIMES = {"released": po_released, "refit_2021": po_refit, "lso_v10": po_lso}

rows = []
consistency_report = []
for regime, per_origin in REGIMES.items():
    pooled = pool_metrics(per_origin)
    log(f"Regime {regime}: pooled AUCs computed. Running {N_BOOT}-resample per-method bootstrap "
        f"(seed {BOOT_SEED})...")
    boot = run_bootstrap_all_methods(per_origin, FIG1_METHODS, n_boot=N_BOOT, seed=BOOT_SEED)
    for m in FIG1_METHODS:
        auc = pooled[m]["auc"]
        ci_lo, ci_hi = boot[m]["auc_ci95"]
        rows.append({"regime": regime, "method": m, "auc": auc, "ci_lo": ci_lo, "ci_hi": ci_hi})
        sealed_val = SEALED_POOLED[regime][m]
        match = abs(auc - sealed_val) < 1e-9
        consistency_report.append({
            "regime": regime, "method": m, "computed_auc": auc, "sealed_auc": sealed_val,
            "match": match,
        })
        if not match:
            log(f"  MISMATCH regime={regime} method={m}: computed={auc!r} sealed={sealed_val!r}")
    log(f"Regime {regime} done.")

fig1_df = pd.DataFrame(rows, columns=["regime", "method", "auc", "ci_lo", "ci_hi"])
fig1_path = os.path.join(OUT_DIR, "fig1_data.csv")
fig1_df.to_csv(fig1_path, index=False)
log(f"Wrote {fig1_path} ({len(fig1_df)} rows)")

n_mismatch = sum(1 for r in consistency_report if not r["match"])
log(f"Consistency check: Figure 1 AUCs vs sealed pooled AUCs (lso_v10 now checked against "
    f"results_repo_f7_lso_v10.json) -- {len(consistency_report)-n_mismatch}/"
    f"{len(consistency_report)} match exactly, {n_mismatch} mismatched.")
pd.DataFrame(consistency_report).to_csv(os.path.join(OUT_DIR, "fig1_consistency_check.csv"), index=False)
assert n_mismatch == 0, "Figure 1 consistency check failed"

# ---- plot (group_stride = group_gap + 1.3; long first-line tick labels are wrapped) ----
import matplotlib.pyplot as plt  # noqa: E402
apply_rcparams()
fig, ax = plt.subplots(figsize=(ONE_COL_W, 3.4))
ax.tick_params(axis="both", labelsize=8)

regimes_order = ["released", "refit_2021", "lso_v10"]
n_methods = len(FIG1_METHODS)
group_gap = 1.0
group_stride = group_gap + 1.3
within_gap = group_gap / (n_methods + 1)

for gi, regime in enumerate(regimes_order):
    if REGIME_POSTHOC_V7[regime]:
        base_x = gi * group_stride
        ax.axvspan(base_x - 0.05, base_x + group_gap + 0.05, color="gray", alpha=0.08, zorder=-1)

for gi, regime in enumerate(regimes_order):
    base_x = gi * group_stride
    for mi, m in enumerate(FIG1_METHODS):
        row = fig1_df[(fig1_df.regime == regime) & (fig1_df.method == m)].iloc[0]
        x = base_x + (mi + 1) * within_gap
        yerr_lo = max(0.0, row["auc"] - row["ci_lo"])
        yerr_hi = max(0.0, row["ci_hi"] - row["auc"])
        ax.errorbar(x, row["auc"], yerr=[[yerr_lo], [yerr_hi]], fmt="o",
                    color=METHOD_COLOR[m], markersize=3.5, capsize=2, linewidth=1,
                    label=METHOD_LABEL[m] if gi == 0 else None)

ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, zorder=0)

tick_x = [gi * group_stride + group_gap / 2 for gi in range(len(regimes_order))]


def wrap_label(lbl):
    lbl = lbl.replace(" (", "\n(", 1)
    first_line = lbl.split("\n")[0]
    if len(first_line) > 15:
        words = first_line.split(" ")
        mid = len(words) // 2
        lbl = " ".join(words[:mid]) + "\n" + " ".join(words[mid:]) + lbl[len(first_line):]
    return lbl


tick_labels = [wrap_label(REGIME_LABEL_V7[r]) for r in regimes_order]
ax.set_xticks(tick_x)
ax.set_xticklabels(tick_labels, fontsize=8)
ax.set_ylabel("Pooled AUC", fontsize=8)
ax.set_ylim(0, 1.02)

ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False,
          fontsize=8, handletextpad=0.4, columnspacing=1.0)

save_fig(fig, os.path.join(OUT_DIR, "fig1"))
plt.close(fig)
log("Wrote fig1.pdf / fig1.png")

log(f"DONE. Total runtime {time.time()-t0:.1f}s")
