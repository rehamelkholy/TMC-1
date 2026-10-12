"""
Gates G1/G2 with protocol amendment 6 (L1): anchors embedded from table_a2_anchors.csv's
`embedding_smiles` column (table_a2_lee2021.TABLE_A2's SMILES) instead of the canonical `smiles`
column. Everything else (K2-K6) is as specified in protocol amendment 5; table_a2_anchors.csv remains
the identity reference.
"""
import os
import sys
import time
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src", "umda_pipeline"))

import numpy as np
import pandas as pd
from scipy import stats

from paths import BUILD, EXTERNAL
from feature_recompute import (
    load_anchors_v2, load_purged_raw, embed_raw_and_anchors, fit_gpr_and_score,
)

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


anchors = load_anchors_v2()
log(f"Loaded {len(anchors)} authoritative anchors (v2, with embedding_smiles).")

raw = load_purged_raw()
log(f"Loaded purged raw universe: {len(raw)} rows.")

X_raw, ok_raw, X_anc, ok_anc, embedder = embed_raw_and_anchors(raw, anchors, anchor_smiles_col="embedding_smiles")
log(f"Embedded {ok_raw.sum()}/{len(raw)} purged raw entries, {ok_anc.sum()}/{len(anchors)} anchors "
    f"(from embedding_smiles).")
assert ok_anc.all(), "expected all 87 anchors to embed from embedding_smiles"

mean_ok, std_ok, dist_ok, lml, kernel_str, X_train0 = fit_gpr_and_score(
    X_query=X_raw[ok_raw], K3_anchor_coords=X_anc[ok_anc], extra_train_X=None, extra_train_y=None,
)
log(f"Fitted kernel_: {kernel_str}")
log(f"log_marginal_likelihood: {lml:.6f}")
log(f"Training set size: {X_train0.shape[0]} (expect 643, no synthetic additions)")

REFERENCE_LML = -632.1523
lml_match = abs(lml - REFERENCE_LML) < 5e-3
log(f"G1a (LML match to reference -632.1523): {'PASS' if lml_match else 'FAIL'} (got {lml:.4f})")

ref = pd.read_csv(os.path.join(EXTERNAL, "lee_list_with_our_predictions.csv"))
ref = ref.set_index("id")

raw_ok = raw.loc[ok_raw].copy().reset_index(drop=True)
raw_ok["our_mean"] = mean_ok
raw_ok["our_sigma"] = std_ok
raw_ok["our_dist"] = dist_ok

merged = raw_ok.merge(ref[["our_mean", "our_sigma", "our_dist"]], left_on="id", right_index=True,
                       suffixes=("_new", "_ref"))
log(f"Matched {len(merged)}/{len(raw_ok)} purged entries against the reference file by id.")

diffs = {}
for col in ["our_mean", "our_sigma", "our_dist"]:
    d = (merged[f"{col}_new"] - merged[f"{col}_ref"]).abs()
    diffs[col] = d.max()
    log(f"  max |{col} diff| = {d.max():.3e}")

features_match = all(v < 1e-6 for v in diffs.values())
log(f"G1b (per-entry features match to 1e-6): {'PASS' if features_match else 'FAIL'}")

G1 = lml_match and features_match
log(f"G1 overall: {'PASS' if G1 else 'FAIL'}")

rho_sigma, p_sigma = stats.spearmanr(raw_ok["our_sigma"], raw_ok["sigma"])
rho_dist, p_dist = stats.spearmanr(raw_ok["our_dist"], raw_ok["dist"])
log(f"Spearman rho (sigma, ours vs Lee): {rho_sigma:.4f}")
log(f"Spearman rho (dist, ours vs Lee):  {rho_dist:.4f}")
G2 = (rho_sigma >= 0.95) and (rho_dist >= 0.95)
log(f"G2 overall: {'PASS' if G2 else 'FAIL'}")

gate_report = {
    "amendment": "6 (L1 applied: anchors embedded from table_a2_authoritative_v2.csv's embedding_smiles)",
    "lml": lml,
    "lml_reference": REFERENCE_LML,
    "lml_match": bool(lml_match),
    "kernel_": kernel_str,
    "training_set_size": int(X_train0.shape[0]),
    "max_abs_diff_mean": float(diffs["our_mean"]),
    "max_abs_diff_sigma": float(diffs["our_sigma"]),
    "max_abs_diff_dist": float(diffs["our_dist"]),
    "features_match_1e6": bool(features_match),
    "G1": bool(G1),
    "spearman_rho_sigma": float(rho_sigma),
    "spearman_rho_dist": float(rho_dist),
    "G2": bool(G2),
    "gates_passed": bool(G1 and G2),
}
import json
with open(os.path.join(BUILD, "gate_report_v2.json"), "w") as f:
    json.dump(gate_report, f, indent=2)
log("Wrote gate_report_v2.json")
print(json.dumps(gate_report, indent=2))

if not (G1 and G2):
    log("GATES FAILED -- stopping per K6. No evaluation will be run.")
    sys.exit(1)
log("Gates passed.")
