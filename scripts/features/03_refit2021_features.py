"""
features_refit_2021.csv (protocol amendment 8, N1) -- the G1 anchors-only fit (87 anchors,
embedding_smiles per amendment 6/L1, 643 shipped precutoff training observations, NO post-2021
detections added; i.e. exactly gate G1's own K3/K4/K5 special case, reusing
feature_recompute.fit_gpr_and_score unmodified) applied identically at all 4 origins (held fixed, per
N1: isolate the refit-vs-Lee's-released-values effect from the effect of updating features with new
detections).

Check: this script's own fit is run TWICE independently (self-consistency/determinism check on the
per-raw-entry mean/sigma/dist arrays), and its LML + fitted kernel_ are compared against the sealed G1
values in gate_report_v2.json (written by 01_anchor_gate_check.py); the same computation, run again
here, must reproduce those recorded gate numbers.
"""
import os
import sys
import time
import json
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src", "umda_pipeline"))

import numpy as np
import pandas as pd

import importlib.util

from paths import ROOT, BUILD
from feature_recompute import (
    load_anchors_v2, load_purged_raw, embed_raw_and_anchors, fit_gpr_and_score, ORIGINS,
)
# build_candidates_purged() is defined in provenance/scripts/13_candidates_purged.py (not a valid module name)
_spec = importlib.util.spec_from_file_location(
    "build_candidates_purged", os.path.join(ROOT, "provenance", "scripts", "13_candidates_purged.py"))
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
build_candidates_purged = _mod.build_candidates_purged

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


anchors = load_anchors_v2()
raw = load_purged_raw()
log(f"Loaded {len(anchors)} anchors (v2), {len(raw)} purged raw entries.")

X_raw, ok_raw, X_anc, ok_anc, embedder = embed_raw_and_anchors(raw, anchors, anchor_smiles_col="embedding_smiles")
log(f"Embedded {ok_raw.sum()}/{len(raw)} raw entries, {ok_anc.sum()}/{len(anchors)} anchors.")
assert ok_anc.all()

# --- run 1 ---
mean1, std1, dist1, lml1, kernel1, X_train1 = fit_gpr_and_score(
    X_query=X_raw[ok_raw], K3_anchor_coords=X_anc[ok_anc], extra_train_X=None, extra_train_y=None,
)
log(f"Run 1: training set {X_train1.shape[0]}, LML {lml1:.6f}, kernel_: {kernel1}")

# --- run 2 (self-consistency / determinism check) ---
mean2, std2, dist2, lml2, kernel2, X_train2 = fit_gpr_and_score(
    X_query=X_raw[ok_raw], K3_anchor_coords=X_anc[ok_anc], extra_train_X=None, extra_train_y=None,
)
log(f"Run 2: training set {X_train2.shape[0]}, LML {lml2:.6f}, kernel_: {kernel2}")

max_diff_mean = float(np.abs(mean1 - mean2).max())
max_diff_std = float(np.abs(std1 - std2).max())
max_diff_dist = float(np.abs(dist1 - dist2).max())
self_consistent = (max_diff_mean < 1e-9) and (max_diff_std < 1e-9) and (max_diff_dist < 1e-9) \
    and (kernel1 == kernel2) and (abs(lml1 - lml2) < 1e-9)
log(f"Self-consistency (run1 vs run2): max|diff| mean={max_diff_mean:.3e}, sigma={max_diff_std:.3e}, "
    f"dist={max_diff_dist:.3e}; kernel_ identical: {kernel1 == kernel2}; "
    f"{'PASS' if self_consistent else 'FAIL'}")

# --- compare against the already-sealed G1 gate numbers (gate_report_v2.json) ---
with open(os.path.join(BUILD, "gate_report_v2.json")) as f:
    gate_report_v2 = json.load(f)
diff_lml_vs_gate = abs(lml1 - gate_report_v2["lml"])
kernel_matches_gate = (kernel1 == gate_report_v2["kernel_"])
gate_match = (diff_lml_vs_gate < 1e-9) and kernel_matches_gate
log(f"vs. gate_report_v2.json: |LML diff| = {diff_lml_vs_gate:.3e}; kernel_ matches: {kernel_matches_gate}; "
    f"{'PASS' if gate_match else 'FAIL'}")

check_passed = bool(self_consistent and gate_match)
check_report = {
    "description": "features_refit_2021.csv's underlying fit vs itself (determinism) and vs the sealed G1 gate "
                    "(gate_report_v2.json, Task F7b)",
    "run1_lml": lml1, "run2_lml": lml2, "gate_report_v2_lml": gate_report_v2["lml"],
    "run1_kernel_": kernel1, "run2_kernel_": kernel2, "gate_report_v2_kernel_": gate_report_v2["kernel_"],
    "self_consistency_max_abs_diff_mean": max_diff_mean,
    "self_consistency_max_abs_diff_sigma": max_diff_std,
    "self_consistency_max_abs_diff_dist": max_diff_dist,
    "self_consistency_kernel_identical": bool(kernel1 == kernel2),
    "self_consistency_passed": bool(self_consistent),
    "vs_gate_report_v2_abs_lml_diff": diff_lml_vs_gate,
    "vs_gate_report_v2_kernel_matches": bool(kernel_matches_gate),
    "vs_gate_report_v2_passed": bool(gate_match),
    "tolerance": 1e-9,
    "check_passed": check_passed,
}
with open(os.path.join(BUILD, "g1_features_check.json"), "w") as f:
    json.dump(check_report, f, indent=2)
print(json.dumps(check_report, indent=2))

if not check_passed:
    log("STEP 2 CHECK FAILED -- stopping. No features_refit_2021.csv sealed evaluation will be run.")
    sys.exit(1)

# --- build candidate-level features (J1-collapse), replicated across all 4 origins ---
cand_skel, lee_from_skel = build_candidates_purged()
log(f"F6/F7 candidate skeleton: {len(cand_skel)} candidates.")

raw_ids_ok = raw.loc[ok_raw, "id"].astype(int).to_numpy()
dist_by_id = dict(zip(raw_ids_ok, dist1))
mean_by_id = dict(zip(raw_ids_ok, mean1))
sigma_by_id = dict(zip(raw_ids_ok, std1))

base_rows = []
for _, c in cand_skel.iterrows():
    member_ids = [int(x) for x in str(c["member_ids"]).split(";")]
    dists = [dist_by_id[i] for i in member_ids if i in dist_by_id]
    means = [mean_by_id[i] for i in member_ids if i in mean_by_id]
    sigmas = [sigma_by_id[i] for i in member_ids if i in sigma_by_id]
    if dists:
        d_val, m_val, s_val = float(np.mean(dists)), float(np.mean(means)), float(np.mean(sigmas))
    else:
        d_val = m_val = s_val = float("nan")
    base_rows.append({
        "species_key": c["species_key"],
        "candidate_id": int(c["candidate_id"]),
        "is_positive": bool(c["positive"]),
        "date": c["date"] if pd.notna(c["date"]) else "",
        "dist": d_val, "sigma": s_val, "logcd": m_val,
        "mu": float(c["mu"]), "mw": float(c["mw"]),
        "charged": bool(c["charged"]), "radical": bool(c["radical"]),
    })
base_df = pd.DataFrame(base_rows)

ORIGIN_STRS = [pd.Timestamp(o).strftime("%Y-%m-%d") for o in ORIGINS]
all_rows = []
for ostr in ORIGIN_STRS:
    for _, r in base_df.iterrows():
        row = r.to_dict()
        row["origin"] = ostr
        all_rows.append(row)

out = pd.DataFrame(all_rows)[["origin", "species_key", "candidate_id", "is_positive", "date",
                               "dist", "sigma", "logcd", "mu", "mw", "charged", "radical"]]
out_path = os.path.join(BUILD, "features_refit2021.csv")
out.to_csv(out_path, index=False)
log(f"Wrote {out_path}: {len(out)} rows ({len(cand_skel)} candidates x {len(ORIGIN_STRS)} origins, identical "
    f"dist/logcd/sigma at every origin by construction).")

n_nan = out["dist"].isna().sum()
log(f"NaN feature rows: {n_nan}")

import hashlib
h = hashlib.sha256(open(out_path, "rb").read()).hexdigest()
with open(out_path + ".sha256", "w") as f:
    f.write(f"{h}  features_refit_2021.csv\n")
log(f"SHA-256: {h}")
log(f"\nDONE. Total runtime {time.time()-t0:.1f}s")
