# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Build features_recomputed.csv --
one row per (origin, candidate), dist/logcd/sigma recomputed per Amendment
5 (K3-K5) + Amendment 6 (L1 anchor embeddings, L2 detection embeddings).
mu/mw/charged/radical/positive/date are NOT origin-varying (K5: "stay as
released") and are taken unchanged from the original (v1) candidate skeleton
(build_candidates_purged(), imported unmodified) so candidate_id/species_key
match candidates_purged_v1.csv exactly.
"""
import os
import sys
import time
import warnings
warnings.filterwarnings("ignore")

ROOT = "<PROJECT_ROOT>"
sys.path.insert(0, os.path.join(ROOT, "active_learning"))
sys.path.insert(0, os.path.join(ROOT, "active_learning", "candidate_pool_sources"))
sys.path.insert(0, os.path.join(ROOT, "f4a_implementation"))
sys.path.insert(0, os.path.join(ROOT, "f4c_implementation"))
sys.path.insert(0, os.path.join(ROOT, "f6_implementation"))
sys.path.insert(0, os.path.join(ROOT, "f7_implementation"))

import numpy as np
import pandas as pd

from f7_common import (
    load_anchors_v2, load_purged_raw, load_post2021, embed_raw_and_anchors,
    build_detection_records, build_K3_for_origin, fit_gpr_and_score, ORIGINS,
)
from build_candidates_purged import build_candidates_purged

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


anchors = load_anchors_v2()
raw = load_purged_raw()
det = load_post2021()
log(f"Loaded {len(anchors)} anchors (v2), {len(raw)} purged raw entries, {len(det)} post-2021 detections.")

X_raw, ok_raw, X_anc, ok_anc, embedder = embed_raw_and_anchors(raw, anchors, anchor_smiles_col="embedding_smiles")
log(f"Embedded {ok_raw.sum()}/{len(raw)} raw entries, {ok_anc.sum()}/{len(anchors)} anchors.")
assert ok_anc.all()
n_raw_embed_fail = (~ok_raw).sum()
if n_raw_embed_fail:
    log(f"WARNING: {n_raw_embed_fail} purged raw entries failed to embed -- "
        f"excluded from this recompute (dist/logcd/sigma left NaN); ids: "
        f"{raw.loc[~ok_raw, 'id'].tolist()}")

detection_records = build_detection_records(raw, det, embedder)
n_det_embed_fail = sum(1 for r in detection_records if r["train_coord"] is None)
log(f"Built {len(detection_records)} detection records "
    f"({sum(1 for r in detection_records if r['source']=='matched')} matched, "
    f"{sum(1 for r in detection_records if r['source']=='unmatched')} unmatched); "
    f"{n_det_embed_fail} failed to embed.")

cand_skel, lee_from_skel = build_candidates_purged()
log(f"F6 candidate skeleton: {len(cand_skel)} candidates (species_key/member_ids/positive/date/mu/mw/"
    f"charged/radical taken unchanged from here).")

id_to_pos = {int(rid): i for i, rid in enumerate(raw["id"].astype(int))}

ORIGIN_TS = [pd.Timestamp(o) for o in ORIGINS]

all_rows = []
distance_only_by_origin = {}
for oi, origin in enumerate(ORIGIN_TS):
    log(f"--- Origin {oi} ({origin.date()}) ---")
    dist_coords, extra_X, extra_y, distance_only, n_known = build_K3_for_origin(
        X_anc[ok_anc], detection_records, origin, oi,
    )
    n_train_extra = 0 if extra_X is None else extra_X.shape[0]
    log(f"  K3 known detections: {n_known}; distance-only (no CD): {len(distance_only)} {distance_only}; "
        f"GPR synthetic training points added: {n_train_extra}")

    mean_all, std_all, dist_all, lml, kernel_str, X_train = fit_gpr_and_score(
        X_query=X_raw[ok_raw], K3_anchor_coords=dist_coords, extra_train_X=extra_X, extra_train_y=extra_y,
    )
    log(f"  GPR training set size: {X_train.shape[0]}; LML: {lml:.4f}; kernel_: {kernel_str}")

    # map raw id -> (dist, logcd, sigma) for this origin
    raw_ids_ok = raw.loc[ok_raw, "id"].astype(int).to_numpy()
    dist_by_id = dict(zip(raw_ids_ok, dist_all))
    mean_by_id = dict(zip(raw_ids_ok, mean_all))
    sigma_by_id = dict(zip(raw_ids_ok, std_all))

    for _, c in cand_skel.iterrows():
        member_ids = [int(x) for x in str(c["member_ids"]).split(";")]
        dists = [dist_by_id[i] for i in member_ids if i in dist_by_id]
        means = [mean_by_id[i] for i in member_ids if i in mean_by_id]
        sigmas = [sigma_by_id[i] for i in member_ids if i in sigma_by_id]
        if not dists:  # every member failed to embed -- leave NaN, flagged
            d_val = m_val = s_val = float("nan")
        else:
            d_val = float(np.mean(dists))
            m_val = float(np.mean(means))
            s_val = float(np.mean(sigmas))

        all_rows.append({
            "origin": origin.strftime("%Y-%m-%d"),
            "species_key": c["species_key"],
            "candidate_id": int(c["candidate_id"]),
            "is_positive": bool(c["positive"]),
            "date": c["date"] if pd.notna(c["date"]) else "",
            "dist": d_val,
            "sigma": s_val,
            "logcd": m_val,
            "mu": float(c["mu"]),
            "mw": float(c["mw"]),
            "charged": bool(c["charged"]),
            "radical": bool(c["radical"]),
        })

    distance_only_by_origin[origin.strftime("%Y-%m-%d")] = distance_only

out = pd.DataFrame(all_rows)
out_path = os.path.join(ROOT, "f7_implementation", "features_recomputed.csv")
out.to_csv(out_path, index=False)
log(f"Wrote {out_path}: {len(out)} rows ({len(cand_skel)} candidates x {len(ORIGIN_TS)} origins).")

n_nan = out["dist"].isna().sum()
log(f"NaN dist/logcd/sigma rows (candidate had zero embeddable members this recompute): {n_nan}")

import hashlib
h = hashlib.sha256(open(out_path, "rb").read()).hexdigest()
with open(out_path + ".sha256", "w") as f:
    f.write(f"{h}  features_recomputed.csv\n")
log(f"SHA-256: {h}")

import json
with open(os.path.join(ROOT, "f7_implementation", "distance_only_species_by_origin.json"), "w") as f:
    json.dump(distance_only_by_origin, f, indent=2)
log("Wrote distance_only_species_by_origin.json")
print(json.dumps(distance_only_by_origin, indent=2))
log(f"\nDONE. Total runtime {time.time()-t0:.1f}s")
