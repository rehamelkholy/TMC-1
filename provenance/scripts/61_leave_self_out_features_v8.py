# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Fix the Amendment 7 M1 self-fit defect re-introduced for the 10
species added in the v6 update (candidate_id 268/272/282/298/383/484/493/524/535/765), and recompute
the leave-self-out features as v8 (features_recomputed_lso_v8.csv -- NOT sealed, like the
v7 features file).

CAUSE (see provenance/README.md): the v6 update scripts set new positive/date
labels in candidates_purged_v6.csv and performed anchor-dup removal on
lee_list_entry_labels_purged_v6.csv, but never populated that raw-list file's pre-existing
`matched_detection_formula` column for the 10 new species' own raw ids. src/umda_pipeline/leave_self_out.py's
build_detection_records_tracked() uses `matched_detection_formula` to tag a raw id's
distance/training point ("raw", raw_id) when it's the source of a detection record; an
unmatched raw id instead produces an ("det", formula)-tagged point sourced from the
detection row's OWN smiles. loo_filter() only ever strips ("raw", rid)-tagged points
(matching Amendment 7 M2's documented scope), so the 10 new species -- whose raw ids all
had `matched_detection_formula` = NaN -- fell into the "unmatched" branch and were never
removable from their own reference/training set during their own leave-self-out refit,
reproducing the pre-M2 self-fit (dist=0) defect. The 47 pre-existing positives were
unaffected because their raw ids already had `matched_detection_formula` populated by the
original matching pipeline.

FIX: this script patches a LOCAL in-memory copy of `raw` (loaded from
lee_list_entry_labels_purged_v6.csv, on-disk file untouched) so that exactly the 10 raw ids
below get `matched_detection_formula` set to their own candidate's detection formula string
(verified character-identical to the corresponding post2021_detections_v6.csv row's
`formula` value), before calling the UNCHANGED leave_self_out.py functions. This is everything
else identical to the version-7 leave-self-out build (feature_recompute.py / leave_self_out.py
imported unmodified; same anchors, same post2021_detections_v6.csv, same
candidates_purged_v6.csv skeleton).

Raw ids patched (candidate_id -> raw id(s) -> formula):
  268 (NO)    -> 275  -> "NO"
  272 (NS)    -> 279  -> "NS"
  282 (CH2N+) -> 290  -> "CH2N+"
  298 (H2S)   -> 307  -> "H2S"
  383 (CHO+)  -> 396  -> "CHO+"
  484 (C2H)   -> 503  -> "C2H"
  493 (CHO)   -> 512  -> "CHO"
  524 (C2HNO) -> 545  -> "C2HNO"
  535 (CHNS)  -> 556  -> "CHNS"
  765 (CO)    -> 800  -> "CO"   (NOT 801/13CO: post2021_detections_v6.csv's single CO
                                  detection row's smiles "[C-]#[O+]" matches raw id 800's
                                  smiles exactly, not raw id 801's "[13C-]#[O+]". This
                                  mirrors the existing, pre-v6 convention already used for
                                  candidate 1143/C4H2N2, where only ONE of two member raw
                                  ids [1193, not 1194] carries matched_detection_formula --
                                  the other member simply never enters the known/anchor set
                                  as a distance point in the first place, so it needs no
                                  separate removal. Candidate 765's own leave-self-out
                                  removal set is still {800, 801} [both members, from
                                  member_ids], so if 801 ever did enter the anchor set this
                                  fix would still strip it; it simply never does, exactly as
                                  for 1194.)

Must be run with the project's dedicated Python 3.10 venv (umda embeddings):
  <PROJECT_ROOT>/venv/bin/python3.10 step5_build_lso_features_v8.py
"""
import os
import sys
import time
import json
import hashlib
import warnings
warnings.filterwarnings("ignore")

ROOT = "<PROJECT_ROOT>"
sys.path.insert(0, os.path.join(ROOT, "active_learning"))
sys.path.insert(0, os.path.join(ROOT, "active_learning", "candidate_pool_sources"))
sys.path.insert(0, os.path.join(ROOT, "f4a_implementation"))
sys.path.insert(0, os.path.join(ROOT, "f4c_implementation"))
sys.path.insert(0, os.path.join(ROOT, "f6_implementation"))
sys.path.insert(0, os.path.join(ROOT, "f7_implementation"))
sys.path.insert(0, os.path.join(ROOT, "f7d_implementation"))

import numpy as np
import pandas as pd

from f7_common import load_anchors_v2, embed_raw_and_anchors, fit_gpr_and_score, ORIGINS  # unchanged
from f7d_common import build_detection_records_tracked, build_K3_for_origin_tracked, loo_filter  # unchanged

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


def load_purged_raw_v6():
    return pd.read_csv(os.path.join(ROOT, "f13_implementation", "lee_list_entry_labels_purged_v6.csv"))


def load_post2021_v6():
    return pd.read_csv(os.path.join(ROOT, "f13_implementation", "post2021_detections_v6.csv"))


# raw id -> matched_detection_formula patch (item 1 fix; see module docstring for the cause
# and why 801/13CO is deliberately NOT included here).
RAW_ID_FORMULA_PATCH = {
    275: "NO",
    279: "NS",
    290: "CH2N+",
    307: "H2S",
    396: "CHO+",
    503: "C2H",
    512: "CHO",
    545: "C2HNO",
    556: "CHNS",
    800: "CO",
}

anchors = load_anchors_v2()
raw = load_purged_raw_v6()
det = load_post2021_v6()
log(f"Loaded {len(anchors)} anchors (v2, unchanged), {len(raw)} purged raw entries (v6), "
    f"{len(det)} post-2021 detections (v6, 94 rows).")

# --- item 1 fix: patch the LOCAL in-memory copy only; on-disk CSV is untouched. ---
n_before = int(raw["matched_detection_formula"].notna().sum())
for rid, formula in RAW_ID_FORMULA_PATCH.items():
    mask = raw["id"].astype(int) == rid
    assert mask.sum() == 1, f"expected exactly 1 raw row for id {rid}, got {mask.sum()}"
    before_val = raw.loc[mask, "matched_detection_formula"].iloc[0]
    assert pd.isna(before_val), f"raw id {rid} unexpectedly already has matched_detection_formula={before_val!r}"
    raw.loc[mask, "matched_detection_formula"] = formula
n_after = int(raw["matched_detection_formula"].notna().sum())
log(f"Patched matched_detection_formula (in-memory only) for {len(RAW_ID_FORMULA_PATCH)} raw ids: "
    f"{n_before} -> {n_after} non-null (delta {n_after - n_before}).")
assert n_after - n_before == len(RAW_ID_FORMULA_PATCH)

X_raw, ok_raw, X_anc, ok_anc, embedder = embed_raw_and_anchors(raw, anchors, anchor_smiles_col="embedding_smiles")
log(f"Embedded {ok_raw.sum()}/{len(raw)} raw entries, {ok_anc.sum()}/{len(anchors)} anchors.")
assert ok_anc.all()

detection_records = build_detection_records_tracked(raw, det, embedder)
n_matched = sum(1 for r in detection_records if r["source"] == "matched")
n_unmatched = sum(1 for r in detection_records if r["source"] == "unmatched")
log(f"Built {len(detection_records)} detection records ({n_matched} matched, {n_unmatched} unmatched).")

cand_skel = pd.read_csv(os.path.join(ROOT, "f13_implementation", "candidates_purged_v6.csv"))
cand_skel["date"] = pd.to_datetime(cand_skel["date"])
log(f"F13 candidate skeleton (v6): {len(cand_skel)} candidates, {int(cand_skel['positive'].sum())} positive.")

raw_id_to_cand = {}
cand_member_ids = {}
for _, c in cand_skel.iterrows():
    members = [int(x) for x in str(c["member_ids"]).split(";")]
    cand_member_ids[int(c["candidate_id"])] = members
    for m in members:
        raw_id_to_cand[m] = int(c["candidate_id"])

raw_ids_ok = raw.loc[ok_raw, "id"].astype(int).to_numpy()
raw_pos_by_id = {rid: pos for pos, rid in enumerate(raw_ids_ok)}

ORIGIN_TS = [pd.Timestamp(o) for o in ORIGINS]

all_rows = []
per_origin_loo_log = []

for oi, origin in enumerate(ORIGIN_TS):
    log(f"=== Origin {oi} ({origin.date()}) ===")
    dist_coords, dist_ids, extra_X, extra_y, train_ids, distance_only, n_known = build_K3_for_origin_tracked(
        X_anc[ok_anc], detection_records, origin, oi,
    )
    n_train_extra = 0 if extra_X is None else extra_X.shape[0]
    log(f"  K3 known detections: {n_known}; distance-only: {distance_only}; "
        f"synthetic training points: {n_train_extra}; total dist anchors: {dist_coords.shape[0]}")

    mean_all, std_all, dist_all, lml, kernel_str, X_train = fit_gpr_and_score(
        X_query=X_raw[ok_raw], K3_anchor_coords=dist_coords, extra_train_X=extra_X, extra_train_y=extra_y,
    )
    log(f"  Baseline GPR: training set size {X_train.shape[0]}; LML {lml:.4f}")
    dist_by_id = dict(zip(raw_ids_ok, dist_all))
    mean_by_id = dict(zip(raw_ids_ok, mean_all))
    sigma_by_id = dict(zip(raw_ids_ok, std_all))

    origin_str = origin.strftime("%Y-%m-%d")
    known_mask = (cand_skel["positive"].to_numpy()) & (cand_skel["date"] < origin).to_numpy()
    known_cand_ids = cand_skel.loc[known_mask, "candidate_id"].astype(int).tolist()
    log(f"  Known (positive & date<origin) candidates this origin: {len(known_cand_ids)}")

    t_loo_start = time.time()
    n_loo_done = 0
    for _, c in cand_skel.iterrows():
        cid = int(c["candidate_id"])
        member_ids = cand_member_ids[cid]
        is_known = cid in known_cand_ids

        if not is_known:
            dists = [dist_by_id[i] for i in member_ids if i in dist_by_id]
            means = [mean_by_id[i] for i in member_ids if i in mean_by_id]
            sigmas = [sigma_by_id[i] for i in member_ids if i in sigma_by_id]
            if dists:
                d_val, m_val, s_val = float(np.mean(dists)), float(np.mean(means)), float(np.mean(sigmas))
            else:
                d_val = m_val = s_val = float("nan")
        else:
            remove_set = set(member_ids)
            dc_loo, ex_loo, ey_loo, n_rm_dist, n_rm_train = loo_filter(
                dist_coords, dist_ids, extra_X, extra_y, train_ids, remove_set,
            )
            member_pos = [raw_pos_by_id[i] for i in member_ids if i in raw_pos_by_id]
            Xq = X_raw[ok_raw][member_pos]
            mean_c, std_c, dist_c, lml_c, kernel_c, Xtr_c = fit_gpr_and_score(
                X_query=Xq, K3_anchor_coords=dc_loo, extra_train_X=ex_loo, extra_train_y=ey_loo,
            )
            d_val, m_val, s_val = float(np.mean(dist_c)), float(np.mean(mean_c)), float(np.mean(std_c))
            n_loo_done += 1
            per_origin_loo_log.append({
                "origin": origin_str, "candidate_id": cid, "species_key": c["species_key"],
                "n_removed_dist_points": n_rm_dist, "n_removed_train_points": n_rm_train,
                "loo_training_set_size": int(Xtr_c.shape[0]), "loo_lml": float(lml_c),
            })
            if n_loo_done % 10 == 0:
                log(f"    LOO refit {n_loo_done}/{len(known_cand_ids)} (elapsed {time.time()-t_loo_start:.1f}s)")

        all_rows.append({
            "origin": origin_str, "species_key": c["species_key"], "candidate_id": cid,
            "is_positive": bool(c["positive"]),
            "date": c["date"].strftime("%Y-%m-%d") if pd.notna(c["date"]) else "",
            "dist": d_val, "sigma": s_val, "logcd": m_val,
            "mu": float(c["mu"]), "mw": float(c["mw"]),
            "charged": bool(c["charged"]), "radical": bool(c["radical"]),
        })
    log(f"  LOO refits this origin: {n_loo_done} (took {time.time()-t_loo_start:.1f}s)")

out = pd.DataFrame(all_rows)
out_path = os.path.join(ROOT, "f13_implementation", "features_recomputed_lso_v8.csv")
out.to_csv(out_path, index=False)
log(f"Wrote {out_path}: {len(out)} rows")

h = hashlib.sha256(open(out_path, "rb").read()).hexdigest()
# NOT sealed / no .sha256 written for this
# specific output -- hash still logged to stdout/log entry for provenance.
log(f"SHA-256 (not sealed, per task instruction): {h}")

with open(os.path.join(ROOT, "f13_implementation", "loo_refit_log_v8.json"), "w") as f:
    json.dump(per_origin_loo_log, f, indent=2)
log(f"Wrote loo_refit_log_v8.json ({len(per_origin_loo_log)} refits total)")

log(f"\nDONE. Total runtime {time.time()-t0:.1f}s")
