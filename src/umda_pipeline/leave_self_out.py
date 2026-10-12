"""
Leave-self-out (LSO) variant of the K3-K5 recompute (protocol amendment 7, M2). Reuses
feature_recompute.py's embedding/GPR machinery unmodified (load_anchors_v2, load_purged_raw,
load_post2021, embed_raw_and_anchors, embed70, fit_gpr_and_score, ORIGINS, RANDOM_STATE,
NOISE_SEED_BASE, NOISE_SIGMA, N_SYNTH_PER_SPECIES) and adds ID-TRACKING variants of
build_detection_records/build_K3_for_origin so that, given an already-built K3(o), exactly one
known candidate's own contribution (all of its raw ids' distance-anchor points, and its GPR
training points if it supplied them) can be identified and removed without perturbing any other
species' points -- in particular without perturbing the other known species' 9 fixed noise draws,
which stay bit-identical because the already-generated arrays are filtered rather than re-drawn.

The M2 fix concerns only species that are POSITIVE candidates (matched to a post-2021 detection via
rules A/B/C; "positive = any member in_list_by_v3_1" in candidate_table.py) -- the 87 Lee Table A2
anchors were purged out of the raw universe by rule A and so never correspond to a raw id /
candidate member_id; they cannot be a "known positive candidate" and are never a removal target.
Every "known" (positive & date<origin) candidate maps back to >=1 raw id that appears in some
detection record's matched_ids.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np
import pandas as pd

from feature_recompute import (  # noqa: F401 -- re-exported unmodified for callers
    load_anchors_v2, load_purged_raw, load_post2021, embed_raw_and_anchors,
    embed70, fit_gpr_and_score, ORIGINS, RANDOM_STATE, NOISE_SEED_BASE,
    NOISE_SIGMA, N_SYNTH_PER_SPECIES,
)

from candidate_table import parse_one_date  # noqa: E402 -- I2: month-only -> 1st of month


def build_detection_records_tracked(raw, det, embedder):
    """Identical to feature_recompute.build_detection_records, except every distance-anchor point and
    the training point carry an explicit id key so a later LOO filter can identify and remove exactly
    one candidate's own contribution:
      dist_ids[i]  = ("raw", raw_id) for a matched entry's i-th distance
                     point, or ("det", formula) for an unmatched detection's
                     (single) distance point.
      train_id     = the id key of whichever SMILES supplied train_coord
                     (same convention), or None if train_coord is None.
    """
    raw_matched = raw[raw["matched_detection_formula"].notna()].copy()
    formula_to_ids = {}
    for _, r in raw_matched.iterrows():
        for f in str(r["matched_detection_formula"]).split(";"):
            f = f.strip()
            formula_to_ids.setdefault(f, []).append(int(r["id"]))
    raw_smiles_by_id = dict(zip(raw["id"].astype(int), raw["smiles"]))

    records = []
    for _, d in det.iterrows():
        formula = d["formula"]
        parsed_date, precision = parse_one_date(str(d["first_public_date"]))
        cd = None if pd.isna(d["log10_column_density_cm2"]) else float(d["log10_column_density_cm2"])

        matched_ids = sorted(formula_to_ids.get(formula, []))
        if matched_ids:
            dist_smiles = [raw_smiles_by_id[i] for i in matched_ids]
            dist_id_candidates = [("raw", i) for i in matched_ids]
            train_smiles = raw_smiles_by_id[matched_ids[0]]
            train_id_candidate = ("raw", matched_ids[0])
            source = "matched"
        else:
            dist_smiles = [d["smiles"]]
            dist_id_candidates = [("det", formula)]
            train_smiles = d["smiles"]
            train_id_candidate = ("det", formula)
            source = "unmatched"

        dist_coords = []
        dist_ids = []
        for smi, idk in zip(dist_smiles, dist_id_candidates):
            c = embed70(smi, embedder)
            if c is not None:
                dist_coords.append(c)
                dist_ids.append(idk)
        train_coord = embed70(train_smiles, embedder)
        train_id = train_id_candidate if train_coord is not None else None

        records.append({
            "formula": formula,
            "parsed_date": parsed_date,
            "date_precision": precision,
            "cd": cd,
            "matched_ids": matched_ids,
            "source": source,
            "dist_coords": dist_coords,
            "dist_ids": dist_ids,
            "train_coord": train_coord,
            "train_id": train_id,
        })
    return records


def build_K3_for_origin_tracked(anchors_X, detection_records, origin, origin_index):
    """Identical to feature_recompute.build_K3_for_origin, except it also returns
    parallel id-arrays for dist_coords and extra_train_X so a caller can
    filter out one species' own contribution post hoc, leaving every other
    species' points -- including their exact fixed noise draws -- untouched.

    anchor points get id key ("anchor", None) (never a removal target: no
    anchor is ever a "known positive candidate", see module docstring).

    Returns: dist_coords (m,70), dist_ids (list len m), extra_train_X
    (n,70) or None, extra_train_y (n,) or None, train_ids (list len n) or
    None, distance_only (list of formulas), n_known_detections (int).
    """
    known = [r for r in detection_records if r["parsed_date"] is not None and r["parsed_date"] < origin]
    known = sorted(known, key=lambda r: r["formula"])

    dist_list = [anchors_X]
    dist_ids = [("anchor", None)] * anchors_X.shape[0]
    for r in known:
        for c, idk in zip(r["dist_coords"], r["dist_ids"]):
            dist_list.append(c.reshape(1, -1))
            dist_ids.append(idk)
    dist_coords = np.vstack(dist_list)
    assert len(dist_ids) == dist_coords.shape[0]

    rng = np.random.default_rng(NOISE_SEED_BASE + origin_index)
    train_X_rows = []
    train_y_rows = []
    train_ids = []
    distance_only = []
    for r in known:
        if r["cd"] is None or r["train_coord"] is None:
            distance_only.append(r["formula"])
            continue
        noise = rng.normal(loc=0.0, scale=NOISE_SIGMA, size=N_SYNTH_PER_SPECIES)
        for eps in noise:
            train_X_rows.append(r["train_coord"])
            train_y_rows.append(r["cd"] + eps)
            train_ids.append(r["train_id"])

    extra_train_X = np.vstack(train_X_rows) if train_X_rows else None
    extra_train_y = np.array(train_y_rows) if train_y_rows else None

    return dist_coords, dist_ids, extra_train_X, extra_train_y, train_ids, distance_only, len(known)


def loo_filter(dist_coords, dist_ids, extra_train_X, extra_train_y, train_ids, remove_raw_ids):
    """Filter out every point whose id key is ("raw", rid) with rid in
    remove_raw_ids (a candidate's own member_ids). Anchor points (id key
    ("anchor", None)) and unmatched-detection points (("det", formula)) are
    never removed by this filter, matching Amendment 7 M2's scope (LSO
    applies only to a known candidate's OWN raw-id contributions).

    Returns the same 4-tuple shape (dist_coords, extra_train_X,
    extra_train_y) filtered, plus n_removed_dist / n_removed_train counts
    for logging/verification.
    """
    keep_dist = np.array([not (idk[0] == "raw" and idk[1] in remove_raw_ids) for idk in dist_ids])
    dist_coords_loo = dist_coords[keep_dist]
    n_removed_dist = int((~keep_dist).sum())

    if extra_train_X is None:
        return dist_coords_loo, None, None, n_removed_dist, 0

    keep_train = np.array([
        not (idk is not None and idk[0] == "raw" and idk[1] in remove_raw_ids) for idk in train_ids
    ])
    n_removed_train = int((~keep_train).sum())
    if keep_train.all():
        return dist_coords_loo, extra_train_X, extra_train_y, n_removed_dist, n_removed_train
    extra_train_X_loo = extra_train_X[keep_train] if keep_train.any() else None
    extra_train_y_loo = extra_train_y[keep_train] if keep_train.any() else None
    return dist_coords_loo, extra_train_X_loo, extra_train_y_loo, n_removed_dist, n_removed_train
