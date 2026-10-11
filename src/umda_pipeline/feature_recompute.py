"""
Shared machinery for the feature recompute (protocol amendment 5): embedding, GPR fit/score, and
K3/K4 construction. Used by scripts/features/01_anchor_gate_check.py (K3 = 87 anchors only) and by
the per-origin recompute (K3 = anchors + known post-2021 detections per origin); the gate validates
this code once, in the simple case, before it is trusted for the per-origin runs.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from paths import DATA_RAW, PROV_DATA  # noqa: E402

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi as rd_inchi
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RationalQuadratic, DotProduct, WhiteKernel

RDLogger.DisableLog("rdApp.*")

from embeddings import load_embedder
from loop import load_precutoff_training_set

from candidate_table import parse_one_date  # noqa: E402 -- I2: month-only -> 1st of month

ANCHORS_PATH = os.path.join(PROV_DATA, "table_a2_authoritative_v1.csv")
ANCHORS_V2_PATH = os.path.join(DATA_RAW, "table_a2_anchors.csv")
PURGED_RAW_PATH = os.path.join(PROV_DATA, "lee_list_entry_labels_purged_v1.csv")
POST2021_PATH = os.path.join(PROV_DATA, "post2021_detections_v3_2.csv")

RANDOM_STATE = 1215677  # same as the within-list reproduction run
NOISE_SEED_BASE = 20260921
NOISE_SIGMA = 0.5
N_SYNTH_PER_SPECIES = 9

ORIGINS = ["2022-01-01", "2023-01-01", "2024-01-01", "2025-01-01"]


def load_anchors():
    return pd.read_csv(ANCHORS_PATH)


def load_anchors_v2():
    """Amendment 6, L1: the 87 anchors with an added embedding_smiles
    column (table_a2_lee2021.TABLE_A2's SMILES -- the representation
    already validated to reproduce Lee's own published Table A2
    embeddings to floating-point precision). table_a2_authoritative.csv's
    own `smiles` column remains the identity/citation reference; embedding
    uses `embedding_smiles` exclusively from here on."""
    return pd.read_csv(ANCHORS_V2_PATH)


def load_purged_raw():
    return pd.read_csv(PURGED_RAW_PATH)


def load_post2021():
    return pd.read_csv(POST2021_PATH)


def inchikey_of(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_KEKULIZE,
            catchErrors=True,
        )
    try:
        return rd_inchi.InchiToInchiKey(rd_inchi.MolToInchi(mol, treatWarningAsError=False))
    except Exception:
        return None


def embed70(smi, embedder):
    try:
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            return None
        m.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(m)
        v = embedder.vectorize(smi)
        if v is None or np.all(v == 0):
            return None
        return np.asarray(v).reshape(-1)
    except Exception:
        return None


def embed_raw_and_anchors(raw, anchors, anchor_smiles_col="smiles"):
    """Embed every raw-universe SMILES and every anchor SMILES into the
    raw 70-d PCA space. Same embed70() logic as Lee et al.'s within-list reproduction.

    anchor_smiles_col: which column of `anchors` to embed from. Amendment
    6 (L1) requires this to be "embedding_smiles" (table_a2_lee2021's
    SMILES) rather than table_a2_authoritative.csv's own "smiles" column
    (RDKit-canonical, chemically identical by InChIKey but not always the
    same literal structure/resonance-form mol2vec embeds)."""
    embedder = load_embedder()

    X_raw = np.full((len(raw), 70), np.nan)
    ok_raw = np.zeros(len(raw), dtype=bool)
    for i, smi in enumerate(raw["smiles"]):
        v = embed70(smi, embedder)
        if v is not None:
            X_raw[i] = v
            ok_raw[i] = True

    X_anc = np.full((len(anchors), 70), np.nan)
    ok_anc = np.zeros(len(anchors), dtype=bool)
    for i, smi in enumerate(anchors[anchor_smiles_col]):
        v = embed70(smi, embedder)
        if v is not None:
            X_anc[i] = v
            ok_anc[i] = True

    return X_raw, ok_raw, X_anc, ok_anc, embedder


def embed_species_list(smiles_list, embedder):
    """Embed an arbitrary list of SMILES (e.g. post2021 known-chemistry
    species); returns (X (n,70) with NaN rows for failures, ok mask)."""
    X = np.full((len(smiles_list), 70), np.nan)
    ok = np.zeros(len(smiles_list), dtype=bool)
    for i, smi in enumerate(smiles_list):
        v = embed70(smi, embedder)
        if v is not None:
            X[i] = v
            ok[i] = True
    return X, ok


def build_detection_records(raw, det, embedder):
    """Amendment 6, L2: for each of the post2021_detections_v3_2.csv rows,
    work out (a) its distance-anchor embedding(s), (b) its single GPR-
    training embedding, (c) its column density (None if not recorded), and
    (d) its parsed first_public_date (I2 month->1st-of-month convention,
    same parse_one_date this project uses everywhere else).

    L2: "a detection matched to Lee-list entries by rules A/B/C is embedded
    from the SMILES of its matched entries: all matched entries serve as
    distance anchors, and the entry with the lowest Lee id supplies the GPR
    training embedding. An unmatched detection is embedded from its SMILES
    in post2021_detections_v3_2.csv."

    raw: the purged raw universe (with its own `matched_detection_formula`
    column, semicolon-joined when one raw row ties to >1 detection formula,
    e.g. the cis-/trans-CH3CHCHCN case).
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
        cd = d["log10_column_density_cm2"]
        cd = None if pd.isna(cd) else float(cd)

        matched_ids = sorted(formula_to_ids.get(formula, []))
        if matched_ids:
            dist_smiles = [raw_smiles_by_id[i] for i in matched_ids]
            train_smiles = raw_smiles_by_id[matched_ids[0]]  # lowest Lee id
            source = "matched"
        else:
            dist_smiles = [d["smiles"]]
            train_smiles = d["smiles"]
            source = "unmatched"

        dist_coords = [embed70(s, embedder) for s in dist_smiles]
        dist_coords = [c for c in dist_coords if c is not None]
        train_coord = embed70(train_smiles, embedder)

        records.append({
            "formula": formula,
            "parsed_date": parsed_date,
            "date_precision": precision,
            "cd": cd,
            "matched_ids": matched_ids,
            "source": source,
            "dist_coords": dist_coords,   # list of (70,) arrays, len>=1 if any embedded
            "train_coord": train_coord,   # (70,) array or None if embedding failed
        })
    return records


def build_K3_for_origin(anchors_X, detection_records, origin, origin_index):
    """K3(o) = 87 anchors + every detection with parsed_date < o (Amendment
    5 K3, using L2's embeddings). Returns:
      dist_coords: (m,70) array -- every anchor + every distance-anchor
        point for every known detection (matched detections contribute one
        point per matched Lee entry; unmatched contribute one).
      extra_train_X, extra_train_y: GPR training additions -- 9 synthetic
        points per known non-anchor species THAT HAS a recorded column
        density (K4), embedding = that species' single train_coord, y =
        cd + Gaussian noise (sigma=0.5), seed = 20260921 + origin_index,
        drawn from one RNG per origin in a fixed (formula-sorted) order.
      distance_only: list of formulas known at this origin but excluded
        from GPR training (no recorded column density).
      n_known_detections: count of detections known at this origin.
    """
    known = [r for r in detection_records if r["parsed_date"] is not None and r["parsed_date"] < origin]
    known = sorted(known, key=lambda r: r["formula"])  # fixed order for reproducible noise draws

    dist_list = [anchors_X]
    for r in known:
        for c in r["dist_coords"]:
            dist_list.append(c.reshape(1, -1))
    dist_coords = np.vstack(dist_list)

    rng = np.random.default_rng(NOISE_SEED_BASE + origin_index)
    train_X_rows = []
    train_y_rows = []
    distance_only = []
    for r in known:
        if r["cd"] is None or r["train_coord"] is None:
            distance_only.append(r["formula"])
            continue
        noise = rng.normal(loc=0.0, scale=NOISE_SIGMA, size=N_SYNTH_PER_SPECIES)
        for eps in noise:
            train_X_rows.append(r["train_coord"])
            train_y_rows.append(r["cd"] + eps)

    extra_train_X = np.vstack(train_X_rows) if train_X_rows else None
    extra_train_y = np.array(train_y_rows) if train_y_rows else None

    return dist_coords, extra_train_X, extra_train_y, distance_only, len(known)


def fit_gpr_and_score(X_query, K3_anchor_coords, extra_train_X=None, extra_train_y=None):
    """K4 + K5's GPR/distance step, for one origin (or the G1/G2 anchors-
    only special case, called with extra_train_X=None).

    X_query: (n, 70) raw-PCA-70 coords to score (dist/mean/sigma).
    K3_anchor_coords: (m, 70) coords of EVERY species in K3(o) (anchors +
        known post-2021 species, including any that are distance-only --
        i.e. lack a column density and so are excluded from extra_train_*
        but still count for "nearest species in K3").
    extra_train_X, extra_train_y: additional GPR training points (the 9x
        synthetic observations per non-anchor K3 species with a recorded
        column density), concatenated onto the 643 shipped points. None
        for the anchors-only gate case, where there are none by
        construction.

    Returns: mean, std, dist (each (n,)), lml, kernel_str, X_train (full
        GPR training matrix actually used, for provenance/size checks).
    """
    X_train0, y_train0 = load_precutoff_training_set()
    if extra_train_X is not None and len(extra_train_X) > 0:
        X_train = np.vstack([X_train0, extra_train_X])
        y_train = np.concatenate([y_train0, extra_train_y])
    else:
        X_train = X_train0
        y_train = y_train0

    kernel_init = RationalQuadratic(alpha=1, length_scale=1) + DotProduct(sigma_0=1) + WhiteKernel(noise_level=1)
    gpr = GaussianProcessRegressor(kernel=kernel_init, random_state=RANDOM_STATE)
    gpr.fit(X_train, y_train)
    lml = gpr.log_marginal_likelihood(gpr.kernel_.theta)

    mean, std = gpr.predict(X_query, return_std=True)

    d2 = ((X_query[:, None, :] - K3_anchor_coords[None, :, :]) ** 2).sum(axis=2)
    dist = np.sqrt(np.maximum(d2, 0.0)).min(axis=1)

    return mean, std, dist, float(lml), str(gpr.kernel_), X_train
