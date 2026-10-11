"""
Embed every molecule in the candidate pool with the mol2vec -> PCA(70) pipeline of the umda
EmbeddingModel.pkl, so the resulting feature vectors live in exactly the same 70-dim space as the
fitted regressors' X_train_.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER1_ROOT = os.path.dirname(HERE)
sys.path.insert(0, PAPER1_ROOT)

from paths import EXTERNAL  # noqa: E402

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from joblib import load
from rdkit import Chem, RDLogger
from mol2vec import features

# mol2vec's mol2alt_sentence() (see features.py) calls the deprecated
# AllChem.GetMorganFingerprint() once per molecule. With a recent RDKit (2026.3.6) that logs a
# "DEPRECATION WARNING: please use MorganGenerator" line through RDKit's OWN logger (the
# Boost-log-based "rdApp.*" channel, not Python's `warnings` module, so
# `warnings.filterwarnings("ignore")` does not affect it) on every call, flooding the console for
# thousands of candidates. Silence RDKit's own logger, not just Python's.
RDLogger.DisableLog("rdApp.*")


def upgrade_legacy_keyedvectors(wv):
    d = wv.__dict__
    if "index_to_key" not in d:
        d["index_to_key"] = list(d["index2word"])
    if "key_to_index" not in d:
        d["key_to_index"] = {w: i for i, w in enumerate(d["index_to_key"])}
    return wv


def _find_umda_file(*parts):
    """Locate a file inside the umda checkout. The clone is expected at external/umda
    (see data/raw/FETCH.md); duplicated in loop.py, which needs the same resolution for regressors.pkl."""
    candidates = [
        os.path.join(EXTERNAL, "umda", *parts),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    tried = "\n  ".join(candidates)
    raise FileNotFoundError(
        f"Can't find {os.path.join(*parts)} in the umda checkout "
        f"-- tried:\n  {tried}\nClone "
        f"https://github.com/laserkelvin/umda.git (pinned commit in data/raw/FETCH.md) into external/umda "
        f"-- these are the umda repository's own shipped "
        f"pretrained artifacts (models/EmbeddingModel.pkl, "
        f"models/regressors.pkl), not something this project generates."
    )


def load_embedder(path=None):
    if path is None:
        path = _find_umda_file("models", "EmbeddingModel.pkl")
    embedder = load(path)
    upgrade_legacy_keyedvectors(embedder.model.wv)
    return embedder


def embed_smiles_list(smiles_list, embedder=None):
    """Returns (X, ok_mask) where X is (n_valid, 70) and ok_mask marks which
    input rows produced a usable (non-null, non-all-zero) embedding."""
    if embedder is None:
        embedder = load_embedder()
    vecs = []
    ok = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            if mol is None:
                vecs.append(None)
                ok.append(False)
                continue
            mol.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(mol)
            vec = embedder.vectorize(smi)
            if vec is None or np.all(vec == 0):
                vecs.append(None)
                ok.append(False)
            else:
                vecs.append(vec)
                ok.append(True)
        except Exception:
            vecs.append(None)
            ok.append(False)
    dim = next((v.shape[0] for v in vecs if v is not None), 70)
    X = np.zeros((len(vecs), dim))
    for i, v in enumerate(vecs):
        if v is not None:
            X[i] = v
    return X, np.array(ok)


def build_embedded_pool(candidate_csv=None, out_path=None):
    if candidate_csv is None:
        candidate_csv = os.path.join(EXTERNAL, "candidate_pool.csv")
    if out_path is None:
        out_path = os.path.join(EXTERNAL, "candidate_pool_embedded.npz")
    df = pd.read_csv(candidate_csv)
    embedder = load_embedder()
    X, ok = embed_smiles_list(df["smiles"].tolist(), embedder)
    df["embeddable"] = ok
    df.to_csv(candidate_csv, index=False)
    np.savez_compressed(out_path, X=X, ok=ok)
    print(f"Embedded {ok.sum()} / {len(df)} candidates ({ok.mean()*100:.1f}%).")
    print("Failures by label:")
    print(df.loc[~ok, "label"].value_counts())
    return df, X, ok


if __name__ == "__main__":
    build_embedded_pool()
