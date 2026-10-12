"""Table 3 (recommended TMC-1 targets), built from the reproducible final ranking.

Written for the release. The original builders of the frozen Table 3 (a chain of scripts, the last of
which is not available) are not shipped; this script rebuilds the same table from the released inputs:

  1. Fit the frozen model (7 standardized features, L2 logistic regression, C = 1.0, lbfgs) on all 1,434
     candidates, y = 1 for the 63 known positives (data/frozen/candidates.csv).
  2. Score the 1,371 not-yet-known candidates and rank them (score descending, candidate_id ascending).
  3. Drop the two structures that are not physically meaningful (original ranks 5 and 71: a triple bond
     forced into a three-membered ring, and five cumulated double bonds in a five-membered ring), keeping
     the original rank numbers, and keep the first 100 remaining rows.
  4. Attach the Reasons text (the two largest positive contributions to the logistic score) and the
     hand-curated annotation columns Species, Elsewhere, LabSpec, Note and ReasonsFlag (an audit flag
     prefixed to Reasons for four rows), which are read from
     data/frozen/table3_annotations.csv (keyed by LeeID). Those four columns come from catalog lookups and a
     literature cross-check and cannot be recomputed from the model.

Outputs: tables/table3_recommended_targets.csv and a journal-format (MRT) copy with a clean header.
The MRT data block is identical to the frozen submission file; only the header lines differ.
"""
import hashlib
import os
import sys
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import numpy as np
import pandas as pd
from astropy.table import Table
from sklearn.linear_model import LogisticRegression

from paths import DATA_FROZEN, TABLES
from evaluation_engine import standardize_features, FEATURE_COLS

NON_PHYSICAL_RANKS = [5, 71]   # original ranks removed (structures that cannot be embedded in 3D)
N_ROWS = 100

REASON_LABEL = {
    "dist": "high similarity to known TMC-1 species",
    "logcd": "high predicted column density",
    "sigma": "low predicted line-confusion/blend risk",
    "mu": "large dipole moment",
    "mw": "low molecular weight",
    "charged": "ionic species",
    "radical": "open-shell radical",
}

cand = pd.read_csv(os.path.join(DATA_FROZEN, "candidates.csv"))
assert len(cand) == 1434 and int(cand["positive"].sum()) == 63
Xs, mean, std = standardize_features(cand)
y = cand["positive"].to_numpy().astype(int)
model = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=1000)
model.fit(Xs, y)
coef = dict(zip(FEATURE_COLS, model.coef_[0]))


def reasons_for(row):
    contribs = {}
    for i, f in enumerate(FEATURE_COLS):
        val = float(row[f]) if f not in ("charged", "radical") else float(bool(row[f]))
        z = (val - mean[i]) / std[i]
        contribs[f] = coef[f] * z
    ranked = sorted(contribs.items(), key=lambda kv: -kv[1])
    top2 = [f for f, v in ranked[:2] if v > 0] or [ranked[0][0]]
    return "; ".join(REASON_LABEL[f] for f in top2)


idx_pool = np.where(~cand["positive"].to_numpy())[0]
assert len(idx_pool) == 1371
pool = cand.iloc[idx_pool].copy()
pool["score"] = model.predict_proba(Xs[idx_pool])[:, 1]
pool = pool.sort_values(["score", "candidate_id"], ascending=[False, True]).reset_index(drop=True)
pool["rank"] = np.arange(1, len(pool) + 1)

kept = pool[~pool["rank"].isin(NON_PHYSICAL_RANKS)].head(N_ROWS).copy()
assert len(kept) == N_ROWS

ann = pd.read_csv(os.path.join(DATA_FROZEN, "table3_annotations.csv"), dtype=str, keep_default_na=False)
assert ann["LeeID"].is_unique and len(ann) == N_ROWS


def klass(r):
    if bool(r["charged"]):
        return "ion"
    if bool(r["radical"]):
        return "radical"
    return "closed-shell"


out = pd.DataFrame({
    "Rank": kept["rank"].to_numpy(),
    "LeeID": kept["member_ids"].astype(str).to_numpy(),
    "Formula": kept["formulas"].to_numpy(),
    "SMILES": kept["smiles_members"].to_numpy(),
    "Class": kept.apply(klass, axis=1).to_numpy(),
    "Score": [round(float(v), 4) for v in kept["score"]],
    "logN": [round(float(v), 3) for v in kept["logcd"]],
    "sigma": [round(float(v), 3) for v in kept["sigma"]],
    "dist": [round(float(v), 3) for v in kept["dist"]],
    "mu": [round(float(v), 3) for v in kept["mu"]],
    "MW": [round(float(v), 3) for v in kept["mw"]],
    "Reasons": kept.apply(reasons_for, axis=1).to_numpy(),
})
out = out.merge(ann, on="LeeID", how="left", validate="one_to_one")
assert out["Species"].notna().all(), "annotation missing for some LeeID"
# audit flags (hand-curated) are prefixed to the computed reasons
out["Reasons"] = ["; ".join([f, r]) if f else r for f, r in zip(out["ReasonsFlag"], out["Reasons"])]
COLS = ["Rank", "Species", "LeeID", "Formula", "SMILES", "Class", "Score", "logN", "sigma", "dist", "mu",
        "MW", "Reasons", "Elsewhere", "LabSpec", "Note"]
out = out[COLS]

csv_path = os.path.join(TABLES, "table3_recommended_targets.csv")
out.to_csv(csv_path, index=False)


def sha256_of(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def ascii_only(s):
    s = str(s) if s is not None else s
    assert s.isascii(), f"non-ASCII content found: {s!r}"
    return s


outascii = out.copy()
for c in outascii.columns:
    outascii[c] = outascii[c].apply(ascii_only)

at = Table.from_pandas(outascii)
DESCR = {
    "Rank": "Rank in the frozen-model score; gaps at 5 and 71 (non-physical structures excluded, see table note)",
    "Species": "Plain conventional species name (formula notation; see Note for any identity caveat)",
    "LeeID": "Lee et al. (2021) Zenodo raw-entry id(s)",
    "Formula": "Molecular formula (ASCII plain formula, as used in the candidate table)",
    "SMILES": "SMILES string used for candidate matching (RDKit canonical, as parsed)",
    "Class": "ion / radical / closed-shell",
    "Score": "Frozen 7-feature logistic-regression score (fitted probability)",
    "logN": "Log10 predicted column density (standardized feature: logcd)",
    "sigma": "Predicted line-confusion/blend-risk proxy (standardized feature: sigma)",
    "dist": "Embedding distance to known TMC-1 species (standardized feature: dist)",
    "mu": "Dipole moment, Debye (standardized feature: mu)",
    "MW": "Molecular weight, amu",
    "Reasons": "Top model-score contributors",
    "Elsewhere": "'-' if never detected toward any astronomical source; else '<source(s)> (<first author> <year>)'",
    "LabSpec": "'CDMS <tag>', 'JPL <tag>', or 'none' -- resolved against the full downloaded catalog indexes; "
               "a column, never a filter",
    "Note": "Controlled-vocabulary identity/observability caveat, or blank",
}
for c in at.colnames:
    at[c].unit = "---"
    at[c].description = DESCR[c]

mrt_path = os.path.join(TABLES, "table3_recommended_targets.mrt")
at.write(mrt_path, format="ascii.mrt", overwrite=True)

HEADER = (
    "Active Learning for Molecular Discovery: A Retrospective Test toward TMC-1",
    "El-Kholy R.I., Hayman Z.M.",
    "Highest-ranked candidates for detection toward TMC-1",
)
lines = open(mrt_path).read().split("\n")
filled = 0
for i, line in enumerate(lines):
    if line == "Title:":
        lines[i] = "Title: " + HEADER[0]; filled += 1
    elif line == "Authors:":
        lines[i] = "Authors: " + HEADER[1]; filled += 1
    elif line == "Table:":
        lines[i] = "Table: " + HEADER[2]; filled += 1
assert filled == 3
with open(mrt_path, "w") as f:
    f.write("\n".join(lines))

for path in [csv_path, mrt_path]:
    assert open(path, "rb").read().decode("ascii")
    print(f"{path} SHA-256: {sha256_of(path)}")
