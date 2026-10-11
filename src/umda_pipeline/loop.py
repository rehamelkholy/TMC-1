"""
The acquisition-loop driver.

Loop structure: train on the pre-cutoff inventory -> rank ALL candidate molecules (not just true
detections) by acquisition score -> reveal the true column density for any of the top-k-per-round
candidates that are real post-2021 detections (else leave as still-unlabeled) -> retrain -> repeat.

Evaluation: at every round, before any revealing happens, record the FULL acquisition-ranking
position (0 = top) of every true detection that has not been revealed yet. This gives a
rank-vs-round trajectory for each of the 64 real detections under each method, from which
"recommendations needed to recover N% of detections" curves can be built at any top-N threshold in
post-processing, not just the threshold matching the per-round reveal budget `k_per_round`.

k_per_round=5 is the headline value; it stands in for a plausible per-cycle follow-up-observation
budget, and the top-N "counted as found" cutoff used at reporting time remains a free
post-processing choice, since full rank trajectories are saved every round.

Epsilon-greedy exploration: a deterministic method (GPR above all) that finds no new true detection
in a round leaves its labeled set, and therefore its own predictions, unchanged, so a pure
top-k-by-acquisition-score policy can freeze onto one static ranking indefinitely. Each round, only
`k_per_round * (1 - exploration_frac)` slots go to the top-ranked-by-acquisition-score candidates
(exploitation); the remaining `k_per_round * exploration_frac` slots are drawn uniformly at random
from the unrevealed pool (exploration), regardless of acquisition score. The default
`exploration_frac=0.2` (1 of 5 slots at k_per_round=5) is a standard epsilon-greedy value. It applies
only to this loop (not to `static_nn_baseline` or `random_rank` in baselines.py) and identically to
all three acquisition methods compared.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER1_ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, PAPER1_ROOT)

from paths import EXTERNAL  # noqa: E402

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, Ridge, BayesianRidge
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RationalQuadratic, DotProduct, WhiteKernel
from sklearn.preprocessing import StandardScaler
from joblib import load

from acquisition import (
    gpr_uncertainty, bayesian_ridge_uncertainty, model_zoo_variance, bootstrap_variance,
)
from baselines import static_nn_rank, random_rank


def _find_umda_file(*parts):
    """Locate a file inside the umda checkout. The clone is expected at external/umda
    (see data/raw/FETCH.md); duplicated in embeddings.py."""
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


def fresh_model_zoo(random_state=1215677):
    """Same model types/hyperparameter families as the umda checkout's models/regressors.pkl, as
    fresh, unfitted instances -- refit from scratch each round on the current labeled set."""
    kernel = RationalQuadratic(alpha=1, length_scale=1) + DotProduct(sigma_0=1) + WhiteKernel(noise_level=1)
    return {
        "linear_regression": LinearRegression(),
        "ridge": Ridge(solver="svd"),  # sklearn 1.0.2's cholesky path uses scipy's removed `sym_pos` kwarg with scipy>=1.11 -- svd sidesteps it
        "br": BayesianRidge(),
        "svr": SVR(),
        "knn": KNeighborsRegressor(n_neighbors=5),
        "rfr": RandomForestRegressor(random_state=random_state),
        "gbr": GradientBoostingRegressor(random_state=random_state),
        "gpr": GaussianProcessRegressor(kernel=kernel, random_state=random_state),
    }


def load_precutoff_training_set(regressors_path=None):
    """Pull the pre-cutoff training set directly out of the fitted GPR in the frozen baseline model: a
    fitted sklearn GaussianProcessRegressor retains its own X_train_/y_train_.

    `regs["gpr"]` is a Pipeline [('scaler', StandardScaler()), ('regressor', GaussianProcessRegressor)],
    so the fitted GPR's own `.X_train_` is the SCALER'S OUTPUT (z-scored, per-dimension mean~0/std~1),
    NOT the raw mol2vec->PCA(70) space every candidate/detection embedding lives in (embeddings.py never
    applies this scaler). This function inverse-transforms back to the same raw space before returning
    (X_train_raw = X_train_ * scaler.scale_ + scaler.mean_); 63/68 distinct rows then match the
    published Table A2 embeddings to floating-point precision."""
    if regressors_path is None:
        regressors_path = _find_umda_file("models", "regressors.pkl")
    regs = load(regressors_path)
    scaler = regs["gpr"].steps[0][1]
    gpr = regs["gpr"].steps[-1][1]
    X_train_raw = gpr.X_train_ * scaler.scale_ + scaler.mean_
    return X_train_raw, gpr.y_train_.copy()


def run_acquisition_loop(X_train0, y_train0, X_pool, pool_df, acquisition="model_zoo_variance",
                          k_per_round=5, max_rounds=40, exploration_frac=0.2, random_state=0,
                          zoo_random_state=1215677, verbose=True):
    """
    X_train0, y_train0: starting (pre-cutoff) labeled set.
    X_pool: (n_candidates, 70) embeddings for the full candidate pool.
    pool_df: DataFrame aligned with X_pool rows; needs "label"
             ("true_detection"/"background_undetected") and, for true
             detections, "log10_column_density_cm2" and "formula".
    exploration_frac: fraction of each round's k_per_round budget reserved
              for uniform-random exploration rather than top-acquisition-
              score exploitation (see module docstring). 0.0 recovers pure
              exploitation.
    random_state: seeds the exploration draws only (independent of
             `bootstrap_variance`'s own per-round `random_state=rnd`, so the
             two sources of randomness don't get confounded).
    zoo_random_state: seeds the model zoo's own internal randomness (RF,
              GBR, GPR -- see fresh_model_zoo). Exposed so that a run's results
              (e.g. an early "burst" of true-detection hits) can be checked
              across different model instantiations.

    Returns: rank_log (DataFrame: formula, round, rank, acquisition_score,
             revealed_this_round) and revealed_order (list of (round, rank,
             formula) for every true detection actually revealed).
    """
    n = len(pool_df)
    is_true_det = (pool_df["label"] == "true_detection").to_numpy()
    true_formulas = pool_df.loc[is_true_det, "formula"].to_numpy()
    true_cd = pool_df.loc[is_true_det, "log10_column_density_cm2"].to_numpy()
    formula_to_cd = dict(zip(true_formulas, true_cd))
    formula_to_poolidx = {f: i for i, f in zip(np.where(is_true_det)[0], true_formulas)}

    revealed = np.zeros(n, dtype=bool)  # candidates whose label has been folded into training
    X_labeled = X_train0.copy()
    y_labeled = y_train0.copy()

    rank_rows = []
    revealed_order = []
    remaining_true = set(true_formulas)

    explore_rng = np.random.RandomState(random_state)
    n_explore = int(round(k_per_round * exploration_frac))
    n_exploit = k_per_round - n_explore

    for rnd in range(max_rounds):
        if not remaining_true:
            break
        zoo = fresh_model_zoo(random_state=zoo_random_state)
        for m in zoo.values():
            m.fit(X_labeled, y_labeled)

        unrevealed_idx = np.where(~revealed)[0]
        X_query = X_pool[unrevealed_idx]

        non_bayesian = {k: v for k, v in zoo.items() if k not in ("gpr", "br")}
        if acquisition == "gpr_variance":
            scores = gpr_uncertainty(zoo["gpr"], X_query)
        elif acquisition == "bayesian_ridge_variance":
            scores = bayesian_ridge_uncertainty(zoo["br"], X_query)
        elif acquisition == "model_zoo_variance":
            scores = model_zoo_variance(non_bayesian, X_query)
        elif acquisition == "bootstrap_variance":
            # Ensemble-based non-Bayesian alternative, run side by side with
            # model_zoo_variance -- see acquisition.py's docstring for the
            # cost/design tradeoff. random_state=rnd so each round's bootstrap
            # resample is reproducible but not identical round-to-round.
            scores = bootstrap_variance(non_bayesian, X_labeled, y_labeled, X_query,
                                         n_bootstrap=5, random_state=rnd)
        else:
            raise ValueError(
                f"Unknown acquisition function {acquisition!r} -- expected one of "
                "'gpr_variance', 'bayesian_ridge_variance', 'model_zoo_variance', "
                "'bootstrap_variance'."
            )

        order = np.argsort(-scores)  # descending: highest uncertainty first
        ranked_pool_idx = unrevealed_idx[order]

        # log the rank of every still-unrevealed true detection this round
        rank_of = {pool_idx: r for r, pool_idx in enumerate(ranked_pool_idx)}
        for f in list(remaining_true):
            pidx = formula_to_poolidx[f]
            rank_rows.append({
                "round": rnd, "formula": f, "rank": rank_of[pidx],
                "acquisition_score": float(scores[np.where(unrevealed_idx == pidx)[0][0]]),
            })

        # reveal top-k: n_exploit slots go to the top of the acquisition
        # ranking, n_explore slots are drawn uniformly at random from
        # whatever is left (see module docstring). Exploration keeps the
        # labeled set moving even when a method's ranking has gone stale
        # (e.g. GPR after it stops finding new true detections).
        exploit_idx = ranked_pool_idx[:n_exploit]
        if n_explore > 0:
            explore_pool = np.setdiff1d(unrevealed_idx, exploit_idx, assume_unique=True)
            n_draw = min(n_explore, len(explore_pool))
            if n_draw > 0:
                explore_idx = explore_rng.choice(explore_pool, size=n_draw, replace=False)
            else:
                explore_idx = np.array([], dtype=int)
        else:
            explore_idx = np.array([], dtype=int)
        top_k_pool_idx = np.concatenate([exploit_idx, explore_idx]).astype(int)
        newly_revealed_true = []
        for pidx in top_k_pool_idx:
            if is_true_det[pidx] and not revealed[pidx]:
                f = pool_df.iloc[pidx]["formula"]
                if f in remaining_true:
                    newly_revealed_true.append((pidx, f))
            revealed[pidx] = True  # queried this round either way (background: no label to add)

        for pidx, f in newly_revealed_true:
            cd = formula_to_cd[f]
            X_labeled = np.vstack([X_labeled, X_pool[pidx:pidx+1]])
            y_labeled = np.append(y_labeled, cd)
            rnk = rank_of[pidx]
            revealed_order.append({"round": rnd, "rank": rnk, "formula": f})
            remaining_true.discard(f)

        if verbose:
            print(f"round {rnd:2d}: revealed {len(newly_revealed_true)} true detection(s), "
                  f"{len(remaining_true)} remaining, labeled set size={len(y_labeled)}")

    rank_log_df = pd.DataFrame(rank_rows, columns=["round", "formula", "rank", "acquisition_score"])
    revealed_df = pd.DataFrame(revealed_order, columns=["round", "rank", "formula"])
    return rank_log_df, revealed_df
