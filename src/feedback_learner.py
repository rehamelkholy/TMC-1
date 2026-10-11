"""
Feedback learner V3 (protocol item I7): cold-start active learning over the whole
candidate pool, plus the fixed-ordering comparator it is benchmarked against.
"""
import math
import numpy as np
from sklearn.linear_model import LogisticRegression

N_SEEDS = 200
N_TIEBREAK_TRIALS = 200
BASE_SEED = 20260921


def obs_to_target_avg(score, y_pos, target, n_trials, seed):
    """Fixed-ordering comparator: expected number of observations (top-down
    in a ranking by `score`) needed to reveal >= target positives, averaged
    over n_trials random tie-break realizations (I4's tie-break procedure,
    reused here per I7's 'ties averaged as in I4')."""
    n = len(y_pos)
    y_pos = np.asarray(y_pos, dtype=bool)
    rng = np.random.default_rng(seed)
    totals = []
    for t in range(n_trials):
        jitter = rng.random(n)
        order = np.lexsort((jitter, -np.asarray(score)))
        y_ord = y_pos[order]
        cum = np.cumsum(y_ord)
        # first index (1-based count) where cumulative positives reaches target
        pos_needed = np.argmax(cum >= target) + 1 if (cum >= target).any() else n
        totals.append(pos_needed)
    return float(np.mean(totals)), totals


def run_v3_feedback(Xs, y_pos, n_seeds=N_SEEDS):
    n = len(y_pos)
    n_positives_total = int(y_pos.sum())
    target = math.ceil(0.5 * n_positives_total)

    obs_results = []
    for s in range(n_seeds):
        rng = np.random.default_rng(s)
        revealed = np.zeros(n, dtype=bool)
        n_revealed_pos = 0
        round_idx = 0
        obs_to_50 = None
        max_rounds = math.ceil(n / 5) + 1
        while round_idx < max_rounds and revealed.sum() < n:
            unrevealed_idx = np.where(~revealed)[0]
            batch_size = min(5, len(unrevealed_idx))
            if round_idx < 2 or n_revealed_pos == 0:
                picks = rng.choice(unrevealed_idx, size=batch_size, replace=False)
            else:
                X_rev = Xs[revealed]
                y_rev = y_pos[revealed].astype(int)
                if len(np.unique(y_rev)) < 2:
                    picks = rng.choice(unrevealed_idx, size=batch_size, replace=False)
                else:
                    model = LogisticRegression(C=0.3, class_weight="balanced",
                                                solver="lbfgs", max_iter=1000)
                    model.fit(X_rev, y_rev)
                    scores_unrev = model.predict_proba(Xs[unrevealed_idx])[:, 1]
                    order = np.argsort(-scores_unrev)
                    picks = unrevealed_idx[order[:batch_size]].copy()
                    for i in range(len(picks)):
                        if rng.random() < 0.1:
                            remaining = np.setdiff1d(unrevealed_idx, picks)
                            if len(remaining) > 0:
                                picks[i] = rng.choice(remaining)
            revealed[picks] = True
            n_revealed_pos += int(y_pos[picks].sum())
            round_idx += 1
            if obs_to_50 is None and n_revealed_pos >= target:
                obs_to_50 = int(revealed.sum())
                break
        if obs_to_50 is None:
            obs_to_50 = int(revealed.sum())  # safety: exhausted pool without reaching target (shouldn't happen)
        obs_results.append(obs_to_50)

    return {
        "n_positives_total": n_positives_total,
        "target_positives": target,
        "obs_to_50pct_per_seed": obs_results,
        "median": float(np.median(obs_results)),
        "mean": float(np.mean(obs_results)),
    }
