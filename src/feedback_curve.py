"""
Figure 3: full learner trajectory (not just "observations to target").
feedback_learner.run_v3_feedback() implements the adaptive policy (2-round cold start, then a
per-round LogisticRegression re-fit on revealed data, 10% random-swap exploration) but only
records the FIRST observation count at which the target is reached, then breaks.

This module repeats that per-round logic (same per-seed literal integer seed 0..n_seeds-1, same
batch size 5, same cold-start/explore-exploit rule, same 10% swap) but removes the early break and
instead records the cumulative positives-revealed count at every checkpoint (every batch of 5
observations) out to max_obs -- so the trajectory up to the original stopping point is identical to
run_v3_feedback()'s own (same rng draws consumed in the same order), and the curve continues beyond
it purely for plotting.
"""
import math
import numpy as np
from sklearn.linear_model import LogisticRegression


def run_v3_feedback_curve(Xs, y_pos, n_seeds=200, max_obs=800, batch_size=5):
    n = len(y_pos)
    y_pos = np.asarray(y_pos, dtype=bool)
    n_checkpoints = max_obs // batch_size
    curves = np.zeros((n_seeds, n_checkpoints + 1), dtype=int)  # [:,0] = n_obs=0

    for s in range(n_seeds):
        rng = np.random.default_rng(s)
        revealed = np.zeros(n, dtype=bool)
        n_revealed_pos = 0
        round_idx = 0
        max_rounds = math.ceil(n / batch_size) + 1
        curves[s, 0] = 0
        while round_idx < n_checkpoints and round_idx < max_rounds and revealed.sum() < n:
            unrevealed_idx = np.where(~revealed)[0]
            bs = min(batch_size, len(unrevealed_idx))
            if round_idx < 2 or n_revealed_pos == 0:
                picks = rng.choice(unrevealed_idx, size=bs, replace=False)
            else:
                X_rev = Xs[revealed]
                y_rev = y_pos[revealed].astype(int)
                if len(np.unique(y_rev)) < 2:
                    picks = rng.choice(unrevealed_idx, size=bs, replace=False)
                else:
                    model = LogisticRegression(C=0.3, class_weight="balanced",
                                                solver="lbfgs", max_iter=1000)
                    model.fit(X_rev, y_rev)
                    scores_unrev = model.predict_proba(Xs[unrevealed_idx])[:, 1]
                    order = np.argsort(-scores_unrev)
                    picks = unrevealed_idx[order[:bs]].copy()
                    for i in range(len(picks)):
                        if rng.random() < 0.1:
                            remaining = np.setdiff1d(unrevealed_idx, picks)
                            if len(remaining) > 0:
                                picks[i] = rng.choice(remaining)
            revealed[picks] = True
            n_revealed_pos += int(y_pos[picks].sum())
            round_idx += 1
            curves[s, round_idx] = n_revealed_pos
        # if the loop ended early (pool exhausted before n_checkpoints), hold the last value
        if round_idx < n_checkpoints:
            curves[s, round_idx + 1:] = n_revealed_pos

    n_obs_grid = np.arange(0, max_obs + 1, batch_size)
    return n_obs_grid, curves  # curves shape (n_seeds, n_checkpoints+1)


def obs_to_target_from_curve(n_obs_grid, curves, target):
    """First n_obs at which each seed's curve reaches >=target; matches
    run_v3_feedback()'s own obs_to_50pct_per_seed definition exactly."""
    out = []
    for row in curves:
        idx = np.argmax(row >= target) if (row >= target).any() else len(row) - 1
        out.append(int(n_obs_grid[idx]))
    return np.array(out)
