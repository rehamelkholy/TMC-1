"""
Baselines.

(a) Static nearest-neighbor-in-embedding-space ranking -- Lee et al. 2021's own method: no
    uncertainty, no retraining, one single-shot ranking pass. For each unlabeled candidate,
    score = -(distance to nearest labeled/training point in the 70-dim PCA embedding space), i.e.
    candidates most SIMILAR to already-known molecules rank highest (the "expand the known chemical
    space outward" logic of the original paper, as opposed to the acquisition loop's "go where the
    model is least certain").
(b) Pure random selection (many draws).
(c) A "most abundant chemical family" heuristic is NOT implemented: no such published heuristic for
    TMC-1 candidate prioritization was found in the literature check, so this baseline is left out
    rather than invented.
"""
import numpy as np
from scipy.spatial import cKDTree


def static_nn_rank(X_train, X_candidates):
    """One single-shot ranking pass: score = -distance to the nearest point
    already in the (pre-cutoff, never-updated) training set. Higher score =
    more similar to something already known = ranked earlier, matching Lee
    et al. 2021's nearest-neighbor recommendation logic exactly (no
    uncertainty, no iteration)."""
    tree = cKDTree(X_train)
    dist, _ = tree.query(X_candidates, k=1)
    return -dist


def random_rank(n_candidates, n_draws=5, random_state=0):
    """Returns n_draws independent random permutations (as rank arrays) over the candidate pool."""
    rng = np.random.RandomState(random_state)
    return [rng.permutation(n_candidates) for _ in range(n_draws)]
