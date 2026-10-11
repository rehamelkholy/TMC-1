"""
Acquisition functions based on regressor uncertainty: GPR posterior variance, Bayesian-ridge
credible-interval width, and two ensemble-based alternatives for the non-Bayesian models:

  - `model_zoo_variance`: cheap, no refitting -- variance of predictions across the 6
    non-Bayesian point-estimate models in the zoo (ridge, linear, SVR, kNN, RF, GBR) at each
    candidate point.
  - `bootstrap_variance`: refits the same 6-model zoo on bootstrap resamples of the current
    labeled set each round and takes the prediction variance across (model x resample); about 30x
    more expensive per round at the default n_bootstrap=5. Variance across resamples is computed
    separately per model type and only then averaged across model types, so the score isolates
    within-model resampling (epistemic) uncertainty rather than between-model-type disagreement.
    A version that pooled all len(zoo_models) * n_bootstrap raw predictions into one std was
    dominated by between-model-type disagreement, i.e. nearly the signal `model_zoo_variance`
    already gives (round-0 rank correlation 0.9999996 on the real candidate pool).

The acquisition-loop driver (loop.py) takes an `acquisition` argument so any of these can be selected.
"""
import numpy as np
from sklearn.base import clone


def gpr_uncertainty(gpr_model, X):
    """GPR posterior standard deviation at each candidate point."""
    _, std = gpr_model.predict(X, return_std=True)
    return std


def bayesian_ridge_uncertainty(br_model, X):
    """BayesianRidge's own predictive std (credible-interval width proxy)."""
    _, std = br_model.predict(X, return_std=True)
    return std


def model_zoo_variance(zoo_models, X):
    """Cheap ensemble-based uncertainty proxy: std across the predictions of
    every non-Bayesian point-estimate model already in the zoo (no extra
    fitting cost -- just predict() with whatever is already fitted this
    round)."""
    preds = np.stack([m.predict(X) for m in zoo_models.values()], axis=0)
    return preds.std(axis=0)


def bootstrap_variance(zoo_models, X_train, y_train, X_query, n_bootstrap=5, random_state=0):
    """Ensemble-based alternative to model_zoo_variance.

    Refits EVERY non-Bayesian model in the zoo (linear, ridge, svr, knn, rfr, gbr) on `n_bootstrap`
    bootstrap resamples of the current labeled set each round. For each candidate point it isolates
    genuine resampling/epistemic uncertainty -- how much ONE model's prediction moves as the training
    data is resampled -- rather than between-model-type disagreement, which model_zoo_variance
    already measures:

      1. For each model type, compute the variance of its `n_bootstrap` resampled predictions at
         each query point (within-model variance).
      2. Average that per-model variance across the `len(zoo_models)` model types (so no single
         model type's scale dominates).
      3. Return the sqrt of that mean variance, so the result is on the same "std of predicted
         log10 column density" scale as model_zoo_variance and the other acquisition scores.

    Cost: len(zoo_models) * n_bootstrap refits per round (6*5=30 by default). Bootstrapping the full
    non-Bayesian zoo (not a single model type) keeps every model exercised; if runtime on a large
    candidate pool is a problem, bootstrapping a single representative model (GBR is a reasonable
    pick: unlike RF it has no internal bootstrap aggregation in its point predictions) and/or
    lowering n_bootstrap is the cheaper alternative.
    """
    rng = np.random.RandomState(random_state)
    n = len(y_train)
    per_model_variances = []
    for base_estimator in zoo_models.values():
        preds = []
        for b in range(n_bootstrap):
            idx = rng.randint(0, n, size=n)
            est = clone(base_estimator)
            est.fit(X_train[idx], y_train[idx])
            preds.append(est.predict(X_query))
        preds = np.stack(preds, axis=0)  # (n_bootstrap, n_query)
        per_model_variances.append(preds.var(axis=0))  # within-model variance across resamples
    per_model_variances = np.stack(per_model_variances, axis=0)  # (n_models, n_query)
    mean_variance = per_model_variances.mean(axis=0)  # averaged across model types
    return np.sqrt(mean_variance)


ACQUISITION_FUNCTIONS = {
    "gpr_variance": gpr_uncertainty,
    "bayesian_ridge_variance": bayesian_ridge_uncertainty,
    "model_zoo_variance": model_zoo_variance,
    "bootstrap_variance": bootstrap_variance,
}
