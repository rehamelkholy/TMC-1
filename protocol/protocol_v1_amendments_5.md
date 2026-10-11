AMENDMENTS 5 to protocol_v1 (recorded before any run): the recomputed-feature sensitivity named in protocol_v1.
K1 Scope: temporal evaluation only (V1, V2 and all per-origin metrics). V3 is not re-run (its features are the
   2021 state). Primary verdicts stand as reported; this is a sensitivity. Outcome categories per A8.
K2 Universe, species key, labels, origins, seeds, logistic model, baselines, bootstrap: identical to the F6 run
   (purged 1,489-entry universe, J1 key, 1,440 candidates, 40 positive candidates).
K3 Known chemistry at origin o = Lee Table A2 anchors (authoritative 87) + every detection in
   post2021_detections_v3_2.csv with first_public_date < o (in-list AND out-of-list) + pre-2021 species known
   under A3 (e.g. C3N-, C5N-).
K4 GPR training set at origin o = the 643 observations stored in the distributed regressors.pkl (inverse-transformed
   to raw PCA-70 coordinates) + for each species in K3 not among the 87 anchors, 9 observations: its raw PCA-70
   embedding with log10 N + Gaussian noise (sigma = 0.5), seed 20260921 + origin index. Species without a recorded
   column density are used as distance anchors only, never as GPR training points; list them.
   GPR: the regressors.pkl pipeline structure (StandardScaler + GaussianProcessRegressor with the rational-quadratic,
   dot-product and white-noise kernel), hyperparameters by maximum likelihood, same optimizer settings as the
   within-list reproduction run.
K5 Recomputed per entry at each origin: dist = Euclidean distance in raw PCA-70 space to the nearest species in K3;
   logcd = GPR posterior mean; sigma = GPR posterior standard deviation. mu, mw, charged, radical stay as released.
   Entries are then collapsed by the J1 key exactly as in I1 (numeric features averaged).
K6 Gates, before any evaluation: G1 with K3 restricted to the 87 anchors, the fitted kernel and log-marginal
   likelihood equal the within-list reproduction run (LML -632.15) and per-entry features match that run to
   1e-6. G2 in that state, Spearman rho vs Lee's released values >= 0.95 for sigma and for dist.
   If either gate fails, stop; no evaluation.
