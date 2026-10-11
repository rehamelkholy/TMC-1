AMENDMENTS 3 to protocol_v1: implementation details (2026-09-21; recorded before any run)
I1 Data: lee_list_entry_labels_v3_1.csv (1,501 entries) plus the A7 date corrections. Positive entries =
   in_list_by_v3_1 True. Species key = (RDKit formula incl. charge, connectivity/H-count skeleton with bond
   orders, charges and stereo removed). Entries with equal key collapse into one candidate: numeric features =
   mean over members; date = earliest first_public_date among members' matched detections; charged = any member
   charge != 0; radical = any member mult > 1. Missing dipole imputed with the pool median before any use.
I2 Dates: month-only dates set to the 1st of the month; window = [origin, next origin); last window ends
   2026-09-21.
I3 Origins 2022-01-01, 2023-01-01, 2024-01-01, 2025-01-01. At origin o: known = candidates dated < o
   (including pre-2021 ones); training labels y=1 for known, 0 for all others; features standardized on all
   candidates; sklearn LogisticRegression L2, C=1.0, lbfgs, max_iter=1000, no class weights. Test pool =
   candidates not known at o; test positives = candidates dated in the window; candidates dated after the
   window stay in the pool as negatives.
I4 Baseline scores: high predicted log CD; detectability = logCD + 2*log10(dipole), dipole 0 -> -99; low
   nearest-anchor distance; random. Non-gating: high sigma; UCB = logCD + sigma. Ties broken randomly (200
   tie-breaks averaged, seed 20260921); AUC uses average ranks.
I5 Metrics: per origin AUC-ROC and detections in the top 25/50/100/200 of the test pool (fixed counts).
   Pooled counts = sums over origins; pooled AUC = positive-weighted mean of per-origin AUCs. Best baseline =
   the gating baseline (not random) with the highest pooled top-100 count.
I6 Bootstrap: 2,000 paired resamples, seed 20260921: within each origin resample positives and negatives
   separately with replacement, identical indices for both methods; recompute pooled AUC and top-100 counts;
   95% percentile intervals of (combined minus best baseline).
I7 Feedback learner (V3): pool = all candidates; positives = all candidates matched to an in-list detection
   (any date). Cold start; rounds of 5; first 2 rounds uniformly random; afterward logistic (C=0.3,
   class_weight balanced, standardized features) fit on all outcomes so far; each of the 5 slots is replaced
   by a random unrevealed candidate with probability 0.1; if no positive has been revealed, pick randomly.
   200 seeds (0..199). Observations to 50% = smallest number revealed with positives >= ceil(0.5*n_pos).
   Best fixed ordering = the gating baseline (not random) with the fewest observations to 50% (ties averaged
   as in I4). V3 PASS if the learner's median over seeds is lower; also report the share of seeds that beat it.
I8 Sensitivities: A4(a) approximate in chat-side runs: drop from the test pool test positives whose date string
   does not indicate an arXiv-derived date ('journal', 'not located', 'year only'); the repo run uses the exact
   basis column. A4(b): move each test positive within 60 days of an origin to the other side of it.
I9 Agreement (replaces "three decimals" in A9): deterministic quantities (origin table, logistic fits, AUCs,
   top-k counts) must match to 3 decimals; stochastic quantities (bootstrap endpoints, learner medians, random
   baseline) within +/-0.01 AUC, +/-1 detection, and +/-10% respectively. Each implementation writes results.json
   and logs its SHA-256 BEFORE anything is displayed; results are then revealed and compared together.
