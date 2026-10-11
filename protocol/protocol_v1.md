PROTOCOL v1 (recorded 2026-09-21, before any evaluation is run; amendments go in a separate file)
Question: within a similarity-generated candidate list, do active-learning rankings find later-detected
TMC-1 molecules faster than the rankings observers would use anyway?
Universe: Lee et al. (2021) released list (1,510 entries; already-known anchors removed).
Labels: TMC-1 detections (validation set v3) matched to list entries by rule: A identical canonical SMILES;
B identical after removing stereo, only if the Lee entry carries no stereo marks; C same formula (incl.
charge) and same connectivity/H-count skeleton, only if the Lee entry carries no stereo marks.
Time: each detection has a first-public-date (earliest arXiv posting, else journal online date).
Origins: annual boundaries 2022-01-01, 2023-01-01, 2024-01-01, 2025-01-01, chosen from counts only; drop the
first origin if fewer than 12 in-list detections are known then. Test window k = [origin k, origin k+1);
last window ends 2026-09-21. At an origin, "known" = anchors plus detections dated before it. Test pool =
list entries not yet known; test positives = in-list detections dated inside the window.
Features (primary): Lee's released values as of 2021-07 (nearest-anchor distance, GPR posterior sigma,
GPR predicted log column density, dipole magnitude, MW, charged, radical). Sensitivity: recomputed with our
pipeline using the anchors known at each origin.
Ranking model: standardized logistic regression (L2, C=1.0), the 7 features above, fit on labels known at
the origin. No tuning.
Feedback learner: cold start; batch 5; first 10 random; logistic (C=0.3, class-balanced) retrained after each
batch on outcomes so far; 10% random exploration; 200 seeds.
Baselines: random; high predicted log column density; detectability heuristic (predicted log column
density + 2*log10 dipole); nearest-anchor distance low. Best baseline = highest pooled top-100 count.
Additional non-gating comparators: highest-sigma (uncertainty sampling), upper confidence bound.
Metrics: detections found within budgets 25/50/100/200; molecules searched per detection; budget to reach
50% of eventual in-list detections; AUC-ROC.
Statistics: paired bootstrap, 2,000 resamples, 95% intervals; counts pooled over origins.
Pass/fail: V1 combined-rule AUC beats the best baseline, interval on the difference excludes 0. V2 at budget
100 it finds >=20% more detections than the best baseline, interval excludes 0. V3 the feedback learner
reaches 50% of eventual detections in fewer observations than the best fixed ordering.
If V1 and V2 do not clear, report exactly that and adjust the claim; no post-hoc feature or threshold changes.
Usability deliverables: U1 timestamped, DOI-archived ranked list of still-undetected entries with plain-
language reasons and a literature check; U2 open-source tool reproducing the paper's numbers; U3 budget-vs-
yield figure; U4 "would we have caught it?" rank table.
