AMENDMENTS 9 (recorded after final results were sealed; definitional; no verdict depends on it). R1 The random
ordering's pooled AUC is reported as its exact expectation, 0.5, and its top-k count as the exact expectation
sum over origins of k x (test positives / test-pool size). R2 Top-k counts of methods with tied scores are computed
as exact expectations under uniform random tie-breaking: candidates strictly above the tie are counted fully, and the
tied block at the cutoff contributes (tied positives) x (remaining slots) / (tied candidates). This replaces the
200-draw average of I4. R3 All remaining random-number seeds are derived from a stable function of their inputs,
not from Python's hash().
