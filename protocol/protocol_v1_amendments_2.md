AMENDMENTS 2 to protocol_v1 (2026-09-21; recorded before any evaluation is run)
A6 Unit: candidates are list entries collapsed by species (same formula incl. charge and same
   connectivity/H-count skeleton; bond order, charge placement and stereo ignored). A test positive is a
   species matched by >= 1 detection dated inside the window; multiple detections or entries count once.
A7 Dates: these supersede v3.1: C3H+ and C5H+ = arXiv:2201.02434 (2022-01); H2CCCH+ = arXiv:2307.01733
   (2023-07). Month-level arXiv-ID dates suffice for annual origins. Sensitivity A4(a) is redefined: drop
   test positives whose date is not an arXiv v1 date. A4(b), the 60-day boundary shift, is unchanged.
A8 Outcomes, pre-specified: PASS = criterion met with the interval excluding 0; INCONCLUSIVE = point
   estimate favors the combined rule but the interval includes 0; FAIL = point estimate at or below the
   baseline. With 24 pooled test positives the design has about 93% power for V1 if the true AUCs are 0.90
   vs 0.75 and about 52% if 0.85 vs 0.75.
A9 Features: Lee's released values per entry as in lee_list_entry_labels_v3_1.csv; missing dipole values
   imputed with the pool median. Results are reported only after two independent implementations agree
   to three decimals.
