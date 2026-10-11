AMENDMENTS 7 to protocol_v1 (recorded AFTER the F7 results were unsealed; post hoc)
M1 Defect in Amendments 5 (K3-K5): species known at origin o are themselves members of K3(o) and of the GPR
   training set, so every training positive had dist = 0 and posterior sigma near its minimum (median 0.54),
   while test positives did not (median dist 2.7-3.1). The ranking model at each origin was trained on
   features no test candidate can have. F7 results are reported as pre-registered, with this defect stated.
M2 Leave-self-out variant: for each candidate known at origin o, dist, logcd and sigma are recomputed with that
   species removed from K3(o) (all of its distance-anchor points) and from the GPR training set (its 9
   observations), refitting the GPR by maximum likelihood with the same settings for each such candidate.
   Features of candidates not known at o are unchanged from F7. Everything else as in Amendments 5-6.
M3 Status: post hoc, labelled as such wherever reported; outcomes classified per A8 but never used to replace a
   pre-registered verdict.
