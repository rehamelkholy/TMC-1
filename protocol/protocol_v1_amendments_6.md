AMENDMENTS 6 to protocol_v1 (recorded after gate G1 failed with the F6 anchor SMILES, max dist diff 7.99, and
before any evaluation): embedding representations for the recomputed-feature sensitivity (Amendments 5).
L1 Anchors: each of the 87 anchors is embedded from its SMILES in table_a2_lee2021.py, the representation that
   reproduces Lee's training embeddings. table_a2_authoritative.csv remains the identity reference; the two
   tables must map one-to-one by InChIKey.
L2 Detections (K3): a detection matched to Lee-list entries by rules A/B/C is embedded from the SMILES of its
   matched entries: all matched entries serve as distance anchors, and the entry with the lowest Lee id supplies
   the GPR training embedding. An unmatched detection is embedded from its SMILES in post2021_detections_v3_2.csv.
L3 Gates K6 and everything else in Amendments 5 are unchanged.
