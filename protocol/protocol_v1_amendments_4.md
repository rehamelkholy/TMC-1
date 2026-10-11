AMENDMENTS 4 to protocol_v1 (2026-09-21; recorded after the first joint comparison; definitional only)
J1 Species key (resolves I1): key = CalcMolFormula(mol) + '|' + InChI formula layer + '|' + InChI /c layer + '|'
   + InChI /h layer, from RDKit MolToInchi (default options), all other layers dropped. If MolToInchi returns
   empty or the SMILES does not parse, the entry is its own singleton candidate. The repo implementation is
   authoritative; the chat-side implementation is conformed to it.
J2 A4(b): "within 60 days of an origin" includes the first origin (2022-01-01). A test positive dated on or
   after an origin and within 60 days moves to origin minus 1 day; one dated before an origin and within 60
   days moves to origin plus 1 day. Moves apply simultaneously to positives with baseline test-window dates
   (>= 2022-01-01).
J3 A4(a): the repo's exact first_public_date_basis is authoritative; the chat-side string rule is retired.
J4 The first joint comparison found identical origin tables, best baseline and top-100 counts, and outcomes
   V1 PASS, V2 INCONCLUSIVE, V3 PASS under both implementations, but per-origin AUCs differing by up to
   0.006 because of J1. Both implementations are re-run under J1-J3; the re-run values are final.
