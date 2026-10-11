# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
candidates_purged_v5.csv (1439 rows, 47 positive)
-> candidates_purged_v6.csv (1434 rows, 57 positive), per the approved
item 2/3 disposition (checkpoint-approved additions A/B), PLUS one further
correction (below) surfaced during item 5's new-top-100-entrant audit and
confirmed before incorporating, per the user's explicit "report and stop
if any differ" instruction:

REMOVE (5): anchor-duplicate candidates, connectivity/charge-placement
variants of 5 of the 87 Lee anchors (rule C, item 2):
  293 (H2COH+/CH3O+), 376 (HNC3/C3HN), 396 (HOCO+/CHO2+),
  490 (HC3NH+/C3H2N+), 523 (NCCNH+/C2HN2+)

RELABEL positive=True (10): TMC-1 inventory matches (rule C, item 2 + the
item-B 2013-2021 ADS screening), all dated before every origin
(2022-01-01), so known at every origin, never a test positive. References
below are POST-CHECKPOINT CORRECTED per item A's "verified reference"
requirement -- two of the review checkpoint's original citations (HCO+, CO) did
not survive independent re-verification and were replaced (see notes).
A tenth candidate (493, HCO) was added after the checkpoint: it surfaced
at rank 86 of item 5's v6 U1 ranking (a top-100 entrant requiring a
TMC-1-detection check per item 5's own instruction), was absent from both
Lee's 87 anchors (confirmed: no HCO entry in table_a2_authoritative_v2.csv,
only HCOOH) and the TMC-1 inventory of 40_match_rule_c.py (a genuine
omission -- neutral HCO was never added to its INVENTORY dict), and was confirmed as a firm, non-tentative,
already-detected TMC-1 species via the 2026 TMC-1 review (arXiv:2604.26439,
Table 1: N(HCO)=1.1e12 cm^-2, ref 43 = Cernicharo, Marcelino, Agundez, Endo,
Cabezas et al. 2020, A&A 642, L17). An initial research pass had misattributed
this to a non-existent "Agundez et al. 2021" paper; re-verified directly
against the review's own reference list and the actual A&A page for the
2020 paper (DOI 10.1051/0004-6361/202039351) before incorporating.
  383 HCO+    -- Pratap, Dickens, Snell, Miralles, Bergin, Irvine & Schloerb
                1997, ApJ 486, 862, "A Study of the Physics and Chemistry
                of TMC-1". CORRECTION: the review checkpoint's original citation
                (Langer, Wilson, Henry & Guelin 1978, ApJ 225, L139) was
                re-checked against its own abstract and found to report
                ONLY H13CO+ and HC18O+ toward TMC-1, not the main HCO+
                line -- disqualified under rule A1 (non-isotopologue).
                Pratap et al. 1997 is this project's own pre-existing
                standing TMC-1/HCO+ citation (already used in the v5 audit's
                table3_v5_audit_flags.csv rank-31 entry, and independently
                the field-standard citation per a 2026 TMC-1 review,
                arXiv:2604.26439). "Earliest TMC-1 report located" per
                item A -- no earlier verifiable firm main-HCO+/TMC-1 report
                found despite a dedicated research pass (candidates
                checked and ruled out: Buhl & Snyder 1970 [W3(OH)/Orion/
                SgrA/W51/L134, not TMC-1]; Guelin/Langer/Snell/Wootten 1977
                DCO+ discovery [L63/L134/L134N, not TMC-1]; Guelin/Langer/
                Wilson 1982 [D13CO+ in TMC-1, an isotopologue]; Loren &
                Wootten 1980 HCO+ self-reversals [source list unconfirmed]).
  298 H2S     -- Minh, Irvine & Ziurys 1989, ApJ 345, L63
  272 NS      -- McGonagle, Irvine & Ohishi 1994, ApJ 422, 621
  535 HNCS    -- Adande, Halfen, Ziurys, Quan & Herbst 2010, ApJ 725, 561,
                "Observations of the [HNCS]/[HSCN] ratio in Sgr B2 and
                TMC-1: Evidence for low-temperature gas-phase chemistry" --
                paper's own text confirms "the J=7(0,7)->6(0,6) and
                8(0,8)->7(0,7) transitions of both HNCS and HSCN were
                detected in TMC-1, the first identification of either
                molecule in a cold, dark cloud." (A 2024 Cernicharo paper,
                A&A 682, L4, discovers a DIFFERENT isomer, HCNS/thiofulminic
                acid, and explicitly calls HNCS "well-known" by then --
                checked and ruled out as the HNCS/TMC-1 reference.)
  524 HC(O)CN -- Cernicharo et al. 2021, A&A 650, L14 (arXiv:2105.12996),
                "The Sulphur Saga in TMC-1: Discovery of HCSCN and
                HCSCCH" -- confirmed via the CDMS ISM:HCOCN reference page
                as this species' own listed TMC-1 reference (formyl
                cyanide/HCOCN reported in the same survey letter alongside
                the two headline sulfur species); also this project's own
                pre-existing v5-audit citation (table3_v5_audit_flags.csv rank
                45) and the 2026 TMC-1 review's citation -- triple
                convergence.
  268 NO      -- Gerin, Viala & Casoli 1993, A&A 268, 212, "The abundance
                of nitric oxide in TMC-1" (title confirmed via the 2026
                TMC-1 review, arXiv:2604.26439, which cites exactly this
                paper for NO/TMC-1). No DOI could be verified for this
                pre-DOI-era A&A Letter despite a genuine search (ADS/
                CrossRef); Ref uses the ADS bibcode. (McGonagle, Ziurys,
                Irvine & Minh 1990, ApJ 359, 121 explicitly reports a
                NON-detection of NO toward TMC-1 at that time -- checked
                and ruled out.)
  282 HCNH+   -- Schilke, Walmsley, Millar & Henkel 1991, A&A 247, 487,
                "Protonated HCN in molecular clouds" -- confirmed via CDMS
                (IRAM 30m, TMC-1) and the 2026 TMC-1 review's citation. No
                DOI could be verified despite a genuine search; Ref uses
                the ADS bibcode.
  765 CO      -- Pratap et al. 1997, ApJ 486, 862 (same paper as HCO+
                above). CORRECTION: the review checkpoint's original citation
                (Penzias, Solomon, Jefferts & Wilson 1972, ApJ 174, L43)
                was re-checked and found to report CO toward the BN/KL
                region of the Orion Nebula, NOT "Heiles Cloud 2"/TMC-1 --
                the "Heiles Cloud 2" identification in the review checkpoint note
                was an error. "Earliest TMC-1 report located" per item A.
  484 CCH     -- Wootten, Bozyan, Garrett, Loren & Snell 1980, ApJ 239, 844,
                "Detection of C2H in cold dark clouds" -- confirmed via its
                own ADS abstract: surveyed 30 dark-cloud regions including
                "L1534 (also called TMC 1)" with "extensive maps," and
                reported "definite detections of C2H ... in 19 of the 30
                regions examined." Item-B finding; candidate 484/C2H was
                not among the 87 anchors or Table 1 v7.
  493 HCO     -- Cernicharo, Marcelino, Agundez, Endo, Cabezas, Bermudez,
                Tercero & de Vicente 2020, A&A 642, L17, "Discovery of
                HC3O+ in space: The chemistry of O-bearing species in
                TMC-1" (DOI 10.1051/0004-6361/202039351; received
                2020-09-06, accepted 2020-09-30, online 2020-10-23).
                Surfaced during item 5's audit (see note above); confirmed
                via the 2026 TMC-1 review's Table 1/reference list, which
                cites exactly this paper for HCO's TMC-1 column density
                (firm, non-tentative). This project's general ISM/CDMS
                first-detection reference for the HCO formula (Snyder,
                Hollis & Ulich 1976, ApJ 208, L91, toward W3/W51/NGC2024C/
                K3-50) is NOT toward TMC-1, so the 2020 paper is used as
                "earliest TMC-1 report located" per item A's convention.

NOT relabeled (checked, rejected after verification -- logged, not silent):
  91, 11  -- candidate SMILES (C#CC#CC=O / C#CC#CC#CC=O, both mult=1 in
             Lee's own raw list) are the closed-shell aldehyde tautomers
             (2,4-pentadiynal / 2,4,6-heptatriynal, formula C5H2O/C7H2O).
             The astronomically-detected "HC5O"/"HC7O" radicals have
             formula C5HO/C7HO (one fewer H) -- that formula is ABSENT
             from the whole 1,439-candidate universe. Initial rule-C
             "matches" against my own hand-built SMILES for HC5O/HC7O were
             a self-inflicted construction error (wrong H count), not a
             real match; retracted after checking Lee's mult column and
             the candidates' actual formulas.
  856     -- SMILES N#CC1C=CC=C1 decodes to 5-cyanocyclopentadiene (CN on
             the sp3 ring carbon), a distinct, UNDETECTED isomer from the
             detected "2-cyanocyclopentadiene" (Lee et al. 2021, ApJL 908,
             L11; CN on an sp2 ring carbon not adjacent to the CH2). Per
             the user's explicit confirmation, this is NOT an uncertain
             identity and NOT a STOP condition -- it is a settled, firm
             negative (a real, distinct, undetected isomer), left
             unlabeled accordingly.
"""
import pandas as pd

v5 = pd.read_csv("<PROJECT_ROOT>/f12_implementation/candidates_purged_v5.csv")
assert len(v5) == 1439
assert int(v5["positive"].sum()) == 47

REMOVE_ANCHOR_DUP = [293, 376, 396, 490, 523]
ADD_POSITIVE = {
    383: ("1997-09-10", "HCO+", "Pratap et al. 1997, ApJ 486, 862; DOI 10.1086/304553"),
    298: ("1989-10-01", "H2S", "Minh, Irvine & Ziurys 1989, ApJ 345, L63; DOI 10.1086/185553"),
    272: ("1994-01-01", "NS", "McGonagle, Irvine & Ohishi 1994, ApJ 422, 621; DOI 10.1086/173755"),
    535: ("2010-12-01", "HNCS", "Adande, Halfen, Ziurys, Quan & Herbst 2010, ApJ 725, 561; "
          "DOI 10.1088/0004-637X/725/1/561"),
    # candidates_purged's date column uses day-padded full dates for datetime arithmetic,
    # per this project's established month-only convention (step2_build_candidates_purged_v5.py:
    # "month-only rule: day set to 1"); Table 1's Date column will separately display "2021-05".
    524: ("2021-05-01", "HC(O)CN", "Cernicharo et al. 2021, A&A 650, L14 (arXiv:2105.12996); "
          "DOI 10.1051/0004-6361/202141297"),
    268: ("1993-01-01", "NO", "Gerin, Viala & Casoli 1993, A&A 268, 212; "
          "ADS bibcode 1993A&A...268..212G (no DOI located)"),
    282: ("1991-07-01", "HCNH+", "Schilke, Walmsley, Millar & Henkel 1991, A&A 247, 487; "
          "ADS bibcode 1991A&A...247..487S (no DOI located)"),
    765: ("1997-09-10", "CO", "Pratap et al. 1997, ApJ 486, 862; DOI 10.1086/304553"),
    484: ("1980-01-01", "CCH", "Wootten, Bozyan, Garrett, Loren & Snell 1980, ApJ 239, 844; "
          "DOI 10.1086/158168"),
    493: ("2020-10-23", "HCO", "Cernicharo, Marcelino, Agundez, Endo, Cabezas, Bermudez, Tercero "
          "& de Vicente 2020, A&A 642, L17; DOI 10.1051/0004-6361/202039351"),
}

# Sanity: confirm pre-change state matches what item 2/B found.
for cid, expected_smiles in [
    (293, "[CH2+]O"), (376, "[C-]#C[C+]=N"), (396, "O=[C+]O"),
    (490, "C#C[C+]=N"), (523, "N#C[C+]=N"),
    (383, "C#[O+]"), (298, "S"), (272, "[N]=S"), (535, "N=C=S"),
    (524, "N#CC=O"), (268, "[N]=O"), (282, "[CH+]=N"),
    (765, "[C-]#[O+];[13C-]#[O+]"),  # multi-member candidate: CO + 13CO isotopologue;
                                      # rule-C match is on the CO member (n_members=2)
    (484, "[C]#C"), (493, "[CH]=O"),
]:
    row = v5[v5["candidate_id"] == cid].iloc[0]
    assert row["smiles_members"] == expected_smiles, (cid, row["smiles_members"])
    assert row["positive"] == False, (cid, "expected positive=False pre-patch")

v6 = v5[~v5["candidate_id"].isin(REMOVE_ANCHOR_DUP)].copy()
assert len(v6) == 1439 - len(REMOVE_ANCHOR_DUP) == 1434

for cid, (date_str, species, ref) in ADD_POSITIVE.items():
    mask = v6["candidate_id"] == cid
    assert mask.sum() == 1, cid
    v6.loc[mask, "positive"] = True
    v6.loc[mask, "date"] = date_str

n_positive = int(v6["positive"].sum())
assert n_positive == 47 + len(ADD_POSITIVE) == 57
print(f"candidates_purged_v6.csv: {len(v6)} rows (was {len(v5)}, removed {len(REMOVE_ANCHOR_DUP)} anchor dups)")
print(f"Positive candidates: {n_positive} (was {int(v5['positive'].sum())} in v5, +{len(ADD_POSITIVE)} new)")

# Explicitly confirm candidate 856 is untouched (left False, flagged not relabeled)
row856 = v6[v6["candidate_id"] == 856].iloc[0]
assert row856["positive"] == False
print(f"candidate 856 (settled negative -- distinct, undetected isomer, not a STOP case): "
      f"positive={row856['positive']}")

out_path = "<PROJECT_ROOT>/f13_implementation/candidates_purged_v6.csv"
v6.to_csv(out_path, index=False)
print(f"Wrote {out_path}")

print("\nNew positive rows:")
print(v6[v6["candidate_id"].isin(ADD_POSITIVE.keys())][
    ["candidate_id", "species_key", "member_ids", "smiles_members", "positive", "date"]
].to_string(index=False))
