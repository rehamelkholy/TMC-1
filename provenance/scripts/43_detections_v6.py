# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
post2021_detections_v5.csv (85 rows) ->
post2021_detections_v6.csv (94 rows).

Adds the 9 TMC-1-inventory matches approved at the review checkpoint + item A/B
additions (see 42_candidates_v6.py for
the full match/disposition rationale and reference-verification notes --
duplicated here only at the summary level). All 9 are dated before every
origin (2022-01-01), known at every origin per A3, never a test positive.

References were independently re-verified for this build (NASA ADS, CDMS
molecule pages, CrossRef, and a 2026 TMC-1 review, arXiv:2604.26439). Two of
the review checkpoint's original candidate references did NOT survive
re-verification and were replaced:
  - HCO+: Langer, Wilson, Henry & Guelin 1978 (ApJ 225, L139) reports ONLY
    H13CO+/HC18O+ toward TMC-1 per its own abstract -- disqualified under
    rule A1 (non-isotopologue). Replaced with Pratap et al. 1997 (ApJ 486,
    862), this project's own pre-existing standing HCO+/TMC-1 citation.
  - CO: Penzias, Solomon, Jefferts & Wilson 1972 (ApJ 174, L43) reports CO
    toward the Orion Nebula (BN/KL), NOT "Heiles Cloud 2"/TMC-1 -- the
    "Heiles Cloud 2" identification in the review checkpoint note was an error.
    Replaced with Pratap et al. 1997 (same paper as HCO+ -- "earliest
    verifiable TMC-1 report" per item A).

log10_column_density_cm2 / cd_source are left BLANK for all 9 new rows:
column-density values were not independently re-extracted/verified from
the primary sources during this step (the research focus of this step was
first-detection dating and reference verification, not column densities),
and fabricating a number here would be worse than leaving it blank -- a
future work should fill these from the cited primary sources directly.

level_A/level_B/level_C/in_lee_list are placeholder True here (all 9 are
confirmed matches to the 1,501-row raw universe, verified directly via
lee_list_entry_labels_v5.csv's own ids: 275, 279, 290, 307, 396, 503, 545,
556, 800 -- see candidates_purged_v6.csv's member_ids for the same 9 rows)
-- consistent with this project's established convention (step1_build_
post2021_detections_v5.py), these placeholders would be recomputed from
scratch by a lee_list_entry_labels refresh step if/when one is run; not
done in this step (out of scope -- the v6 pipeline re-run consumes
candidates_purged_v6.csv and post2021_detections_v6.csv directly, plus a
purged-raw-universe file that only needs the 5 anchor-dup raw ids removed,
not a full rule A/B/C relabel of the 1,501-row universe).

Does not touch post2021_detections_v5.csv.
"""
import pandas as pd

v5 = pd.read_csv("<PROJECT_ROOT>/f12_implementation/post2021_detections_v5.csv", dtype={"arxiv_id": str})
assert len(v5) == 85

PRE2022_NOTE = (
    "pre-2021-08 (F13 item-2/B addition -- known at every origin per A3, never a test "
    "positive; NOT the '2021-08_to_2023-09' query-window tag used for in-window rows, "
    "matching the CH3CO+ precedent in post2021_detections_v4.csv)"
)

new_rows = [
    dict(
        formula="NO", name="nitric oxide", smiles="[N]=O",
        log10_column_density_cm2=None,
        cd_source=None,
        first_public_date="1993-01", first_public_date_basis="journal_month_only",
        arxiv_v1_date=None,
        journal_received=None, journal_accepted=None, journal_online="1993-01",
        arxiv_id=None, doi=None,
        journal_ref="A&A 268, 212 (1993) [ADS bibcode 1993A&A...268..212G; no DOI located]",
        authors="M. Gerin, Y. Viala, F. Casoli",
        telescope_survey="unconfirmed this pass",
        formal_charge=0, n_radical_electrons=1, class_="radical",
        date_confidence="medium",
        flags=(
            "F13 ITEM 2/3 ADDITION (rule-C inventory match, candidate_id 268). Title "
            "'The abundance of nitric oxide in TMC-1' confirmed via a 2026 TMC-1 review "
            "(arXiv:2604.26439), which cites exactly this paper for NO/TMC-1. No DOI "
            "could be verified for this pre-DOI-era A&A Letter despite a genuine search "
            "(ADS/CrossRef robots-blocked / no match); Ref uses the ADS bibcode. "
            "McGonagle, Ziurys, Irvine & Minh 1990 (ApJ 359, 121) explicitly reports a "
            "NON-detection of NO toward TMC-1 at that earlier time -- checked and ruled "
            "out as an earlier TMC-1 report."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="1993-01-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="NS", name="nitrogen sulfide (NS radical)", smiles="[N]=S",
        log10_column_density_cm2=None,
        cd_source=None,
        first_public_date="1994-01", first_public_date_basis="journal_month_only",
        arxiv_v1_date=None,
        journal_received=None, journal_accepted=None, journal_online="1994-01",
        arxiv_id=None, doi="10.1086/173755",
        journal_ref="ApJ 422, 621 (1994)",
        authors="D. McGonagle, W. M. Irvine, M. Ohishi",
        telescope_survey="unconfirmed this pass",
        formal_charge=0, n_radical_electrons=1, class_="radical",
        date_confidence="medium",
        flags=(
            "F13 ITEM 2/3 ADDITION (rule-C inventory match, candidate_id 272). Title "
            "'Nitrogen sulfide in quiescent dark clouds' confirmed; paper reports NS "
            "'in TMC-1 and L134', the first TMC-1 observation per CDMS's own NS "
            "reference page. DOI confirmed via CrossRef."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="1994-01-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="CH2N+", name="protonated hydrogen cyanide (HCNH+)", smiles="[CH+]=N",
        log10_column_density_cm2=None,
        cd_source=None,
        first_public_date="1991-07", first_public_date_basis="journal_month_only",
        arxiv_v1_date=None,
        journal_received=None, journal_accepted=None, journal_online="1991-07",
        arxiv_id=None, doi=None,
        journal_ref="A&A 247, 487 (1991) [ADS bibcode 1991A&A...247..487S; no DOI located]",
        authors="P. Schilke, C. M. Walmsley, T. J. Millar, C. Henkel",
        telescope_survey="IRAM 30m",
        formal_charge=1, n_radical_electrons=0, class_="ion",
        date_confidence="medium",
        flags=(
            "F13 ITEM 2/3 ADDITION (rule-C inventory match, candidate_id 282). Title "
            "'Protonated HCN in molecular clouds'. TMC-1 detection (IRAM 30m) confirmed "
            "via CDMS's own HCNH+ reference page and a 2026 TMC-1 review "
            "(arXiv:2604.26439), which cites exactly this paper for HCNH+/TMC-1. No DOI "
            "could be verified despite a genuine search; Ref uses the ADS bibcode."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="1991-07-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="H2S", name="hydrogen sulfide", smiles="S",
        log10_column_density_cm2=None,
        cd_source=None,
        first_public_date="1989-10", first_public_date_basis="journal_month_only",
        arxiv_v1_date=None,
        journal_received=None, journal_accepted=None, journal_online="1989-10",
        arxiv_id=None, doi="10.1086/185553",
        journal_ref="ApJ 345, L63 (1989)",
        authors="Y. C. Minh, W. M. Irvine, L. M. Ziurys",
        telescope_survey="unconfirmed this pass",
        formal_charge=0, n_radical_electrons=0, class_="closed-shell",
        date_confidence="medium",
        flags=(
            "F13 ITEM 2/3 ADDITION (rule-C inventory match, candidate_id 298). Title "
            "'Detection of Interstellar Hydrogen Sulfide in Cold, Dark Clouds'. TMC-1 "
            "among the surveyed cold dark clouds (CDMS H2S reference page, DOI "
            "confirmed). This project's own F12f audit (table3_v5_audit_flags.csv rank "
            "38) independently cites the same paper."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="1989-10-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="CHO+", name="HCO+ cation (formylium)", smiles="C#[O+]",
        log10_column_density_cm2=None,
        cd_source=None,
        first_public_date="1997-09-10", first_public_date_basis="journal_online_day",
        arxiv_v1_date=None,
        journal_received=None, journal_accepted=None, journal_online="1997-09-10",
        arxiv_id=None, doi="10.1086/304553",
        journal_ref="ApJ 486, 862 (1997)",
        authors="P. Pratap, J. E. Dickens, R. L. Snell, M. P. Miralles, E. A. Bergin, "
                 "W. M. Irvine, F. P. Schloerb",
        telescope_survey="FCRAO 14m (TMC-1 spectral line survey)",
        formal_charge=1, n_radical_electrons=0, class_="ion",
        date_confidence="medium",
        flags=(
            "F13 ITEM 2/3 ADDITION (rule-C inventory match, candidate_id 383). "
            "CORRECTED REFERENCE: the F13 checkpoint's original citation (Langer, "
            "Wilson, Henry & Guelin 1978, ApJ 225, L139) was re-checked against its own "
            "abstract and found to report ONLY H13CO+ and HC18O+ toward TMC-1/L134N, "
            "not the main HCO+ line -- disqualified under rule A1 (non-isotopologue). "
            "Pratap et al. 1997, 'A Study of the Physics and Chemistry of TMC-1', is "
            "this project's own pre-existing standing TMC-1/HCO+ citation (already used "
            "in F12f's table3_v5_audit_flags.csv rank-31 entry) and the field-standard "
            "citation per a 2026 TMC-1 review (arXiv:2604.26439). 'Earliest TMC-1 report "
            "located' per item A -- no earlier verifiable firm main-HCO+/TMC-1 report "
            "found despite a dedicated research pass (candidates checked and ruled out: "
            "Buhl & Snyder 1970 [W3(OH)/Orion/SgrA/W51/L134, not TMC-1]; Guelin/Langer/"
            "Snell/Wootten 1977 DCO+ discovery [L63/L134/L134N, not TMC-1]; Guelin/"
            "Langer/Wilson 1982 [D13CO+ in TMC-1, an isotopologue]; Loren & Wootten 1980 "
            "HCO+ self-reversals [source list unconfirmed])."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="1997-09-10", date_precision="day",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="C2H", name="ethynyl radical (CCH)", smiles="[C]#C",
        log10_column_density_cm2=None,
        cd_source=None,
        first_public_date="1980", first_public_date_basis="journal_year_only",
        arxiv_v1_date=None,
        journal_received=None, journal_accepted=None, journal_online="1980",
        arxiv_id=None, doi="10.1086/158168",
        journal_ref="ApJ 239, 844 (1980)",
        authors="A. Wootten, E. P. Bozyan, D. B. Garrett, R. B. Loren, R. L. Snell",
        telescope_survey="unconfirmed this pass (dark-cloud survey, 30 regions)",
        formal_charge=0, n_radical_electrons=1, class_="radical",
        date_confidence="medium",
        flags=(
            "F13 ITEM-B ADDITION (2013-2021 ADS screening cross-check surfaced CCH as a "
            "rule-C candidate match not previously in Table 1 v7 or the 87 anchors; "
            "candidate_id 484). Title 'Detection of C2H in cold dark clouds'; own "
            "abstract confirms 'L1534 (also called TMC 1)' among 30 surveyed dark-cloud "
            "regions with 'extensive maps', and 'definite detections of C2H ... in 19 of "
            "the 30 regions examined'. DOI confirmed via CrossRef. Month/day not "
            "independently confirmed this pass; year-level date used."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="1980-06-15", date_precision="year",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="C2HNO", name="formyl cyanide (cyanoformaldehyde, HC(O)CN)", smiles="N#CC=O",
        log10_column_density_cm2=None,
        cd_source=None,
        first_public_date="2021-05", first_public_date_basis="arxiv_id_month_only",
        arxiv_v1_date="2021-05 (month only, from arXiv ID prefix)",
        journal_received=None, journal_accepted=None, journal_online=None,
        arxiv_id="2105.12996", doi="10.1051/0004-6361/202141297",
        journal_ref="A&A 650, L14 (2021)",
        authors="J. Cernicharo, C. Cabezas, Y. Endo, M. Agundez, B. Tercero, "
                 "J. R. Pardo, N. Marcelino, P. de Vicente",
        telescope_survey="QUIJOTE (Yebes 40m)",
        formal_charge=0, n_radical_electrons=0, class_="closed-shell",
        date_confidence="medium",
        flags=(
            "F13 ITEM 2/3 ADDITION (rule-C inventory match, candidate_id 524). Title "
            "'The Sulphur Saga in TMC-1: Discovery of HCSCN and HCSCCH' -- formyl "
            "cyanide (HCOCN) is reported in the same survey letter alongside the two "
            "headline sulfur species (its abundance is used as a direct comparison "
            "point in the paper's own abstract). Confirmed via the CDMS ISM:HCOCN "
            "reference page as this species' own listed TMC-1 reference; also this "
            "project's own pre-existing F12f citation (table3_v5_audit_flags.csv rank "
            "45) and the 2026 TMC-1 review's citation -- triple convergence. arXiv v1 "
            "day not independently confirmed this pass; month-level date used, per this "
            "project's established 'arxiv_id_month_only' convention."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="2021-05-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="CHNS", name="isothiocyanic acid (HNCS)", smiles="N=C=S",
        log10_column_density_cm2=None,
        cd_source=None,
        first_public_date="2010-12", first_public_date_basis="journal_month_only",
        arxiv_v1_date=None,
        journal_received=None, journal_accepted=None, journal_online="2010-12",
        arxiv_id=None, doi="10.1088/0004-637X/725/1/561",
        journal_ref="ApJ 725, 561 (2010)",
        authors="G. R. Adande, D. T. Halfen, L. M. Ziurys, D. Quan, E. Herbst",
        telescope_survey="ARO 12m / KP 10m (unconfirmed exact telescope this pass)",
        formal_charge=0, n_radical_electrons=0, class_="closed-shell",
        date_confidence="medium",
        flags=(
            "F13 ITEM 2/3 ADDITION (rule-C inventory match, candidate_id 535). "
            "CORRECTED TITLE from the F13 checkpoint: actual title is 'Observations of "
            "the [HNCS]/[HSCN] ratio in Sgr B2 and TMC-1: Evidence for low-temperature "
            "gas-phase chemistry'. Paper's own text confirms 'the J=7(0,7)->6(0,6) and "
            "8(0,8)->7(0,7) transitions of both HNCS and HSCN were detected in TMC-1, "
            "the first identification of either molecule in a cold, dark cloud.' "
            "Checked and ruled out: a 2024 Cernicharo paper (A&A 682, L4, DOI "
            "10.1051/0004-6361/202349105) discovers a DIFFERENT isomer (HCNS, "
            "thiofulminic acid) and explicitly calls HNCS 'well-known' in interstellar "
            "clouds by that point -- NOT the first HNCS/TMC-1 report, despite appearing "
            "in this project's own prior F12f audit table under that assumption; that "
            "F12f-table entry should be corrected in a future task."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="2010-12-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="CO", name="carbon monoxide", smiles="[C-]#[O+]",
        log10_column_density_cm2=None,
        cd_source=None,
        first_public_date="1997-09-10", first_public_date_basis="journal_online_day",
        arxiv_v1_date=None,
        journal_received=None, journal_accepted=None, journal_online="1997-09-10",
        arxiv_id=None, doi="10.1086/304553",
        journal_ref="ApJ 486, 862 (1997)",
        authors="P. Pratap, J. E. Dickens, R. L. Snell, M. P. Miralles, E. A. Bergin, "
                 "W. M. Irvine, F. P. Schloerb",
        telescope_survey="FCRAO 14m (TMC-1 spectral line survey)",
        formal_charge=0, n_radical_electrons=0, class_="closed-shell",
        date_confidence="medium",
        flags=(
            "F13 ITEM 2/3 ADDITION (rule-C inventory match, candidate_id 765; this "
            "candidate also carries a 13CO isotopologue member, raw id 801, unaffected "
            "by this relabel -- the candidate-level label applies to the group as a "
            "whole, consistent with this project's existing multi-member-candidate "
            "convention). CORRECTED REFERENCE: the F13 checkpoint's original citation "
            "(Penzias, Solomon, Jefferts & Wilson 1972, ApJ 174, L43) was re-checked and "
            "found to report CO toward the BN/KL region of the Orion Nebula, NOT "
            "'Heiles Cloud 2'/TMC-1 -- the 'Heiles Cloud 2' identification in the "
            "checkpoint note was an error (that paper's own discovery narrative "
            "explicitly names Orion). Replaced with Pratap et al. 1997 (same paper as "
            "HCO+ above) as the 'earliest TMC-1 report located' per item A."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="1997-09-10", date_precision="day",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="CHO", name="formyl radical", smiles="[CH]=O",
        log10_column_density_cm2=12.041,
        cd_source="Cernicharo et al. 2020, A&A 642, L17, Table 1 of the 2026 TMC-1 "
                   "review (arXiv:2604.26439), N(HCO)=1.1e12 cm^-2",
        first_public_date="2020-10-23", first_public_date_basis="journal_online_day",
        arxiv_v1_date=None,
        journal_received="2020-09-06", journal_accepted="2020-09-30", journal_online="2020-10-23",
        arxiv_id=None, doi="10.1051/0004-6361/202039351",
        journal_ref="A&A 642, L17 (2020)",
        authors="J. Cernicharo, N. Marcelino, M. Agundez, Y. Endo, C. Cabezas, "
                 "C. Bermudez, B. Tercero, P. de Vicente",
        telescope_survey="Yebes 40m (QUIJOTE line survey)",
        formal_charge=0, n_radical_electrons=1, class_="radical",
        date_confidence="medium",
        flags=(
            "F13 ITEM 5 ADDITION (surfaced during the top-100-entrant TMC-1-detection "
            "check, candidate_id 493, rank 86 of the v6 U1 ranking; NOT part of the "
            "original item 2/3 rule-C pass -- a genuine omission from this task's own "
            "item-1 TMC-1 inventory, since neutral HCO was never added to the "
            "INVENTORY dict in step2_match_rule_c.py, and it is also absent from Lee's "
            "87 anchors (table_a2_authoritative_v2.csv has only HCOOH, not HCO)). "
            "Confirmed via the 2026 TMC-1 review (arXiv:2604.26439), Table 1: HCO "
            "listed with N=1.1e12 cm^-2, not flagged tentative, citing this exact "
            "paper (ref 43) as the source of the adopted TMC-1 value. An initial "
            "research pass misattributed this finding to a non-existent 'Agundez et "
            "al. 2021' paper (confused with the unrelated A&A 656, L21 Cernicharo "
            "2021 letter, aa42634-21); corrected after directly checking the review's "
            "own reference list and the actual A&A page for this 2020 paper "
            "(DOI 10.1051/0004-6361/202039351). The general ISM first-detection "
            "reference for the HCO formula (Snyder, Hollis & Ulich 1976, ApJ 208, "
            "L91, toward W3/W51/NGC2024C/K3-50) is not toward TMC-1, so this 2020 "
            "paper is used as 'earliest TMC-1 report located' per item A's "
            "convention; log10_column_density_cm2 above is this project's own "
            "log10 of the review's adopted linear value (not independently "
            "re-derived from the primary paper's own tables)."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="2020-10-23", date_precision="day",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
]

new_df = pd.DataFrame(new_rows).rename(columns={"class_": "class"})
new_df = new_df[v5.columns]

v6 = pd.concat([v5, new_df], ignore_index=True)
assert len(v6) == 95

out_path = "<PROJECT_ROOT>/f13_implementation/post2021_detections_v6.csv"
v6.to_csv(out_path, index=False)
print(f"Wrote {out_path}: {len(v6)} rows ({len(v5)} v5 + {len(new_df)} new)")
print(new_df[["formula", "name", "smiles", "first_public_date", "first_public_date_basis", "doi"]].to_string(index=False))
