# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
post2021_detections_v4.csv (81 rows) -> post2021_detections_v5.csv (85 rows).

Adds 4 rows for the anchor/label-error candidates verified as detections: C5N (candidate_id 171), l-HC4N
(candidate_id 651), NS+ (candidate_id 263), HC3S+ (candidate_id 113). All four pass A1 (firm, TMC-1,
non-isotopologue, SMILES matches the discovery paper's stated connectivity/charge) and are dated before
2022-01-01, so they are labelled positive and known at every origin (A3), never a test positive.

Dates follow A2: arXiv v1 date where an arXiv version exists, at month level if the day is unconfirmed (all 4
have an arXiv id but the v1 posting DAY was not independently confirmed, so all 4 use the
"arxiv_id_month_only" basis label, as for CH3CO+ in 20_universe_v4.py).

detection_window_v2 follows the CH3CO+ precedent: these are pre-2021-08 (they predate Lee et al. 2021's own
anchor-table cutoff and all 4 test origins), so known at every origin per A3, removed from every origin's
test pool, never a test positive -- NOT the '2021-08_to_2023-09' query-window tag used for in-window rows.

level_A/level_B/level_C/in_lee_list are placeholder True here (all 4 are known matches to the 1,501-row raw
universe) and are unconditionally recomputed against the TRUE 1,501-row raw universe by 31_labels_v5.py, as
21_labels_v4.py did for v4.

Does not touch post2021_detections_v4.csv.
"""
import pandas as pd

v4 = pd.read_csv("<PROJECT_ROOT>/f9_implementation/post2021_detections_v4.csv", dtype={"arxiv_id": str})
assert len(v4) == 81

PRE2022_NOTE = (
    "pre-2021-08 (F12c addition -- known at every origin per A3, never a test positive; "
    "NOT the '2021-08_to_2023-09' query-window tag used for in-window rows, matching the "
    "CH3CO+ precedent in post2021_detections_v4.csv)"
)

new_rows = [
    dict(
        formula="C5N", name="C5N radical", smiles="[C]#CC#CC#N",
        log10_column_density_cm2=11.49,
        cd_source="Guelin, Neininger & Cernicharo (1998), A&A 335, L1: N=3.1e11 cm-2",
        first_public_date="1998-05", first_public_date_basis="arxiv_id_month_only",
        arxiv_v1_date="1998-05 (month only, from arXiv ID prefix)",
        journal_received=None, journal_accepted=None, journal_online=None,
        arxiv_id="astro-ph/9805105", doi=None,
        journal_ref="A&A 335, L1 (1998)",
        authors="M. Guelin, N. Neininger, J. Cernicharo",
        telescope_survey="IRAM 30m",
        formal_charge=0, n_radical_electrons=1, class_="radical",
        date_confidence="medium",
        flags=(
            "F12c ADDITION (Xue et al. 2025 TMC-1 cross-check, F12 item 1 -> F12b -> F12c "
            "verification chain). Firm detection, TMC-1, matches candidate_id 171 by Rule A "
            "(exact canonical SMILES). No DOI (pre-DOI-era A&A Letter). Labelled positive, "
            "known at every origin per A3 (pre-2021-08, well before all 4 test origins)."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="1998-05-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="C4HN", name="l-HC4N radical (linear)", smiles="C#C[C+]=C=[N-]",
        log10_column_density_cm2=11.57,
        cd_source="Cernicharo et al. (2021), A&A 647, L2: N=(3.7+/-0.4)e11 cm-2",
        first_public_date="2021-02", first_public_date_basis="arxiv_id_month_only",
        arxiv_v1_date="2021-02 (month only, from arXiv ID prefix)",
        journal_received=None, journal_accepted="2021-02-15", journal_online=None,
        arxiv_id="2102.10920", doi="10.1051/0004-6361/202140434",
        journal_ref="A&A 647, L2 (2021)",
        authors=(
            "J. Cernicharo, M. Agundez, C. Cabezas, N. Marcelino, B. Tercero, J. R. Pardo, "
            "J. D. Gallego, F. Tercero, J. A. Lopez-Perez, P. de Vicente"
        ),
        telescope_survey="QUIJOTE (Yebes 40m)",
        formal_charge=0, n_radical_electrons=1, class_="radical",
        date_confidence="medium",
        flags=(
            "F12c ADDITION (F12 item 1 -> F12b -> F12c verification chain). Firm detection, "
            "TMC-1, matches candidate_id 651 by Rule C (InChI-skeleton match). Labelled "
            "positive, known at every origin per A3 (pre-2021-08, well before all 4 test "
            "origins). Journal-accepted day (2021-02-15) is known even though the arXiv v1 "
            "posting day was not independently reconfirmed this session; A2 basis is "
            "arxiv_id_month_only per Reham's explicit approval-message rule."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="2021-02-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="NS+", name="NS+ cation", smiles="N#[S+]",
        log10_column_density_cm2=10.72,
        cd_source="Cernicharo et al. (2018), ApJL 853, L22: N=5.2e10 cm-2",
        first_public_date="2018-01", first_public_date_basis="arxiv_id_month_only",
        arxiv_v1_date="2018-01 (month only, from arXiv ID prefix)",
        journal_received=None, journal_accepted=None, journal_online=None,
        arxiv_id="1801.05559", doi="10.3847/2041-8213/aaa83a",
        journal_ref="ApJL 853, L22 (2018)",
        authors="J. Cernicharo et al.",
        telescope_survey="IRAM 30m / Yebes 40m",
        formal_charge=1, n_radical_electrons=0, class_="ion",
        date_confidence="medium",
        flags=(
            "F12c ADDITION (anchor charge-audit item 1, resolving the HCCN charge-encoding "
            "discrepancy raised in F12b). Firm detection, TMC-1, matches candidate_id 263. "
            "Labelled positive, known at every origin per A3 (pre-2021-08, well before all 4 "
            "test origins)."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="2018-01-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
    dict(
        formula="C3HS+", name="HC3S+ cation", smiles="C#CC#[S+]",
        log10_column_density_cm2=11.60,
        cd_source="Cernicharo et al. (2021), A&A 646, L3: N=(4.00+/-1.50)e11 cm-2",
        first_public_date="2021-01", first_public_date_basis="arxiv_id_month_only",
        arxiv_v1_date="2021-01 (month only, from arXiv ID prefix)",
        journal_received=None, journal_accepted="2021-01-09", journal_online=None,
        arxiv_id="2101.05163", doi="10.1051/0004-6361/202040013",
        journal_ref="A&A 646, L3 (2021)",
        authors=(
            "J. Cernicharo, C. Cabezas, M. Agundez, B. Tercero, J. R. Pardo, Y. Endo, "
            "N. Marcelino, P. de Vicente"
        ),
        telescope_survey="QUIJOTE (Yebes 40m)",
        formal_charge=1, n_radical_electrons=0, class_="ion",
        date_confidence="medium",
        flags=(
            "F12c ADDITION (McGuire/astromol census check, item 3). Firm detection, TMC-1, "
            "matches candidate_id 113. Labelled positive, known at every origin per A3 "
            "(pre-2021-08, well before all 4 test origins). Journal-accepted day "
            "(2021-01-09) is known even though the arXiv v1 posting day was not "
            "independently reconfirmed this session; A2 basis is arxiv_id_month_only per "
            "Reham's explicit approval-message rule."
        ),
        detection_window_v2=PRE2022_NOTE,
        parsed_date="2021-01-15", date_precision="month",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
    ),
]

new_df = pd.DataFrame(new_rows).rename(columns={"class_": "class"})
new_df = new_df[v4.columns]

v5 = pd.concat([v4, new_df], ignore_index=True)
assert len(v5) == 85

out_path = "<PROJECT_ROOT>/f12_implementation/post2021_detections_v5.csv"
v5.to_csv(out_path, index=False)
print(f"Wrote {out_path}: {len(v5)} rows ({len(v4)} v4 + {len(new_df)} new)")
print(new_df[["formula", "smiles", "first_public_date", "first_public_date_basis"]].to_string(index=False))
