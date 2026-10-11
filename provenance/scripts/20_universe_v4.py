# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
post2021_detections_v3_2.csv (77 rows) -> post2021_detections_v4.csv; v3.2 is not overwritten. All additions and
corrections were verified against primary sources.

ADD (completeness gaps confirmed absent from v3.2 and not a Lee Table A2 anchor):
  - CH3NC   (Tennis et al. 2023, MNRAS 525, 2154)
  - ethynylbenzene (Loru et al. 2023, A&A 677, A166)
  - HS2     (Esplugues et al. 2025, A&A 699, L2) -- first FIRM detection TOWARD TMC-1/a cold dark cloud
    specifically; NOT the first-ever ISM detection of HS2 (Horsehead PDR, Fuente et al. 2017/2025) -- flagged
    in `flags`, included per A1's TMC-1-scoped reading.

ADD (A3 known-at-origin only, never a positive/test candidate -- same treatment as the C3N-/C5N- rows):
  - CH3CO+  (Cernicharo et al. 2021c, A&A 646, L7) -- first-public-date 2021-01, before v3.2/v4's own
    Aug-2021 coverage start AND before all 4 test origins (earliest 2022-01-01) -- included for A3 "known at
    all origins" purposes only.

DO NOT ADD:
  - l-H2C3O  -- already a Lee Table A2 anchor (table_a2_authoritative.csv row 76, SMILES C=C=C=O, cd 11.0414,
    Cernicharo et al. 2020), not a new species.
  - C10H-    -- already in v3.2 (row 8, arxiv_id 2301.0776, a float-truncated match to 2301.07760 / Remijan
    et al. 2023).
  - HNC5     -- already in v3.2 (row 68); a metadata gap in v3.2, not a completeness gap.
  - cis-HCOOH -- identical SMILES to the anchor HCOOH (table_a2_authoritative.csv row 72, SMILES O=CO), so
    the same species key under the A1(iv) machinery, whatever the true spectroscopic conformer distinction.

CORRECT: 1-cyanopyrene's arxiv_id, 2410.0067 (float-truncated 2410.00670, itself wrong) -> 2410.00657 (the
  paper is arXiv:2410.00657, not arXiv:2410.00670, which is the separate 2-/4-cyanopyrene paper).
"""
import pandas as pd

v3_2 = pd.read_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3_2.csv")
df = v3_2.copy()
for c in df.columns:
    if df[c].dtype != object:
        df[c] = df[c].astype("object")

change_log_rows = []

# ---- 1. Fix 1-cyanopyrene arxiv_id (float truncation -> wrong paper entirely) ----
idx = df.index[df["formula"] == "1-cyanopyrene"]
assert len(idx) == 1
i = idx[0]
old_arxiv = df.loc[i, "arxiv_id"]
old_ref = df.loc[i, "cd_source"]
df.loc[i, "arxiv_id"] = "2410.00657"
df.loc[i, "cd_source"] = "Wenzel et al. 2024, Nature Astronomy (arXiv:2410.00657)"
df.loc[i, "flags"] = (
    "F9 CORRECTION: arxiv_id was 2410.0067 (float-truncated 2410.00670), itself the WRONG "
    "paper -- 2410.00670 is the separate 2-/4-cyanopyrene paper. Correct id for THIS "
    "species (1-cyanopyrene) is 2410.00657, confirmed via F8 Step 2 primary-source read."
)
change_log_rows.append(dict(
    formula="1-cyanopyrene", action="CORRECT_ARXIV_ID",
    old_arxiv_id=old_arxiv, new_arxiv_id="2410.00657",
    old_cd_source=old_ref, new_cd_source=df.loc[i, "cd_source"],
    evidence="F8 Task Step 2 (C08/R0097) primary-source read, this session's build_v4.py",
))

# ---- 2. ADD new genuine completeness-gap species ----
new_rows = [
    dict(
        formula="CH3NC", name="methyl isocyanide", smiles="C[N+]#[C-]",
        log10_column_density_cm2=None,
        cd_source="Tennis et al. 2023 (value not extracted this session -- flagged)",
        first_public_date="2023-08-03", first_public_date_basis="journal_accepted",
        arxiv_v1_date=None,  # re-verified: no arXiv preprint exists
        journal_received="2023-01-18", journal_accepted="2023-08-03",
        journal_online="2023-08-08", arxiv_id=None, doi="10.1093/mnras/stad2398",
        journal_ref="MNRAS 525(2), 2154-2171 (2023)",
        authors="J. D. Tennis, C. Xue, D. Talbi, P. B. Changala, M. L. Sita, B. A. McGuire, E. Herbst",
        telescope_survey="GBT (GOTHAM)", formal_charge=0, n_radical_electrons=0,
        class_="closed_shell", date_confidence="high",
        flags=(
            "F9 ADDITION. TASK TEXT ASSUMED an arXiv v1 date exists ('record arXiv v1 date') -- "
            "RE-VERIFIED TWICE (independently in F8 and again this session, both via direct arXiv "
            "title/author search and Semantic Scholar's externalIds field) that NO arXiv preprint "
            "exists for this paper. first_public_date_basis is therefore journal_accepted (A2's own "
            "fallback), not arxiv_v1 as the task text assumed -- flagged for Reham, not silently "
            "substituted. 13.4-sigma detection ('report detection of CH3NC in the cold prestellar "
            "core TMC-1 using the Green Bank Telescope'). Distinct structural isomer of the "
            "already-known TMC-1 species CH3CN (methyl cyanide), not an isotopologue of it. NOT the "
            "first-ever interstellar detection of CH3NC (prior: Sgr B2(OH) 1988 tentative, Sgr B2(N) "
            "2005, Horsehead PDR 2013, protostar IRAS 16293-2422 2018) -- this is the first detection "
            "in a cold/dark-cloud (prestellar) environment specifically, consistent with A1's "
            "TMC-1-scoped reading (same treatment as CH3CO+/HS2's prior-elsewhere detections)."
        ),
        detection_window_v2="2021-08_to_2023-09",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
        parsed_date="2023-08-03", date_precision="day",
    ),
    dict(
        formula="C6H5CCH", name="ethynylbenzene (phenylacetylene)", smiles="C#Cc1ccccc1",
        log10_column_density_cm2=12.48, cd_source="Loru et al. 2023 ((3.0+/-0.5)e12 cm^-2)",
        first_public_date="2023-07-17", first_public_date_basis="journal_accepted",
        arxiv_v1_date=None,  # re-verified: no arXiv preprint exists (3 independent checks)
        journal_received="2023-05-26", journal_accepted="2023-07-17",
        journal_online="2023-09-22", arxiv_id=None, doi="10.1051/0004-6361/202347023",
        journal_ref="A&A 677, A166 (2023)",
        authors="D. Loru, C. Cabezas, J. Cernicharo, M. Schnell, A. L. Steber",
        telescope_survey="QUIJOTE (Yebes 40m), presumed -- not independently confirmed this session",
        formal_charge=0, n_radical_electrons=0, class_="closed_shell", date_confidence="high",
        flags=(
            "F9 ADDITION. Re-verified: no arXiv preprint exists (3 independent checks: arXiv "
            "title/author search, Semantic Scholar externalIds, A&A article page). "
            "first_public_date_basis = journal_accepted (A2 fallback). Paper: 'We confirm the "
            "detection of ethynylbenzene in TMC-1' via 25 rotational lines, 'definitively been "
            "detected' -- upgrades a 2021 TENTATIVE report (Cernicharo et al. 2021c) to FIRM; this "
            "2023 paper is correctly the first FIRM detection under A1(iii). Companion species "
            "1,2-diethynylbenzene in the SAME paper is a non-detection/3-sigma-upper-limit only -- "
            "confirmed NOT added. No earlier astronomical detection of ethynylbenzene found anywhere "
            "(comets, circumstellar envelopes, protoplanetary disks) -- this is the first-ever firm "
            "astronomical detection of this species, not just first-toward-TMC-1."
        ),
        detection_window_v2="2021-08_to_2023-09",
        level_A=False, level_B=False, level_C=False, in_lee_list=False,
        parsed_date="2023-07-17", date_precision="day",
    ),
    dict(
        formula="HS2", name="disulfanyl radical", smiles="S[S]",
        log10_column_density_cm2=None,
        cd_source="Esplugues et al. 2025 (value not extracted this session -- flagged)",
        first_public_date="2025-06-15", first_public_date_basis="arxiv_v1",
        arxiv_v1_date="2025-06-15", journal_received=None, journal_accepted=None,
        journal_online="2025-07 (month only, from A&A URL path /2025/07/aa55384-25/)",
        arxiv_id="2506.12974", doi="10.1051/0004-6361/202555384",
        journal_ref="A&A 699, L2 (2025)",
        authors="G. Esplugues, M. Agundez, G. Molpeceres, B. Tercero, C. Cabezas, N. Marcelino, R. Fuentetaja, J. Cernicharo",
        telescope_survey="QUIJOTE (Yebes 40m), presumed -- cyanopolyyne-peak position of TMC-1, not independently confirmed this session",
        formal_charge=0, n_radical_electrons=1, class_="radical", date_confidence="high",
        flags=(
            "F9 ADDITION. 'We report the first detection of HS2 towards the cold dark cloud TMC-1.' "
            "4 lines >=4-sigma, 1 line at 3.5-sigma, confirmed via fine/hyperfine components of 2 "
            "rotational transitions (paper's own quoted significance). NOT the first-ever ISM "
            "detection of HS2 -- paper itself states 'HS2 has only been detected in the Horsehead "
            "photon-dominated region (Fuente et al. 2017, 2025)' -- a different source TYPE (PDR, not "
            "cold dark cloud); this TMC-1 paper is scoped as first-toward-a-cold-dark-cloud, "
            "consistent with A1's TMC-1-scoped reading (same treatment as CH3NC/CH3CO+). SMILES "
            "radical-site placement (terminal S vs the H-bearing S) is an analogy to HO2, not a "
            "quoted electronic-state label from the paper -- documented per A1(iv), not paper-verified. "
            "SMILES corrected from the verifying subagent's invalid RDKit syntax '[S]SH' to the valid, "
            "canonicalization-equivalent 'S[S]' (same structure: HS2, one bracket-radical terminal S)."
        ),
        detection_window_v2="post-2023-09",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
        parsed_date="2025-06-15", date_precision="day",
    ),
    dict(
        formula="CH3CO+", name="acetyl cation", smiles="CC#[O+]",
        log10_column_density_cm2=None, cd_source="Cernicharo et al. 2021c (value not re-extracted this session)",
        first_public_date="2021-01", first_public_date_basis="arxiv_id_month_only",
        arxiv_v1_date="2021-01 (month only, from arXiv ID prefix)",
        journal_received=None, journal_accepted=None, journal_online=None,
        arxiv_id="2101.04603", doi="10.1051/0004-6361/202140368",
        journal_ref="A&A 646, L7 (2021)", authors="J. Cernicharo et al. (full author list not re-verified this session)",
        telescope_survey="QUIJOTE (Yebes 40m), presumed",
        formal_charge=1, n_radical_electrons=0, class_="ion", date_confidence="medium",
        flags=(
            "F9 ADDITION, A3 KNOWN-AT-ORIGIN ONLY -- NEVER A POSITIVE/TEST CANDIDATE. First-public-date "
            "(Jan 2021, arXiv:2101.04603) predates v3.2/v4's own post-2021-08 target coverage window "
            "AND all 4 test origins (earliest 2022-01-01) -- known & date < every origin, so removed "
            "from every origin's test pool per protocol_v1_amendments A3, exactly as C3N-/C5N- (also "
            "pre-Aug-2021) are already treated in this same file. Carried from F8's Step 2 "
            "primary-source read (C17b); NOT independently re-verified from the primary source again "
            "this session beyond that prior read -- exact day, received/accepted/online dates, and "
            "column density not confirmed. R0210 (a 2022 conference proceeding) is NOT the discovery "
            "source and was excluded; true discovery is Cernicharo et al. 2021c."
        ),
        detection_window_v2="pre-2021-08 (F9/A3 addition -- NOT the '2021-08_to_2023-09' label used "
                             "for the existing pre-cutoff C3N-/C5N- rows; that label looks like a "
                             "provenance/query-window tag, not a true date bucket, and is not "
                             "reproduced here to avoid implying this row was found in that query)",
        level_A=True, level_B=True, level_C=True, in_lee_list=True,
        parsed_date="2021-01-15", date_precision="month",
    ),
]

new_df = pd.DataFrame(new_rows).rename(columns={"class_": "class"})
new_df = new_df[df.columns]  # enforce same column order

df = pd.concat([df, new_df], ignore_index=True)

for r in new_rows:
    change_log_rows.append(dict(
        formula=r["formula"], action="ADD",
        old_arxiv_id=None, new_arxiv_id=r["arxiv_id"],
        old_cd_source=None, new_cd_source=r["cd_source"],
        evidence=f"F9 primary-source verification this session (or carried from F8 Step 2 for CH3CO+); "
                 f"first_public_date={r['first_public_date']} ({r['first_public_date_basis']})",
    ))

out_path = "<PROJECT_ROOT>/f9_implementation/post2021_detections_v4.csv"
df.to_csv(out_path, index=False)
cl_path = "<PROJECT_ROOT>/f9_implementation/v4_change_log.csv"
pd.DataFrame(change_log_rows).to_csv(cl_path, index=False)

print(f"v4 written: {len(df)} rows ({len(v3_2)} in v3.2 + {len(new_rows)} added, 1 corrected)")
print(pd.DataFrame(change_log_rows)[["formula", "action"]].to_string(index=False))
