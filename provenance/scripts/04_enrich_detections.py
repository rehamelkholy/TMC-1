# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Enrich the 76-row base table with name, verified SMILES, log10 CD + source, first-public-date (+evidence),
DOI/arXiv id, telescope/survey, class and date_confidence, from a primary-source literature-research pass.

Sourcing legend for date_confidence:
  high        = arXiv v1 date read directly off arXiv's own metadata (primary)
  medium      = arXiv month/ID known + journal received/accepted/online dates bracket it, but the exact
                v1 day is from a third-party mirror or inferred, not read directly off arxiv.org itself
  low         = only journal accepted/received dates available (no arXiv), or a pre-existing citation
                whose arXiv id/DOI was carried over without independent re-verification
  unresolved  = genuine open question (see flags column) -- do not treat the date as settled
"""
import pandas as pd
from rdkit import Chem
import hashlib

df = pd.read_csv("<PROJECT_ROOT>/f3_validation_freeze/_v3_base_76.csv")

# citation registry: short-form reference key (as substring-matched against the 'reference' column)
# -> dict of resolved fields. All data below is transcribed from primary-source fetches of arXiv and
# A&A/ApJ pages, consolidated by hand.
CIT = {
 "Agundez et al. 2021": dict(
    authors="M. Agundez, N. Marcelino, B. Tercero, C. Cabezas, P. de Vicente, J. Cernicharo",
    journal="A&A 649, L4 (2021)", doi="10.1051/0004-6361/202140978", arxiv="2104.11506",
    arxiv_date="2021-04-23", jrecv="2021-04-01", jacc="2021-04-21", jonline="2021-05-05",
    telescope="QUIJOTE (Yebes 40m)", date_conf="high", flags=""),
 "Siebert et al. 2022": dict(
    authors="Siebert et al. (not independently re-verified this session)",
    journal="A&A (2022)", doi="", arxiv="2110.13152", arxiv_date="",
    jrecv="", jacc="", jonline="", telescope="QUIJOTE (presumed, not re-verified)",
    date_conf="low", flags="arXiv id carried over from prior project pass (excluded_v2.csv); "
      "full author list/DOI/exact date not independently re-verified this session"),
 "Sita et al. 2022": dict(
    authors="M. L. Sita et al.", journal="ApJL 938, L12 (2022)", doi="",
    arxiv="2209.06851", arxiv_date="", jrecv="", jacc="", jonline="",
    telescope="GOTHAM (GBT)", date_conf="low",
    flags="full journal ref resolved in a prior session (references_and_sources.md SC item16); "
      "exact arXiv v1 date not independently pulled this session"),
 "Remijan et al. 2023": dict(
    authors="A. Remijan, H. N. Scolati, A. M. Burkhardt, et al.",
    journal="ApJL 944, L45 (2023)", doi="10.3847/2041-8213/acb648", arxiv="2301.07760",
    arxiv_date="", jrecv="", jacc="", jonline="", telescope="GOTHAM (GBT)",
    date_conf="low", flags="resolved in a prior session (SC item 25); exact arXiv v1 day not "
      "independently re-pulled this session"),
 "Cooke et al. 2023": dict(
    authors="Cooke et al. (not independently re-verified this session)",
    journal="A&A (2023)", doi="", arxiv="2303.12221", arxiv_date="", jrecv="", jacc="",
    jonline="", telescope="GOTHAM (presumed, not re-verified)", date_conf="low",
    flags="arXiv id carried over from prior project pass (excluded_v2.csv); full bib not "
      "independently re-verified this session"),
 "Cernicharo et al. 2021d": dict(
    authors="J. Cernicharo, M. Agundez, R. I. Kaiser, C. Cabezas, B. Tercero, N. Marcelino, "
            "J. R. Pardo, P. de Vicente",
    journal="A&A 652, L9 (2021)", doi="10.1051/0004-6361/202141660", arxiv="2108.02308",
    arxiv_date="2021-08-04", jrecv="2021-06-28", jacc="2021-07-30", jonline="2021-08-18",
    telescope="QUIJOTE (Yebes 40m)", date_conf="medium",
    flags="arXiv v1 date from third-party mirror (alphaxiv), not read directly off arxiv.org "
      "(that fetch was blocked repeatedly this session); brackets consistently with journal dates"),
 "Cernicharo et al. 2021c": dict(
    authors="J. Cernicharo, M. Agundez, C. Cabezas, B. Tercero, N. Marcelino, R. Fuentetaja, "
            "J. R. Pardo, P. de Vicente",
    journal="A&A 656, L21 (2021)", doi="10.1051/0004-6361/202142634", arxiv="2112.01130",
    arxiv_date="2021-12-02", jrecv="2021-11-10", jacc="2021-11-26", jonline="2021-12-17",
    telescope="QUIJOTE (Yebes 40m)", date_conf="medium",
    flags="arXiv v1 date from third-party mirror, not primary-confirmed; C5O's electronic "
      "ground state (radical vs closed-shell) is NOT stated anywhere in this paper's text -- "
      "abstract calls companion species 'the radical HCCCO' but never applies 'radical' to C5O "
      "-- current SMILES [C]=C=C=C=C=O (radical) is UNVERIFIED, not confirmed or contradicted"),
 "Cernicharo et al. 2022a": dict(
    authors="J. Cernicharo, M. Agundez, C. Cabezas, R. Fuentetaja, B. Tercero, N. Marcelino, "
            "Y. Endo, J. R. Pardo, P. de Vicente",
    journal="A&A 657, L16 (2022)", doi="10.1051/0004-6361/202142992", arxiv="",
    arxiv_date="", jrecv="", jacc="", jonline="2022 (year only, from citation label "
      "'2022a'; exact date not located, no journal dates or arXiv id found)",
    telescope="QUIJOTE (Yebes 40m)", date_conf="unresolved",
    flags="resolved in a prior session (SC item 26); arXiv id not captured there and not "
      "independently re-searched this session; no exact date located at all -- half-year "
      "bucket for C5H+/C3H+ is UNCERTAIN (could be either half of 2022)"),
 "Fuentetaja et al. 2022b": dict(
    authors="Fuentetaja et al. (not independently re-verified this session)",
    journal="A&A (2022)", doi="", arxiv="2206.08992", arxiv_date="", jrecv="", jacc="",
    jonline="", telescope="QUIJOTE (presumed, not re-verified)", date_conf="low",
    flags="arXiv id carried over from prior project pass (excluded_v2.csv)"),
 "Cernicharo et al. 2022c": dict(
    authors="Cernicharo et al. (not independently re-verified this session)",
    journal="A&A, Sept 2022", doi="", arxiv="2207.09369", arxiv_date="",
    jrecv="", jacc="", jonline="", telescope="QUIJOTE",
    date_conf="low", flags="confirmed via full-text fetch in a prior session (excluded_v2.csv): "
      "'we report the detection of fulvenallene (c-C5H4CCH2)...with the QUIJOTE line survey'"),
 "Fuentetaja et al. 2022a": dict(
    authors="R. Fuentetaja, M. Agundez, C. Cabezas, B. Tercero, N. Marcelino, J. R. Pardo, "
            "P. de Vicente, J. Cernicharo",
    journal="A&A 667, L4 (2022)", doi="10.1051/0004-6361/202245018", arxiv="2210.11379",
    arxiv_date="2022-10-20", jrecv="2022-09-19", jacc="2022-10-11", jonline="2022-11-04",
    telescope="QUIJOTE (Yebes 40m)", date_conf="high",
    flags="HCCCHCCC: paper confirms singlet carbene ('1A ground state', terminal-carbon "
      "carbene) -- consistent with current SMILES, no correction needed. HCCCCS: paper confirms "
      "radical (2Pi ground state) but does not pin down exact radical-center atom vs current "
      "SMILES's terminal-CH placement -- flagged unverified, not contradicted"),
 "Agundez et al. 2023a": dict(
    authors="M. Agundez, C. Cabezas, N. Marcelino, R. Fuentetaja, B. Tercero, P. de Vicente, "
            "J. Cernicharo",
    journal="A&A 669, L1 (2023)", doi="10.1051/0004-6361/202245492", arxiv="2212.06903",
    arxiv_date="2022-12-13", jrecv="2022-11-17", jacc="2022-12-13", jonline="2022-12-22",
    telescope="QUIJOTE (Yebes 40m)", date_conf="high",
    flags="SMILES LIKELY WRONG: paper describes NC4NH+ as 'the protonated form of NC4N' "
      "(dicyanopolyyne) formed by proton transfer -- by direct analogy to its own explicitly "
      "named sibling NCCNH+ (N-protonated) this should be a linear, terminal-N-protonated "
      "cation, not the current SMILES's carbon-centered-charge structure C(#C[C+]=N)C#N. "
      "NOT corrected here (no direct quote naming the atom) -- recommend Reham review before "
      "using in the manuscript; candidate corrected form ~ N#CC#CC#C[NH+] (linear, N-protonated)."),
 "Cernicharo et al. 2023": dict(
    authors="J. Cernicharo, J. R. Pardo, C. Cabezas, M. Agundez, B. Tercero, N. Marcelino, "
            "R. Fuentetaja, M. Guelin, P. de Vicente",
    journal="A&A 670, L19 (2023)", doi="10.1051/0004-6361/202245816", arxiv="2301.12819",
    arxiv_date="", jrecv="2022-12-29", jacc="2023-01-30", jonline="2023-02-15",
    telescope="QUIJOTE (Yebes 40m)", date_conf="medium",
    flags="arXiv v1 exact day not locatable this session (ID prefix fixes Jan 2023); SMILES "
      "correction to anion form [C-]#CC#CC#CC#N CONFIRMED directly: paper states 'the radical "
      "can be discarded, as the observed lines belong to a singlet species' -- applied in Step 1"),
 "Fuentetaja et al. 2023": dict(
    authors="R. Fuentetaja, C. Bermudez, C. Cabezas, M. Agundez, B. Tercero, N. Marcelino, "
            "J. R. Pardo, L. Margules, R. A. Motiyenko, J.-C. Guillemin, P. de Vicente, "
            "J. Cernicharo",
    journal="A&A 671, L6 (2023)", doi="10.1051/0004-6361/202245732", arxiv="2302.08825",
    arxiv_date="2023-02-17", jrecv="2022-12-20", jacc="2023-02-02", jonline="2023-03-03",
    telescope="QUIJOTE (Yebes 40m)", date_conf="high",
    flags="closed-shell confirmed directly ('planar molecular skeleton, Cs frame')"),
 "Agundez et al. 2023b": dict(
    authors="M. Agundez, J.-C. Loison, K. M. Hickson, V. Wakelam, R. Fuentetaja, C. Cabezas, "
            "N. Marcelino, B. Tercero, P. de Vicente, J. Cernicharo",
    journal="A&A 673, A34 (2023)", doi="10.1051/0004-6361/202346076", arxiv="2303.16121",
    arxiv_date="2023-03-28", jrecv="2023-02-03", jacc="2023-03-27", jonline="",
    telescope="QUIJOTE (Yebes 40m)", date_conf="high",
    flags="published-online journal date not located (received/accepted dates and arXiv date "
      "are confirmed primary-source); all 3 species (ethanol, acetone, propanal) closed-shell"),
 "Cernicharo et al. 2022d": dict(
    authors="J. Cernicharo, R. Fuentetaja, C. Cabezas, M. Agundez, N. Marcelino, B. Tercero, "
            "J. R. Pardo, P. de Vicente",
    journal="A&A 663, L5 (2022)", doi="10.1051/0004-6361/202244255", arxiv="2207.03903",
    arxiv_date="2022-07-08", jrecv="2022-06-12", jacc="2022-07-07", jonline="2022-07-20",
    telescope="QUIJOTE (Yebes 40m)", date_conf="medium",
    flags="arXiv v1 date is a plausible secondary indicator (Semantic Scholar publicationDate), "
      "not read directly off arxiv.org; paper actually reports 5 species (title: 'five cyano "
      "derivatives of propene') -- allyl cyanide (CH2CHCH2CN) specifically has 2 detected "
      "conformers (gauche, cis) in this paper; our single generic SMILES C=CCC#N doesn't "
      "specify which -- flagged, not resolved"),
 "Cernicharo et al. 2021f": dict(
    authors="J. Cernicharo, C. Cabezas, M. Agundez, B. Tercero, J. R. Pardo, N. Marcelino, "
            "J. D. Gallego, F. Tercero, J. A. Lopez-Perez, P. de Vicente",
    journal="A&A 648, L3 (2021)", doi="10.1051/0004-6361/202140642", arxiv="2103.12431",
    arxiv_date="2021-03-23", jrecv="2021-02-23", jacc="2021-03-17", jonline="2021-04-07",
    telescope="Yebes 40m Q-band survey (pre-dates the 'QUIJOTE' acronym in this paper's own "
      "text, but is the same program later branded QUIJOTE)", date_conf="medium",
    flags="arXiv v1 date is circumstantial (astro-ph mailing-list archive), not primary-"
      "confirmed. Per-species: NCS confirmed radical (2Pi3/2, explicit term symbol). HCCS: "
      "paper never uses the word 'radical' or gives a term symbol for it -- chemically expected "
      "open-shell (odd valence-e-) but UNVERIFIED by this source. H2CCS, H2CCCS: confirmed "
      "closed-shell (symmetric ortho/para transition pattern, no radical language). C4S: paper "
      "explicitly calls it 'the C4S radical' but gives no singlet/triplet term symbol -- C4S has "
      "an even valence-electron count, so 'radical' language plus even-electron count suggests "
      "diradical/triplet character that the current single-radical-electron SMILES may not "
      "correctly represent -- UNRESOLVED. C5S: paper gives NO radical language and NO term "
      "symbol at all for C5S -- current SMILES's 2 radical electrons are UNVERIFIED by this "
      "source in either direction."),
 "Cernicharo et al. 2021b": dict(
    authors="J. Cernicharo, M. Agundez, C. Cabezas, B. Tercero, N. Marcelino, J. R. Pardo, "
            "P. de Vicente",
    journal="A&A 649, L15 (2021)", doi="10.1051/0004-6361/202141156", arxiv="2104.13991",
    arxiv_date="2021-04-28", jrecv="2021-04-22", jacc="2021-04-28", jonline="2021-05-13",
    telescope="QUIJOTE (Yebes 40m) + IRAM 30m", date_conf="high",
    flags="c-C5H6 confirmed ordinary closed-shell hydrocarbon, consistent with current SMILES. "
      "c-C3HCCH POSSIBLE MISCHARACTERIZATION: paper's stated formation route (CCH radical + "
      "closed-shell c-C3H2 carbene -> closed-shell product + H atom) implies c-C3HCCH should be "
      "closed-shell, but current SMILES C#CC1=C[C]1 encodes an open-shell/biradical ring carbon "
      "-- NOT corrected here, recommend Reham review / expert cheminformatics check."),
 "Cernicharo et al. 2021e": dict(
    authors="J. Cernicharo, M. Agundez, R. I. Kaiser, C. Cabezas, B. Tercero, N. Marcelino, "
            "J. R. Pardo, P. de Vicente",
    journal="A&A 655, L1 (2021)", doi="10.1051/0004-6361/202142226", arxiv="2110.09105",
    arxiv_date="2021-10-18", jrecv="2021-09-15", jacc="2021-10-14", jonline="2021-11-09",
    telescope="QUIJOTE (Yebes 40m)", date_conf="medium",
    flags="arXiv v1 date from third-party mirror, corroborated by the PDF's own manuscript "
      "compile-date footer (Oct 19, one day after) -- moderately confident but not a direct "
      "arxiv.org metadata read. Both isomers confirmed closed-shell (simple H-to-ethynyl ring "
      "substitutions on cyclopentadiene); exact regiochemistry (which ring carbon) not "
      "independently verified beyond the 1-/2- isomer naming."),
 "Cernicharo et al. 2021g": dict(
    authors="J. Cernicharo, C. Cabezas, Y. Endo, M. Agundez, B. Tercero, J. R. Pardo, "
            "N. Marcelino, P. de Vicente",
    journal="A&A 650, L14 (2021)", doi="10.1051/0004-6361/202141297", arxiv="2105.12996",
    arxiv_date="2021-05-27", jrecv="2021-05-12", jacc="2021-05-21", jonline="2021-06-23",
    telescope="QUIJOTE (Yebes 40m)", date_conf="high",
    flags="both species (HCSCN 'cyano thioformaldehyde', HCSCCH 'propynethial') confirmed "
      "closed-shell, no radical/ionic language in paper; E/Z geometry not discussed"),
 "Cabezas et al. 2021a": dict(
    authors="C. Cabezas, M. Agundez, N. Marcelino, B. Tercero, J. R. Pardo, P. de Vicente, "
            "J. Cernicharo",
    journal="A&A 654, L9 (2021)", doi="10.1051/0004-6361/202142156", arxiv="2110.03779",
    arxiv_date="2021-10-07", jrecv="2021-09-04", jacc="2021-10-06", jonline="2021-10-19",
    telescope="QUIJOTE (Yebes 40m)", date_conf="high",
    flags="radical confirmed (2B1 doublet ground state), consistent with current SMILES, no "
      "correction needed"),
 "Cabezas et al. 2021c": dict(
    authors="C. Cabezas, B. Tercero, M. Agundez, N. Marcelino, J. R. Pardo, P. de Vicente, "
            "J. Cernicharo",
    journal="A&A 650, L9 (2021)", doi="10.1051/0004-6361/202141274", arxiv="2106.00635",
    arxiv_date="2021-06 (exact day not located)", jrecv="2021-05-08", jacc="2021-05-31",
    jonline="2021-06-15", telescope="QUIJOTE (Yebes 40m)", date_conf="medium",
    flags="exact arXiv v1 day not extractable this session (month/ID confirmed); cumulene "
      "carbene, 'two non-bonded electrons on terminal C' language conventionally implies "
      "closed-shell singlet but no explicit term symbol given -- current SMILES not "
      "contradicted, not fully confirmed either"),
 "Cabezas et al. 2022b": dict(
    authors="C. Cabezas, M. Agundez, N. Marcelino, B. Tercero, Y. Endo, R. Fuentetaja, "
            "J. R. Pardo, P. de Vicente, J. Cernicharo",
    journal="A&A 657, L4 (2022)", doi="10.1051/0004-6361/202142815", arxiv="2112.11855",
    arxiv_date="2021-12-22", jrecv="2021-12-02", jacc="2021-12-21", jonline="2022-01-07",
    telescope="QUIJOTE (Yebes 40m)", date_conf="high",
    flags="SMILES LIKELY WRONG: paper explicitly gives ground state as X~3Sigma- (a LINEAR "
      "TRIPLET/diradical cation, analogous to neutral CCS's own triplet ground state), not the "
      "closed-shell singlet implied by current SMILES C#C[S+]. NOT corrected here (paper gives "
      "no bond-length/charge-localization detail beyond the term symbol) -- recommend Reham "
      "review / expert cheminformatics check before using in the manuscript."),
 "Cabezas et al. 2022c": dict(
    authors="C. Cabezas, M. Agundez, N. Marcelino, B. Tercero, R. Fuentetaja, P. de Vicente, "
            "J. Cernicharo",
    journal="A&A 659, L8 (2022)", doi="10.1051/0004-6361/202142972", arxiv="2203.05617",
    arxiv_date="2022-03 (exact day not located)", jrecv="2021-12-21", jacc="2022-03-06",
    jonline="2022-03-23", telescope="QUIJOTE (Yebes 40m)", date_conf="medium",
    flags="exact arXiv v1 day not extractable this session (month/ID confirmed); 1Sigma "
      "closed-shell ground state confirmed, consistent with current SMILES, no correction "
      "needed (N-protonation site inferred from naming/analogy, not an explicit quote)"),
 "Cabezas et al. 2023": dict(
    authors="C. Cabezas, J. Tang, M. Agundez, K. Seiki, Y. Sumiyoshi, Y. Ohshima, B. Tercero, "
            "N. Marcelino, R. Fuentetaja, P. de Vicente, Y. Endo, J. Cernicharo",
    journal="A&A 676, L5 (2023)", doi="10.1051/0004-6361/202347385", arxiv="2307.09127",
    arxiv_date="2023-07-18", jrecv="2023-07-06", jacc="2023-07-18", jonline="2023-07-31",
    telescope="QUIJOTE (Yebes 40m)", date_conf="high",
    flags="H2CCCN SMILES correction CONFIRMED directly: paper states alpha-CVR has a bent "
      "(147 deg), '2A prime' doublet ground state with the radical on the central carbon "
      "(alpha-H removed) -- applied in Step 1 as C=[C]C#N"),
 "Cabezas et al. 2022a": dict(
    authors="C. Cabezas et al.", journal="A&A (2022)", doi="", arxiv="2206.00955",
    arxiv_date="", jrecv="", jacc="", jonline="", telescope="QUIJOTE (Yebes 40m)",
    date_conf="low", flags="resolved in a prior session (SC item 17); full journal "
      "volume/page not independently re-confirmed there or this session"),
 "Silva et al. 2023": dict(
    authors="W. G. D. P. Silva, J. Cernicharo, S. Schlemmer, et al.",
    journal="A&A 676, L1 (2023)", doi="10.1051/0004-6361/202347174", arxiv="",
    arxiv_date="", jrecv="", jacc="", jonline="2023 (year only, from citation label; "
      "exact date not located, no journal dates or arXiv id found)",
    telescope="QUIJOTE (Yebes 40m)", date_conf="unresolved",
    flags="resolved in a prior session (SC item 27); arXiv id not captured there and not "
      "independently re-searched this session; no exact date located at all -- half-year "
      "bucket for H2CCCH+ is UNCERTAIN (could be either half of 2023)"),
 "Agundez et al. 2022a": dict(
    authors="M. Agundez, C. Cabezas, N. Marcelino, R. Fuentetaja, B. Tercero, P. de Vicente, "
            "J. Cernicharo",
    journal="A&A 659, L9 (2022)", doi="10.1051/0004-6361/202243396", arxiv="2203.06999",
    arxiv_date="2022-03-14", jrecv="2022-02-22", jacc="2022-03-11", jonline="2022-03-24",
    telescope="QUIJOTE (Yebes 40m) + IRAM 30m", date_conf="high",
    flags="linear cation structure confirmed, consistent with current SMILES, no correction "
      "needed"),
 "Agundez et al. 2022b": dict(
    authors="M. Agundez, N. Marcelino, C. Cabezas, R. Fuentetaja, B. Tercero, P. de Vicente, "
            "J. Cernicharo",
    journal="A&A 657, A96 (2022)", doi="10.1051/0004-6361/202142678", arxiv="2112.05370",
    arxiv_date="2021-12-10", jrecv="2021-11-16", jacc="2021-12-09", jonline="2022-01-17",
    telescope="IRAM 30m (3mm) + Yebes 40m (QUIJOTE, 8mm)", date_conf="high",
    flags="ARXIV ID CORRECTED this session: prior record risked conflating this paper with the "
      "companion discovery paper Agundez et al. 2021 (A&A 647 L10, arXiv:2103.03807, the "
      "original 8mm detection) -- 2112.05370 is the correct id for the '2022b' (3mm "
      "confirmation) paper cited here. Radical (doublet, two equivalent H nuclei) confirmed; "
      "exact resonance-form (propargyl vs allenyl radical) not settled by this paper's text, "
      "current SMILES C=C=[CH] (allenyl form) is a standard, reasonable representation"),
 "Cernicharo et al. 2020b": dict(
    authors="J. Cernicharo, N. Marcelino, J. R. Pardo, M. Agundez, B. Tercero, "
            "P. de Vicente, C. Cabezas, C. Bermudez",
    journal="A&A 641, L9 (2020)", doi="10.1051/0004-6361/202039231", arxiv="2009.03240",
    arxiv_date="2020-09-07", jrecv="2020-08-21", jacc="2020-08-30", jonline="2020-09-23",
    telescope="QUIJOTE (Yebes 40m) + IRAM 30m", date_conf="high",
    flags="C3N-/C5N- both confirmed closed-shell singlet anions (electron attachment to "
      "open-shell doublet neutral radicals C3N/C5N yields closed-shell singlet anions) -- "
      "consistent with current SMILES, no correction needed. Paper's own anion-specific "
      "electronic-state language not directly quoted (inferred from attachment mechanism)."),
 "Marcelino et al. 2023": dict(
    authors="N. Marcelino, C. Puzzarini, M. Agundez, R. Fuentetaja, B. Tercero, "
            "P. de Vicente, J. Cernicharo",
    journal="A&A 674, L13 (2023)", doi="10.1051/0004-6361/202346935", arxiv="",
    arxiv_date="", jrecv="2023-05-18", jacc="2023-06-05", jonline="",
    telescope="QUIJOTE (Yebes 40m)", date_conf="unresolved",
    flags="No arXiv preprint locatable for this paper -- checked independently twice this "
      "session (direct search + Semantic Scholar metadata, which shows no arXiv field at all "
      "for this DOI); genuinely appears A&A-only. Using journal 'accepted' date as the best "
      "available first-public-date proxy (no 'published online' date located); well within the "
      "Aug2021-Sep2023 core window regardless, so half-year bucket assignment is not in doubt "
      "even though the exact day is not pinned down. Species electronic structure: paper gives "
      "only the term symbol 2A'' with no statement of which atom (terminal O per our S[O] vs. "
      "some other placement) carries the unpaired electron -- SMILES not independently "
      "confirmed or contradicted."),
 "Fatima et al. 2023": dict(
    authors="M. Fatima, H. S. P. Muller, O. Zingsheim, F. Lewen, V. M. Rivilla, "
            "I. Jimenez-Serra, J. Martin-Pintado, S. Schlemmer",
    journal="A&A 680, A25 (2023)", doi="10.1051/0004-6361/202347112", arxiv="2309.17236",
    arxiv_date="2023-09-29", jrecv="2023-06-05", jacc="2023-09-14", jonline="2023-12-04",
    telescope="G+0.693 (Galactic Center cloud) -- NOT TMC-1", date_conf="unresolved",
    flags="CRITICAL: this paper does NOT report a TMC-1 detection of C4H8. Verified directly "
      "from its own text: isobutene was detected only in G+0.693; TMC-1 is mentioned exactly "
      "3 times, all as a suggested FUTURE target, not a result ('isobutene is a promising "
      "candidate to be searched for in space'). Shay et al. (2025a) Table C1 appears to have "
      "mis-cited this row. NO valid TMC-1 source for C4H8 was found in this pass. This row's "
      "first-public-date/telescope are left UNRESOLVED pending Reham's decision: drop C4H8 "
      "entirely, or find the correct citation."),
}

# rows whose 'reference' field in v2 already carries enough information -- parsed directly rather than
# re-listed above.
ALREADY_INLINE = {
    "CH3CHS": dict(authors="M. Agundez et al.", journal="A&A (2025)", doi="",
        arxiv="2501.05125", arxiv_date="", jrecv="", jacc="", jonline="",
        telescope="QUIJOTE (presumed)", date_conf="low",
        flags="arXiv id from project's existing record; full bib not independently "
              "re-verified this session"),
    "1-cyanoacenaphthylene": dict(authors="J. Cernicharo et al.", journal="A&A (2024)",
        doi="10.1051/0004-6361/202452196", arxiv="", arxiv_date="", jrecv="", jacc="",
        jonline="2024 (year only, from citation label; exact date not located)",
        telescope="QUIJOTE (presumed)", date_conf="unresolved",
        flags="DOI from project's existing record; arXiv id/full bib not independently "
              "re-verified this session; no exact date located"),
    "5-cyanoacenaphthylene": dict(authors="J. Cernicharo et al.", journal="A&A (2024)",
        doi="10.1051/0004-6361/202452196", arxiv="", arxiv_date="", jrecv="", jacc="",
        jonline="2024 (year only, from citation label; exact date not located)",
        telescope="QUIJOTE (presumed)", date_conf="unresolved",
        flags="same paper as 1-cyanoacenaphthylene"),
    "3-cyanoacenaphthylene": dict(authors="J. Cernicharo et al.", journal="A&A (2026)",
        doi="10.1051/0004-6361/202557893", arxiv="", arxiv_date="", jrecv="", jacc="",
        jonline="2026 (year only, from citation label; exact date not located)",
        telescope="QUIJOTE (presumed)", date_conf="unresolved",
        flags="follow-up paper completing the isomer set; DOI from project's existing "
              "record, arXiv id not captured, no exact date located"),
    "4-cyanoacenaphthylene": dict(authors="J. Cernicharo et al.", journal="A&A (2026)",
        doi="10.1051/0004-6361/202557893", arxiv="", arxiv_date="", jrecv="", jacc="",
        jonline="2026 (year only, from citation label; exact date not located)",
        telescope="QUIJOTE (presumed)", date_conf="unresolved",
        flags="same paper as 3-cyanoacenaphthylene"),
    "1-cyanopyrene": dict(authors="G. Wenzel et al.", journal="Nature Astronomy (2024)",
        doi="", arxiv="2410.00670", arxiv_date="", jrecv="", jacc="", jonline="",
        telescope="GOTHAM (GBT)", date_conf="low",
        flags="2:1:2 abundance ratio reported across the 3 cyanopyrene isomers in one paper"),
    "2-cyanopyrene": dict(authors="G. Wenzel et al.", journal="Nature Astronomy (2024)",
        doi="", arxiv="2410.00670", arxiv_date="", jrecv="", jacc="", jonline="",
        telescope="GOTHAM (GBT)", date_conf="low", flags="same paper as 1-cyanopyrene"),
    "4-cyanopyrene": dict(authors="G. Wenzel et al.", journal="Nature Astronomy (2024)",
        doi="", arxiv="2410.00670", arxiv_date="", jrecv="", jacc="", jonline="",
        telescope="GOTHAM (GBT)", date_conf="low", flags="same paper as 1-cyanopyrene"),
    "cyanocoronene (C24H11CN)": dict(
        authors="G. Wenzel, S. Gong, C. Xue, P. B. Changala, M. S. Holdren, T. H. Speak, "
                "D. A. Stewart, Z. T. P. Fried, R. H. J. Willis, E. A. Bergin, "
                "A. M. Burkhardt, A. N. Byrne, S. B. Charnley, A. Lipnicky, R. A. Loomis, "
                "C. N. Shingledecker, I. R. Cooke, A. J. Remijan, M. C. McCarthy, "
                "A. E. Wendlandt, B. A. McGuire",
        journal="ApJL 984, L36 (2025)", doi="10.3847/2041-8213/adc911", arxiv="2504.05232",
        arxiv_date="", jrecv="", jacc="", jonline="", telescope="GOTHAM (GBT)",
        date_conf="low", flags="fully resolved in a prior session (SC item 23, ADS bicode "
              "2025ApJ...984L..36W confirmed); exact arXiv v1 day not independently re-pulled "
              "this session"),
    "phenalene (c-C13H10)": dict(
        authors="C. Cabezas, M. Agundez, C. Perez, D. Villar-Castro, G. Molpeceres, "
                "D. Perez, A. L. Steber, R. Fuentetaja, B. Tercero, N. Marcelino, "
                "A. Lesarri, P. de Vicente, J. Cernicharo",
        journal="A&A 701, L8 (2025)", doi="10.1051/0004-6361/202556687", arxiv="2508.13857",
        arxiv_date="2025-08-19", jrecv="", jacc="", jonline="", telescope="QUIJOTE (Yebes 40m)",
        date_conf="high", flags="author list + telescope resolved this session (previously "
              "flagged 'authors not independently re-extracted'); note this is QUIJOTE, not "
              "GOTHAM -- double-check if previously filed otherwise"),
    "1-cyano propargyl radical (HCCCHCN)": dict(
        authors="C. Cabezas, M. Agundez, N. Marcelino, C. H. Chang, R. Fuentetaja, "
                "B. Tercero, M. Nakajima, Y. Endo, P. de Vicente, J. Cernicharo",
        journal="A&A 693, L14 (2025)", doi="10.1051/0004-6361/202453419", arxiv="2501.01938",
        arxiv_date="2025-01-03", jrecv="", jacc="", jonline="", telescope="QUIJOTE (Yebes 40m)",
        date_conf="high", flags="author list + telescope resolved this session (previously "
              "flagged 'authors not independently re-extracted'); note this is QUIJOTE, not "
              "GOTHAM -- double-check if previously filed otherwise"),
    # 14 census-v2 additions (13 kept after the l-H2C3O drop) -- the reference field in v2 already carries
    # most of what is needed; carried over as-is (all originally primary-source-verified).
    "C5H4": dict(authors="Fuentetaja et al. (from prior session's census)",
        journal="A&A (2024)", doi="", arxiv="2407.07466", arxiv_date="", jrecv="", jacc="",
        jonline="", telescope="QUIJOTE (presumed)", date_conf="low",
        flags="carried over from prior session's completeness-census pass, not re-verified "
              "this session"),
    "CH2CHCHS": dict(authors="(from prior session's census)",
        journal="A&A, June 2025 (aa54670-25)", doi="", arxiv="", arxiv_date="", jrecv="",
        jacc="", jonline="2025-06 (month only, from journal issue code aa54670-25)",
        telescope="QUIJOTE (presumed)", date_conf="low",
        flags="'Discovery of propenethial in TMC-1'; carried over from prior session, arXiv "
              "id not captured, not re-verified this session"),
    "HC5N+": dict(authors="(from prior session's census)",
        journal="A&A, June 2024 (aa50801-24)", doi="", arxiv="2406.06204", arxiv_date="",
        jrecv="", jacc="", jonline="", telescope="QUIJOTE", date_conf="low",
        flags="2Pi doublet ground state confirmed in prior session's HC3N+/HC5N+/HC7N+ "
              "verification pass; SMILES correction applied in Step 1"),
    "HC7N+": dict(authors="(from prior session's census)",
        journal="A&A, June 2024 (aa50801-24)", doi="", arxiv="2406.06204", arxiv_date="",
        jrecv="", jacc="", jonline="", telescope="QUIJOTE", date_conf="low",
        flags="same paper as HC5N+; SMILES correction applied in Step 1"),
    "HC3N+": dict(authors="C. Cabezas et al. (arXiv:2407.02121)",
        journal="A&A, July 2024 (aa51081-24)", doi="", arxiv="2407.02121", arxiv_date="",
        jrecv="", jacc="", jonline="", telescope="QUIJOTE", date_conf="low",
        flags="2Pi doublet, TMC-1 unambiguously confirmed in prior session ('We report the "
              "discovery of the radical cation HC3N+ in TMC-1'); SMILES correction applied "
              "in Step 1"),
    "HNC5": dict(authors="Fuentetaja et al. (from prior session's census)",
        journal="A&A, Aug 2024 (aa51319-24)", doi="", arxiv="", arxiv_date="", jrecv="",
        jacc="", jonline="2024-08 (month only, from journal issue code aa51319-24)",
        telescope="QUIJOTE (presumed)", date_conf="low",
        flags="carried over from prior session's completeness-census pass, not re-verified "
              "this session"),
    "c-C3H2S": dict(authors="A. Remijan et al.", journal="ApJ 982, 191 (2025)",
        doi="10.3847/1538-4357/adb84c", arxiv="2501.06343", arxiv_date="", jrecv="",
        jacc="", jonline="", telescope="GOTHAM", date_conf="medium",
        flags="quote-verified in prior session: 'we report the first interstellar detection "
              "of c-C3H2S toward TMC-1'"),
    "NC3S": dict(authors="J. Cernicharo et al.", journal="A&A 688, L13 (2024)",
        doi="10.1051/0004-6361/202451256", arxiv="2407.15275", arxiv_date="", jrecv="",
        jacc="", jonline="", telescope="QUIJOTE", date_conf="medium",
        flags="quote-verified in prior session: 'We present the detection of the free "
              "radicals NC3S and HC3S towards TMC-1 with the QUIJOTE line survey'"),
    "HC3S": dict(authors="J. Cernicharo et al.", journal="A&A 688, L13 (2024)",
        doi="10.1051/0004-6361/202451256", arxiv="2407.15275", arxiv_date="", jrecv="",
        jacc="", jonline="", telescope="QUIJOTE", date_conf="medium",
        flags="same paper as NC3S"),
    "NCCHCS": dict(authors="C. Cabezas et al.", journal="A&A 686, A121 (2024)",
        doi="10.1051/0004-6361/202450323", arxiv="2405.02033", arxiv_date="", jrecv="",
        jacc="", jonline="", telescope="QUIJOTE", date_conf="medium",
        flags="quote-verified in prior session ('Laboratory and astronomical discovery of "
              "cyanothioketene, NCCHCS, in the cold starless core TMC-1')"),
    "CH2(CN)2": dict(authors="M. Agundez et al.", journal="A&A 688, L31 (2024)",
        doi="10.1051/0004-6361/202451525", arxiv="2408.02843", arxiv_date="", jrecv="",
        jacc="", jonline="", telescope="QUIJOTE", date_conf="medium",
        flags="quote-verified in prior session ('rich interstellar reservoir of dinitriles: "
              "malononitrile and maleonitrile in TMC-1')"),
    "NC-CHCH-CN(Z)": dict(authors="M. Agundez et al.", journal="A&A 688, L31 (2024)",
        doi="10.1051/0004-6361/202451525", arxiv="2408.02843", arxiv_date="", jrecv="",
        jacc="", jonline="", telescope="QUIJOTE", date_conf="medium",
        flags="same paper as CH2(CN)2 (maleonitrile, the Z/cis dinitrile isomer)"),
    "HCNS": dict(authors="J. Cernicharo et al.", journal="A&A 682, L4 (2024)",
        doi="10.1051/0004-6361/202349105", arxiv="2401.11785", arxiv_date="", jrecv="",
        jacc="", jonline="", telescope="QUIJOTE", date_conf="medium",
        flags="quote-verified in prior session ('first detection of HCNS (thiofulminic acid) "
              "in space with the QUIJOTE line survey')"),
}

CIT.update(ALREADY_INLINE)

def lookup(ref):
    for k, v in CIT.items():
        if k in str(ref):
            return v
    return None

rows = []
for _, r in df.iterrows():
    f = r["formula"]
    cit = CIT.get(f) or lookup(r["reference"])
    if cit is None:
        cit = dict(authors="", journal="", doi="", arxiv="", arxiv_date="", jrecv="",
                    jacc="", jonline="", telescope="", date_conf="unresolved",
                    flags="NOT MATCHED to a citation-registry entry -- needs manual check")
    arxiv_id_month = ""
    if cit["arxiv"] and not cit["arxiv_date"]:
        aid = cit["arxiv"].split(".")[0]
        if len(aid) == 4 and aid.isdigit():
            yy, mm = aid[:2], aid[2:]
            arxiv_id_month = f"20{yy}-{mm} (month only, from arXiv ID prefix)"
    fpd = (cit["arxiv_date"] or arxiv_id_month or cit["jonline"] or cit["jacc"]
           or cit["jrecv"])
    fpd_basis = ("arxiv_v1" if cit["arxiv_date"] else
                 "arxiv_id_month_only" if arxiv_id_month else
                 "journal_online" if cit["jonline"] else
                 "journal_accepted" if cit["jacc"] else
                 "journal_received" if cit["jrecv"] else "none")
    mol = Chem.MolFromSmiles(r["smiles"])
    if mol is not None:
        charge = Chem.GetFormalCharge(mol)
        nrad = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        klass = "ion" if charge != 0 else ("radical" if nrad > 0 else "closed_shell")
    else:
        charge, nrad, klass = None, None, "PARSE_FAILED"
    rows.append({
        "formula": f, "name": r.get("note", "") or f, "smiles": r["smiles"],
        "log10_column_density_cm2": r["log10_column_density_cm2"],
        "cd_source": r["reference"], "first_public_date": fpd,
        "first_public_date_basis": fpd_basis,
        "arxiv_v1_date": cit["arxiv_date"], "journal_received": cit["jrecv"],
        "journal_accepted": cit["jacc"], "journal_online": cit["jonline"],
        "arxiv_id": cit["arxiv"], "doi": cit["doi"], "journal_ref": cit["journal"],
        "authors": cit["authors"], "telescope_survey": cit["telescope"],
        "formal_charge": charge, "n_radical_electrons": nrad, "class": klass,
        "date_confidence": cit["date_conf"], "flags": cit["flags"],
        "detection_window_v2": r["detection_window"],
    })

out = pd.DataFrame(rows)
out.to_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3.csv", index=False)
print(f"v3 written: {len(out)} rows")
print(out["class"].value_counts())
print("\ndate_confidence breakdown:")
print(out["date_confidence"].value_counts())
print("\nrows with no first_public_date:")
print(out.loc[out["first_public_date"].isna() | (out["first_public_date"] == ""), "formula"].tolist())
