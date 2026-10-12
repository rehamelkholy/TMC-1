# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Detections v3 (76 rows) -> v3.1; v3 is not overwritten.
DROP C4H8 (not toward TMC-1, per A1(i): Fatima et al. 2023 report only G+0.693, never TMC-1).
CORRECT NC4NH+ -> N#CC#CC#[NH+] (linear, N-protonated cation, matching the paper's own "protonated form
    of NC4N via proton transfer" framing and its sibling NCCNH+'s N-protonation).
CORRECT c-C3HCCH -> closed-shell (per its discovery paper's formation mechanism: CCH radical + closed-shell
    c-C3H2 carbene -> closed-shell product + H atom). New SMILES C#CC1=C=C1 follows the Table A2 convention
    for c-C3H2 itself (anchor SMILES "C1=C=C1", a 0-radical-electron ring form) with an ethynyl substitution at
    one ring position; verified via RDKit to give the correct formula (C5H2) and 0 radical electrons. Per
    amendment A1(iv): bond-order/spin representation is a documented convention choice, not a literal
    quantum-chemical claim -- connectivity, H-placement and charge (0) are what is matched.
LEAVE UNCHANGED (flagged, not corrected): HCCS+, C4S, C5S, HCCS, C5O, HCCCCS, CH2CHCH2CN, with notes carried over.
ADD CH3CH2CCH (qualifying event: Cernicharo et al. 2024, arXiv:2404.07585, "Confirmation of the detection of
    CH3CH2CCH in TMC-1"; CD (6.2+-0.2)e11 cm^-2, log10=11.79, DOI 10.1051/0004-6361/202449531, A&A 686 A139,
    QUIJOTE).
ADD c-C11H8 (Fuentetaja et al. 2026, arXiv:2601.13403; SMILES cross-validated by two methods -- OPSIN
    name-to-structure on the IUPAC name "1H-cyclopent[cd]indene" and the CACTUS resolver -- which agree on
    canonical SMILES C1=Cc2cccc3c2C1=CC3, formula C11H8; CD (6.0+-0.5)e12 = log10 12.78).
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

v3 = pd.read_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3.csv")
# drop v3's own Step-5 match/bucket columns -- stale for v3.1 (recomputed fresh below)
v3 = v3.drop(columns=["level_A", "level_B", "level_C", "in_lee_list", "half_year"],
             errors="ignore")

df = v3[v3["formula"] != "C4H8"].copy()
assert len(df) == 75

df.loc[df["formula"] == "NC4NH+", "smiles"] = "N#CC#CC#[NH+]"
df.loc[df["formula"] == "c-C3HCCH", "smiles"] = "C#CC1=C=C1"

CH3CH2CCH_ROW = {
    "formula": "CH3CH2CCH", "name": "ethyl acetylene (1-butyne); qualifying event = "
        "Cernicharo et al. 2024 CONFIRMATION (2021b was tentative-only)",
    "smiles": "CCC#C", "log10_column_density_cm2": 11.79,
    "cd_source": "Cernicharo et al. 2024, A&A 686, A139, arXiv:2404.07585 "
        "((6.2+-0.2)e11 cm^-2, Trot=9K)",
    "first_public_date": "2024-04 (month only, from arXiv ID prefix; exact "
        "day not located despite direct primary-source attempts)",
    "first_public_date_basis": "arxiv_id_month_only",
    "arxiv_v1_date": "", "journal_received": "2024-02-07",
    "journal_accepted": "2024-03-08", "journal_online": "2024-06-07",
    "arxiv_id": "2404.07585", "doi": "10.1051/0004-6361/202449531",
    "journal_ref": "A&A 686, A139 (2024)",
    "authors": "J. Cernicharo et al. (not independently re-extracted this session)",
    "telescope_survey": "QUIJOTE (Yebes 40m)",
    "date_confidence": "medium",
    "flags": "ADDED in F3b per amendments -- qualifying event is the 2024 CONFIRMATION "
        "paper, not the 2021b tentative report (direct quote confirmed: 'We confirm the "
        "presence of CH3CH2CCH in TMC-1, which was previously reported as tentative by "
        "our group'). Exact arXiv v1 day not located (arxiv.org metadata inaccessible to "
        "the fetch tool this session); month-level estimate used.",
    "detection_window_v2": "NEW (F3b addition)",
}

C11H8_ROW = {
    "formula": "c-C11H8", "name": "1H-cyclopent[cd]indene",
    "smiles": "C1=Cc2cccc3c2C1=CC3", "log10_column_density_cm2": 12.78,
    "cd_source": "Fuentetaja et al. 2026, arXiv:2601.13403 ((6.0+-0.5)e12 cm^-2, "
        "log10 form confirmed directly in the paper)",
    "first_public_date": "2026-01-19",
    "first_public_date_basis": "arxiv_v1",
    "arxiv_v1_date": "2026-01-19", "journal_received": "2025-12-04",
    "journal_accepted": "2025-12-31", "journal_online": "",
    "arxiv_id": "2601.13403", "doi": "10.1051/0004-6361/202558398",
    "journal_ref": "A&A, Feb 2026 issue",
    "authors": "R. Fuentetaja, C. Cabezas, M. Agundez, B. Tercero, N. Marcelino, "
        "P. de Vicente, J. Cernicharo",
    "telescope_survey": "QUIJOTE (Yebes 40m)",
    "date_confidence": "high",
    "flags": "ADDED in F3b per amendments. SMILES independently cross-validated this "
        "session via TWO independent methods: OPSIN name-to-structure resolution of the "
        "IUPAC name '1H-cyclopent[cd]indene' (-> C1C=C2C=3C(=CC=CC13)C=C2) and the prior "
        "session's CACTUS/NCI resolver (-> C1C=C2C=Cc3cccc1c23) -- both canonicalize to "
        "the IDENTICAL SMILES C1=Cc2cccc3c2C1=CC3 (RDKit-confirmed formula C11H8), a much "
        "higher-confidence structure determination than either method alone.",
    "detection_window_v2": "NEW (F3b addition)",
}

for r in (CH3CH2CCH_ROW, C11H8_ROW):
    mol = Chem.MolFromSmiles(r["smiles"])
    r["formal_charge"] = Chem.GetFormalCharge(mol)
    r["n_radical_electrons"] = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    r["class"] = ("ion" if r["formal_charge"] != 0 else
                  ("radical" if r["n_radical_electrons"] > 0 else "closed_shell"))

df = pd.concat([df, pd.DataFrame([CH3CH2CCH_ROW, C11H8_ROW])], ignore_index=True)
assert len(df) == 77, f"expected 77 rows, got {len(df)}"

# recompute class/charge/radelec for the two corrected rows to keep them consistent
for f in ["NC4NH+", "c-C3HCCH"]:
    smi = df.loc[df["formula"] == f, "smiles"].iloc[0]
    mol = Chem.MolFromSmiles(smi)
    charge = Chem.GetFormalCharge(mol)
    nrad = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    klass = "ion" if charge != 0 else ("radical" if nrad > 0 else "closed_shell")
    df.loc[df["formula"] == f, "formal_charge"] = charge
    df.loc[df["formula"] == f, "n_radical_electrons"] = nrad
    df.loc[df["formula"] == f, "class"] = klass

# add/extend flags for the two corrected + the 7 left-unchanged-with-notes rows
LEAVE_UNCHANGED_NOTE = ("F3b amendments left this unchanged (flag stands from F3): ")
for f in ["HCCS+", "C4S", "C5S", "HCCS", "C5O", "HCCCCS", "CH2CHCH2CN"]:
    idx = df.index[df["formula"] == f]
    if len(idx):
        cur = df.loc[idx[0], "flags"]
        df.loc[idx[0], "flags"] = (cur if isinstance(cur, str) else "") + \
            " [F3b: left unchanged per task instruction, flag carried over from F3.]"

df.loc[df["formula"] == "NC4NH+", "flags"] = (
    "F3b CORRECTED: SMILES changed to N#CC#CC#[NH+] (linear, N-protonated cation) "
    "per task instruction -- matches the discovery paper's 'protonated form of NC4N "
    "via proton transfer' framing and direct analogy to its explicitly named sibling "
    "NCCNH+ (N-protonated). Formula/charge verified: C4HN2+, no radical electrons.")
df.loc[df["formula"] == "c-C3HCCH", "flags"] = (
    "F3b CORRECTED: SMILES changed to C#CC1=C=C1 (closed-shell) per task instruction -- "
    "the discovery paper's own formation mechanism (CCH radical + closed-shell c-C3H2 "
    "carbene -> closed-shell product + H atom) implies a closed-shell product, "
    "inconsistent with the prior open-shell/biradical SMILES. New form verified via "
    "RDKit: formula C5H2 (matches), 0 radical electrons (was 2), matching this project's "
    "own Table A2 convention for the parent c-C3H2 ring ('C1=C=C1', 0-radical-electron "
    "closed-shell form). Bond-order representation is a documented convention choice per "
    "amendment A1(iv), not a literal quantum-chemical claim.")

df.to_csv("<PROJECT_ROOT>/f3_validation_freeze/post2021_detections_v3_1.csv", index=False)
print(f"v3.1 written: {len(df)} rows")
print(df[["formula","class"]].tail(10).to_string())
