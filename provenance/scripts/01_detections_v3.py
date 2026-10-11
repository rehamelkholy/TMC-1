# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Build post2021_detections_v3.csv from v2 plus verified corrections and a primary-source literature pass.
Open questions are flagged in the 'flags' column rather than silently resolved.
"""
import pandas as pd
from rdkit import Chem

v2 = pd.read_csv("<PROJECT_ROOT>/validation_freeze/post2021_detections_v2.csv")

# ---------------------------------------------------------------------------
# STEP 1: drops, SMILES corrections, addition
# ---------------------------------------------------------------------------
DROP_FORMULAS = ["l-H2C3O", "HCCCO", "c-C9H8"]
# Each canonical SMILES is identical to a Table A2 (Lee 2021) pre-cutoff anchor in every case
# (l-H2C3O==H2C3O/Cernicharo2020, HCCCO==HC3O/Cernicharo2020, c-C9H8==C9H8/Burkhardt2021 [no "b"]).

SMILES_CORRECTIONS = {
    # formula: (new_smiles, reason)
    "H2CCCN": ("C=[C]C#N", "cyanovinyl (alpha-CVR) radical, not vinyl cyanide; "
               "confirmed this session against Cabezas et al. 2023 (A&A 676, L5): "
               "'2A prime ground state', central-carbon radical from alpha-H removal"),
    "C7N-": ("[C-]#CC#CC#CC#N", "anion (closed-shell singlet), not neutral radical; "
             "confirmed this session against Cernicharo et al. 2023 (A&A 670, L19): "
             "'the radical can be discarded, as the observed lines belong to a singlet "
             "species' -- neutral C7N is 2Pi radical, detected carrier is the anion"),
    "HC3N+": ("C#CC#[N+]", "acetylenic-skeleton cation, doublet (2Pi); pre-existing "
              "correction, re-confirmed via arXiv:2407.02121"),
    "HC5N+": ("C#CC#CC#[N+]", "acetylenic-skeleton cation, doublet (2Pi); pre-existing "
              "correction, re-confirmed via arXiv:2406.06204"),
    "HC7N+": ("C#CC#CC#CC#[N+]", "acetylenic-skeleton cation, doublet (2Pi); pre-existing "
              "correction, re-confirmed via arXiv:2406.06204"),
}

HSO_ROW = {
    "formula": "HSO", "smiles": "S[O]", "log10_column_density_cm2": 11.04,
    "reference": "Marcelino et al. 2023, A&A 674, L13, doi:10.1051/0004-6361/202346935",
    "detection_window": "2021-08_to_2023-09", "source": "gotham_quijote_census_v2 2026-09-19",
    "confidence": "high", "n_fragments": None, "n_oov": None, "coverage": None,
    "vectorize_ok": True,
    "note": "radical; QUIJOTE; absent from Shay et al. 2025a Table C1 despite falling "
            "inside its stated coverage window -- flagged, recommend include (see log).",
}

df = v2[~v2["formula"].isin(DROP_FORMULAS)].copy()
df = pd.concat([df, pd.DataFrame([HSO_ROW])], ignore_index=True)

change_log_rows = []
for f in DROP_FORMULAS:
    change_log_rows.append({
        "formula": f, "action": "DROP",
        "reason": "canonical SMILES identical to a Table A2 (Lee et al. 2021) "
                   "pre-cutoff anchor -- verified this session via RDKit canonical-SMILES "
                   "comparison against active_learning/candidate_pool_sources/table_a2_lee2021.py",
        "old_smiles": v2.loc[v2["formula"] == f, "smiles"].iloc[0], "new_smiles": "",
    })
for f, (new_smi, reason) in SMILES_CORRECTIONS.items():
    old_smi = v2.loc[v2["formula"] == f, "smiles"].iloc[0]
    df.loc[df["formula"] == f, "smiles"] = new_smi
    change_log_rows.append({
        "formula": f, "action": "CORRECT_SMILES", "reason": reason,
        "old_smiles": old_smi, "new_smiles": new_smi,
    })
change_log_rows.append({
    "formula": "HSO", "action": "ADD", "reason": HSO_ROW["note"],
    "old_smiles": "", "new_smiles": HSO_ROW["smiles"],
})

assert len(df) == 76, f"expected 76 rows, got {len(df)}"
df = df.reset_index(drop=True)

pd.DataFrame(change_log_rows).to_csv(
    "<PROJECT_ROOT>/f3_validation_freeze/v3_change_log.csv", index=False)

df.to_csv("<PROJECT_ROOT>/f3_validation_freeze/_v3_base_76.csv", index=False)
print(f"Base v3 table: {len(df)} rows. Change log: {len(change_log_rows)} entries.")
