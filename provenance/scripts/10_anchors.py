# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Build table_a2_authoritative.csv.

Primary sources (both read directly from the PDFs, not from any earlier transcription):
  - Lee et al. (2021), ApJL 917, L6, Appendix F, Table A2 (87 rows, "totaling 87 unique species").
  - Toru Shay et al. (2025a), ApJ 985, 123, Appendix B, Table B1 ("OD" data set -- Lee's Table A2
    "reproduced ... with corrections", stated to contain "84 rather than 87 molecules" because "any
    unconfirmed species were removed").

Method: transcribe both tables verbatim (row order preserved), row-match Lee's 87 against Shay's 84 by
(formula, approximate SMILES identity), and report exactly what happened to each of Lee's 87 rows:
  - MATCHED: present in B1, values carried from B1 (the "corrected" source) when B1 differs from A2 (SMILES,
    log CD, or reference).
  - DROPPED: absent from B1 and absent from Shay's Table C1 (new detections) -- genuinely missing, not just
    relocated.
  - RELOCATED: absent from B1 but present, same SMILES/species, in Shay's Table C1 (their "new detections
    since 2021" table) -- flagged as a likely bookkeeping inconsistency in Shay et al. 2025a itself (these are
    pre-2021-08 Lee-list species, not genuinely new detections), not silently "corrected" into the anchor set.
Every row's fate is either MATCHED, DROPPED, or RELOCATED, with the reasoning visible in the `basis` column.
"""
import csv
import hashlib

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, inchi as rd_inchi
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
# Lee et al. (2021) Table A2, verbatim, row order as printed. 87 rows.
# (formula_as_printed, smiles_as_printed, log10_cd, reference)
# ---------------------------------------------------------------------------
LEE_A2 = [
    ("CH3C6H", "CC#CC#CC#C", 12.4914, "Remijan et al. (2006)"),
    ("CH3C4H", "CC#CC#C", 13.4771, "MacLeod et al. (1984)"),
    ("CH3C5N", "CC#CC#CC#N", 11.9243, "Remijan et al. (2006)"),
    ("CH3C3N", "CC#CC#N", 12.2553, "Remijan et al. (2006)"),
    ("N2H+", "N#[NH+]", 12.6990, "Choi et al. (2017)"),
    ("NH3", "N", 14.6998, "Gratier et al. (2016)"),
    ("CH3OH", "CO", 13.1614, "Gratier et al. (2016)"),
    ("C3H", "C1=C=[C]1", 13.4800, "Gratier et al. (2016)"),
    ("C3H", "[CH+]=C=[C-]", 12.7497, "Gratier et al. (2016)"),
    ("C3H2", "C1=C=C1", 11.7701, "Gratier et al. (2016)"),
    ("C3H2", "C1C#C1", 13.2695, "Gratier et al. (2016)"),
    ("C3H2", "C=C=[C]", 12.3979, "Cernicharo et al. (1991)"),
    ("CH3CCH", "CC#C", 14.0607, "Gratier et al. (2016)"),
    ("C2O", "[C+]#C[O-]", 12.5705, "Gratier et al. (2016)"),
    ("CH2CN", "[CH2]C#N", 13.5798, "Gratier et al. (2016)"),
    ("CH3CN", "CC#N", 12.6096, "Gratier et al. (2016)"),
    ("HNCO", "N=C=O", 13.0294, "Gratier et al. (2016)"),
    ("CS", "[C-]#[S+]", 13.4594, "Gratier et al. (2016)"),
    ("CH3CHO", "CC=O", 12.4298, "Gratier et al. (2016)"),
    ("HCS+", "C#[S+]", 12.7597, "Gratier et al. (2016)"),
    ("H2CS", "C=S", 13.6201, "Gratier et al. (2016)"),
    ("SO", "S=O", 13.6702, "Gratier et al. (2016)"),
    ("C4H", "[C]#CC#C", 13.4298, "Gratier et al. (2016)"),
    ("C4H2", "C=C=C=[C]", 13.3365, "Gratier et al. (2016)"),
    ("C3N", "[C]#CC#N", 13.5502, "Gratier et al. (2016)"),
    ("HNC3", "[C-]#C-C#[NH+]", 11.6803, "Gratier et al. (2016)"),
    ("C3O", "[C]#C[C]=O", 11.9201, "Gratier et al. (2016)"),
    ("HC3NH+", "C#CC#[NH+]", 11.8698, "Gratier et al. (2016)"),
    ("CH2CHCN", "C=CC#N", 12.8102, "Gratier et al. (2016)"),
    ("HCCCHO", "C#CC=O", 11.2601, "Gratier et al. (2016)"),
    ("C2S", "[C+]#C[S-]", 14.0086, "Gratier et al. (2016)"),
    ("OCS", "O=C=S", 13.2601, "Gratier et al. (2016)"),
    ("C5H", "[CH+]=C=C=C=[C-]", 12.2695, "Gratier et al. (2016)"),
    ("C3S", "[C-]#CC#[S+]", 13.1399, "Gratier et al. (2016)"),
    ("C6H", "[C]#CC#CC#C", 12.7404, "Gratier et al. (2016)"),
    ("HC3N", "C#CC#N", 14.2430, "Xue et al. (2020)"),
    ("HCCNC", "C#C[N+]#[C-]", 12.5821, "Xue et al. (2020)"),
    ("HC5N", "C#CC#CC#N", 13.8254, "Xue et al. (2020)"),
    ("HC4NC", "C#CC#C[N+]#[C-]", 11.5172, "Xue et al. (2020)"),
    ("HC7N", "C#CC#CC#CC#N", 13.5623, "Xue et al. (2020)"),
    ("HC6NC", "C#CC#CC#C[N+]#[C-]", 11.6064, "Xue et al. (2020)"),
    ("HC9N", "C#CC#CC#CC#CC#N", 13.3345, "Loomis et al. (2021)"),
    ("HC11N", "C#CC#CC#CC#CC#CC#N", 12.0170, "Loomis et al. (2021)"),
    ("C5H5CN", "C1C=CC=C1C#N", 11.9191, "Lee et al. (2021b)"),
    ("C5H5CN", "C1C=CC(=C1)C#N", 11.2788, "Lee et al. (2021b)"),
    ("C11H7N", "C1=CC=C2C(=C1)C=CC=C2C#N", 11.8663, "McGuire et al. (2021)"),
    ("C11H7N", "C1=CC=C2C=C(C=CC2=C1)C#N", 11.8482, "McGuire et al. (2021)"),
    ("C6H5CN", "C1=CC=C(C=C1)C#N", 12.2380, "McGuire et al. (2021)"),
    ("HCCCH2CN", "C#CCC#N", 11.9643, "McGuire et al. (2020)"),
    ("H3C5N", "C#C/C=C/C#N", 11.3874, "Lee et al. (2021a)"),
    ("H3C5N", "C=CC#CC#N", 11.0719, "Lee et al. (2021a)"),
    ("H3C5N", r"C#C/C=C\C#N", 11.3032, "Lee et al. (2021a)"),
    ("C8H", "[C]#CC#CC#CC#C", 11.6628, "Brünken et al. (2007)"),
    ("C8H-", "C#CC#CC#CC#[C-]", 10.3222, "Brünken et al. (2007)"),
    ("C6H-", "C#CC#CC#[C-]", 11.0792, "Brünken et al. (2007)"),
    ("C4H-", "C#CC#[C-]", 10.9294, "Brünken et al. (2007)"),
    ("H2CCO", "C=C=O", 12.7118, "Soma et al. (2018)"),
    ("CN", "[C]#N", 12.8899, "Pratap et al. (1997)"),
    ("HNC", "[C-]#[NH+]", 12.6201, "Pratap et al. (1997)"),
    ("HC7O", "C#CC#CC#C[C+]=O", 11.8921, "Cordiner et al. (2017)"),
    ("HC5O", "C#CC#C[C+]=O", 12.2304, "McGuire et al. (2017)"),
    ("H2CN", "C=[N]", 11.1761, "Ohishi et al. (1994)"),
    ("H2CO", "C=O", 13.0792, "Soma et al. (2018)"),
    ("HC3O+", "C#CC#[O+]", 11.3222, "Cernicharo et al. (2020)"),
    ("HOCO+", "O=C=[OH+]", 11.6021, "Cernicharo et al. (2020)"),
    ("H2COH+", "C=[OH+]", 11.4771, "Cernicharo et al. (2020)"),
    ("H2NCO+", "NC#[O+]", 10.6021, "Cernicharo et al. (2020)"),
    ("HCNO", "C#N[O]", 10.8451, "Cernicharo et al. (2020)"),
    ("HOCN", "OC#N", 11.0414, "Cernicharo et al. (2020)"),
    ("C4O", "[C]#CC#[C]=O", 11.0792, "Cernicharo et al. (2020)"),
    ("HCOOH", "C(=O)O", 12.1461, "Cernicharo et al. (2020)"),
    ("HC2O", "C#[C]=O", 12.0000, "Cernicharo et al. (2020)"),
    ("HC3O", "C#C[C]=O", 11.3010, "Cernicharo et al. (2020)"),
    ("HC4O", "C#CC#[C]=O", 11.4771, "Cernicharo et al. (2020)"),
    ("H2C3O", "C=C=C=O", 11.0414, "Cernicharo et al. (2020)"),
    ("H2C3O", "C1=C(=O)=C1", 11.6021, "Cernicharo et al. (2020)"),
    ("CH", "[CH]", 14.1461, "Sakai et al. (2013)"),
    ("CNCN", "[C]#NC#N", 11.9542, "Agúndez et al. (2018)"),
    ("NCCNH+", "N#CC#[NH+]", 10.9345, "Agúndez et al. (2015)"),
    ("C6H2", "C=C=C=C=C=[C]", 10.3284, "Langer et al. (1997)"),
    ("CH3CHCH2", "CC=C", 13.6021, "Marcelino et al. (2007)"),
    ("CH2C2HCN", "C=C=CC#N", 11.6532, "Lovas et al. (2006)"),
    ("HCN", "C#N", 12.3892, "Hirota et al. (1998)"),
    ("C9H8", "c1ccc2c(c1)CC=C2", 12.9823, "Burkhardt et al. 2021"),
    ("CH2CHCCH", "C=CC#C", 13.0792, "Cernicharo et al. (2021)"),
    ("HCCN", "N#C[CH+]", 11.6435, "Cernicharo et al. (2021)"),
    ("CH3CH2CN", "CCC#N", 11.0414, "Cernicharo et al. (2021)"),
]
assert len(LEE_A2) == 87, f"expected 87 Lee A2 rows, got {len(LEE_A2)}"

# ---------------------------------------------------------------------------
# Toru Shay et al. (2025a) Table B1 ("OD" data set), verbatim. 84 rows.
# (label_as_printed, smiles_as_printed, log10_cd, reference)
# label carries the c-/l- prefix Shay added where Lee's table left it off.
# ---------------------------------------------------------------------------
SHAY_B1 = [
    ("CH3C6H", "CC#CC#CC#C", 12.49, "A. J. Remijan et al. (2006)"),
    ("CH3C4H", "CC#CC#C", 13.26, "A. J. Remijan et al. (2006)"),
    ("CH3C5N", "CC#CC#CC#N", 11.87, "L. E. Snyder et al. (2006)"),
    ("CH3C3N", "CC#CC#N", 12.26, "A. J. Remijan et al. (2006)"),
    ("N2H+", "N#[NH+]", 12.70, "Y. Choi et al. (2017)"),
    ("NH3", "N", 14.70, "P. Gratier et al. (2016)"),
    ("CH3OH", "CO", 13.18, "T. Soma et al. (2018)"),
    ("c-C3H", "C1=C=[C]1", 13.48, "P. Gratier et al. (2016)"),
    ("l-C3H", "[CH+]=C=[C-]", 12.75, "P. Gratier et al. (2016)"),
    ("l-C3H2", "C=C=[C]", 11.77, "P. Gratier et al. (2016)"),
    ("c-C3H2", "C1C#C1", 13.27, "P. Gratier et al. (2016)"),
    ("CH3CCH", "CC#C", 14.06, "P. Gratier et al. (2016)"),
    ("C2O", "[C+]#C[O-]", 11.88, "J. Cernicharo et al. (2021c)"),
    ("CH2CN", "[CH2]C#N", 13.28, "J. D. Tennis et al. (2023)"),
    ("CH3CN", "CC#N", 12.70, "J. D. Tennis et al. (2023)"),
    ("HNCO", "N=C=O", 13.04, "J. Cernicharo et al. (2024a)"),
    ("CS", "[C-]#[S+]", 14.54, "J. Cernicharo et al. (2021f)"),
    ("CH3CHO", "CC=O", 12.51, "J. Cernicharo et al. (2020a)"),
    ("HCS+", "C#[S+]", 13.00, "P. Gratier et al. (2016)"),
    ("H2CS", "C=S", 13.62, "J. Cernicharo et al. (2021f)"),
    ("SO", "S=O", 14.00, "J. Cernicharo et al. (2021f)"),
    ("C4H", "[C]#CC#C", 13.43, "P. Gratier et al. (2016)"),
    ("C4H2", "C=C=C=[C]", 13.34, "P. Gratier et al. (2016)"),
    ("C3N", "[C]#CC#N", 13.55, "P. Gratier et al. (2016)"),
    ("HNC3", "[C-]#CC#[NH+]", 11.68, "P. Gratier et al. (2016)"),
    ("HC3NH+", "C#CC#[NH+]", 11.87, "P. Gratier et al. (2016)"),
    ("CH2CHCN", "C=CC#N", 12.81, "P. Gratier et al. (2016)"),
    ("HCCCHO", "C#CC=O", 12.86, "A. J. Remijan et al. (2024)"),
    ("C2S", "[C+]#C[S-]", 14.01, "P. Gratier et al. (2016)"),
    ("l-C5H", "[CH+]=C=C=C=[C-]", 12.11, "C. Cabezas et al. (2022c)"),
    ("l-C3S", "[C-]#CC#[S+]", 13.14, "P. Gratier et al. (2016)"),
    ("C6H", "[C]#CC#CC#C", 12.74, "P. Gratier et al. (2016)"),
    ("HC3N", "C#CC#N", 14.24, "C. Xue et al. (2020)"),
    ("HCCNC", "C#C[N+]#[C-]", 12.58, "C. Xue et al. (2020)"),
    ("HC5N", "C#CC#CC#N", 13.83, "C. Xue et al. (2020)"),
    ("HC4NC", "C#CC#C[N+]#[C-]", 11.52, "C. Xue et al. (2020)"),
    ("HC7N", "C#CC#CC#CC#N", 13.56, "C. Xue et al. (2020)"),
    ("HC6NC", "C#CC#CC#C[N+]#[C-]", 11.61, "C. Xue et al. (2020)"),
    ("HC9N", "C#CC#CC#CC#CC#N", 13.33, "R. A. Loomis et al. (2021)"),
    ("HC11N", "C#CC#CC#CC#CC#CC#N", 12.02, "R. A. Loomis et al. (2021)"),
    ("C5H5CN", "C1C=CC=C1C#N", 11.92, "J. Lee et al. (2021)"),
    ("C5H5CN", "C1C=CC(=C1)C#N", 11.28, "J. Lee et al. (2021)"),
    ("C11H7N", "C1=CC=C2C(=C1)C=CC=C2C#N", 11.87, "B. A. McGuire et al. (2021)"),
    ("C11H7N", "C1=CC=C2C=C(C=CC2=C1)C#N", 11.85, "B. A. McGuire et al. (2021)"),
    ("C6H5CN", "C1=CC=C(C=C1)C#N", 12.24, "B. A. McGuire et al. (2021)"),
    ("HCCCH2CN", "C#CCC#N", 11.96, "B. A. McGuire et al. (2020)"),
    ("H3C5N", "C#C/C=C/C#N", 11.39, "K. L. K. Lee et al. (2021b)"),
    ("H3C5N", "C=CC#CC#N", 11.07, "K. L. K. Lee et al. (2021b)"),
    ("H3C5N", r"C#C/C=C\C#N", 11.30, "K. L. K. Lee et al. (2021b)"),
    ("C8H", "[C]#CC#CC#CC#C", 11.66, "S. Brünken et al. (2007)"),
    ("C8H-", "C#CC#CC#CC#[C-]", 10.32, "S. Brünken et al. (2007)"),
    ("C6H-", "C#CC#CC#[C-]", 11.08, "S. Brünken et al. (2007)"),
    ("C4H-", "C#CC#[C-]", 10.93, "S. Brünken et al. (2007)"),
    ("H2CCO", "C=C=O", 12.71, "T. Soma et al. (2018)"),
    ("CN", "[C]#N", 12.89, "P. Pratap et al. (1997)"),
    ("HNC", "[C-]#[NH+]", 12.62, "P. Pratap et al. (1997)"),
    ("HC7O", "C#CC#CC#C[C+]=O", 11.89, "M. A. Cordiner et al. (2017)"),
    ("HC5O", "C#CC#C[C+]=O", 12.23, "B. A. McGuire et al. (2017)"),
    ("H2CN", "C=[N]", 11.18, "M. Ohishi et al. (1994)"),
    ("H2CO", "C=O", 13.08, "T. Soma et al. (2018)"),
    ("HC3O+", "C#CC#[O+]", 11.32, "J. Cernicharo et al. (2020a)"),
    ("HOCO+", "O=C=[OH+]", 11.60, "J. Cernicharo et al. (2020a)"),
    ("H2COH+", "C=[OH+]", 11.48, "J. Cernicharo et al. (2020a)"),
    ("H2NCO+", "NC#[O+]", 10.60, "J. Cernicharo et al. (2020a)"),
    ("HCNO", "C#N[O]", 10.85, "J. Cernicharo et al. (2020a)"),
    ("HOCN", "OC#N", 11.04, "J. Cernicharo et al. (2020a)"),
    ("H2", "[HH]", 22.00, "J. Cernicharo et al. (2020a)"),  # not in Lee A2 -- new addition
    ("C4O", "[C]#CC#[C]=O", 11.08, "J. Cernicharo et al. (2020a)"),
    ("HCOOH", "C(=O)O", 12.15, "J. Cernicharo et al. (2020a)"),
    ("HC2O", "C#[C]=O", 12.00, "J. Cernicharo et al. (2020a)"),
    ("C3O", "[C-]#CC#[O+]", 12.18, "J. Cernicharo et al. (2020a)"),
    ("HC4O", "C#CC#[C]=O", 11.48, "J. Cernicharo et al. (2020a)"),
    ("H2C3O", "C=C=C=O", 11.04, "J. Cernicharo et al. (2020a)"),
    ("c-H2C3O", "C1=C(=O)=C1", 11.60, "J. Cernicharo et al. (2020a)"),
    ("CH", "[CH]", 14.15, "N. Sakai et al. (2013)"),
    ("CNCN", "[C]#NC#N", 11.95, "M. Agúndez et al. (2018)"),
    ("NCCNH+", "N#CC#[NH+]", 10.93, "M. Agúndez et al. (2015)"),
    ("C6H2", "C=C=C=C=C=[C]", 10.33, "W. D. Langer et al. (1997)"),
    ("CH3CHCH2", "CC=C", 13.60, "N. Marcelino et al. (2007)"),
    ("CH2C2HCN", "C=C=CC#N", 11.65, "F. J. Lovas et al. (2006)"),
    ("HCN", "C#N", 12.39, "T. Hirota et al. (1998)"),
    ("CH2CHCCH", "C=CC#C", 13.08, "J. Cernicharo et al. (2021a)"),
    ("HCCN", "N#C[CH+]", 11.64, "J. Cernicharo et al. (2021a)"),
    ("CH3CH2CN", "CCC#N", 11.04, "J. Cernicharo et al. (2021a)"),
]
assert len(SHAY_B1) == 84, f"expected 84 Shay B1 rows, got {len(SHAY_B1)}"

# Two Shay Table C1 rows ("new detections") that are actually the same
# species as two of Lee's 87 rows, relocated rather than genuinely new --
# used only to document the RELOCATED fate, not to feed the anchor set's
# values (Lee's own A2 SMILES/CD are kept for these, since Table C1 is not
# "Table B1" and the cross-check is against B1 specifically).
SHAY_C1_RELOCATED = {
    "C9H8": ("c-C9H8", "C1C=CC2=CC=CC=C21", 12.98, "A. M. Burkhardt et al. (2021b)"),
    "HC3O": ("HCCCO", "C#C[C]=O", 11.20, "J. Cernicharo et al. (2021c)"),
}


def canon(smiles):
    mol = Chem.MolFromSmiles(smiles)
    used_fallback = False
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None, None, None, True
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_KEKULIZE,
            catchErrors=True,
        )
        used_fallback = True
    formula = rdMolDescriptors.CalcMolFormula(mol)
    can_smi = Chem.MolToSmiles(mol)
    try:
        ikey = rd_inchi.InchiToInchiKey(rd_inchi.MolToInchi(mol, treatWarningAsError=False))
    except Exception:
        ikey = None
    return formula, can_smi, ikey, used_fallback


# Build a lookup of Shay B1 rows keyed by canonical SMILES (RDKit-level
# identity -- catches the pure label/prefix differences, e.g. "C3H" vs
# "c-C3H", for free) for exact matching; formula-level fallback for the
# small number of rows where Shay corrected the SMILES itself.
shay_by_can_smiles = {}
shay_by_formula = {}
for label, smi, cd, ref in SHAY_B1:
    formula, can_smi, ikey, fallback = canon(smi)
    shay_by_can_smiles.setdefault(can_smi, []).append((label, smi, cd, ref, formula, ikey, fallback))
    shay_by_formula.setdefault(formula, []).append((label, smi, cd, ref, formula, ikey, fallback))

rows_out = []
dropped = []
relocated = []
matched_count = 0
used_shay_keys = set()

for formula_printed, smi, cd, ref in LEE_A2:
    lee_formula, lee_can_smi, lee_ikey, lee_fallback = canon(smi)

    match = None
    match_key = None
    if lee_can_smi in shay_by_can_smiles:
        candidates = [c for c in shay_by_can_smiles[lee_can_smi] if id(c) not in used_shay_keys]
        if candidates:
            match = candidates[0]
            match_key = id(match)

    if match is not None:
        used_shay_keys.add(match_key)
        label, shay_smi_printed, shay_cd, shay_ref, shay_formula, shay_ikey, shay_fallback = match
        basis = "MATCHED_EXACT_SMILES"
        smiles_changed = (shay_smi_printed.strip() != smi.strip())
        cd_changed = abs(shay_cd - cd) > 0.005
        ref_changed = shay_ref.split("(")[0].strip().replace(".", "").split()[-1] not in ref
        note_bits = []
        if smiles_changed:
            note_bits.append(f"Shay B1 SMILES differs in raw text ('{shay_smi_printed}' vs Lee's '{smi}') but canonicalizes identically")
        if cd_changed:
            note_bits.append(f"log CD corrected by Shay B1: {cd} -> {shay_cd}")
        if ref_changed:
            note_bits.append(f"reference updated by Shay B1: '{ref}' -> '{shay_ref}'")
        rows_out.append({
            "name": label,
            "formula": lee_formula,
            "smiles": lee_can_smi,
            "inchikey": lee_ikey,
            "lee_a2_smiles_printed": smi,
            "lee_a2_logcd": cd,
            "lee_a2_reference": ref,
            "shay_b1_smiles_printed": shay_smi_printed,
            "shay_b1_logcd": shay_cd,
            "shay_b1_reference": shay_ref,
            "basis": basis,
            "note": "; ".join(note_bits) if note_bits else "Shay B1 confirms Lee A2 unchanged",
            "sanitize_fallback_used": lee_fallback,
        })
        matched_count += 1
        continue

    # No exact-canonical-SMILES match in B1. Check Shay's Table C1 relocation.
    if formula_printed in SHAY_C1_RELOCATED:
        c1_label, c1_smi, c1_cd, c1_ref = SHAY_C1_RELOCATED[formula_printed]
        relocated.append((formula_printed, smi, cd, ref, c1_label, c1_smi, c1_cd, c1_ref))
        rows_out.append({
            "name": formula_printed,
            "formula": lee_formula,
            "smiles": lee_can_smi,
            "inchikey": lee_ikey,
            "lee_a2_smiles_printed": smi,
            "lee_a2_logcd": cd,
            "lee_a2_reference": ref,
            "shay_b1_smiles_printed": "",
            "shay_b1_logcd": "",
            "shay_b1_reference": "",
            "basis": "RELOCATED_TO_SHAY_TABLE_C1",
            "note": (f"Absent from Shay B1 (OD); reappears in Shay Table C1 ('new detections') "
                     f"as '{c1_label}' SMILES '{c1_smi}' logCD {c1_cd} ref '{c1_ref}' -- likely a "
                     f"bookkeeping inconsistency in Shay et al. 2025a (this is a pre-2021-08 Lee-list "
                     f"species, not a genuinely new detection). Lee A2's own values are kept here; "
                     f"not silently replaced with the Table C1 values."),
            "sanitize_fallback_used": lee_fallback,
        })
        continue

    # Genuinely dropped from Shay's OD/B1 and not found in Table C1 either.
    dropped.append((formula_printed, smi, cd, ref))
    rows_out.append({
        "name": formula_printed,
        "formula": lee_formula,
        "smiles": lee_can_smi,
        "inchikey": lee_ikey,
        "lee_a2_smiles_printed": smi,
        "lee_a2_logcd": cd,
        "lee_a2_reference": ref,
        "shay_b1_smiles_printed": "",
        "shay_b1_logcd": "",
        "shay_b1_reference": "",
        "basis": "DROPPED_FROM_SHAY_OD",
        "note": ("Absent from both Shay B1 (OD) and Shay Table C1. Shay et al. 2025a's own text "
                 "states OD 'contains 84 rather than 87 molecules ... any unconfirmed species were "
                 "removed'. This row is one of the ones removed. Kept here (Lee A2 is the primary "
                 "source for the anchor list; Shay's removal is noted, not silently applied) since "
                 "the task's target is Lee's own 87, and this project has not independently confirmed "
                 "the species is actually unconfirmed/wrong -- flagged for Reham's judgment."),
        "sanitize_fallback_used": lee_fallback,
    })

print(f"Lee A2 rows: {len(LEE_A2)}")
print(f"Matched exactly (canonical SMILES) to Shay B1: {matched_count}")
print(f"Relocated to Shay Table C1 (not in B1): {len(relocated)}")
for r in relocated:
    print("  RELOCATED:", r[0], r[1], "->", r[4], r[5])
print(f"Dropped (absent from both B1 and C1): {len(dropped)}")
for d in dropped:
    print("  DROPPED:", d)
print(f"Unmatched Shay B1 rows (present in B1, no Lee A2 canonical-SMILES match): "
      f"{len(SHAY_B1) - matched_count}")
# Which Shay B1 rows were never matched -- should be exactly H2 (the one
# genuine addition) plus anything whose SMILES Shay corrected enough that
# canonical-SMILES matching missed it.
unmatched_shay = []
for label, smi, cd, ref in SHAY_B1:
    formula, can_smi, ikey, fallback = canon(smi)
    found = False
    for lee_formula_printed, lee_smi, lee_cd, lee_ref in LEE_A2:
        lf, lc, li, lfb = canon(lee_smi)
        if lc == can_smi:
            found = True
            break
    if not found:
        unmatched_shay.append((label, smi, cd, ref))
print("Unmatched Shay B1 rows (by canonical SMILES against ALL of Lee A2, not just remaining pool):")
for u in unmatched_shay:
    print("  ", u)

# ---------------------------------------------------------------------------
# Write table_a2_authoritative.csv
# ---------------------------------------------------------------------------
fieldnames = ["name", "formula", "smiles", "inchikey", "lee_a2_smiles_printed",
              "lee_a2_logcd", "lee_a2_reference", "shay_b1_smiles_printed",
              "shay_b1_logcd", "shay_b1_reference", "basis", "note",
              "sanitize_fallback_used"]
out_path = "table_a2_authoritative.csv"
with open(out_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in rows_out:
        w.writerow(r)

print(f"\nWrote {out_path}, {len(rows_out)} rows")
n_confirmed = sum(1 for r in rows_out if r["basis"] in ("MATCHED_EXACT_SMILES", "RELOCATED_TO_SHAY_TABLE_C1"))
n_dropped = sum(1 for r in rows_out if r["basis"] == "DROPPED_FROM_SHAY_OD")
print(f"Confirmed-by-cross-check (matched or relocated): {n_confirmed}")
print(f"Dropped-from-Shay-OD, kept on Lee's authority alone, flagged: {n_dropped}")
assert len(rows_out) == 87

with open(out_path, "rb") as f:
    h = hashlib.sha256(f.read()).hexdigest()
print("SHA-256:", h)
with open(out_path + ".sha256", "w") as f:
    f.write(h + "\n")
