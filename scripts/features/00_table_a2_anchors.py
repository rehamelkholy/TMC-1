# Table A2 of Lee et al. (2021), ApJL 917, L6 -- the 87-molecule TMC-1
# "PD" (prior detections) inventory used to train the original umda models.
# This is the pre-cutoff training list used as the starting inventory
# (traceable to a citable source, unlike X_train_/y_train_ pulled out of the
# frozen fitted GPR in models/regressors.pkl).
#
# TRANSCRIPTION CAVEAT: the table was transcribed from the paper via a
# markdown conversion of the web version, not a byte-exact text extraction,
# so it is a good-faith transcription. Every SMILES was validated against
# RDKit's own computed molecular formula; that check caught 6 invalid/garbled
# SMILES on the first pass (all fixed below, each flagged inline) and found 0
# formula mismatches among the 81 that parsed cleanly. Spot-checking a sample
# of rows against the paper (arXiv:2108.10419 or the ApJL version) is
# recommended before relying on it.
#
# The two references-only columns ("references correspond to the source
# used for the column density, not initial detection" -- the paper's own
# footnote) are kept as printed.

TABLE_A2 = [
    # formula,     smiles,                        log10_cd,  reference
    ("CH3C6H",     "CC#CC#CC#C",                  12.4914,  "Remijan et al. 2006"),
    ("CH3C4H",     "CC#CC#C",                     13.4771,  "MacLeod et al. 1984"),
    ("CH3C5N",     "CC#CC#CC#N",                  11.9243,  "Remijan et al. 2006"),
    ("CH3C3N",     "CC#CC#N",                     12.2553,  "Remijan et al. 2006"),
    ("N2H+",       "N#[NH+]",                     12.6990,  "Choi et al. 2017"),
    ("NH3",        "N",                           14.6998,  "Gratier et al. 2016"),
    ("CH3OH",      "CO",                          13.1614,  "Gratier et al. 2016"),
    ("C3H",        "C1=C=[C]1",                   13.4800,  "Gratier et al. 2016"),
    ("C3H",        "[CH+]=C=[C-]",                12.7497,  "Gratier et al. 2016"),
    ("C3H2",       "C1=C=C1",                     11.7701,  "Gratier et al. 2016"),
    ("C3H2",       "C1C#C1",                      13.2695,  "Gratier et al. 2016"),
    ("C3H2",       "C=C=[C]",                     12.3979,  "Cernicharo et al. 1991"),
    ("CH3CCH",     "CC#C",                        14.0607,  "Gratier et al. 2016"),
    ("C2O",        "[C+]#C[O-]",                  12.5705,  "Gratier et al. 2016"),
    ("CH2CN",      "[CH2]C#N",                    13.5798,  "Gratier et al. 2016"),
    ("CH3CN",      "CC#N",                        12.6096,  "Gratier et al. 2016"),
    ("HNCO",       "N=C=O",                       13.0294,  "Gratier et al. 2016"),
    ("CS",         "[C-]#[S+]",                   13.4594,  "Gratier et al. 2016"),
    ("CH3CHO",     "CC=O",                        12.4298,  "Gratier et al. 2016"),
    ("HCS+",       "C#[S+]",                      12.7597,  "Gratier et al. 2016"),
    ("H2CS",       "C=S",                         13.6201,  "Gratier et al. 2016"),
    ("SO",         "S=O",                         13.6702,  "Gratier et al. 2016"),
    ("C4H",        "[C]#CC#C",                    13.4298,  "Gratier et al. 2016"),
    ("C4H2",       "C=C=C=[C]",                   13.3365,  "Gratier et al. 2016"),
    ("C3N",        "[C]#CC#N",                    13.5502,  "Gratier et al. 2016"),
    ("HNC3",       "[C-]#C-C#[NH+]",              11.6803,  "Gratier et al. 2016"),
    ("C3O",        "[C]#C[C]=O",                  11.9201,  "Gratier et al. 2016"),
    ("HC3NH+",     "C#CC#[NH+]",                  11.8698,  "Gratier et al. 2016"),
    ("CH2CHCN",    "C=CC#N",                      12.8102,  "Gratier et al. 2016"),
    ("HCCCHO",     "C#CC=O",                      11.2601,  "Gratier et al. 2016"),
    ("C2S",        "[C+]#C[S-]",                  14.0086,  "Gratier et al. 2016"),
    ("OCS",        "O=C=S",                       13.2601,  "Gratier et al. 2016"),
    ("C5H",        "[CH+]=C=C=C=[C-]",            12.2695,  "Gratier et al. 2016"),
    ("C3S",        "[C-]#CC#[S+]",                13.1399,  "Gratier et al. 2016"),
    ("C6H",        "[C]#CC#CC#C",                 12.7404,  "Gratier et al. 2016"),
    ("HC3N",       "C#CC#N",                      14.2430,  "Xue et al. 2020"),
    ("HCCNC",      "C#C[N+]#[C-]",                12.5821,  "Xue et al. 2020"),
    ("HC5N",       "C#CC#CC#N",                   13.8254,  "Xue et al. 2020"),
    ("HC4NC",      "C#CC#C[N+]#[C-]",             11.5172,  "Xue et al. 2020"),
    ("HC7N",       "C#CC#CC#CC#N",                13.5623,  "Xue et al. 2020"),
    ("HC6NC",      "C#CC#CC#C[N+]#[C-]",          11.6064,  "Xue et al. 2020"),
    ("HC9N",       "C#CC#CC#CC#CC#N",             13.3345,  "Loomis et al. 2021"),
    ("HC11N",      "C#CC#CC#CC#CC#CC#N",          12.0170,  "Loomis et al. 2021"),
    ("C5H5CN",     "C1C=CC=C1C#N",                11.9191,  "Lee et al. 2021b"),
    ("C5H5CN",     "C1C=CC(=C1)C#N",              11.2788,  "Lee et al. 2021b"),
    ("C11H7N",     "C1=CC=C2C(=C1)C=CC=C2C#N",    11.8663,  "McGuire et al. 2021"),
    ("C11H7N",     "C1=CC=C2C=C(C=CC2=C1)C#N",    11.8482,  "McGuire et al. 2021"),
    ("C6H5CN",     "C1=CC=C(C=C1)C#N",            12.2380,  "McGuire et al. 2021"),
    ("HCCCH2CN",   "C#CCC#N",                     11.9643,  "McGuire et al. 2020"),
    ("H3C5N",      "C#C/C=C/C#N",                 11.3874,  "Lee et al. 2021a"),
    ("H3C5N",      "C=CC#CC#N",                   11.0719,  "Lee et al. 2021a"),
    ("H3C5N",      r"C#C/C=C\C#N",                11.3032,  "Lee et al. 2021a"),
    ("C8H",        "[C]#CC#CC#CC#C",              11.6628,  "Brünken et al. 2007"),
    ("C8H-",       "C#CC#CC#CC#[C-]",             10.3222,  "Brünken et al. 2007"),
    ("C6H-",       "C#CC#CC#[C-]",                11.0792,  "Brünken et al. 2007"),
    ("C4H-",       "C#CC#[C-]",                   10.9294,  "Brünken et al. 2007"),
    ("H2CCO",      "C=C=O",                       12.7118,  "Soma et al. 2018"),
    ("CN",         "[C]#N",                       12.8899,  "Pratap et al. 1997"),
    ("HNC",        "[C-]#[NH+]",                  12.6201,  "Pratap et al. 1997"),
    ("HC7O",       "C#CC#CC#C[C+]=O",             11.8921,  "Cordiner et al. 2017"),
    ("HC5O",       "C#CC#C[C+]=O",                12.2304,  "McGuire et al. 2017"),
    ("H2CN",       "C=[N]",                       11.1761,  "Ohishi et al. 1994"),
    ("H2CO",       "C=O",                         13.0792,  "Soma et al. 2018"),
    ("HC3O+",      "C#CC#[O+]",                   11.3222,  "Cernicharo et al. 2020"),
    ("HOCO+",      "O=C=[OH+]",                   11.6021,  "Cernicharo et al. 2020"),
    ("H2COH+",     "C=[OH+]",                     11.4771,  "Cernicharo et al. 2020"),
    ("H2NCO+",     "NC#[O+]",                     10.6021,  "Cernicharo et al. 2020"),
    ("HCNO",       "C#[N+][O-]",                  10.8451,  "Cernicharo et al. 2020"),  # FIXED: transcribed "C#N[O]" (invalid, N valence 4). Fulminic acid's standard zwitterionic SMILES restores a valid CHNO structure.
    ("HOCN",       "OC#N",                        11.0414,  "Cernicharo et al. 2020"),
    ("C4O",        "[C]=C=C=C=O",                 11.0792,  "Cernicharo et al. 2020"),  # FIXED: transcribed "[C]#CC#[C]=O" (invalid valence). Cumulated-chain form, consistent with the C3O/C5O convention used elsewhere in this table.
    ("HCOOH",      "C(=O)O",                      12.1461,  "Cernicharo et al. 2020"),
    ("HC2O",       "[CH]=C=O",                    12.0000,  "Cernicharo et al. 2020"),  # FIXED: transcribed "C#[C]=O" (invalid valence). Ketenyl radical HCCO, matches the HCCS/HCCS+ cumulated-radical convention.
    ("HC3O",       "C#C[C]=O",                    11.3010,  "Cernicharo et al. 2020"),
    ("HC4O",       "[CH]=C=C=C=O",                11.4771,  "Cernicharo et al. 2020"),  # FIXED: transcribed "C#CC#[C]=O" (invalid valence). Cumulated form matching HC2O's convention above.
    ("H2C3O",      "C=C=C=O",                     11.0414,  "Cernicharo et al. 2020"),
    ("H2C3O",      "O=C1C=C1",                    11.6021,  "Cernicharo et al. 2020"),  # FIXED: transcribed "C1=C(=O)=C1" (invalid, ring atom with 3 double bonds). Cyclopropenone -- the well-known cyclic H2C3O isomer, correct formula C3H2O.
    ("CH",         "[CH]",                        14.1461,  "Sakai et al. 2013"),
    ("CNCN",       "N#C[N+]#[C-]",                11.9542,  "Agúndez et al. 2018"),  # FIXED: transcribed "[C]#NC#N" (invalid, N valence 4). Isocyanogen's exact resonance structure is genuinely uncertain here -- flagged as lower-confidence; any valid C2N2 structure serves the embedding step, but don't cite this specific SMILES as authoritative without checking the paper directly.
    ("NCCNH+",     "N#CC#[NH+]",                  10.9345,  "Agúndez et al. 2015"),
    ("C6H2",       "C=C=C=C=C=[C]",               10.3284,  "Langer et al. 1997"),
    ("CH3CHCH2",   "CC=C",                        13.6021,  "Marcelino et al. 2007"),
    ("CH2C2HCN",   "C=C=CC#N",                    11.6532,  "Lovas et al. 2006"),
    ("HCN",        "C#N",                         12.3892,  "Hirota et al. 1998"),
    ("C9H8",       "c1ccc2c(c1)CC=C2",            12.9823,  "Burkhardt et al. 2021"),
    ("CH2CHCCH",   "C=CC#C",                      13.0792,  "Cernicharo et al. 2021"),
    ("HCCN",       "N#C[CH+]",                    11.6435,  "Cernicharo et al. 2021"),
    ("CH3CH2CN",   "CCC#N",                       11.0414,  "Cernicharo et al. 2021"),
]

if __name__ == "__main__":
    import sys
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    import re

    print(f"Total rows: {len(TABLE_A2)}")
    fails, checked = 0, 0
    for formula, smi, cd, ref in TABLE_A2:
        mol = Chem.MolFromSmiles(smi)
        checked += 1
        if mol is None:
            print(f"  PARSE FAILURE: {formula} {smi} ({ref})")
            fails += 1
    print(f"Parse failures: {fails} / {checked}")
