"""
Fetch PCBA / PubChem BioAssay (Wang et al. 2012, Nucleic Acids Res. 40, D400) -- 437,929 molecules per
Lee et al. 2021's own citation.

This is the standard MoleculeNet "pcba" benchmark dataset -- a SMILES-plus-128-bioassay-outcome-columns
table. Only the SMILES column is used here.

Output: external/pool/raw/pcba_smiles.smi
"""
import os

import http_common

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "external", "pool", "raw")
os.makedirs(OUT_DIR, exist_ok=True)
CSV_GZ_PATH = os.path.join(OUT_DIR, "pcba.csv.gz")
OUT_PATH = os.path.join(OUT_DIR, "pcba_smiles.smi")

# If this URL has moved, locate the standard MoleculeNet "pcba.csv.gz" benchmark file and place it at
# external/pool/raw/pcba.csv.gz, then re-run.
PCBA_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/pcba.csv.gz"


def main():
    if not os.path.exists(CSV_GZ_PATH):
        print(f"Downloading PCBA from {PCBA_URL} (this is a few hundred MB) ...")
        http_common.download_to_file(PCBA_URL, CSV_GZ_PATH)
    else:
        print(f"Found existing {CSV_GZ_PATH}, skipping download.")

    import pandas as pd
    from rdkit import Chem

    df = pd.read_csv(CSV_GZ_PATH, compression="gzip", usecols=lambda c: c.lower() == "smiles" or c == "smiles")
    smi_col = df.columns[0]
    n_written = 0
    with open(OUT_PATH, "w") as f:
        for smi in df[smi_col].dropna():
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                f.write(Chem.MolToSmiles(mol, canonical=True) + "\n")
                n_written += 1
    print(f"Wrote {n_written} canonical SMILES to {OUT_PATH} "
          f"(paper reports 437,929 -- some loss to parse failures/dedup is normal).")


if __name__ == "__main__":
    main()
