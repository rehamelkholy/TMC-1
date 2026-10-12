"""
Fetch QM9 (Ramakrishnan et al. 2014, Sci. Data 1, 1) -- 133,885 molecules, per Lee et al. 2021's own
citation.

The original umda paper worked from the raw QM9 .xyz files (Figshare collection:
https://doi.org/10.6084/m9.figshare.978904), which requires a bond-perception step to get SMILES from
Cartesian coordinates. MoleculeNet redistributes the same dataset with SMILES already computed (same
molecules, same source paper); this script uses that route.

Output: external/pool/raw/qm9_smiles.smi (one canonical SMILES per line)
"""
import os

import http_common

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "external", "pool", "raw")
os.makedirs(OUT_DIR, exist_ok=True)
CSV_PATH = os.path.join(OUT_DIR, "qm9.csv")
OUT_PATH = os.path.join(OUT_DIR, "qm9_smiles.smi")

# MoleculeNet's hosted copy (same underlying Ramakrishnan et al. 2014 QM9 set). If this URL has moved,
# download qm9.csv manually from https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.csv (or
# locate the MoleculeNet "qm9.csv" file) and place it at external/pool/raw/qm9.csv, then re-run this
# script -- it skips the download if the file exists.
QM9_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.csv"


def main():
    if not os.path.exists(CSV_PATH):
        print(f"Downloading QM9 from {QM9_URL} ...")
        http_common.download_to_file(QM9_URL, CSV_PATH)
    else:
        print(f"Found existing {CSV_PATH}, skipping download.")

    import pandas as pd
    from rdkit import Chem

    df = pd.read_csv(CSV_PATH)
    smi_col = "smiles" if "smiles" in df.columns else df.columns[df.columns.str.lower() == "smiles"][0]
    n_written = 0
    with open(OUT_PATH, "w") as f:
        for smi in df[smi_col].dropna():
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                f.write(Chem.MolToSmiles(mol, canonical=True) + "\n")
                n_written += 1
    print(f"Wrote {n_written} canonical SMILES to {OUT_PATH} "
          f"(paper reports 133,885 -- some loss to parse failures is normal).")


if __name__ == "__main__":
    main()
