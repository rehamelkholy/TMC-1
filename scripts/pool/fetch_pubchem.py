"""
Fetch PubChem (Kim et al. 2021, Nucleic Acids Res. 49, D1388) -- Lee et al. 2021 report 2,444,441 entries
from this source.

The source is PubChem's official bulk compound-SMILES dump: umda's own companion script
(umda/scripts/pubchem_preprocess.py) reads a local file named "CID-SMILES" and takes the first 3,000,000
lines before deduplicating, which is why the final count differs slightly from the paper's reported
number (duplicates and unparsable entries fall out along the way). This script reproduces that logic.

WARNING: the official file (https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz) is very
large (multiple GB compressed, >100M rows). This script streams and decompresses it without loading it
into memory and stops reading once it has enough candidate lines, but the download itself can take a
while.

Output: external/pool/raw/pubchem_smiles.smi
"""
import gzip
import os

import http_common

MAX_NUM = 3_000_000  # matches umda's own scripts/pubchem_preprocess.py exactly

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "external", "pool", "raw")
os.makedirs(OUT_DIR, exist_ok=True)
GZ_PATH = os.path.join(OUT_DIR, "CID-SMILES.gz")
OUT_PATH = os.path.join(OUT_DIR, "pubchem_smiles.smi")

PUBCHEM_URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz"


def main():
    if not os.path.exists(GZ_PATH):
        print(f"Downloading PubChem CID-SMILES from {PUBCHEM_URL} ...")
        print("(This is a large file -- expect this to take a while.)")
        http_common.download_to_file(PUBCHEM_URL, GZ_PATH)
    else:
        print(f"Found existing {GZ_PATH}, skipping download.")

    from rdkit import Chem

    seen = set()
    n_read = 0
    with gzip.open(GZ_PATH, "rt", errors="ignore") as read_file:
        for line in read_file:
            if n_read >= MAX_NUM:
                break
            parts = line.split()
            if len(parts) < 2:
                continue
            smi = parts[1]
            n_read += 1
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                seen.add(Chem.MolToSmiles(mol, canonical=True))
            if n_read % 200_000 == 0:
                print(f"  processed {n_read:,} raw lines, {len(seen):,} unique canonical SMILES so far")

    with open(OUT_PATH, "w") as f:
        for smi in seen:
            f.write(smi + "\n")
    print(f"Wrote {len(seen):,} unique canonical SMILES to {OUT_PATH} "
          f"(paper reports 2,444,441 from this source).")


if __name__ == "__main__":
    main()
