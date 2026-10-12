"""Parse the Lee et al. (2021) recommendation list (PDF) into data/raw/lee_recommendations_parsed.csv.

Written for the release. The parser that originally produced data/raw/lee_recommendations_parsed.csv is not
available, so this script re-implements it from the layout of the PDF. Run on version 2 of the deposit
(Zenodo 10.5281/zenodo.5146276, tmc1_recommendations.pdf, MD5 9d3664ca507b5542e2bf3a279599f891; fetched by
scripts/00_fetch_external.sh) it reproduces the shipped CSV byte for byte.

Each of the 1,510 entries of the PDF is a block that starts with "geom<N>" and contains the lines
    SMILES: <smiles>
    Nearest TMC-1 molecule (distance): <smiles> (<distance>)
    Formula <formula>
    Charge, Multiplicity <charge>, <multiplicity>
    Predicted log column density <logcd>±<sigma>
The remaining fields of the PDF (molecular weight, IUPAC name, dipole moment, rotational constants,
electronic energy) are not used by the analysis and are not written.

Output columns, in order: id (= the N of geomN), smiles, canonical_smiles, formula, charge, mult, nearest,
nearest_canonical, dist, logcd, sigma. The two canonical columns are RDKit's sanitised canonical SMILES
(Chem.MolToSmiles(Chem.MolFromSmiles(s))) of the smiles and nearest columns; they are empty where RDKit cannot
sanitise the SMILES (7 and 56 entries). The CSV uses CRLF line endings, as the shipped file does.

Needs the poppler `pdftotext` command (apt install poppler-utils / brew install poppler) and RDKit.

Usage: python scripts/01_parse_lee_recommendations.py [pdf] [output.csv]
Defaults: external/tmc1_recommendations.pdf -> build/lee_recommendations_parsed.csv (the shipped file in data/raw
is never overwritten). The script then compares the result with data/raw/lee_recommendations_parsed.csv.
"""
import hashlib
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import pandas as pd
from rdkit import Chem, RDLogger

from paths import BUILD, DATA_RAW, EXTERNAL

RDLogger.DisableLog("rdApp.*")

N_ENTRIES = 1510
COLUMNS = ["id", "smiles", "canonical_smiles", "formula", "charge", "mult",
           "nearest", "nearest_canonical", "dist", "logcd", "sigma"]


def pdf_text(pdf_path):
    out = subprocess.run(["pdftotext", "-layout", pdf_path, "-"], check=True,
                         capture_output=True, text=True).stdout
    # the PDF uses non-breaking hyphens (U+2011) and minus signs (U+2212)
    return out.replace("‑", "-").replace("−", "-")


def parse_entries(text):
    parts = re.split(r"\n\s*geom(\d+)\s*\n", "\n" + text)
    rows = []
    for i in range(1, len(parts), 2):
        gid, block = int(parts[i]), parts[i + 1]

        def grab(pattern):
            m = re.search(pattern, block)
            if m is None:
                raise ValueError(f"geom{gid}: pattern not found: {pattern}")
            return m.groups()

        (smiles,) = grab(r"SMILES:\s*(\S+)")
        nearest, dist = grab(r"Nearest TMC-1 molecule \(distance\):\s*(\S+)\s*\(([\d.]+)\)")
        (formula,) = grab(r"Formula[ \t]+(\S+)")
        charge, mult = grab(r"Charge, Multiplicity[ \t]+(-?\d+),\s*(\d+)")
        logcd, sigma = grab(r"Predicted log column density[ \t]+(-?[\d.]+)\s*±\s*([\d.]+)")
        rows.append({"id": gid, "smiles": smiles, "formula": formula, "charge": int(charge), "mult": int(mult),
                     "nearest": nearest, "dist": float(dist), "logcd": float(logcd), "sigma": float(sigma)})
    return pd.DataFrame(rows)


def canonical(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m is not None else ""


def parse_pdf(pdf_path):
    df = parse_entries(pdf_text(pdf_path))
    assert len(df) == N_ENTRIES, f"expected {N_ENTRIES} entries, found {len(df)}"
    assert df["id"].tolist() == list(range(1, N_ENTRIES + 1)), "entry ids are not 1..1510 in order"
    df["canonical_smiles"] = df["smiles"].map(canonical)
    df["nearest_canonical"] = df["nearest"].map(canonical)
    return df[COLUMNS]


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


if __name__ == "__main__":
    pdf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(EXTERNAL, "tmc1_recommendations.pdf")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(BUILD, "lee_recommendations_parsed.csv")
    df = parse_pdf(pdf)
    df.to_csv(out, index=False, lineterminator="\r\n")
    print(f"{len(df)} entries -> {out}")
    ref = os.path.join(DATA_RAW, "lee_recommendations_parsed.csv")
    if os.path.exists(ref) and os.path.abspath(out) != os.path.abspath(ref):
        same = open(out, "rb").read() == open(ref, "rb").read()
        print(f"identical to data/raw/lee_recommendations_parsed.csv: {same}")
        print(f"  parsed {sha256(out)}\n  shipped {sha256(ref)}")
        sys.exit(0 if same else 1)
