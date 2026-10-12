"""
Per-chunk worker for canonicalizing the combined candidate pool. Applies Lee-style composition
filters and computes charge/radical/canonical-SMILES annotations for one chunk of
combined_candidate_pool.csv (columns: smiles, source, mw; no header).

Filters:
  - elements_ok: every atom is one of H, C, N, O, Si, P, S
  - no_dot: raw SMILES has no '.' (no disconnected fragments)
  - no_isotope: no atom has a non-zero isotope label
  - mw_ok: mw <= 500 (applied upstream in combine_pool.py; verified, not re-filtered)
An entry "passes filters" iff elements_ok AND no_dot AND no_isotope AND mw_ok AND it parses at all
under the embedder-compatible convention (sanitize=False + UpdatePropertyCache(strict=False)
+ GetSymmSSSR).

Charge is read directly off the sanitize=False parse (formal charge perception does not require
sanitization). The radical-electron count needs FULL sanitization to be reliable (sanitize=False
silently drops radical-electron assignment for under-valent bracket atoms); when full sanitization
fails, radical status is reported as "unknown_sanitize_failed" rather than guessed. Canonical SMILES
(for the cross-list join) is likewise only taken from the full-sanitize parse, so one consistent
convention applies across the pool and the reference lists.
"""
import sys
import time

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ALLOWED_ELEMENTS = {"H", "C", "N", "O", "Si", "P", "S"}


def process(in_path, out_path):
    t0 = time.time()
    df = pd.read_csv(in_path, header=None, names=["smiles", "source", "mw"])
    n = len(df)

    elements_ok = [False] * n
    no_dot = [False] * n
    no_isotope = [False] * n
    embed_parse_ok = [False] * n
    charged = [None] * n
    n_atoms_heavy = [None] * n
    radical_status = ["unknown_sanitize_failed"] * n
    n_radical_electrons = [None] * n
    canonical_smiles_full = [None] * n

    for i, smi in enumerate(df["smiles"].to_numpy()):
        no_dot[i] = "." not in smi

        # sanitize=False parse: embedder-compatible convention, used for
        # element/isotope/charge checks and the "does this even count as
        # a pool entry the embedder could touch" gate.
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is not None:
            try:
                m.UpdatePropertyCache(strict=False)
                Chem.GetSymmSSSR(m)
                embed_parse_ok[i] = True
                syms = {a.GetSymbol() for a in m.GetAtoms()}
                elements_ok[i] = syms.issubset(ALLOWED_ELEMENTS)
                no_isotope[i] = all(a.GetIsotope() == 0 for a in m.GetAtoms())
                charged[i] = any(a.GetFormalCharge() != 0 for a in m.GetAtoms())
                n_atoms_heavy[i] = m.GetNumAtoms()
            except Exception:
                pass

        # full-sanitize parse: reliable radical-electron perception +
        # canonical SMILES for the cross-list join.
        try:
            mf = Chem.MolFromSmiles(smi)  # default sanitize=True
            if mf is not None:
                n_re = sum(a.GetNumRadicalElectrons() for a in mf.GetAtoms())
                n_radical_electrons[i] = n_re
                radical_status[i] = "radical" if n_re > 0 else "closed_shell"
                canonical_smiles_full[i] = Chem.MolToSmiles(mf)
        except Exception:
            pass

    df["elements_ok"] = elements_ok
    df["no_dot"] = no_dot
    df["no_isotope"] = no_isotope
    df["mw_ok"] = df["mw"] <= 500.0
    df["embed_parse_ok"] = embed_parse_ok
    df["passes_filters"] = (
        df["elements_ok"] & df["no_dot"] & df["no_isotope"] & df["mw_ok"] & df["embed_parse_ok"]
    )
    df["charged"] = charged
    df["n_atoms_heavy"] = n_atoms_heavy
    df["radical_status"] = radical_status
    df["n_radical_electrons"] = n_radical_electrons
    df["canonical_smiles_full"] = canonical_smiles_full

    df.to_parquet(out_path, index=False)
    elapsed = time.time() - t0
    print(f"{in_path}: {n} rows in {elapsed:.1f}s ({n/elapsed:.0f} rows/s) -> {out_path}", flush=True)
    return elapsed, n


if __name__ == "__main__":
    in_path, out_path = sys.argv[1], sys.argv[2]
    process(in_path, out_path)
