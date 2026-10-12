"""
Rule C key (protocol_v1): "formula incl. net charge + heavy-atom connectivity
with per-atom H counts; bond orders and charge placement ignored."

Algorithm:
  1. Parse SMILES with RDKit (sanitize=True; on failure, partial-sanitize
     fallback, as in provenance/scripts/11_purge.py).
  2. net_charge = Chem.GetFormalCharge(mol) (sum over atoms -- invariant to
     which atom carries the charge).
  3. AddHs to make every H an explicit atom (so per-atom H count becomes
     part of the graph itself).
  4. Build a flattened copy: every bond forced to SINGLE, every atom's
     formal charge zeroed, no radicals, no stereo (isotopes ARE kept, since
     isotopologues are treated as distinct species -- e.g. HDCO is kept
     separate from H2CO).
  5. Canonical SMILES of the flattened, all-explicit-H graph = the
     connectivity + per-atom-H-count signature.
  6. Key = "<heavy-atom formula>|<net_charge:+d>|<flattened canonical SMILES>"

This key is charge-placement-invariant and bond-order-invariant by construction
(step 4), distinguishes any difference in which heavy atom an H is attached to
(steps 3-5), and distinguishes true connectivity differences (isomers / different
skeletons), as rule C requires.
"""
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def _parse(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    mol.UpdatePropertyCache(strict=False)
    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_KEKULIZE,
            catchErrors=True,
        )
    except Exception:
        pass
    return mol


def rule_c_key(smiles):
    """Returns (key, net_charge, heavy_formula) or (None, None, None) if
    the SMILES does not parse at all."""
    mol = _parse(smiles)
    if mol is None:
        return None, None, None

    try:
        net_charge = Chem.GetFormalCharge(mol)
    except Exception:
        net_charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())

    heavy_formula = rdMolDescriptors.CalcMolFormula(mol)
    # CalcMolFormula includes the net charge suffix (e.g. "CH3O+"); strip it
    # for a charge-free heavy/H formula, since charge is reported separately.
    import re
    heavy_formula_nocharge = re.sub(r"[+-]\d*$", "", heavy_formula)

    try:
        molH = Chem.AddHs(mol)
    except Exception:
        return None, None, None

    # Build flattened copy: same atoms (incl. isotopes), all bonds -> SINGLE,
    # all formal charges -> 0, all radical electrons -> 0.
    rw = Chem.RWMol()
    idx_map = {}
    for atom in molH.GetAtoms():
        a = Chem.Atom(atom.GetAtomicNum())
        if atom.GetIsotope():
            a.SetIsotope(atom.GetIsotope())
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        a.SetNumRadicalElectrons(0)
        new_idx = rw.AddAtom(a)
        idx_map[atom.GetIdx()] = new_idx
    seen_bonds = set()
    for bond in molH.GetBonds():
        i, j = idx_map[bond.GetBeginAtomIdx()], idx_map[bond.GetEndAtomIdx()]
        key_ij = (min(i, j), max(i, j))
        if key_ij in seen_bonds:
            continue
        seen_bonds.add(key_ij)
        rw.AddBond(i, j, Chem.BondType.SINGLE)

    flat = rw.GetMol()
    try:
        Chem.SanitizeMol(
            flat,
            sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_KEKULIZE
            ^ Chem.SANITIZE_SETAROMATICITY,
            catchErrors=True,
        )
    except Exception:
        pass

    try:
        canon = Chem.MolToSmiles(flat, canonical=True)
    except Exception:
        return None, None, None

    charge_str = f"{net_charge:+d}" if net_charge != 0 else "0"
    key = f"{heavy_formula_nocharge}|{charge_str}|{canon}"
    return key, net_charge, heavy_formula_nocharge


if __name__ == "__main__":
    pairs = [
        ("293 H2COH+", "[CH2+]O", "C=[OH+]"),
        ("376 HNC3", "[C-]#C[C+]=N", "[C-]#CC#[NH+]"),
        ("396 HOCO+", "O=[C+]O", "O=C=[OH+]"),
        ("490 HC3NH+", "C#C[C+]=N", "C#CC#[NH+]"),
        ("523 NCCNH+", "N#C[C+]=N", "N#CC#[NH+]"),
    ]
    all_match = True
    for label, cand_smi, anchor_smi in pairs:
        k1, c1, f1 = rule_c_key(cand_smi)
        k2, c2, f2 = rule_c_key(anchor_smi)
        match = (k1 == k2)
        all_match &= match
        print(f"{label}: match={match}")
        print(f"  candidate {cand_smi!r:20s} -> {k1}")
        print(f"  anchor    {anchor_smi!r:20s} -> {k2}")
    print("\nALL MATCH:", all_match)
