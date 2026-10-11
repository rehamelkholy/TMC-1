"""
Fetch the NASA Ames PAH IR Spectroscopic Database (Boersma et al. 2014, ApJS 211, 8; Bauschlicher et al.
2018, ApJS 234, 32; Mattioda et al. 2020, ApJS 251, 22) -- Lee et al. 2021 report 3,139 entries,
converted from the database's Cartesian-coordinate structures to SMILES via OpenBabel.

PAHdb is gated behind a free email-registration form, not a direct download link, so one manual step is
needed (the downloaded file is subject to the database's own terms of use and is not redistributed here):
  1. Go to https://www.astrochemistry.org/pahdb/theoretical/4.00/download/view
  2. Fill in the short form (email, name, affiliation) and submit.
  3. A download link for the XML library file arrives by email.
  4. Save the downloaded file as external/pool/raw/pahdb.xml (it may be gzip-compressed under that
     name; this script sniffs the first two bytes and decompresses transparently).
Then run this script.

Parsing notes: the document declares a default XML namespace (xmlns="http://www.astrochemistry.org/
pahdb/theoretical"), so tag names are matched by their local part; each atom's element is given as an
atomic NUMBER in a <type> tag (an atomic-number-to-symbol table is used); <position> next to each atom's
<x>/<y>/<z> is that atom's 1-based INDEX in the list, with x/y/z as flat siblings. The file is parsed with
ET.iterparse and per-<specie> clearing to keep memory bounded (~10,749 species with vibrational-mode data).

Formula cross-check: for every species, the formula of the molecule built from atomic coordinates is
independently recomputed by RDKit and compared against the file's own stated <formula> field.
Mismatches are counted and reported, not silently dropped (the stated formula may reflect isotope
labeling or other conventions a plain formula-string comparison cannot fully account for).

Coordinates are formatted as fixed-point (`f"{x:.8f}"`): Python's default float-to-string conversion
switches to bare scientific notation for small magnitudes (e.g. "1e-06"), which RDKit's
MolFromXYZBlock cannot parse.

Output: external/pool/raw/pahdb_smiles.smi
"""
import gzip
import os
import xml.etree.ElementTree as ET

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "external", "pool", "raw")
XML_PATH = os.path.join(OUT_DIR, "pahdb.xml")
OUT_PATH = os.path.join(OUT_DIR, "pahdb_smiles.smi")

# Full atomic-number -> symbol table (Z=1..118). PAHdb species are
# overwhelmingly C/H/N/O (plus occasional Mg/Si/Fe for grain/metal-PAH
# studies), but there's no reason to arbitrarily truncate this.
_SYMBOLS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
ATOMIC_SYMBOLS = {i + 1: sym for i, sym in enumerate(_SYMBOLS)}


def _open_xml_source(path):
    """Return a readable binary file object for path, transparently
    decompressing if it turns out to be gzip-compressed under the hood
    (the PAHdb download link hands back gzip bytes even though the file
    you save is conventionally named with a plain .xml extension)."""
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return open(path, "rb")


def _local(tag):
    """Strip a namespace prefix ('{uri}tag' -> 'tag'); a no-op if there
    isn't one. Used everywhere instead of raw tag/find comparisons so this
    script doesn't care whether the document declares a default xmlns."""
    return tag.rsplit("}", 1)[-1]


def _child(elem, name):
    """First direct child of elem whose local tag name == name, namespace-
    agnostic (see _local)."""
    for child in elem:
        if _local(child.tag) == name:
            return child
    return None


def _child_text(elem, name):
    child = _child(elem, name)
    return child.text.strip() if child is not None and child.text else None


def mol_from_geometry(symbols, coords, charge=0):
    """Build an RDKit mol from element symbols + XYZ coordinates via
    distance-based bond perception (RDKit's rdDetermineBonds, the same kind
    of step OpenBabel does implicitly)."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    # Fixed-point formatting is deliberate: Python's default float-to-string switches to bare scientific
    # notation for small magnitudes (e.g. str(1e-06) == "1e-06"), which RDKit's XYZ-block parser cannot
    # parse ("Cannot convert '1e-06' to double"). PAHdb's optimized geometries have many near-zero
    # coordinate components. Plain decimal notation always parses.
    xyz_block = f"{len(symbols)}\n\n"
    for sym, (x, y, z) in zip(symbols, coords):
        xyz_block += f"{sym} {x:.8f} {y:.8f} {z:.8f}\n"
    mol = Chem.MolFromXYZBlock(xyz_block)
    if mol is None:
        return None
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
    except Exception:
        return None
    return mol


def _parse_specie(specie_elem):
    """Return (symbols, coords, stated_formula, charge) for one <specie>
    element, or None if it has no usable geometry."""
    geometry = _child(specie_elem, "geometry")
    if geometry is None:
        return None

    symbols, coords = [], []
    for atom_el in geometry:
        if _local(atom_el.tag) != "atom":
            continue
        type_text = _child_text(atom_el, "type")
        x_text = _child_text(atom_el, "x")
        y_text = _child_text(atom_el, "y")
        z_text = _child_text(atom_el, "z")
        if type_text is None or x_text is None or y_text is None or z_text is None:
            continue
        try:
            z_num = int(type_text)
            x, y, z = float(x_text), float(y_text), float(z_text)
        except ValueError:
            continue
        symbol = ATOMIC_SYMBOLS.get(z_num)
        if symbol is None:
            continue
        symbols.append(symbol)
        coords.append((x, y, z))

    if len(symbols) < 3:
        return None

    stated_formula = _child_text(specie_elem, "formula")
    charge_text = _child_text(specie_elem, "charge")
    try:
        charge = int(charge_text) if charge_text is not None else 0
    except ValueError:
        charge = 0

    return symbols, coords, stated_formula, charge


def main():
    if not os.path.exists(XML_PATH):
        print(f"Missing {XML_PATH} -- see the manual registration step described "
              f"at the top of this script (fetch_nasa_pahdb.py) before running it.")
        return

    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    seen = set()
    n_species = 0
    n_ok = 0
    n_formula_match = 0
    n_formula_mismatch = 0

    with _open_xml_source(XML_PATH) as f:
        for event, elem in ET.iterparse(f, events=("end",)):
            if _local(elem.tag) != "specie":
                continue

            n_species += 1
            parsed = _parse_specie(elem)
            elem.clear()  # free this species' subtree (atoms + vibrational
                          # modes, the bulk of the file's size) once we're
                          # done with it, so memory stays bounded regardless
                          # of how many thousand species the file has.
            if parsed is None:
                continue
            symbols, coords, stated_formula, charge = parsed

            mol = mol_from_geometry(symbols, coords, charge=charge)
            if mol is None:
                continue
            try:
                smi = Chem.MolToSmiles(mol, canonical=True)
            except Exception:
                continue

            seen.add(smi)
            n_ok += 1

            if stated_formula:
                computed = rdMolDescriptors.CalcMolFormula(mol)
                # CalcMolFormula can append a charge suffix (e.g. "C6H6O2+")
                # that PAHdb's plain <formula> field doesn't carry -- strip
                # that before comparing.
                computed_bare = computed.rstrip("+-0123456789") if computed and computed[-1] in "+-" else computed
                if computed_bare == stated_formula or computed == stated_formula:
                    n_formula_match += 1
                else:
                    n_formula_mismatch += 1

            if n_species % 1000 == 0:
                print(f"  processed {n_species:,} species so far, "
                      f"{len(seen):,} unique SMILES, {n_formula_mismatch} formula mismatches")

    print(f"\nScanned {n_species:,} <specie> entries, {n_ok:,} converted to a molecule "
          f"({len(seen):,} unique canonical SMILES after dedup).")
    print(f"Formula cross-check: {n_formula_match:,} matched the file's own <formula> "
          f"field, {n_formula_mismatch:,} didn't (see the formula cross-check note "
          f"in this script's docstring -- not necessarily errors, just flagged).")

    if not seen:
        print("\nNothing converted -- if n_species above is also 0, the schema must "
              "have changed again; if it's nonzero but n_ok stayed 0, open a fresh "
              "issue by re-running inspect_pahdb_specie.py and sending the output.")
        return

    with open(OUT_PATH, "w") as out_f:
        for smi in seen:
            out_f.write(smi + "\n")
    print(f"\nWrote {len(seen):,} unique canonical SMILES to {OUT_PATH} "
          f"(paper reports 3,139 from this source; this file's fuller v4.00 "
          f"database will likely give more).")


if __name__ == "__main__":
    main()
