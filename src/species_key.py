"""
Species key (protocol amendment 4, J1):

    key = CalcMolFormula(mol) + '|' + InChI formula layer + '|' + InChI /c layer
          + '|' + InChI /h layer,

taken from RDKit MolToInchi (default options), all other layers dropped. If MolToInchi
returns empty or the SMILES does not parse, the entry is its own singleton candidate.
No fallback parsing (e.g. sanitize=False retry) is used.
"""
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def species_key_v2(smiles):
    """Returns a string key, or None if the SMILES does not parse or
    MolToInchi returns empty (caller must treat None as 'own singleton
    candidate', per J1)."""
    mol = Chem.MolFromSmiles(smiles)  # default options, sanitize=True
    if mol is None:
        return None
    formula = rdMolDescriptors.CalcMolFormula(mol)
    try:
        inchi = Chem.MolToInchi(mol)  # default options
    except Exception:
        inchi = None
    if not inchi:
        return None
    layers = inchi.split("/")
    if len(layers) < 2:
        return None
    inchi_formula_layer = layers[1]
    c_layer = ""
    h_layer = ""
    for layer in layers[2:]:
        if layer.startswith("c") and not c_layer:
            c_layer = layer
        elif layer.startswith("h") and not h_layer:
            h_layer = layer
        # all other layers (b/t/m/s/q/p/...) dropped, per J1
    return f"{formula}|{inchi_formula_layer}|{c_layer}|{h_layer}"


if __name__ == "__main__":
    import pandas as pd
    from paths import PROV_DATA
    lee = pd.read_csv(f"{PROV_DATA}/lee_list_entry_labels_v3_1.csv")
    keys = lee["smiles"].apply(species_key_v2)
    n_fail = keys.isna().sum()
    print("total:", len(lee), "singleton (parse/InChI failure):", n_fail)
    print(lee.loc[keys.isna(), ["id", "smiles", "formula"]].to_string())
    lee["species_key_v2"] = keys
    print("\ndistinct keys (excluding singletons, which are each unique by construction):",
          lee.loc[keys.notna(), "species_key_v2"].nunique())
