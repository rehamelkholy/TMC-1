"""Re-match every Table 1 (version 13) row against Lee's 1,501 raw entries (rules A/B/C) and compare
with the LeeID/Rule columns already in the table. Read-only on the inputs; writes a check CSV."""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import pandas as pd
from rdkit import Chem

from paths import PROV_DATA, BUILD, RESULTS
from rule_c import rule_c_key

L = pd.read_csv(os.path.join(PROV_DATA, 'lee_list_entry_labels_v3_1.csv'), dtype=str, keep_default_na=False)
t = pd.read_csv(os.path.join(BUILD, 'table1_validation_set_v13.csv'), dtype=str, keep_default_na=False)


def canon(s):
    m = Chem.MolFromSmiles(s)
    if m is None:
        m = Chem.MolFromSmiles(s, sanitize=False)
        if m is None:
            return None
    return Chem.MolToSmiles(m)


def nostereo(s):
    m = Chem.MolFromSmiles(s)
    if m is None:
        return None
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)


L['A'] = L.smiles.map(canon)
L['B'] = L.smiles.map(nostereo)
L['C'] = L.smiles.map(lambda s: rule_c_key(s)[0])
rows = []
for _, r in t.iterrows():
    a, b, c = canon(r.SMILES), nostereo(r.SMILES), rule_c_key(r.SMILES)[0]
    mA = L[L.A == a].id.tolist()
    mB = L[L.B == b].id.tolist()
    mC = L[L.C == c].id.tolist()
    rule = 'A' if mA else ('B' if mB else ('C' if mC else ''))
    ids = mA or mB or mC
    rows.append((r.Species, r.LeeID, r.Rule, rule, '|'.join(ids), len(mC)))
d = pd.DataFrame(rows, columns=['Species', 'LeeID', 'Rule', 'calcRule', 'calcIDs', 'nC'])
d['has'] = d.LeeID != ''
d['calc_has'] = d.calcIDs != ''
print('existing with/without', d.has.sum(), (~d.has).sum(),
      '| recomputed with/without', d.calc_has.sum(), (~d.calc_has).sum())
print('\nrows where existing != recomputed:')
print(d[(d.LeeID != d.calcIDs) | (d.Rule != d.calcRule)].to_string())
d.to_csv(os.path.join(RESULTS, 'table1_lee_rematch_check.csv'), index=False)
