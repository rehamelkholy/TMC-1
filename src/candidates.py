"""
Candidate table under the J1 species key (protocol amendment 4). Everything except the
species-key function is reused from candidate_table.py (A7 date correction, date parsing,
feature aggregation, dipole imputation).
"""
import numpy as np
import pandas as pd
from candidate_table import (  # noqa: E402  -- reused unchanged from candidate_table.py
    LEE_PATH, V32_PATH, apply_a7, parse_date_field,
)
from species_key import species_key_v2


def build_candidates_v2():
    lee = pd.read_csv(LEE_PATH)
    v32 = pd.read_csv(V32_PATH)

    lee, n_a7 = apply_a7(lee, v32)  # reused unchanged from candidate_table.py

    raw_keys = lee["smiles"].apply(species_key_v2)
    # J1: parse/InChI failure -> the entry is its own singleton candidate.
    # Give each such row a key unique to itself so it can never collide with
    # another entry (including another singleton with the same raw SMILES,
    # of which there are none here, but this is the literal, safe reading).
    singleton_ids = lee.loc[raw_keys.isna(), "id"].tolist()
    final_keys = raw_keys.copy()
    for idx in lee.index[raw_keys.isna()]:
        final_keys.loc[idx] = f"__SINGLETON__{lee.loc[idx, 'id']}"
    lee["species_key"] = final_keys

    parsed = lee["first_public_date_of_match"].apply(parse_date_field)  # reused unchanged from candidate_table.py
    lee["parsed_date"] = [p[0] for p in parsed]
    lee["date_precision"] = [p[1] for p in parsed]

    pos_no_date = lee[(lee["in_list_by_v3_1"] == True) & (lee["parsed_date"].isna())]
    assert len(pos_no_date) == 0, f"positive entries with no parseable date: {pos_no_date['id'].tolist()}"

    groups = []
    for key, g in lee.groupby("species_key", sort=False):
        positive = bool(g["in_list_by_v3_1"].any())
        date = g.loc[g["in_list_by_v3_1"] == True, "parsed_date"].min() if positive else pd.NaT
        groups.append({
            "species_key": key,
            "n_members": len(g),
            "member_ids": ";".join(str(x) for x in g["id"]),
            "formulas": ";".join(sorted(set(g["formula"]))),
            "smiles_members": ";".join(g["smiles"]),
            "dist": g["dist"].mean(),
            "logcd": g["logcd"].mean(),
            "sigma": g["sigma"].mean(),
            "mu": g["mu"].mean(),
            "mw": g["mw"].mean(),
            "charged": bool((g["charge"] != 0).any()),
            "radical": bool((g["mult"] > 1).any()),
            "positive": positive,
            "date": date,
        })
    cand = pd.DataFrame(groups).reset_index(drop=True)
    cand["candidate_id"] = np.arange(len(cand))

    n_missing_mu = cand["mu"].isna().sum()
    if n_missing_mu:
        med = cand["mu"].median()
        cand.loc[cand["mu"].isna(), "mu"] = med

    cand.attrs["n_missing_mu_imputed"] = int(n_missing_mu)
    cand.attrs["n_a7_dates_updated"] = int(n_a7)
    cand.attrs["n_raw_entries"] = len(lee)
    cand.attrs["singleton_ids"] = singleton_ids
    cand.attrs["n_singletons"] = len(singleton_ids)
    return cand, lee


if __name__ == "__main__":
    cand, lee = build_candidates_v2()
    print("raw entries:", len(lee), " -> candidates (v2, J1 key):", len(cand))
    print("singleton (parse/InChI-failure) ids:", cand.attrs["singleton_ids"])
    print("n_singletons:", cand.attrs["n_singletons"])
    print("positive candidates:", cand["positive"].sum())
