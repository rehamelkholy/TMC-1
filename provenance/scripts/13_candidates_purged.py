# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Candidate table from the purged raw universe: identical to src/candidates.py (same J1 species key, same A7
date correction, same feature aggregation); the ONLY difference is the raw universe it reads -- the purged
1,489-row lee_list_entry_labels_purged.csv instead of the original 1,501-row
lee_list_entry_labels_v3_1.csv.
"""
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "<PROJECT_ROOT>/f4a_implementation")
from step2_build_candidates import V32_PATH, apply_a7, parse_date_field  # noqa: E402

sys.path.insert(0, "<PROJECT_ROOT>/f4c_implementation")
from species_key_v2 import species_key_v2  # noqa: E402

PURGED_LEE_PATH = "<PROJECT_ROOT>/f6_implementation/lee_list_entry_labels_purged.csv"


def build_candidates_purged():
    lee = pd.read_csv(PURGED_LEE_PATH)
    v32 = pd.read_csv(V32_PATH)

    lee, n_a7 = apply_a7(lee, v32)

    raw_keys = lee["smiles"].apply(species_key_v2)
    singleton_ids = lee.loc[raw_keys.isna(), "id"].tolist()
    final_keys = raw_keys.copy()
    for idx in lee.index[raw_keys.isna()]:
        final_keys.loc[idx] = f"__SINGLETON__{lee.loc[idx, 'id']}"
    lee["species_key"] = final_keys

    parsed = lee["first_public_date_of_match"].apply(parse_date_field)
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
    cand, lee = build_candidates_purged()
    print("purged raw entries:", len(lee), " -> candidates:", len(cand))
    print("singleton ids:", cand.attrs["singleton_ids"])
    print("positive candidates:", cand["positive"].sum())
