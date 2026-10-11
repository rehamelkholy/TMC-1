"""
Build the candidate table from lee_list_entry_labels_v3_1.csv (1,501 entries) plus the A7 date
corrections (already applied in post2021_detections_v3_2.csv), per protocol amendment 3 (I1/I2).

Output: a per-candidate table (species-key collapsed) with columns:
  candidate_id, n_members, member_ids, formulas,
  dist, logcd, sigma, mu, mw (means over members),
  charged (any member charge != 0), radical (any member mult > 1),
  positive (any member in_list_by_v3_1), date (earliest parsed date among
  positive members, else NaT), date_raw (the raw string(s) that produced it),
  date_precision ('day'/'month'/'year'/None)
"""
import re
import pandas as pd
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

from paths import PROV_DATA

RDLogger.DisableLog("rdApp.*")

LEE_PATH = f"{PROV_DATA}/lee_list_entry_labels_v3_1.csv"
V32_PATH = f"{PROV_DATA}/post2021_detections_v3_2.csv"

# ---- A7 corrections (from protocol_v1_amendments_2.md, cross-checked against
# post2021_detections_v3_2.csv, which already carries them) ----
A7_FORMULAS = ["C3H+", "C5H+", "H2CCCH+"]


def species_key(smiles):
    """RDKit formula incl. charge, + connectivity/H-count skeleton with bond
    orders, charges and stereo removed (I1, verbatim)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # fallback for SMILES RDKit's default sanitizer rejects (e.g. CH5+,
        # SiH5+ -- hypervalent protonated cations). Partial sanitize so we can
        # still compute a formula and an InChI-based skeleton.
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_KEKULIZE,
            catchErrors=True,
        )
    formula = rdMolDescriptors.CalcMolFormula(mol)  # includes net charge, e.g. "C3H3N+"
    try:
        inchi = Chem.MolToInchi(mol)
    except Exception:
        inchi = None
    if not inchi:
        skeleton = Chem.MolToSmiles(mol, isomericSmiles=False)
    else:
        layers = inchi.split("/")
        keep = []
        for layer in layers[1:]:
            if not layer:
                continue
            tag = layer[0]
            if tag in ("b", "t", "m", "s", "q", "p"):
                continue  # drop stereo (b/t/m/s) and charge (q/p) layers
            keep.append(layer)
        if keep:
            keep = keep[1:]  # drop InChI's own formula layer (redundant with `formula` above)
        skeleton = "/".join(keep)
    return (formula, skeleton)


def parse_one_date(s):
    """Best-effort parse of a single first_public_date string per I2: day-level
    if present, else month-level with day set to 1, else year-level (fallback,
    not expected to be hit post-A7: uses Jan 1 of that year, flagged via
    precision='year')."""
    if not isinstance(s, str) or not s:
        return None, None
    s = s.strip()
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})", s)
    if m:
        y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
        return pd.Timestamp(year=y, month=mo, day=d), "day"
    m = re.match(r"^(\d{4})-(\d{2})\b", s)
    if m:
        y, mo = int(m.group(1)), int(m.group(2))
        return pd.Timestamp(year=y, month=mo, day=1), "month"  # I2: month-only -> 1st of month
    m = re.match(r"^(\d{4})\b", s)
    if m:
        y = int(m.group(1))
        return pd.Timestamp(year=y, month=1, day=1), "year"  # fallback; not expected post-A7
    return None, None


def parse_date_field(field):
    """A first_public_date_of_match field may be semicolon-joined (multiple
    matching detections, e.g. the cis-/trans-CH3CHCHCN tie). Take the earliest
    parsed date among all segments."""
    if not isinstance(field, str) or not field:
        return None, None
    parsed = [parse_one_date(seg) for seg in field.split(";")]
    parsed = [p for p in parsed if p[0] is not None]
    if not parsed:
        return None, None
    parsed.sort(key=lambda p: p[0])
    return parsed[0]


def apply_a7(lee, v32):
    """Overwrite first_public_date_of_match for lee entries matched to
    C3H+/C5H+/H2CCCH+ with the A7-corrected dates from v3.2 (the lee labels
    file was built before A7 was recorded)."""
    lee = lee.copy()
    lee["first_public_date_of_match"] = lee["first_public_date_of_match"].astype("object")
    n_updated = 0
    for formula in A7_FORMULAS:
        new_date = v32.loc[v32["formula"] == formula, "first_public_date"].iloc[0]
        mask = lee["matched_detection_formula"].astype(str).str.contains(
            re.escape(formula), na=False, regex=True
        )
        # guard against substring collisions (e.g. "C3H+" inside a longer formula string)
        def is_exact_member(field, target=formula):
            if not isinstance(field, str):
                return False
            return target in field.split(";")
        mask = lee["matched_detection_formula"].apply(is_exact_member)
        assert mask.sum() == 1, f"expected exactly 1 lee entry matched to {formula}, got {mask.sum()}"
        lee.loc[mask, "first_public_date_of_match"] = new_date
        n_updated += mask.sum()
    return lee, n_updated


def build_candidates():
    lee = pd.read_csv(LEE_PATH)
    v32 = pd.read_csv(V32_PATH)

    lee, n_a7 = apply_a7(lee, v32)

    lee["species_key"] = lee["smiles"].apply(species_key)
    assert lee["species_key"].notna().all(), "species_key failed to build for some rows"

    parsed = lee["first_public_date_of_match"].apply(parse_date_field)
    lee["parsed_date"] = [p[0] for p in parsed]
    lee["date_precision"] = [p[1] for p in parsed]

    # sanity: every positive entry has a parseable date
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

    # I1: missing dipole imputed with the pool median, before any use
    n_missing_mu = cand["mu"].isna().sum()
    if n_missing_mu:
        med = cand["mu"].median()
        cand.loc[cand["mu"].isna(), "mu"] = med
    cand.attrs["n_missing_mu_imputed"] = int(n_missing_mu)
    cand.attrs["n_a7_dates_updated"] = int(n_a7)
    cand.attrs["n_raw_entries"] = len(lee)
    return cand, lee


if __name__ == "__main__":
    cand, lee = build_candidates()
    print("raw entries:", len(lee), " -> candidates:", len(cand))
    print("a7 dates updated:", cand.attrs["n_a7_dates_updated"])
    print("missing mu imputed:", cand.attrs["n_missing_mu_imputed"])
    print("positive candidates:", cand["positive"].sum())
    print(cand[cand["positive"]][["candidate_id", "formulas", "n_members", "date"]].sort_values("date").to_string())
