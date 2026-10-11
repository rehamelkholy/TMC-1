# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Parse both RIS exports, dedupe, and
match against post2021_detections_v3_2.csv by DOI/arXiv id so the "already
in v3.2" bucket is deterministic (not a text-similarity guess). Writes
ads_records_parsed.csv (one row per unique ADS record, source flags,
v3.2 match) for the next (literature-screening) step to consume.
"""
import csv
import json
import re
import sys

import pandas as pd

ROOT = "<PROJECT_ROOT>"


def parse_ris(path, source_tag):
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    # records separated by ER  - lines
    raw_records = re.split(r"\nER  -\s*\n", text)
    records = []
    for raw in raw_records:
        raw = raw.strip("\n")
        if not raw.strip():
            continue
        lines = raw.split("\n")
        rec = {"AU": [], "KW": [], "N2": [], "TI": [], "C1": []}
        cur_tag = None
        for line in lines:
            m = re.match(r"^([A-Z][A-Z0-9])  - (.*)$", line)
            if m:
                tag, val = m.group(1), m.group(2)
                cur_tag = tag
                if tag in ("AU", "KW", "C1"):
                    rec[tag].append(val)
                elif tag in ("TI", "N2"):
                    rec[tag].append(val)
                else:
                    rec[tag] = val
            else:
                # continuation line (wrapped abstract/title), append to last tag
                if cur_tag in ("TI", "N2") and rec.get(cur_tag):
                    rec[cur_tag][-1] += " " + line.strip()
                elif cur_tag and isinstance(rec.get(cur_tag), str):
                    rec[cur_tag] += " " + line.strip()
        rec["TI"] = " ".join(rec["TI"]).strip()
        rec["N2"] = " ".join(rec["N2"]).strip()
        rec["source_file"] = source_tag
        records.append(rec)
    return records


def extract_arxiv_id(rec):
    # Preprint-type records (TY - Preprint): SP field is "arXiv:2609.24285".
    sp = rec.get("SP", "") or ""
    m = re.search(r"arXiv:(\d{4}\.\d{4,5})", sp)
    if m:
        return m.group(1)
    # Published journal-type records (TY - JOUR): the arXiv id instead
    # lives in a C1 note field, "eprint: arXiv:2501.05125" (SP holds the
    # journal's own article id, e.g. "L20", and is NOT the arXiv id here).
    for c1 in rec.get("C1", []) or []:
        m = re.search(r"arXiv:(\d{4}\.\d{4,5})", c1)
        if m:
            return m.group(1)
    # try DO field for an arxiv-DOI (10.48550/arXiv.XXXX.XXXXX)
    do = rec.get("DO", "") or ""
    m = re.search(r"arXiv\.(\d{4}\.\d{4,5})", do)
    if m:
        return m.group(1)
    # last resort: UR field sometimes carries an arXiv abs URL directly
    ur = rec.get("UR", "") or ""
    m = re.search(r"arxiv\.org/abs/(\d{4}\.\d{4,5})", ur, re.IGNORECASE)
    if m:
        return m.group(1)
    return None


def normalize_arxiv_id(aid):
    """Right-pad the post-dot digit string to 5 digits (the standard width
    for all arXiv ids in our 2020-2026 date range). This repairs a
    pre-existing data-quality issue found in post2021_detections_v3_2.csv
    itself: its arxiv_id column was at some point round-tripped through a
    float64 pandas column, and a float silently drops trailing zeros
    (2410.00670 == 2410.0067 as a *number*, so the on-disk CSV literally
    stores "2410.0067", four digits, for the three cyanopyrene rows).
    Right-padding with zeros exactly undoes this: floats only ever drop
    TRAILING zeros (never leading/interior digits), so padding the
    fractional part back out to 5 digits reconstructs the true id in
    every case, not just this one. Applied to both v3.2's stored ids and
    the ADS-extracted ids so the comparison is on equal footing; the
    v3.2 CSV file itself is never modified."""
    if aid is None:
        return None
    m = re.match(r"^(\d{2,4})\.(\d{1,5})$", aid.strip())
    if not m:
        return aid.strip()
    yymm, frac = m.groups()
    return f"{yymm}.{frac.ljust(5, '0')}"


def normalize_doi(doi):
    if not doi:
        return None
    return doi.strip().lower().rstrip(".")


def normalize_title(t):
    t = re.sub(r"<[^>]+>", "", t or "")  # strip HTML-ish tags like <SUB>
    t = re.sub(r"[^a-z0-9]+", " ", t.lower())
    return t.strip()


tmc1_records = parse_ris(f"{ROOT}/f8_implementation/ads_tmc1_2021-2026.ris", "tmc1_query")
groups_records = parse_ris(f"{ROOT}/f8_implementation/ads_groups_2021-2026.ris", "groups_query")
print(f"Parsed {len(tmc1_records)} records from tmc1 query, {len(groups_records)} from groups query", file=sys.stderr)

all_records = tmc1_records + groups_records
for r in all_records:
    r["arxiv_id"] = normalize_arxiv_id(extract_arxiv_id(r))
    r["doi_norm"] = normalize_doi(r.get("DO"))
    r["title_norm"] = normalize_title(r.get("TI"))

# --- Dedupe: group by (arxiv_id) else (doi_norm) else (title_norm) ---
groups = {}
order = []
for r in all_records:
    if r["arxiv_id"]:
        key = ("arxiv", r["arxiv_id"])
    elif r["doi_norm"]:
        key = ("doi", r["doi_norm"])
    else:
        key = ("title", r["title_norm"])
    if key not in groups:
        groups[key] = []
        order.append(key)
    groups[key].append(r)

print(f"Total raw records: {len(all_records)}; unique after dedupe: {len(groups)}", file=sys.stderr)

unique_rows = []
for i, key in enumerate(order):
    members = groups[key]
    sources = sorted(set(m["source_file"] for m in members))
    primary = members[0]
    unique_rows.append({
        "record_id": f"R{i+1:04d}",
        "dedupe_key_type": key[0],
        "dedupe_key_value": key[1],
        "n_duplicate_members": len(members),
        "source_files": ";".join(sources),
        "TY": primary.get("TY", ""),
        "title": primary.get("TI", ""),
        "authors": "; ".join(primary.get("AU", [])),
        "journal": primary.get("JO", primary.get("JF", "")),
        "year_field": primary.get("Y1", ""),
        "arxiv_id": primary.get("arxiv_id"),
        "doi": primary.get("DO", ""),
        "ads_url": primary.get("UR", ""),
        "abstract": primary.get("N2", ""),
    })

df = pd.DataFrame(unique_rows)
df.to_csv(f"{ROOT}/f8_implementation/ads_records_parsed.csv", index=False)
print(f"Wrote ads_records_parsed.csv: {len(df)} unique records", file=sys.stderr)

# --- Match against post2021_detections_v3_2.csv by arxiv_id / doi ---
# IMPORTANT: build keys only from genuinely-non-null values (checked with
# pd.notna on the RAW column, before any str conversion) -- astype(str)
# on a float64 column with NaNs does not reliably turn every NaN into the
# literal string "nan" (a handful survive as actual float('nan') objects
# after .str.strip()), and float('nan') is truthy and != "nan", so a naive
# `if aid and aid != "nan"` filter lets stray NaN floats slip in as a
# dict key -- which then spuriously "matches" every OTHER NaN-arxiv-id
# row via `aid in dict`, corrupting nearly every unmatched record. Fixed
# by testing definite non-nullness on the original numeric/str column.
v32 = pd.read_csv(f"{ROOT}/f3_validation_freeze/post2021_detections_v3_2.csv")

# Compute cleaned id/doi as plain Python dicts keyed by v32's row position,
# never as a pandas Series -- a Series built from a list containing None
# (via .apply, or even a Python list assigned to a column) silently turns
# some/all of those Nones back into float('nan') on construction, which
# reintroduces exactly the bug this rewrite fixes. Plain dict/list only.
v32_arxiv_clean = {}  # row position -> str or None
v32_doi_clean = {}
for pos in range(len(v32)):
    raw_aid = v32["arxiv_id"].iat[pos]
    v32_arxiv_clean[pos] = None if pd.isna(raw_aid) else normalize_arxiv_id(str(raw_aid).strip())
    raw_doi = v32["doi"].iat[pos]
    v32_doi_clean[pos] = None if pd.isna(raw_doi) else normalize_doi(str(raw_doi))

arxiv_to_species = {}
doi_to_species = {}
for pos in range(len(v32)):
    formula = v32["formula"].iat[pos]
    aid = v32_arxiv_clean[pos]
    if aid is not None:
        arxiv_to_species.setdefault(aid, []).append(formula)
    doi = v32_doi_clean[pos]
    if doi is not None:
        doi_to_species.setdefault(doi, []).append(formula)

assert not any(isinstance(k, float) for k in arxiv_to_species), "nan leaked into arxiv_to_species"
assert not any(isinstance(k, float) for k in doi_to_species), "nan leaked into doi_to_species"

matched_species_col = []
match_basis_col = []
for _, row in df.iterrows():
    aid_raw = row["arxiv_id"]
    aid = aid_raw if isinstance(aid_raw, str) and aid_raw.strip() else None
    doi_raw = row["doi"]
    doi = normalize_doi(doi_raw) if isinstance(doi_raw, str) and doi_raw.strip() else None
    species = None
    basis = None
    if aid is not None and aid in arxiv_to_species:
        species = arxiv_to_species[aid]
        basis = "arxiv_id"
    elif doi is not None and doi in doi_to_species:
        species = doi_to_species[doi]
        basis = "doi"
    matched_species_col.append(";".join(species) if species else "")
    match_basis_col.append(basis or "")

df["v32_matched_species"] = matched_species_col
df["v32_match_basis"] = match_basis_col
df["decision_prelim"] = df["v32_matched_species"].apply(lambda x: "already_in_v3.2" if x else "")

df.to_csv(f"{ROOT}/f8_implementation/ads_records_parsed.csv", index=False)

n_matched = (df["decision_prelim"] == "already_in_v3.2").sum()
n_unmatched = len(df) - n_matched
print(f"Matched to v3.2 by DOI/arXiv id: {n_matched} records", file=sys.stderr)
print(f"Remaining unmatched (need literature screening): {n_unmatched} records", file=sys.stderr)

# --- Step 4 prep: which v3.2 species/entries have NO matching ADS record at all? ---
v32_arxiv_matched = set()
v32_doi_matched = set()
for _, row in df.iterrows():
    if row["v32_match_basis"] == "arxiv_id":
        v32_arxiv_matched.add(row["arxiv_id"])
    elif row["v32_match_basis"] == "doi":
        v32_doi_matched.add(normalize_doi(row["doi"]))

found_flags = []
for pos in range(len(v32)):
    aid = v32_arxiv_clean[pos]
    doi = v32_doi_clean[pos]
    found = (aid is not None and aid in v32_arxiv_matched) or (doi is not None and doi in v32_doi_matched)
    found_flags.append(found)
v32["found_in_ads_export"] = found_flags
missing = v32[~v32["found_in_ads_export"]]
missing_out = missing[["formula", "first_public_date", "arxiv_id", "doi", "journal_ref", "authors"]]
missing_out.to_csv(f"{ROOT}/f8_implementation/v32_species_missing_from_ads_export.csv", index=False)
print(f"v3.2 species/entries with NO matching ADS record (by DOI/arXiv id): {len(missing)} "
      f"(of {len(v32)} total rows / {v32['formula'].nunique()} unique formulas)", file=sys.stderr)
