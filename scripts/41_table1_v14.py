"""
Table 1 version 13 (101 rows) -> version 14 (101 rows, no data changes).
Fills LeeID/Rule for OH, HSCN, HCS, C7H, NCO, CH3SH from the re-match of every Table 1 SMILES
against Lee's 1,501 raw entries (42_table1_lee_rematch_check.py; all six match rule A, exact canonical
SMILES). LeeID = Lee raw-entry id (member_ids), NOT the candidate number. OH Ref -> bibcode only.
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import hashlib
import pandas as pd
from astropy.table import Table

from paths import BUILD, TABLES

LEE_2021_PUBLISHED = "2021-08-12"
LEE_NOTE = "Reported before Lee et al. (2021) was published"
SANCTIONED_NOTE_PHRASES = {
    "Reported before Lee et al. (2021) was published",
    "Not in Table C1 of Toru Shay et al. (2025)",
    "SMILES corrected",
    "First detection toward TMC-1; previously detected in the Horsehead",
    "Earliest TMC-1 report located",
}


def pad_to_day(d):
    if len(d) == 4:
        return d + "-01-01"
    if len(d) == 7:
        return d + "-01"
    return d


def sha256_of(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def note_is_sanctioned(note):
    if note == "":
        return True
    if note == "First detection toward TMC-1; previously detected in the Horsehead":
        return True
    parts = note.split("; ")
    return all(p in SANCTIONED_NOTE_PHRASES for p in parts)


t11 = pd.read_csv(os.path.join(BUILD, "table1_validation_set_v13.csv"),
                   keep_default_na=False, dtype=str)
assert len(t11) == 101

FILL = {"OH": ("1365", "A"), "HSCN": ("543", "A"), "HCS": ("517", "A"),
        "C7H": ("1372", "A"), "NCO": ("533", "A"), "CH3SH": ("328", "A")}
t14 = t11.copy()
assert len(t14) == 101
changed_cells = []
for sp, (lid, rule) in FILL.items():
    idx = t14.index[t14["Species"] == sp]
    assert len(idx) == 1
    i = idx[0]
    assert t14.at[i, "LeeID"] == "" and t14.at[i, "Rule"] == ""
    changed_cells.append((sp, "LeeID", "", lid)); changed_cells.append((sp, "Rule", "", rule))
    t14.at[i, "LeeID"] = lid; t14.at[i, "Rule"] = rule
i = t14.index[t14["Species"] == "OH"][0]
old = t14.at[i, "Ref"]
assert old == "2000A&A...353.1065H (ADS bibcode; no DOI located)"
t14.at[i, "Ref"] = "2000A&A...353.1065H"
changed_cells.append(("OH", "Ref", old, "2000A&A...353.1065H"))
print("Changed cells:")
for c in changed_cells: print(" ", c)
assert len(changed_cells) == 13
n_with = (t14["LeeID"] != "").sum()
print("LeeID with/without:", n_with, len(t14) - n_with)
assert (n_with, len(t14) - n_with) == (64, 37)

note_change_log = []
for i, r in t14.iterrows():
    padded = pad_to_day(r["Date"])
    should_have_note = padded < LEE_2021_PUBLISHED
    old_note = r["Note"]
    parts = [p for p in old_note.split("; ") if p]
    has_note = LEE_NOTE in parts

    if should_have_note and not has_note:
        new_parts = [LEE_NOTE] + parts
    elif not should_have_note and has_note:
        new_parts = [p for p in parts if p != LEE_NOTE]
    else:
        new_parts = parts

    new_note = "; ".join(new_parts)
    if new_note != old_note:
        note_change_log.append({"Species": r["Species"], "Column": "Note",
                                 "old": old_note, "new": new_note})
        t14.at[i, "Note"] = new_note

print(f"Note re-derivation: {len(note_change_log)} cells changed.")
for c in note_change_log:
    print(" -", c)

violations = [r["Species"] for _, r in t14.iterrows() if not note_is_sanctioned(r["Note"])]
assert not violations, f"Note column has non-controlled-vocabulary text: {violations}"

out_dir = TABLES
t14.to_csv(f"{out_dir}/table1_validation_set.csv", index=False)


def ascii_only(s):
    s = str(s) if s is not None else s
    assert s.isascii(), f"non-ASCII content found: {s!r}"
    return s


t14ascii = t14.copy()
for c in t14ascii.columns:
    t14ascii[c] = t14ascii[c].apply(ascii_only)

at = Table.from_pandas(t14ascii)
DESCR = {
    "Species": "Species formula as used in this paper (ASCII plain formula)",
    "SMILES": "SMILES string used for candidate matching (RDKit canonical, as parsed)",
    "Date": "First public date (YYYY-MM-DD if day-level known, else YYYY-MM or YYYY)",
    "DateSrc": "Date source: arXiv (v1 preprint), accepted (journal acceptance), or online (journal online pub)",
    "LeeID": "Lee et al. (2021) Zenodo raw-entry id(s) matching this species; blank if not in that list",
    "Rule": "Strongest match rule to the Lee et al. (2021) list: A=exact canonical SMILES, B=stereo-free canonical SMILES, C=InChI skeleton; blank if not in that list",
    "Ref": "Reference identifier: DOI, else ADS bibcode or arXiv identifier (identifier only, no explanatory text)",
    "Note": "Controlled-vocabulary tag(s) only, joined with '; ' where more than one applies; blank if none",
}
for c in at.colnames:
    at[c].unit = "---"
    at[c].description = DESCR[c]

mrt_path = f"{out_dir}/table1_validation_set.mrt"
at.write(mrt_path, format="ascii.mrt", overwrite=True)

MRT_FILL = (
    "(to be added)",
    "(to be added)",
    "Validation set of post-2021 TMC-1 species detections (101 rows)",
)
lines = open(mrt_path).read().split("\n")
filled = 0
for i, line in enumerate(lines):
    if line == "Title:":
        lines[i] = "Title: " + MRT_FILL[0]; filled += 1
    elif line == "Authors:":
        lines[i] = "Authors: " + MRT_FILL[1]; filled += 1
    elif line == "Table:":
        lines[i] = "Table: " + MRT_FILL[2]; filled += 1
assert filled == 3
with open(mrt_path, "w") as f:
    f.write("\n".join(lines))

for path in [mrt_path, f"{out_dir}/table1_validation_set.csv"]:
    content = open(path, "rb").read()
    assert content.decode("ascii")
    print(f"\n{path} SHA-256: {sha256_of(path)}")

print(f"\nTotal rows: {len(t14)} (unchanged)")
print("\nNew rows:")
print(t14[t14["Species"].isin(["C7H","NCO","CH3SH"])].to_string(index=False))
