# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Split the 329 ADS records that did NOT match a v3.2
detection by DOI/arXiv id into batches for parallel title+abstract
screening (irrelevant vs candidate new TMC-1 species detection). Batches
of ~28 records each -> ~12 batches, sized so each screening pass stays
comfortably bounded.
"""
import json
import math

import pandas as pd

ROOT = "<PROJECT_ROOT>"
df = pd.read_csv(f"{ROOT}/f8_implementation/ads_records_parsed.csv")
unmatched = df[df["decision_prelim"] != "already_in_v3.2"].copy()
unmatched = unmatched.fillna("")
print(f"Records needing title+abstract screening: {len(unmatched)}")
print(f"Records with empty abstract: {(unmatched['abstract'].str.strip() == '').sum()}")

BATCH_SIZE = 28
records = unmatched.to_dict("records")
n_batches = math.ceil(len(records) / BATCH_SIZE)
for i in range(n_batches):
    batch = records[i * BATCH_SIZE:(i + 1) * BATCH_SIZE]
    slim = [
        {
            "record_id": r["record_id"],
            "title": r["title"],
            "authors": r["authors"],
            "journal": r["journal"],
            "year_field": r["year_field"],
            "arxiv_id": r["arxiv_id"],
            "doi": r["doi"],
            "abstract": r["abstract"],
        }
        for r in batch
    ]
    with open(f"{ROOT}/f8_implementation/screen_batch_{i+1:02d}.json", "w") as f:
        json.dump(slim, f, indent=2)
print(f"Wrote {n_batches} batch files.")
