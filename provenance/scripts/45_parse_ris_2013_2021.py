# NON-RUNNABLE PROVENANCE: this script is shipped as a record of how an intermediate frozen file was
# produced. Paths refer to the original working tree (<PROJECT_ROOT>) and inputs not shipped here.
# It is not part of the reproduction pipeline; see provenance/README.md.
"""
Parse ads_tmc1_2013-2021.ris into one row per record, for
literature screening (reusing the parsing logic of 46_parse_ris_2021_2026.py, adapted to a single input file with no v3.2 pre-matching step --
this window's job is to screen EVERY record for ANY species match, not just
new-to-space detections).
"""
import re
import pandas as pd

RIS_PATH = "<PROJECT_ROOT>/ads_tmc1_2013-2021.ris.txt"
OUT_PATH = "<PROJECT_ROOT>/f13_implementation/ads_2013_2021_records_parsed.csv"


def parse_ris(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
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
                if cur_tag in ("TI", "N2") and rec.get(cur_tag):
                    rec[cur_tag][-1] += " " + line.strip()
                elif cur_tag and isinstance(rec.get(cur_tag), str):
                    rec[cur_tag] += " " + line.strip()
        rec["TI"] = " ".join(rec["TI"]).strip()
        rec["N2"] = " ".join(rec["N2"]).strip()
        records.append(rec)
    return records


def extract_arxiv_id(rec):
    sp = rec.get("SP", "") or ""
    m = re.search(r"arXiv:(\d{4}\.\d{4,5})", sp)
    if m:
        return m.group(1)
    for c1 in rec.get("C1", []) or []:
        m = re.search(r"arXiv:(\d{4}\.\d{4,5})", c1)
        if m:
            return m.group(1)
    do = rec.get("DO", "") or ""
    m = re.search(r"arXiv\.(\d{4}\.\d{4,5})", do)
    if m:
        return m.group(1)
    ur = rec.get("UR", "") or ""
    m = re.search(r"arxiv\.org/abs/(\d{4}\.\d{4,5})", ur, re.IGNORECASE)
    if m:
        return m.group(1)
    return None


records = parse_ris(RIS_PATH)
print(f"Parsed {len(records)} records from {RIS_PATH}")

rows = []
for i, r in enumerate(records):
    rows.append({
        "record_id": f"B{i+1:04d}",
        "TY": r.get("TY", ""),
        "title": r.get("TI", ""),
        "authors": "; ".join(r.get("AU", [])),
        "journal": r.get("JO", r.get("JF", "")),
        "year_field": r.get("Y1", ""),
        "arxiv_id": extract_arxiv_id(r),
        "doi": r.get("DO", ""),
        "ads_url": r.get("UR", ""),
        "abstract": r.get("N2", ""),
    })

df = pd.DataFrame(rows)
df.to_csv(OUT_PATH, index=False)
print(f"Wrote {OUT_PATH}: {len(df)} records")
print(df[["record_id", "TY", "year_field"]].groupby("TY").count())
