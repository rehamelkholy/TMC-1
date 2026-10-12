"""Cell-level change logs between consecutive frozen versions of each curated table.

Written for the release. Several intermediate versions were produced by scripts that are no longer
available (see provenance/README.md); the frozen files themselves are all shipped, so every label change
is still traceable by comparing consecutive versions. For each pair of consecutive versions this script
writes provenance/changelogs/<series>_<old>_to_<new>.csv with one row per changed cell:

    row_key, column, old_value, new_value

Rows present in only one version are reported with column "(row)"; columns present in only one version
with column "(column)". Values are compared as strings.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from paths import DATA_FROZEN, PROV_DATA, ROOT, TABLES  # noqa: E402

OUT_DIR = os.path.join(ROOT, "provenance", "changelogs")
os.makedirs(OUT_DIR, exist_ok=True)


def prov(name):
    return os.path.join(PROV_DATA, name)


# series name -> (row-key column, ordered list of (version label, path))
SERIES = {
    "post2021_detections": ("formula", [
        ("v2", prov("post2021_detections_v2.csv")), ("v3", prov("post2021_detections_v3.csv")),
        ("v3_1", prov("post2021_detections_v3_1.csv")), ("v3_2", prov("post2021_detections_v3_2.csv")),
        ("v4", prov("post2021_detections_v4.csv")), ("v5", prov("post2021_detections_v5.csv")),
        ("v6", prov("post2021_detections_v6.csv")),
        ("v7", os.path.join(DATA_FROZEN, "post2021_detections.csv"))]),
    "lee_list_entry_labels": ("id", [
        ("v3_1", prov("lee_list_entry_labels_v3_1.csv")), ("v4", prov("lee_list_entry_labels_v4.csv")),
        ("v5", prov("lee_list_entry_labels_v5.csv"))]),
    "lee_list_entry_labels_purged": ("id", [
        ("v1", prov("lee_list_entry_labels_purged_v1.csv")), ("v4", prov("lee_list_entry_labels_purged_v4.csv")),
        ("v5", prov("lee_list_entry_labels_purged_v5.csv")),
        ("v6", os.path.join(DATA_FROZEN, "lee_list_entry_labels_purged.csv"))]),
    "candidates": ("member_ids", [
        ("v1", prov("candidates_purged_v1.csv")), ("v4", prov("candidates_purged_v4.csv")),
        ("v5", prov("candidates_purged_v5.csv")), ("v6", prov("candidates_purged_v6.csv")),
        ("v7", os.path.join(DATA_FROZEN, "candidates.csv"))]),
    "table1_validation_set": ("Species", [
        ("v5", prov("table1_validation_set_v5.csv")), ("v6", prov("table1_validation_set_v6.csv")),
        ("v7", prov("table1_validation_set_v7.csv")), ("v8", prov("table1_validation_set_v8.csv")),
        ("v9", prov("table1_validation_set_v9.csv")), ("v10", prov("table1_validation_set_v10.csv")),
        ("v12", prov("table1_validation_set_v12.csv")), ("v13", prov("table1_validation_set_v13.csv")),
        ("v14", os.path.join(TABLES, "table1_validation_set.csv"))]),
    "table3_recommended_targets": ("LeeID", [
        ("v5", prov("table3_recommended_targets_v5.csv")), ("v6", prov("table3_recommended_targets_v6.csv")),
        ("v7", prov("table3_recommended_targets_v7.csv")),
        ("v8", os.path.join(TABLES, "table3_recommended_targets.csv"))]),
}


def load(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def diff(a, b, key):
    rows = []
    a = a.set_index(key, drop=False)
    b = b.set_index(key, drop=False)
    cols_a, cols_b = list(a.columns), list(b.columns)
    for c in cols_a:
        if c not in cols_b:
            rows.append(("", "(column)", c, "(absent)"))
    for c in cols_b:
        if c not in cols_a:
            rows.append(("", "(column)", "(absent)", c))
    shared = [c for c in cols_b if c in cols_a and c != key]
    for k in b.index:
        if k not in a.index:
            rows.append((k, "(row)", "(absent)", "(present)"))
            continue
        for c in shared:
            old, new = a.at[k, c], b.at[k, c]
            if old != new:
                rows.append((k, c, old, new))
    for k in a.index:
        if k not in b.index:
            rows.append((k, "(row)", "(present)", "(absent)"))
    return pd.DataFrame(rows, columns=["row_key", "column", "old_value", "new_value"])


summary = []
for name, (key, versions) in SERIES.items():
    for (va, pa), (vb, pb) in zip(versions[:-1], versions[1:]):
        out = diff(load(pa), load(pb), key)
        out_path = os.path.join(OUT_DIR, f"{name}_{va}_to_{vb}.csv")
        out.to_csv(out_path, index=False, lineterminator="\n")
        summary.append((name, va, vb, len(out)))
        print(f"{name}: {va} -> {vb}: {len(out)} changes")
