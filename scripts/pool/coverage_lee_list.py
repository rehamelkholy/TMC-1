"""
Coverage of Lee's 1,510-entry list and the 40 pos_primary in-list detections in the (filtered)
broader candidate pool, using ONE consistent canonicalization convention (full sanitization) on both
sides.

Memory/time note: a full canonical-smiles -> {sources} map for all ~8.1M filtered pool rows would be
prohibitive (groupby().apply(set) runs out of memory); it is only needed for the ~1,500 reference-list
molecules that turn out to match, so it is built AFTER narrowing to matches.
"""
import glob
import os
import time

import pandas as pd

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "external", "pool")
t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


files = sorted(glob.glob(os.path.join(out_dir, "chunk_out", "*.parquet")))
cols = ["canonical_smiles_full", "source", "passes_filters"]
parts = [pd.read_parquet(f, columns=cols) for f in files]
pool = pd.concat(parts, ignore_index=True)
del parts
log(f"Loaded pool (narrow columns): {len(pool)} rows.")

pool_filt = pool[pool["passes_filters"]]
log(f"Filtered pool: {len(pool_filt)} rows.")

pool_canon_series = pool_filt["canonical_smiles_full"].dropna()
pool_canon_set = set(pool_canon_series.unique())
log(f"Distinct canonical SMILES in filtered pool: {len(pool_canon_set)}")


def sources_for(smiles_list):
    """Only look up sources for a SMALL target list of canonical SMILES (the reference-list matches),
    not the whole pool -- keeps this cheap."""
    target = set(smiles_list)
    sub = pool_filt[pool_filt["canonical_smiles_full"].isin(target)]
    return sub.groupby("canonical_smiles_full")["source"].apply(set)


def coverage_report(ref_df, label, id_cols):
    ref = ref_df.copy()
    comparable = ref["canonical_smiles_full"].notna()
    n_comparable = int(comparable.sum())
    n_fail_canon = int((~comparable).sum())
    found = ref["canonical_smiles_full"].isin(pool_canon_set) & comparable
    n_found = int(found.sum())
    log(f"=== {label} ===")
    print(f"Total: {len(ref)}  |  canonicalizable: {n_comparable}  |  "
          f"sanitize-failed (excluded): {n_fail_canon}")
    print(f"Found in pool: {n_found}/{n_comparable} "
          f"({100*n_found/n_comparable:.1f}% of comparable, "
          f"{100*n_found/len(ref):.1f}% of total)")

    print("By class:")
    class_rows = []
    for cls, grp in ref[comparable].groupby("class"):
        f = grp["canonical_smiles_full"].isin(pool_canon_set)
        print(f"  {cls:10s}: {int(f.sum())}/{len(grp)} ({100*f.sum()/len(grp):.1f}%)")
        class_rows.append({"class": cls, "found": int(f.sum()), "total": len(grp)})
    pd.DataFrame(class_rows).to_csv(os.path.join(out_dir, f"{label}_by_class.csv"), index=False)

    canon_to_sources = sources_for(ref.loc[found, "canonical_smiles_full"].tolist())
    src_counts = {}
    for smi in ref.loc[found, "canonical_smiles_full"]:
        for s in canon_to_sources.get(smi, set()):
            src_counts[s] = src_counts.get(s, 0) + 1
    print("Found entries, by contributing pool source (can be >1 per molecule):")
    for s, c in sorted(src_counts.items(), key=lambda x: -x[1]):
        print(f"  {s:12s}: {c}")
    pd.Series(src_counts, name="n_matches").to_csv(os.path.join(out_dir, f"{label}_by_source.csv"))

    not_found_out = ref[~comparable | ~found].copy()
    not_found_out["reason"] = not_found_out.apply(
        lambda r: "sanitize_failed" if pd.isna(r["canonical_smiles_full"]) else "not_in_pool", axis=1)
    keep_cols = [c for c in id_cols + ["smiles", "canonical_smiles_full", "class", "reason"] if c in not_found_out.columns]
    not_found_out[keep_cols].to_csv(os.path.join(out_dir, f"{label}_not_found.csv"), index=False)
    log(f"Wrote {label}_not_found.csv ({len(not_found_out)} rows).")

    return {"label": label, "total": len(ref), "comparable": n_comparable,
            "sanitize_failed": n_fail_canon, "found": n_found,
            "pct_found_of_comparable": round(100 * n_found / n_comparable, 2) if n_comparable else None,
            "pct_found_of_total": round(100 * n_found / len(ref), 2)}


lee1510 = pd.read_csv(os.path.join(out_dir, "lee1510_canonicalized.csv"))
pos40 = pd.read_csv(os.path.join(out_dir, "pos40_canonicalized.csv"))

summary = []
summary.append(coverage_report(lee1510, "lee1510", ["id"]))
summary.append(coverage_report(pos40, "pos40_in_list_detections", ["id", "formula"]))

pd.DataFrame(summary).to_csv(os.path.join(out_dir, "item2_coverage_summary.csv"), index=False)
log("Wrote item2_coverage_summary.csv")

log("DONE.")
