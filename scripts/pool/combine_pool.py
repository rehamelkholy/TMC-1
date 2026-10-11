"""
Combine every source fetched by fetch_qm9.py / fetch_zinc.py / fetch_pubchem.py / fetch_pcba.py /
fetch_kida.py / fetch_nasa_pahdb.py into one candidate pool, mirroring umda's own combination logic
(umda/scripts/pool_smiles.py): concatenate, canonicalize, drop exact duplicate canonical SMILES, drop
anything containing an element on the same exclude list the original authors used, compute molecular
weight.

Not included: umda's small hand-picked `hand_pah.smi` and `yalamanchi_2020.smi` files (present in
umda's repo code but not cited in Lee et al. 2021 and not independently obtainable). To add a
hand-curated addendum, drop a `.smi` file (one SMILES per line) into external/pool/raw/extra_*.smi;
this script picks it up automatically (see the EXTRA_GLOB pattern below).

Run AFTER whichever of the fetch_*.py scripts should be included: any source whose output file does
not exist is skipped with a note, so a smaller pool can be built from a subset of the sources.

MW cutoff: candidates heavier than MW_MAX_DA = 500 are dropped, uniformly across every source. The
value follows the size distribution of the post-2021 TMC-1 detections of the original validation set
(37-325 amu, median 69, 95th percentile 227): 500 is generous headroom above the largest known
detection (cyanocoronene, 325 amu), and excludes an astrochemically implausible mass range.

Performance and robustness (the pool has ~14M raw entries):
  - Canonicalization and the molecular-weight step run across a multiprocessing.Pool
    (N_WORKERS = cpu_count - 1). _canon_smiles() and _mw_of() are module-level functions so they can
    be pickled under the "spawn" start method (Windows). Progress is printed every 500,000 molecules.
  - RDKit's kekulizer can take pathologically long on very large fused aromatic systems (NASA PAHdb
    contains several). Chunks are therefore submitted through a sliding window of pool.apply_async()
    calls (one slot per worker), each with its OWN stall deadline: a chunk that does not return in time
    is written off on its own (its items become None and are dropped by dropna()), while every other
    chunk keeps being processed. A wedged worker cannot be cancelled, so each stall costs one worker of
    parallelism for the rest of the step. The default chunksize (25) bounds the collateral loss per
    pathological entry; if every worker gets wedged, the remaining items are abandoned in one bulk step.

Source attribution: SOURCES are processed in the order qm9, zinc, pubchem, pcba, kida, nasa_pahdb and
the deduplication keeps the FIRST-seen row per canonical SMILES, so a molecule present in several
sources is attributed to the earliest one. Per-source counts therefore understate later sources
(notably KIDA) without any molecule being missing from the pool.

Output: external/pool/combined_candidate_pool.csv (columns: smiles, source, mw)
"""
import glob
import itertools
import multiprocessing as mp
import os
import re
import time

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

# RDKit logs a "Can't kekulize mol" / "not removing hydrogen atom without neighbors" / SMILES-parse-error
# line straight to the terminal (bypassing Python's own stdout) for every candidate it can't fully parse.
# This is expected at the scale of a ~14M-entry pool from public databases (_canon_smiles() returns None
# for these and dropna() drops them a few lines later, so nothing is silently miscounted), and printing
# millions of un-buffered lines can slow a run down. Silenced here; comment out to debug what a specific
# new source is failing on.
RDLogger.DisableLog("rdApp.*")

POOL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "external", "pool")
RAW_DIR = os.path.join(POOL_DIR, "raw")
OUT_PATH = os.path.join(POOL_DIR, "combined_candidate_pool.csv")

MW_MAX_DA = 500.0  # see "MW CUTOFF" in the module docstring above

# leave one core free for the OS/other work
N_WORKERS = max(1, (os.cpu_count() or 4) - 1)

# how long _parallel_map_with_progress waits on ONE chunk with no result before concluding that THAT
# chunk (not the whole run) is stuck rather than just slow, and abandoning it
STALL_TIMEOUT_S = 120

SOURCES = {
    "qm9": os.path.join(RAW_DIR, "qm9_smiles.smi"),
    "zinc": os.path.join(RAW_DIR, "zinc_smiles.smi"),
    "pubchem": os.path.join(RAW_DIR, "pubchem_smiles.smi"),
    "pcba": os.path.join(RAW_DIR, "pcba_smiles.smi"),
    "kida": os.path.join(RAW_DIR, "kida_smiles.smi"),
    "nasa_pahdb": os.path.join(RAW_DIR, "pahdb_smiles.smi"),
}
EXTRA_GLOB = os.path.join(RAW_DIR, "extra_*.smi")

# same exclusion list as umda's scripts/pool_smiles.py, verbatim
EXCLUDE_LIST = [
    "P", "Ho", "As", "Zr", "Sn", "V", "Au", "Br", "F", "Re", "U", "Al", "In",
    "Cl", "Ag", "Ce", "Hg", "Ta", "Te", "Ni", "Be", "Li", "Sb", "Bi", "Mn",
    "Co", "Tm", "Si", "Mo", "Se",
]


def _canon_smiles(smi):
    """Canonicalize one SMILES string, or return None if RDKit can't parse it. Module-level (not a
    closure inside main()) so it can be pickled and sent to worker processes."""
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol, canonical=True) if mol else None


def _mw_of(canonical_smi):
    """Exact molecular weight of an already-canonicalized SMILES string. Module-level for the same
    picklability reason as _canon_smiles."""
    mol = Chem.MolFromSmiles(canonical_smi, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    return Descriptors.ExactMolWt(mol)


def _chunked(iterable, size):
    it = iter(iterable)
    while True:
        batch = list(itertools.islice(it, size))
        if not batch:
            return
        yield batch


def _indexed_call_chunk(fn_and_chunk):
    """Run fn over one chunk of (index, item) pairs, returning a list of (index, result) pairs.
    Module-level (picklable) wrapper, so results can be placed back at the right position in the output
    list regardless of the order in which chunks complete. Chunks are submitted individually via
    apply_async rather than through pool.imap_unordered(fn, items, chunksize=N): with chunksize > 1 that
    call returns a plain generator, and a single stalled chunk's timeout discards every other chunk
    still outstanding. apply_async's per-call AsyncResult objects let each chunk be timed and abandoned
    independently of every other chunk."""
    fn, chunk = fn_and_chunk
    return [(i, fn(item)) for i, item in chunk]


def _parallel_map_with_progress(pool, fn, items, label, chunksize=25,
                                 report_every=500_000,
                                 stall_timeout_s=STALL_TIMEOUT_S):
    """Runs fn over items across pool, printing progress every `report_every` items (a silent
    hour-long step is indistinguishable from a hung one).

    Chunks are submitted through a sliding window of pool.apply_async() calls, one slot per worker
    (N_WORKERS), refilled as each slot's chunk completes or times out:

    1. Each in-flight chunk has its OWN independent deadline (stall_timeout_s). A chunk that does not
       return in time is abandoned on its own -- its items are left as None in the results list (same
       treatment as any other unparseable candidate, dropped later by dropna()) -- while every other
       chunk, in flight or not yet submitted, is unaffected and keeps being processed normally.
    2. A worker wedged on a genuinely stuck item stays wedged for the rest of this step --
       multiprocessing.Pool cannot cancel a running task, only terminate the whole pool. Each stall
       permanently costs one worker's worth of parallelism (not correctness).
    3. chunksize defaults to 25. Since a pathological item can take its whole chunk down with it, the
       collateral loss per pathological entry is bounded by chunksize regardless of the number of
       workers, and a small source (e.g. PAHdb's 3,419 entries) survives largely intact even if it
       contains several pathological entries.
    4. Submission is capped at N_WORKERS - (workers already permanently wedged) concurrent chunks.
       If EVERY worker eventually gets wedged, the remaining not-yet-submitted items are abandoned
       immediately in one bulk step rather than each separately waiting out stall_timeout_s against a
       pool with zero free workers (which would degrade to roughly stall_timeout_s * remaining_chunks).
    """
    n = len(items)
    if n == 0:
        return []
    results = [None] * n
    t0 = time.time()
    next_report = report_every
    done = 0
    n_abandoned = 0
    n_stalled_chunks = 0
    n_wedged_workers = 0
    n_chunks_total = -(-n // chunksize)  # ceil division

    chunks_iter = _chunked(enumerate(items), chunksize)
    in_flight = []  # list of [AsyncResult, chunk, submit_time]

    def _capacity():
        # shrinks by one for every permanently-wedged worker, so no more chunks are kept in flight than
        # there are usable worker processes left.
        return N_WORKERS - n_wedged_workers

    def _submit_next():
        if len(in_flight) >= _capacity():
            return False
        chunk = next(chunks_iter, None)
        if chunk is None:
            return False
        ar = pool.apply_async(_indexed_call_chunk, ((fn, chunk),))
        in_flight.append([ar, chunk, time.time()])
        return True

    while _submit_next():
        pass

    # Poll granularity for the wait below. Blocking on AsyncResult.wait() returns as soon as that result is
    # ready, whereas a flat time.sleep(POLL_INTERVAL_S) would impose a floor of POLL_INTERVAL_S per chunk
    # (hundreds of thousands of chunks at chunksize=25) regardless of how fast the actual work is.
    POLL_INTERVAL_S = 0.05

    while in_flight:
        in_flight[0][0].wait(timeout=POLL_INTERVAL_S)
        for slot in list(in_flight):
            ar, chunk, submit_time = slot
            timed_out = (time.time() - submit_time) > stall_timeout_s
            if not ar.ready() and not timed_out:
                continue

            in_flight.remove(slot)
            if ar.ready():
                try:
                    chunk_results = ar.get()
                    for i, result in chunk_results:
                        results[i] = result
                except Exception as exc:
                    n_abandoned += len(chunk)
                    print(f"  [{label}] WARNING: a chunk of {len(chunk)} "
                          f"item(s) raised {exc!r} -- dropped, same as an "
                          f"unparseable candidate.")
            else:
                n_abandoned += len(chunk)
                n_stalled_chunks += 1
                n_wedged_workers += 1
                print(f"  [{label}] WARNING: no result for one chunk of "
                      f"{len(chunk)} item(s) after {stall_timeout_s}s -- "
                      f"almost certainly a pathological input in it (e.g. "
                      f"an unusually large fused-ring system RDKit's "
                      f"kekulizer struggles with). Abandoning just this "
                      f"chunk (dropped, same as an unparseable candidate) "
                      f"and continuing with the rest of the run -- that "
                      f"chunk's worker is now permanently stuck, so "
                      f"effective parallelism drops by one for the rest of "
                      f"this step ({n_wedged_workers}/{N_WORKERS} worker(s) "
                      f"now permanently wedged).")
            done += len(chunk)
            _submit_next()

            if done >= next_report or done == n:
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0.0
                remaining = (n - done) / rate if rate > 0 else float("nan")
                abandoned_note = f" ({n_abandoned} abandoned so far)" if n_abandoned else ""
                print(f"  [{label}] {done}/{n} ({100 * done / n:.1f}%) -- "
                      f"{elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining "
                      f"at current rate{abandoned_note}")
                next_report += report_every

    if done < n:
        # every worker got permanently wedged before all chunks were even submitted. Abandon the rest in one
        # bulk step (already None in results) instead of waiting out stall_timeout_s per remaining chunk.
        n_remaining = n - done
        n_abandoned += n_remaining
        done = n
        print(f"  [{label}] WARNING: all {N_WORKERS} worker process(es) "
              f"became permanently stuck on pathological inputs -- this "
              f"points at an unusually large number of them in this "
              f"source. Abandoning the remaining {n_remaining} not-yet-"
              f"started item(s) immediately rather than waiting on a pool "
              f"with no free workers left.")

    if n_stalled_chunks:
        print(f"  [{label}] finished: {n_abandoned} item(s) abandoned across "
              f"{n_stalled_chunks} stalled chunk(s) out of {n_chunks_total} "
              f"total chunks ({100 * n_abandoned / n:.3f}% of {n}).")
    return results


def main():
    rows = []
    for source, path in SOURCES.items():
        if not os.path.exists(path):
            print(f"  [skip] {source}: {path} not found -- run fetch_{source}.py first if you want it included.")
            continue
        with open(path) as f:
            smis = [line.strip() for line in f if line.strip()]
        print(f"  [{source}] {len(smis)} SMILES loaded")
        for smi in smis:
            rows.append({"smiles": smi, "source": source})

    for path in glob.glob(EXTRA_GLOB):
        with open(path) as f:
            smis = [line.strip() for line in f if line.strip()]
        label = os.path.basename(path).replace("extra_", "").replace(".smi", "")
        print(f"  [extra:{label}] {len(smis)} SMILES loaded")
        for smi in smis:
            rows.append({"smiles": smi, "source": f"extra:{label}"})

    if not rows:
        print("\nNo source files found yet -- run at least one fetch_*.py script first.")
        return

    df = pd.DataFrame(rows)
    print(f"\nTotal raw entries across all sources: {len(df)}")

    # canonicalize (defensive re-check; fetch scripts already canonicalize, but this matches
    # pool_smiles.py's own belt-and-suspenders re-pass). Parallelized across N_WORKERS processes.
    print(f"\nCanonicalizing {len(df)} SMILES across {N_WORKERS} worker "
          f"processes (CPU-bound -- this is the slow step; progress prints "
          f"every 500k so it never looks silently stuck)...")
    with mp.Pool(N_WORKERS) as pool:
        df["canonical_smiles"] = _parallel_map_with_progress(
            pool, _canon_smiles, df["smiles"].tolist(), "canonicalize"
        )
    df = df.dropna(subset=["canonical_smiles"])
    df = df.drop_duplicates(subset=["canonical_smiles"]).reset_index(drop=True)
    print(f"After canonicalization + dedup: {len(df)}")

    exclude_pattern = "|".join(EXCLUDE_LIST)
    mask = ~df["canonical_smiles"].str.contains(
        r'\b(?:{})\b'.format(exclude_pattern), regex=True, na=False
    )
    df = df.loc[mask].reset_index(drop=True)
    print(f"After excluded-element filter: {len(df)}")

    print(f"\nComputing molecular weights for {len(df)} unique candidates "
          f"across {N_WORKERS} worker processes...")
    with mp.Pool(N_WORKERS) as pool:
        mws = _parallel_map_with_progress(
            pool, _mw_of, df["canonical_smiles"].tolist(), "molecular weight"
        )
    df["mw"] = mws

    n_before_mw = len(df)
    df = df[df["mw"] <= MW_MAX_DA].reset_index(drop=True)
    n_dropped_mw = n_before_mw - len(df)
    print(f"After MW <= {MW_MAX_DA:.0f} Da filter: {len(df)} "
          f"({n_dropped_mw} dropped -- see \"MW CUTOFF\" in this script's "
          f"own docstring for the reasoning)")

    df = df[["canonical_smiles", "source", "mw"]].rename(columns={"canonical_smiles": "smiles"})
    df.to_csv(OUT_PATH, index=False)
    print(f"\nWrote {len(df)} candidates to {OUT_PATH}")
    print("\nPer-source breakdown:")
    print(df["source"].value_counts())
    print(f"\n(Lee et al. 2021 report 3,316,454 unique canonical SMILES from all "
          f"seven paper-cited sources combined, WITHOUT this project's MW <= "
          f"{MW_MAX_DA:.0f} Da cutoff applied -- your total isn't directly "
          f"comparable to that figure for that reason as well as depending on "
          f"which sources you included.)")


if __name__ == "__main__":
    main()
