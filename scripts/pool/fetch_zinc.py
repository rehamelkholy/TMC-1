"""
Fetch ZINC (Sterling & Irwin 2015, J. Chem. Inf. Model. 55, 2324) -- Lee et al. 2021 report 3,862,980
entries, their single largest source.

Unlike the other five sources, ZINC has no fixed, scriptable "download everything" URL: the bulk-tranche
system (zinc15.docking.org/tranches/home/) requires one manual, interactive step to generate the list of
tranche download links (there is no public API endpoint that returns this list programmatically).
Everything else in this folder is a single script run.

ONE-TIME MANUAL STEP (a couple of minutes):
  1. Go to https://zinc15.docking.org/tranches/home/
  2. Under "2D Format", pick "SMILES" (structure only, which keeps the download far smaller).
  3. MW cutoff: restrict the molecular-weight grid to <= 500 Da (deselect every MW row above 500,
     leaving the logP columns as-is). This reduces the number of tranches to fetch, which helps with
     the rate limiting described below. combine_pool.py enforces the same cutoff on the final pool
     regardless of exactly which tranches are selected here.
  4. Click the download/output option that gives a plain list of URLs (e.g. "URLs (default)"; NOT the
     wget-wrapped script), one URL per line.
  5. Save that file as external/pool/raw/zinc_urls.txt

Then run this script -- it downloads every tranche in that list and canonicalizes the SMILES, like the
other fetch_*.py scripts.

Implementation notes:
  - files.docking.org returns 403 Forbidden to requests carrying Python's default urllib User-Agent, so
    every request goes through http_common.py, which sets a browser-like one.
  - The tranche downloads are independent and network-bound, so each tranche is downloaded and parsed in
    a worker thread (ThreadPool). The script is fully resumable: a tranche file already on disk is
    parsed, not re-downloaded.
  - files.docking.org rate-limits sustained load (mid-transfer connection resets, SSL EOFs, HTTP 503).
    Mitigations: low concurrency (MAX_WORKERS = 3), a global pacing gate (MIN_REQUEST_INTERVAL_SECONDS)
    shared across ALL worker threads that caps how often a new download can START, retries with backoff
    (N_DOWNLOAD_ATTEMPTS, RETRY_BACKOFF_SECONDS), and a longer dedicated backoff for HTTP 503 that
    respects the server's Retry-After header when present (30 s otherwise).
  - Tranche URLs that still fail after all retries are written to raw/zinc_failed_urls.txt at the end;
    re-running the script later retries only those. If tranches keep failing, the cheapest remedy is
    time, or narrowing the tranche selection on the ZINC site.

Output: external/pool/raw/zinc_smiles.smi
"""
import concurrent.futures as cf
import os
import threading
import time
import urllib.error

import http_common

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "external", "pool", "raw")
os.makedirs(OUT_DIR, exist_ok=True)
URLS_PATH = os.path.join(OUT_DIR, "zinc_urls.txt")
TRANCHE_DIR = os.path.join(OUT_DIR, "zinc_tranches")
OUT_PATH = os.path.join(OUT_DIR, "zinc_smiles.smi")
FAILED_URLS_PATH = os.path.join(OUT_DIR, "zinc_failed_urls.txt")

MAX_WORKERS = 3  # kept low to stay under the server's rate limit
N_DOWNLOAD_ATTEMPTS = 5
RETRY_BACKOFF_SECONDS = (5, 15, 30, 60)  # wait before attempts 2..5 respectively
                                          # (last value repeats if N_DOWNLOAD_ATTEMPTS grows)
SERVICE_UNAVAILABLE_BACKOFF_SECONDS = 30  # fallback wait for a 503 with no Retry-After header
MIN_REQUEST_INTERVAL_SECONDS = 0.75  # global pacing floor between ANY two download starts,
                                      # across all worker threads combined

_rate_lock = threading.Lock()
_last_request_started_at = [0.0]


def _paced_download(url, fname):
    """download_to_file, but gated so that across ALL worker threads combined, a new download never
    STARTS less than MIN_REQUEST_INTERVAL_SECONDS after the previous one started. This bounds the start
    rate, not the concurrency (MAX_WORKERS already bounds how many downloads can be in flight)."""
    with _rate_lock:
        wait = _last_request_started_at[0] + MIN_REQUEST_INTERVAL_SECONDS - time.time()
        if wait > 0:
            time.sleep(wait)
        _last_request_started_at[0] = time.time()
    http_common.download_to_file(url, fname)


def _fetch_and_parse_one(i, url):
    """Runs in a worker thread: download this one tranche (skip if already
    on disk -- resumability), parse its SMILES, return the set of canonical
    SMILES found (or an error string). Returns a set, not a shared
    collection, specifically so nothing here needs a lock: all the actual
    merging happens back in the main thread."""
    from rdkit import Chem

    fname = os.path.join(TRANCHE_DIR, url.split("/")[-1])
    if not os.path.exists(fname):
        last_err = None
        next_wait = 0.0
        for attempt in range(N_DOWNLOAD_ATTEMPTS):
            if next_wait > 0:
                time.sleep(next_wait)
            backoff_idx = min(attempt, len(RETRY_BACKOFF_SECONDS) - 1)
            next_wait = RETRY_BACKOFF_SECONDS[backoff_idx]
            try:
                _paced_download(url, fname)
                last_err = None
                break
            except urllib.error.HTTPError as e:
                last_err = e
                if e.code == 503:
                    retry_after = e.headers.get("Retry-After") if e.headers else None
                    try:
                        service_wait = float(retry_after) if retry_after else SERVICE_UNAVAILABLE_BACKOFF_SECONDS
                    except ValueError:
                        service_wait = SERVICE_UNAVAILABLE_BACKOFF_SECONDS
                    next_wait = max(next_wait, service_wait)
            except Exception as e:
                last_err = e
        if last_err is not None:
            return i, url, None, f"FAILED to fetch after {N_DOWNLOAD_ATTEMPTS} attempts: {last_err}"

    smis = set()
    try:
        with open(fname) as tf:
            next(tf, None)  # ZINC .smi tranche files have a header row
            for line in tf:
                parts = line.split()
                if not parts:
                    continue
                mol = Chem.MolFromSmiles(parts[0])
                if mol is not None:
                    smis.add(Chem.MolToSmiles(mol, canonical=True))
    except Exception as e:
        return i, url, None, f"FAILED to parse {fname}: {e}"
    return i, url, smis, None


def main():
    if not os.path.exists(URLS_PATH):
        print(f"Missing {URLS_PATH} -- see the manual step described at the top of this "
              f"script's docstring (fetch_zinc.py) before running it.")
        return

    os.makedirs(TRANCHE_DIR, exist_ok=True)
    with open(URLS_PATH) as f:
        urls = [line.strip() for line in f if line.strip()]
    n_total = len(urls)
    print(f"{n_total} tranche URLs to fetch ({MAX_WORKERS} at a time).")

    seen = set()
    n_done = 0
    n_failed = 0
    failed_urls = []
    print_lock = threading.Lock()

    with cf.ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {ex.submit(_fetch_and_parse_one, i, url): url for i, url in enumerate(urls)}
        for fut in cf.as_completed(futures):
            _, url, smis, err = fut.result()
            with print_lock:
                n_done += 1
                if err:
                    n_failed += 1
                    failed_urls.append(url)
                    print(f"  [{n_done}/{n_total}] {err}: {url}")
                else:
                    seen.update(smis)
                if n_done % 50 == 0 or n_done == n_total:
                    print(f"  [{n_done}/{n_total}] tranches processed, "
                          f"{len(seen):,} unique SMILES so far ({n_failed} failed)")

    with open(OUT_PATH, "w") as f:
        for smi in seen:
            f.write(smi + "\n")
    print(f"Wrote {len(seen):,} unique canonical SMILES to {OUT_PATH} "
          f"(paper reports 3,862,980 from this source).")

    if failed_urls:
        with open(FAILED_URLS_PATH, "w") as f:
            f.write("\n".join(failed_urls) + "\n")
        print(f"\n{len(failed_urls)} tranche(s) still failed after {N_DOWNLOAD_ATTEMPTS} attempts each "
              f"-- listed in {FAILED_URLS_PATH} for reference. This script is fully resumable: just "
              f"run it again (now, or later) and it will only retry these, leaving every "
              f"already-downloaded tranche untouched.")
    elif os.path.exists(FAILED_URLS_PATH):
        os.remove(FAILED_URLS_PATH)  # clean up a stale list from a previous run that's now fully resolved


if __name__ == "__main__":
    main()
