"""Canonicalize all pool chunks in parallel (two worker processes), writing one parquet per chunk."""
import multiprocessing as mp
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
from canonicalize_pool_chunk import process

POOL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "external", "pool")
CHUNK_DIR = os.path.join(POOL_DIR, "chunks")
OUT_DIR = os.path.join(POOL_DIR, "chunk_out")
os.makedirs(OUT_DIR, exist_ok=True)

chunks = sorted(f for f in os.listdir(CHUNK_DIR) if f.endswith(".csv"))


def worker(fname):
    in_path = os.path.join(CHUNK_DIR, fname)
    out_path = os.path.join(OUT_DIR, fname.replace(".csv", ".parquet"))
    return process(in_path, out_path)


if __name__ == "__main__":
    t0 = time.time()
    with open(os.path.join(OUT_DIR, "PROGRESS.log"), "w") as f:
        f.write(f"Starting {len(chunks)} chunks at {time.strftime('%H:%M:%S')}\n")
    with mp.Pool(2) as pool:
        for fname, (elapsed, n) in zip(chunks, pool.imap(worker, chunks)):
            with open(os.path.join(OUT_DIR, "PROGRESS.log"), "a") as f:
                f.write(f"DONE {fname}: {n} rows in {elapsed:.1f}s at {time.strftime('%H:%M:%S')}\n")
    with open(os.path.join(OUT_DIR, "PROGRESS.log"), "a") as f:
        f.write(f"ALL DONE in {time.time()-t0:.1f}s\n")
