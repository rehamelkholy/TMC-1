"""Filesystem layout of the release.

Every script and module resolves its input and output locations through this file, so the
repository can be cloned anywhere and run from any working directory.
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_RAW = os.path.join(ROOT, "data", "raw")
DATA_FROZEN = os.path.join(ROOT, "data", "frozen")
DATA_DERIVED = os.path.join(ROOT, "data", "derived")
PROV_DATA = os.path.join(ROOT, "provenance", "data")
RESULTS = os.path.join(ROOT, "results")
SEALED = os.path.join(RESULTS, "sealed")
TABLES = os.path.join(ROOT, "tables")
FIGURES = os.path.join(ROOT, "figures")
BUILD = os.path.join(ROOT, "build")  # intermediate files written by the pipeline (git-ignored)
os.makedirs(BUILD, exist_ok=True)
EXTERNAL = os.path.join(ROOT, "external")  # third-party inputs fetched by data/raw/FETCH.md (git-ignored)
