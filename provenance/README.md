# Provenance (not runnable)

Everything in this folder documents how the frozen inputs in `data/frozen/` and the intermediate label versions were built. It is shipped
for transparency and is **not part of the reproduction pipeline**: nothing here is executed by the scripts in `scripts/`.

| Folder | Contents |
|---|---|
| `scripts/` | The scripts that produced the successive detection sets, Lee-list labels, anchors, purged universes and candidate tables (versions 1 to 6), the ADS export parsers and the screening-batch builder, numbered in the order they were run. They are shipped as they were used, apart from comment and docstring clean-up. They contain absolute paths of the original working tree, written `<PROJECT_ROOT>`, and read files that are not shipped (raw ADS exports, working directories of earlier steps), so they cannot be run from this repository. |
| `data/` | The intermediate versions of the curated tables: detections v2 to v6, Lee-list labels v3.1 to v5 (raw and purged), candidates v1 to v6, Table 1 v5 to v13, Table 3 v5 to v7, and the first anchor table. Together with `data/frozen/` and `tables/` this is every version that exists. |
| `changelogs/` | Cell-level change logs between consecutive versions (`row_key, column, old_value, new_value`), produced by `scripts/70_changelogs.py`. They make every label change traceable even where the script that made it is no longer available. |
| `sealed_superseded/` | The earlier sealed result files (38 files) from the intermediate label versions, kept because the manuscript discusses results from intermediate label versions. They were superseded by `results/sealed/*.json`; do not use them for the headline numbers. |
| `superseded_figures/` | Figure 3 files produced before amendment 9 (R2); replaced by the regenerated files in `figures/`. See `superseded_figures/NOTE.md`. |

## Versions that have no producer script in this release

The files below are frozen outputs. The script that produced them is not available (or, where noted, is shipped but not runnable), so they
are included as data and **no script is presented as their producer**. Where two consecutive versions exist, `changelogs/` records what changed.

| File | Status |
|---|---|
| `data/frozen/candidates.csv` (candidate table, version 7) | Obtained from version 6 by six label edits (candidate ids 1305, 522, 498, 318, 512, 1312); the script that applied them is not available. Change log: `changelogs/candidates_v6_to_v7.csv`. |
| `data/frozen/post2021_detections.csv` (version 7) | Builder not available. Change log: `changelogs/post2021_detections_v6_to_v7.csv`. The file is kept byte-exact (it is hashed by the seal), including its free-text `notes` column. |
| `data/frozen/table1_validation_set_v12.csv` | Builders of Table 1 versions 11 and 12 are not available. Versions 13 and 14 are rebuilt from it by `scripts/40_table1_v13.py` and `scripts/41_table1_v14.py`. Change logs: `changelogs/table1_validation_set_v10_to_v12.csv` (version 11 does not exist as a file), `..._v12_to_v13.csv`, `..._v13_to_v14.csv`. |
| `data/frozen/features_leave_self_out.csv` | Features for the leave-self-out regime (versions 9 and 10 are byte-identical, SHA-256 `f86d3424...97dd`). The version-8 builder is shipped as `scripts/61_leave_self_out_features_v8.py` (non-runnable); the version-9/10 builders are not available. The computation is sensitive to the BLAS build and took about 36 minutes, so it was verified by hash and not re-run. |
| `data/frozen/features_refit2021.csv` | Produced by `scripts/features/03_refit2021_features.py` (embedding environment, needs umda). Not re-run in the verification run as a whole; its SHA-256 matches the recorded value and the fit it contains was re-computed in the embedding environment to within 3e-10 (see README). |
| `data/frozen/lee_list_entry_labels_purged.csv` | Version 6 of the purged Lee-list labels; the producer is `scripts/44_labels_v6.py` (shipped, non-runnable). |
| `data/raw/lee_recommendations_parsed.csv` | Parsed from the Lee et al. (2021) recommendation PDF, Zenodo version 2 (record 5146276). The original parser is not available; `scripts/01_parse_lee_recommendations.py` was written for the release and reproduces this file byte for byte from the version-2 PDF (not an original). Version 1 of the deposit (record 5080543) has different `logcd` and `sigma` values and is not used. |
| `data/raw/lee_list_entry_labels.csv`, `data/raw/detections_labeled_vs_lee.csv` | Inputs of the early lineage steps; their producers are not available. |
| `data/raw/table_a2_anchors.csv` | Authoritative anchor table (87 species, version 2). Version 1 and the script that built the series (`scripts/10_anchors.py`) are in this folder (non-runnable). Kept byte-exact together with its `.sha256` file. |
| `results/coverage_check_37.csv` | Result of an ad hoc check whether 37 species of interest are found in the broader pool (13 found, all by rule A). The check script and the chunked pool outputs it used are not available. |
| `provenance/sealed_superseded/*.json` | The scripts that sealed them are not available (only the final versions in `results/sealed/` are reproducible). |
| `tables/table3_recommended_targets.csv` | Reproduced by the new script `scripts/60_table3_recommended_targets.py`, written for the release (not an original). Its four hand-curated annotation columns (`Species`, `Elsewhere`, `LabSpec`, `Note`) and the audit flags in `Reasons` are read from `data/frozen/table3_annotations.csv`, which was extracted from the frozen Table 3; the original catalog-lookup scripts are not available. |
| `provenance/data/table3_recommended_targets_v5..v7.csv`, `table1_validation_set_v5..v10.csv` | Builders exist for some of them but are not shipped; the cell-level change logs are provided instead. |

## Byte-exact files that contain internal wording

Several files are covered by hashes recorded when they were frozen and therefore were not edited, even though they contain wording about
how they were produced: the three sealed result files (`results/sealed/*.json`, keys `task` and `implementation_notes`),
`data/frozen/post2021_detections.csv` and `provenance/data/post2021_detections_v*.csv` (free-text `notes`), `data/raw/table_a2_anchors.csv` and the protocol and amendment files in
`protocol/` (a pre-registration record, kept verbatim; its dated entries refer to earlier stages of the work by their internal labels). The strings written into the sealed JSONs by the seal scripts
(`src/seal_released.py`, `scripts/21_*`, `scripts/22_*`) are left unchanged for the same reason.
