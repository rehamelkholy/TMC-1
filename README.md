# TMC-1: ranking unseen interstellar-molecule candidates with a frozen model

Code, frozen inputs, results and tables for the TMC-1 detection-ranking study.

* Paper: El-Kholy, R. I. & Hayman, Z. M. 2026, "Active Learning for Molecular Discovery: A Retrospective Test toward TMC-1", submitted to ApJ
* Repository: https://github.com/rehamelkholy/TMC-1
* Contact: R. I. El-Kholy, relkholy@sci.cu.edu.eg
* Zenodo DOI: [10.5281/zenodo.23312984](https://doi.org/10.5281/zenodo.23312984)

The study asks whether a model fitted only on what was known before a cut-off date ranks the later-detected TMC-1 molecules high among
~1,400 candidate species, using the recommendation list of Lee et al. (2021, ApJL 917, L6) as the candidate universe. The analysis plan
(`protocol/`) was fixed before the results were seen; later changes are recorded as dated amendments.

## What is reproduced

Starting from the frozen inputs in `data/frozen/` (candidate table, detection set, Lee-list labels, precomputed features), the pipeline
re-computes

* the three sealed result files `results/sealed/{released,refit2021,leave_self_out}.json` (primary evaluation at four origin dates, the two
  date-handling sensitivities, bootstraps, baselines and the feedback learner),
* Figures 1-3 (data CSVs, consistency checks, PDF/PNG), Tables 1-3 (CSV and journal-format MRT), and the descriptive statistics.

The features themselves (mol2vec embedding, Gaussian process, leave-self-out refits) come from the umda pipeline in a second, older
environment; they are shipped as frozen files and verified by hash (see "Environments" and "What was verified"). The construction of the
frozen inputs from the raw lists is documented in `provenance/` but is not runnable.

## Quick start

```bash
python3.11 -m venv .venv && .venv/bin/pip install -r env/requirements-analysis.txt
bash scripts/00_fetch_external.sh            # only for the embedding / pool tiers and the optional parser check; not needed below
mkdir -p results/sealed
P=.venv/bin/python
$P scripts/20_seal_released.py       results/sealed/released.json
$P scripts/21_seal_refit2021.py      results/sealed/refit2021.json     # also writes data/derived/features_refit2021_relabelled.csv
$P scripts/22_seal_leave_self_out.py results/sealed/leave_self_out.json
$P scripts/30_descriptive_discussion.py data/frozen/candidates.csv results/descriptive.json
$P scripts/01_parse_lee_recommendations.py   # optional, needs external/tmc1_recommendations.pdf and pdftotext; writes build/ and compares with data/raw
$P scripts/40_table1_v13.py; $P scripts/41_table1_v14.py; $P scripts/42_table1_lee_rematch_check.py
$P scripts/60_table3_recommended_targets.py
$P scripts/50_fig1_pooled_auc.py; $P scripts/51_fig2_budget_yield.py; $P scripts/52_fig3_feedback.py
$P scripts/53_table2_ranks.py; $P scripts/70_changelogs.py
```

Run the scripts in this order: Figure 1 and 3 read the sealed files, Figure 1 and Table 2 read the relabelled features and Table 1.
All scripts resolve paths through `src/paths.py`, so they can be started from any directory. Re-running writes into the tracked
files (`results/`, `tables/`, `figures/`); `git diff` shows whether anything changed. Intermediate files go to `build/` (git-ignored).

Run times, measured on a 2-core machine with the analysis environment:

| Script | Writes | Time |
|---|---|---|
| `20_seal_released.py` | `results/sealed/released.json` | about 1 min |
| `21_seal_refit2021.py` | `results/sealed/refit2021.json`, `data/derived/features_refit2021_relabelled.csv` | about 40 s |
| `22_seal_leave_self_out.py` | `results/sealed/leave_self_out.json` | about 40 s |
| `01_parse_lee_recommendations.py` | `build/lee_recommendations_parsed.csv` (compared with `data/raw/lee_recommendations_parsed.csv`) | seconds |
| `30_descriptive_discussion.py` | `results/descriptive.json` | seconds |
| `40`-`42_table1_*.py` | `tables/table1_validation_set.{csv,mrt}`, `results/table1_lee_rematch_check.csv` | seconds |
| `60_table3_recommended_targets.py` | `tables/table3_recommended_targets.{csv,mrt}` | seconds |
| `50_fig1_pooled_auc.py` | `figures/fig1.*`, `fig1_data.csv`, `fig1_consistency_check.csv` | 2.5-3 min |
| `51_fig2_budget_yield.py` | `figures/fig2.*`, `fig2_data.csv`, `fig2_consistency_check.csv`, `fig2_everyk_check.csv` | 15 s |
| `52_fig3_feedback.py` | `figures/fig3.*`, `fig3_data.csv`, `fig3_consistency_check.csv` | about 2.5 min |
| `53_table2_ranks.py` | `tables/table2_ranks.csv` | seconds |
| `70_changelogs.py` | `provenance/changelogs/*.csv` | seconds |

The whole analysis takes about 10 minutes. The embedding tier (`scripts/features/`) and the candidate-pool tier (`scripts/pool/`) are
not part of this and are described below.

## Environments

Two Python environments are needed (details, the patch to `mol2vec` and the dropped pins are in `env/README_env.md`):

* **Analysis, Python 3.11** - `env/requirements-analysis.txt`. scikit-learn 1.8.0, numpy 2.4.4, scipy 1.17.1, pandas 3.0.2, matplotlib,
  rdkit, astropy, joblib. All results above come from this environment.
* **Embedding, Python 3.10** - `env/requirements-embedding.txt`, `bash env/install_embedding_env.sh`. The umda pipeline's pretrained
  regressors unpickle only under scikit-learn 1.0.2, which has no wheel for Python 3.11 or later. Uses the mol2vec fork
  laserkelvin/mol2vec @ 850d944 with the patched `features.py` shipped in `env/`.

## Layout

```
README.md  LICENSE  CITATION.cff
env/          pinned requirements for both environments, install script, patched mol2vec features.py (+ its licence)
protocol/     pre-registered protocol and its 9 dated amendments (with .sha256 files), kept verbatim
data/raw/     Lee-list entries and labels, parsed recommendation list, anchor table (87 species); FETCH.md for third-party inputs
data/frozen/  frozen inputs of the evaluation (candidates, detections, labels, precomputed features, Table 1 v12, Table 3 annotations)
data/derived/ features_refit2021_relabelled.csv (written by scripts/21_seal_refit2021.py)
literature/   ADS queries and ID-only screening tables (bibcode, DOI, arXiv id, title, decision, reason)
src/          evaluation engine, bootstrap, feedback learner, candidate construction, plotting; src/umda_pipeline/ = embedding / GPR code
scripts/      numbered in run order (see above); features/ = embedding tier; pool/ = candidate-pool tier; 00_fetch_external.sh
results/      sealed result files (+ SHA256SUMS), descriptive.json, coverage_check_37.csv, table1_lee_rematch_check.csv
tables/       Table 1 (CSV + MRT), Table 2 (CSV), Table 3 (CSV + MRT)
figures/      Figures 1-3 (PDF, PNG, data CSVs, consistency checks)
provenance/   NON-RUNNABLE lineage: intermediate versions, the scripts that built them, change logs, superseded sealed results
```

## Paper items and where they come from

Section numbers are not given; the items are named as in the manuscript.

| Paper item | Script | Output |
|---|---|---|
| Figure 1 (pooled AUC by regime and method) | `scripts/50_fig1_pooled_auc.py` | `figures/fig1.pdf`, `fig1.png`, `fig1_data.csv`; consistency check against the sealed AUCs: `fig1_consistency_check.csv` |
| Figure 2 (budget vs yield, released regime) | `scripts/51_fig2_budget_yield.py` | `figures/fig2.pdf`, `fig2.png`, `fig2_data.csv`, `fig2_everyk_check.csv` |
| Figure 3 (feedback simulation) | `scripts/52_fig3_feedback.py` | `figures/fig3.pdf`, `fig3.png`, `fig3_data.csv` (amendment 9 tie handling; earlier version in `provenance/superseded_figures/`) |
| Table 1 (validation set of post-2021 detections) | `scripts/40_table1_v13.py` -> `41_table1_v14.py` (from `data/frozen/table1_validation_set_v12.csv`) | `tables/table1_validation_set.csv`, `.mrt`; match check `results/table1_lee_rematch_check.csv` |
| Table 2 (rank of each later detection) | `scripts/53_table2_ranks.py` | `tables/table2_ranks.csv` |
| Table 3 (top-100 recommended targets) | `scripts/60_table3_recommended_targets.py` (written for the release, see below) | `tables/table3_recommended_targets.csv`, `.mrt` |
| Pooled AUC, top-k yield, per-origin results, baselines, bootstrap comparison, sensitivities A4(a) and A4(b), feedback-learner median | `scripts/20_seal_released.py` (released labels), `21_seal_refit2021.py`, `22_seal_leave_self_out.py` | `results/sealed/released.json`, `refit2021.json`, `leave_self_out.json` (keys `primary`, `sensitivity_*`, `v3_feedback_learner`) |
| Descriptive statistics of the discussion section | `scripts/30_descriptive_discussion.py` | `results/descriptive.json` |
| Coverage of the 37 species of interest in the broader pool ("13 of 37") | none (frozen result) | `results/coverage_check_37.csv` |
| Results from intermediate label versions | none (frozen results) | `provenance/sealed_superseded/`, with `provenance/changelogs/` for the label changes |
| Pre-registration | - | `protocol/` |

## Data provenance and licences

| Data | Provenance | Licence / terms |
|---|---|---|
| Authors' derived tables and results (`data/frozen`, `data/derived`, `tables`, `results`, `figures`, `provenance/data`, annotations) | this study | CC BY 4.0 |
| `data/raw/lee_recommendations_parsed.csv` | parsed from the recommendation list of Lee et al. (2021), Zenodo version 2, 10.5281/zenodo.5146276 (`scripts/01_parse_lee_recommendations.py`). Version 1 (record 5080543, 7 July 2021) lists the same candidates with different `logcd` and `sigma` values; the analysis uses version 2. | CC BY 4.0; cite Lee et al. (2021) |
| `data/raw/table_a2_anchors.csv` | the 87 TMC-1 species used to train umda (Table A2 of Lee et al. 2021, ApJL 917, L6), structures checked by the authors | factual data, shipped with citation; the journal's permission for reuse has not been checked |
| Detection dates and references (`data/frozen/post2021_detections.csv`, Table 1) | the primary literature (DOIs / bibcodes / arXiv ids given in the files) | facts with citation |
| `literature/ads_screening_*_ids.csv` | NASA ADS records (identifiers and titles only) | identifiers; ADS abstracts are not redistributed |
| umda (Lee et al. 2021) | fetched at a pinned commit | MIT, Copyright 2020 Kin Long Kelvin Lee |
| mol2vec fork, patched `env/mol2vec_features_patched.py` | laserkelvin/mol2vec @ 850d944 | BSD 3-Clause (`env/LICENSE-mol2vec.txt`) |
| Candidate-pool sources: QM9, ZINC, PubChem, PCBA, KIDA, NASA Ames PAHdb | downloaded by `scripts/pool/fetch_*.py`; the combined pool (8,910,820 rows) is not redistributed | each source keeps its own terms: ZINC restricts redistribution of major portions; PubChem contains third-party material; PCBA derives from PubChem BioAssay (dataset licence not stated); QM9 (figshare licence not retrieved); KIDA (no data licence found); NASA Ames PAHdb (terms of service require consent for reuse). Check each before use. |
| Toru & Shay (2025) and Lee et al. (2021) papers | cited only | not redistributed |

The code is MIT-licensed (`LICENSE`).

## Frozen files without a runnable producer

Several frozen files were produced by scripts that are no longer available. They are shipped as frozen data and no script is presented as
their producer; the table in `provenance/README.md` lists each file and its status, and `provenance/changelogs/` records every
cell that changed between consecutive versions. In short: `candidates.csv` (v6 -> v7 label edits), `post2021_detections.csv` (v7),
`table1_validation_set_v12.csv`, the leave-self-out features, `coverage_check_37.csv` and the superseded
sealed results.

`scripts/60_table3_recommended_targets.py` is new for the release: it rebuilds Table 3 from the released candidate table (model fit,
ranking, removal of the two non-physical structures at ranks 5 and 71, reasons) and reads the hand-curated columns from
`data/frozen/table3_annotations.csv`. It reproduces the frozen table byte for byte.

## Byte-exact files

The sealed JSONs, `data/frozen/post2021_detections.csv`, the anchor table and the protocol files are hashed records and were not edited.
They (and the strings that the seal scripts write into the JSONs) therefore still contain internal wording about how they were produced.
`tables/table1_validation_set.mrt` and `tables/table3_recommended_targets.mrt` are the journal submission files; their headers were
regenerated (only the title, authors and table-title lines and the column descriptions were written) and their data blocks are byte-identical to the frozen files.

## What was verified

A scratch copy of this tree was run from `env/` in fresh environments (Python 3.11.16 and 3.11.17 for the analysis, 3.10.12 for the embedding tier):

* byte-identical to the frozen files: `released.json`, `refit2021.json`, `descriptive.json`, `features_refit2021_relabelled.csv`, Table 1 (v13 and v14 CSV) and its MRT data block, Table 2, Table 3 (CSV and MRT data block), `fig1_data.csv`, `fig2_data.csv` and the consistency / every-k checks of Figures 1 and 2, and the Figure 3 consistency check;
* `leave_self_out.json`: all computed values identical; the file also records the interpreter (`environment.python`, `python_executable`), so it differs byte-wise from the sealed copy when a different interpreter is used;
* `figures/fig3_data.csv`, `fig3.pdf`, `fig3.png` were regenerated by `scripts/52_fig3_feedback.py` under amendment 9 (R2, exact tie expectation). The pre-amendment files differ in 30 of 805 rows (the `high_logcd` and `low_dist` baseline curves, which averaged 200 random tie-breaks) and are kept in `provenance/superseded_figures/`. The feedback-learner curve is unchanged, and the consistency check against `results/sealed/released.json` is identical (learner median 200.0, IQR 185-221.25, target 32 of 63; fixed orderings high_logcd 428.49, detectability 348.0, low_distance 332.01);
* PDF files are not compared byte-wise (they embed creation metadata); `fig2.png` is byte-identical, and the regenerated `fig3.png` differs from the superseded one in 0.09% of its pixels (the redrawn baseline curves above);
* the frozen feature files were verified by SHA-256 (`features_leave_self_out.csv` `f86d3424...`, `features_refit2021.csv` `f0673be8...`) and not rebuilt; in the embedding environment the underlying anchor fit reproduces the reference log marginal likelihood (-632.1523) and the stored `features_refit2021.csv` values to within 3e-10.

SHA-256 of the main outputs: `results/sealed/SHA256SUMS`.

## Citation

See `CITATION.cff`. Please cite the paper (El-Kholy & Hayman 2026, submitted to ApJ) and the Zenodo archive (doi:10.5281/zenodo.23312984), and Lee et al. (2021, ApJL 917, L6) for
the recommendation list and the umda pipeline (doi:10.3847/2041-8213/ac194b).
