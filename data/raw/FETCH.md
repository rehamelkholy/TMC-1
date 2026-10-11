# Third-party inputs that are not redistributed

Run `bash scripts/00_fetch_external.sh` to fetch both items into `external/` (git-ignored), or do it by hand:

| Item | Where from | Licence | Needed by |
|---|---|---|---|
| umda, commit `a95cc1c1eb98a1ff64b37a4c6ed92f9546d0215b` | `git clone https://github.com/laserkelvin/umda.git external/umda && git -C external/umda checkout a95cc1c1eb98a1ff64b37a4c6ed92f9546d0215b` | MIT, Copyright 2020 Kin Long Kelvin Lee | the embedding environment: `src/umda_pipeline/` loads `models/EmbeddingModel.pkl` and `models/regressors.pkl` from this checkout; the `umda` package must be importable (`PYTHONPATH=external/umda`) to unpickle the embedding model |
| Lee et al. (2021) recommendation list, `tmc1_recommendations.pdf` version 2 (10.2 MB, md5 `9d3664ca507b5542e2bf3a279599f891`) | Zenodo, https://doi.org/10.5281/zenodo.5146276 | CC BY 4.0 | only to re-derive `data/raw/lee_recommendations_parsed.csv` with `scripts/01_parse_lee_recommendations.py`; the evaluation does not read it. Version 1 (record 5080543, 7 July 2021) lists the same candidates with different `logcd` and `sigma` values and is not used |

The candidate-pool sources (QM9, ZINC, PubChem, PCBA, KIDA, NASA Ames PAHdb) are fetched by `scripts/pool/fetch_*.py` into
`external/pool/raw/`; their data are not redistributed (see the README for the terms of each source). Only the reference-list canonicalization
and the pool scripts need them; the evaluation does not.
