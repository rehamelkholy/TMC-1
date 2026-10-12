# Environments

The release needs two Python environments. They cannot be merged.

| Environment | Python | Used by | Why separate |
|---|---|---|---|
| `requirements-analysis.txt` | 3.11 | everything in `scripts/` (except `scripts/features/`) and `src/` (except `src/umda_pipeline/`), including the pool scripts | the evaluation results were produced with scikit-learn 1.8.0 / numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.2; the lbfgs fits and the bootstraps depend on those builds |
| `requirements-embedding.txt` + `install_embedding_env.sh` | 3.10 | `src/umda_pipeline/`, `scripts/features/` | the pretrained regressors shipped with umda unpickle only under scikit-learn 1.0.2, which has no wheel for Python 3.11+ |

Create them with

```
python3.11 -m venv .venv && .venv/bin/pip install -r env/requirements-analysis.txt
bash env/install_embedding_env.sh .venv-embed
```

Verified for this release: the analysis environment was built from `requirements-analysis.txt` with CPython 3.11.16
and 3.11.17 (the original results were produced with 3.11.15); the embedding environment was built on CPython 3.10.12
and used to load umda's `EmbeddingModel.pkl`, `mol2vec_model.pkl` and `regressors.pkl` and to embed test molecules.

## Patched `mol2vec/features.py`

`mol2vec_features_patched.py` replaces `site-packages/mol2vec/features.py` of the mol2vec fork
(laserkelvin/mol2vec @ 850d944d5f48a58e26ed0264332b5741f72555aa; BSD 3-Clause, `LICENSE-mol2vec.txt`).
`install_embedding_env.sh` does the copy. Without the patch the embedding silently fails (all test molecules returned "not embeddable" when the
unpatched fork was used). The patch changes:

* RDKit sanitization is controllable. `generate_corpus`, `_parallel_job` and `_read_smi` gain a `sanitize=False` argument that is passed
  to `Chem.MolFromSmiles(..., sanitize=sanitize)`; after parsing without sanitization, `mol.UpdatePropertyCache(strict=False)` and
  `Chem.GetSymmSSSR(mol)` are called so that radicals, charged and valence-unusual species (the interstellar molecules in this study)
  can still be fingerprinted.
* gensim >= 4.0 API: `Word2Vec(..., vector_size=...)` instead of `size=`; `wv.key_to_index` / `wv.get_vector` are used when
  available, with a fallback to `wv.vocab` / `wv.word_vec` for models pickled under gensim 3.x (`sentences2vec`).
* The "molecules successfully processed" message prints the count (the upstream format string had no argument).
* A license/provenance header comment.

## One-line patch to dask-ml

dask-ml 1.8.0 (`_typing.py`) uses the `np.float` / `np.int` aliases that were removed from numpy 1.24+. `install_embedding_env.sh` runs
`sed -i 's/np\.float\b/float/g; s/np\.int\b/int/g' site-packages/dask_ml/_typing.py`. umda needs dask-ml to unpickle its embedding model.

## Pins that were dropped or added

Dropped from the original analysis lock file (`requirements.lock`, 199 packages captured from the general-purpose Python 3.11 system
environment of the project): everything that no script in the release imports. The release imports numpy, scipy, pandas,
scikit-learn, matplotlib, rdkit and astropy (+ joblib, which scikit-learn uses and which the embedding scripts import). Kept in addition: the
indirect dependencies of those packages. Dropped (178 packages): absl-py, annotated-types, anyio, argcomplete, asttokens, attrs, babel, backrefs, beautifulsoup4, blinker, camelot-py, certifi, cffi, chardet, charset-normalizer, claude-agent-sdk, click, cloudpickle, colorama, conan, cryptography, dask, dask-glm, dask-ml, dbus-python, defusedxml, Deprecated, distributed, distro, docopt, et_xmlfile, executing, fasteners, Flask, flatbuffers, fsspec, gensim, ghp-import, greenlet, grip, h11, h5py, httpcore, httplib2, httpx, httpx-sse, idna, ImageIO, imageio-ffmpeg, img2pdf, importlib_metadata, ipython, ipython_pygments_lexers, itsdangerous, jedi, Jinja2, jsonschema, jsonschema-specifications, launchpadlib, lazr.restfulclient, lazr.uri, lazy-loader, llvmlite, locket, loguru, lxml, magika, Markdown, markdown-it-py, markdownify, markitdown, marko, MarkupSafe, matplotlib-inline, mcp, mdurl, mediapipe, mergedeep, mistune, mkdocs, mkdocs-get-deps, mkdocs-material, mkdocs-material-extensions, mol2vec, msgpack, multipledispatch, networkx, numba, oauthlib, odfpy, onnxruntime, opencv-contrib-python, opencv-python, opencv-python-headless, openpyxl, paginate, palettable, parso, partd, patch-ng, path-and-address, pathspec, pdf2image, pdfkit, pdfminer.six, pdfplumber, pexpect, pikepdf, pipx, platformdirs, playwright, prompt_toolkit, protobuf, psutil, ptyprocess, pure_eval, py2opsin, pyarrow, pycparser, pydantic, pydantic-settings, pydantic_core, pyee, Pygments, PyGObject, PyJWT, pymdown-extensions, pyoo, pypdf, pypdfium2, pytesseract, python-apt, python-docx, python-dotenv, python-multipart, python-pptx, pyyaml_env_tag, referencing, reportlab, requests, rpds-py, scikit-image, seaborn, smart_open, sniffio, sortedcontainers, sounddevice, soupsieve, sparse, sse-starlette, stack-data, starlette, tabula-py, tabulate, tblib, tifffile, toml, toolz, tornado, tqdm, traitlets, typing-inspection, typing_extensions, unoserver, urllib3, userpath, uvicorn, wadllib, Wand, watchdog, wcwidth, Werkzeug, wrapt, xlsxwriter, xmltodict, yq, zict, zipp.

Dropped from the original embedding pin list (`requirements_windows_py310.txt`): `palettable` and `click` (neither is imported by the
code or needed to load the umda models; click is still installed as a dependency of dask). Added: `tqdm` (imported by the mol2vec fork)
and `setuptools==80.9.0` (provides `pkg_resources`, which umda's unpickling imports; it was removed from later setuptools). Everything
else is unchanged, including the mol2vec commit pin.
