#!/usr/bin/env bash
# Build the Python 3.10 embedding environment (see env/requirements-embedding.txt and env/README_env.md).
# Usage: bash env/install_embedding_env.sh [venv-dir]      (default: .venv-embed, created in the current directory)
set -euo pipefail

VENV="${1:-.venv-embed}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python3.10}"

"$PY" -m venv "$VENV"
"$VENV/bin/python" -m pip install --quiet --upgrade pip
"$VENV/bin/python" -m pip install -r "$HERE/requirements-embedding.txt"
# mol2vec fork: --no-deps because its setup.py pins gensim==3.8.3
"$VENV/bin/python" -m pip install --no-deps \
  "mol2vec @ git+https://github.com/laserkelvin/mol2vec.git@850d944d5f48a58e26ed0264332b5741f72555aa"

SITE="$("$VENV/bin/python" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
# 1. patched mol2vec/features.py (RDKit sanitize switch, gensim>=4 API); see env/README_env.md for the diff
cp "$HERE/mol2vec_features_patched.py" "$SITE/mol2vec/features.py"
# 2. dask-ml 1.8.0 still uses the np.float / np.int aliases removed from numpy
sed -i 's/np\.float\b/float/g; s/np\.int\b/int/g' "$SITE/dask_ml/_typing.py"

"$VENV/bin/python" -c "import sklearn, gensim, mol2vec.features; print('embedding environment OK, scikit-learn', sklearn.__version__)"
