#!/usr/bin/env bash
# Fetch the third-party inputs that are not redistributed in this repository (see data/raw/FETCH.md).
#   external/umda/                          umda (Lee et al. 2021), pinned commit; MIT licence
#   external/tmc1_recommendations.pdf       Lee et al. 2021 recommendation list, Zenodo version 2 (CC BY 4.0), md5-checked
# Usage: bash scripts/00_fetch_external.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXT="$ROOT/external"
UMDA_COMMIT="a95cc1c1eb98a1ff64b37a4c6ed92f9546d0215b"
LEE_URL="https://zenodo.org/records/5146276/files/tmc1_recommendations.pdf?download=1"
LEE_MD5="9d3664ca507b5542e2bf3a279599f891"
mkdir -p "$EXT"

if [ ! -d "$EXT/umda/.git" ]; then
  git clone https://github.com/laserkelvin/umda.git "$EXT/umda"
fi
git -C "$EXT/umda" checkout --quiet "$UMDA_COMMIT"
echo "umda at $(git -C "$EXT/umda" rev-parse HEAD)"

if [ ! -f "$EXT/tmc1_recommendations.pdf" ]; then
  curl -L --fail -o "$EXT/tmc1_recommendations.pdf" "$LEE_URL" \
    || { echo "Automatic download failed: download tmc1_recommendations.pdf from https://doi.org/10.5281/zenodo.5146276 into $EXT manually and re-run."; exit 1; }
fi
if command -v md5sum >/dev/null; then SUM="$(md5sum "$EXT/tmc1_recommendations.pdf" | cut -d' ' -f1)"; else SUM="$(md5 -q "$EXT/tmc1_recommendations.pdf")"; fi
[ "$SUM" = "$LEE_MD5" ] && echo "tmc1_recommendations.pdf md5 OK" || { echo "md5 mismatch for tmc1_recommendations.pdf ($SUM)"; exit 1; }
