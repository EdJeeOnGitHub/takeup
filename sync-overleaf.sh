#!/usr/bin/env bash
# sync-overleaf.sh
# Copy fit tables and figures into ../overleaf/overleaf-takeup/.
#
# Usage (from project root):
#   bash sync-overleaf.sh
#   bash sync-overleaf.sh 106

set -euo pipefail

VERSION="${1:-105}"
OVERLEAF="../overleaf/overleaf-takeup"

TABLE_SRC="presentations/tables/fit${VERSION}"
TABLE_DST="${OVERLEAF}/tables/fit${VERSION}"

FIG_SRC="presentations/optim-figures/fit${VERSION}"
FIG_DST="${OVERLEAF}/optim-figures/fit${VERSION}"

STRUCT_FIG_SRC="presentations/figures/fit${VERSION}"
STRUCT_FIG_DST="${OVERLEAF}/figures/fit${VERSION}"

echo "Syncing fit${VERSION} tables..."
mkdir -p "${TABLE_DST}"
cp "${TABLE_SRC}"/*.tex "${TABLE_DST}/"
echo "  → ${TABLE_DST}/"

echo "Syncing fit${VERSION} optim figures..."
mkdir -p "${FIG_DST}"
cp "${FIG_SRC}"/*.pdf "${FIG_DST}/"
echo "  → ${FIG_DST}/"

echo "Syncing fit${VERSION} structural figures..."
mkdir -p "${STRUCT_FIG_DST}"
cp "${STRUCT_FIG_SRC}"/*.pdf "${STRUCT_FIG_DST}/"
echo "  → ${STRUCT_FIG_DST}/"

echo "Done."
