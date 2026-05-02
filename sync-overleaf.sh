#!/usr/bin/env bash
# sync-overleaf.sh
# Copy fit105 tables and optim figures into ../overleaf-takeup/.
#
# Usage (from project root):
#   bash sync-overleaf.sh

set -euo pipefail

VERSION=105
OVERLEAF="../overleaf-takeup"

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
