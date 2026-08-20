#!/usr/bin/env bash
# run_optim_tables.sh
# Run locally after syncing outputs from the cluster (archive/code/policy-v1/slurm_run_optim.sh).
# Checks that all cluster outputs are present, then prints the remaining
# manual steps to regenerate the two optimization tables and update the TeX.
#
# Usage:
#   bash run_optim_tables.sh          # check outputs and print next steps
#   bash run_optim_tables.sh --sync   # rsync outputs from cluster first

set -euo pipefail

VERSION=105
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
DATA_DIR="optim/data/${MODEL}/agg-full-many-pots"
TABLE_DIR="presentations/tables/fit${VERSION}"
FIG_DIR="presentations/optim-figures/fit${VERSION}"
CLUSTER="midway"
CLUSTER_ROOT="~/projects/takeup"

# ── Optional sync from cluster ────────────────────────────────────────────────
if [[ "${1:-}" == "--sync" ]]; then
    echo "Syncing outputs from ${CLUSTER}..."

    mkdir -p "${FIG_DIR}" "${DATA_DIR}"

    # Summary CSV
    rsync -avz "${CLUSTER}:${CLUSTER_ROOT}/${DATA_DIR}/posterior-clean-summ-optim.csv" \
        "${DATA_DIR}/"

    # Panel figure (written directly to presentations/optim-figures/ on cluster)
    rsync -avz "${CLUSTER}:${CLUSTER_ROOT}/presentations/optim-figures/panel-scenarios-compare-optimal-allocation-plot-distconstraint-3500.pdf" \
        "${FIG_DIR}/"

    # Demand-vstar figure (misc-optim-plots writes to DATA_DIR on cluster)
    rsync -avz "${CLUSTER}:${CLUSTER_ROOT}/${DATA_DIR}/plot-scaled-${MODEL}-agg-identity-full-many-pots-pred-demand-vstar-comp-all.pdf" \
        "${FIG_DIR}/"

    # Comp-dist figure (create-presentation-plots writes to DATA_DIR with distconstraint in name)
    rsync -avz \
        "${CLUSTER}:${CLUSTER_ROOT}/${DATA_DIR}/comp-dist-plot3-fit${VERSION}-distconstraint-3500-util-identity-${MODEL}.pdf" \
        "${FIG_DIR}/comp-dist-plot3-fit${VERSION}-util-identity-${MODEL}.pdf"

    echo "Sync complete."
    echo ""
fi

# ── Check that cluster outputs were synced ────────────────────────────────────
echo "Checking for synced cluster outputs..."

OK=1
check() {
    if [[ -f "$1" ]]; then
        echo "  [ok]  $1"
    else
        echo "  [MISSING]  $1"
        OK=0
    fi
}

check "${DATA_DIR}/posterior-clean-summ-optim.csv"
check "${FIG_DIR}/comp-dist-plot3-fit${VERSION}-util-identity-${MODEL}.pdf"
check "${FIG_DIR}/panel-scenarios-compare-optimal-allocation-plot-distconstraint-3500.pdf"
check "${FIG_DIR}/plot-scaled-${MODEL}-agg-identity-full-many-pots-pred-demand-vstar-comp-all.pdf"

if [[ $OK -eq 0 ]]; then
    echo ""
    echo "Missing outputs — sync from cluster first:"
    echo "  bash run_optim_tables.sh --sync"
    exit 1
fi

echo ""
echo "All cluster outputs present. Generating tables..."
echo ""

Rscript --no-save --no-restore presentations/create-optim-tables.R \
    --optim-input-path="${DATA_DIR}" \
    --table-output-path="${TABLE_DIR}" \
    --model="${MODEL}"

echo ""
echo "══════════════════════════════════════════════════════════════"
echo " Tables written — remaining manual steps"
echo "══════════════════════════════════════════════════════════════"
echo ""
echo "1. Review/diff the generated fit${VERSION} table against the prior TeX-referenced"
echo "   filename before replacing anything:"
echo "   ${TABLE_DIR}/optim-summ-table.tex"
echo "   presentations/tables/manual-optim-summ-table.tex"
echo ""
echo "2. Mark all optimization rows as 'done' in:"
echo "   docs/legacy-structural-workflows.csv"
echo ""
echo "══════════════════════════════════════════════════════════════"
