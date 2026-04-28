#!/usr/bin/env bash
# run_optim_tables.sh
# Run locally after syncing outputs from the cluster (slurm_run_optim.sh).
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
TABLE_DIR="presentations/tables/fit105"
CLUSTER="midway"
CLUSTER_ROOT="~/projects/takeup"

# ── Optional sync from cluster ────────────────────────────────────────────────
if [[ "${1:-}" == "--sync" ]]; then
    echo "Syncing outputs from ${CLUSTER}..."

    mkdir -p presentations/figures presentations/misc-figures presentations/optim-figures "${DATA_DIR}"

    # Summary CSV
    rsync -avz "${CLUSTER}:${CLUSTER_ROOT}/${DATA_DIR}/posterior-clean-summ-optim.csv" \
        "${DATA_DIR}/"

    # Panel figure (written directly to presentations/optim-figures/ on cluster)
    rsync -avz "${CLUSTER}:${CLUSTER_ROOT}/presentations/optim-figures/" \
        presentations/optim-figures/

    # Demand-vstar figure (misc-optim-plots writes to DATA_DIR, not presentations/figures/)
    rsync -avz "${CLUSTER}:${CLUSTER_ROOT}/${DATA_DIR}/plot-scaled-${MODEL}-agg-identity-full-many-pots-pred-demand-vstar-comp-all.pdf" \
        presentations/figures/

    # Comp-dist figure (create-presentation-plots writes to DATA_DIR with distconstraint in name)
    rsync -avz \
        "${CLUSTER}:${CLUSTER_ROOT}/${DATA_DIR}/comp-dist-plot3-fit${VERSION}-distconstraint-3500-util-identity-${MODEL}.pdf" \
        "presentations/misc-figures/comp-dist-plot3-fit${VERSION}-util-identity-${MODEL}.pdf"

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
check "presentations/misc-figures/comp-dist-plot3-fit${VERSION}-util-identity-${MODEL}.pdf"
check "presentations/optim-figures/panel-scenarios-compare-optimal-allocation-plot-distconstraint-3500.pdf"
check "presentations/figures/plot-scaled-${MODEL}-agg-identity-full-many-pots-pred-demand-vstar-comp-all.pdf"

if [[ $OK -eq 0 ]]; then
    echo ""
    echo "Missing outputs — sync from cluster first:"
    echo "  rsync -avz midway:~/projects/takeup/optim/data/${MODEL}/agg-full-many-pots/posterior-clean-summ-optim.csv \\"
    echo "    ${DATA_DIR}/"
    echo "  rsync -avz midway:~/projects/takeup/presentations/misc-figures/ presentations/misc-figures/"
    echo "  rsync -avz midway:~/projects/takeup/presentations/optim-figures/ presentations/optim-figures/"
    echo "  rsync -avz midway:~/projects/takeup/presentations/figures/ presentations/figures/"
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
echo "1. Review/diff the generated fit105 table against the prior TeX-referenced"
echo "   filename before replacing anything:"
echo "   ${TABLE_DIR}/optim-summ-table.tex"
echo "   presentations/tables/manual-optim-summ-table.tex"
echo ""
echo "2. Update 'ECM ReStud.tex' line 1402 — comp-dist-plot fit version:"
echo "   FROM: {misc-figures/comp-dist-plot3-fit86-util-identity-${MODEL}.pdf}"
echo "   TO:   {misc-figures/comp-dist-plot3-fit${VERSION}-util-identity-${MODEL}.pdf}"
echo ""
echo "3. Mark all optimization rows as 'done' in:"
echo "   recreate_structural_robustness_and_optimization.csv"
echo ""
echo "══════════════════════════════════════════════════════════════"
