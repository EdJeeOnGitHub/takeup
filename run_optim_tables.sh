#!/usr/bin/env bash
# run_optim_tables.sh
# Run locally after syncing outputs from the cluster (slurm_run_optim.sh).
# Checks that all cluster outputs are present, then prints the remaining
# manual steps to regenerate the two optimization tables and update the TeX.

set -euo pipefail

VERSION=105
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
DATA_DIR="optim/data/${MODEL}/agg-full-many-pots"
TABLE_DIR="presentations/tables"

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
echo "All cluster outputs present."
echo ""
echo "══════════════════════════════════════════════════════════════"
echo " Remaining steps to complete the optimization results update"
echo "══════════════════════════════════════════════════════════════"
echo ""
echo "1. Enable the two disabled optim chunks in presentations/tables.Rmd:"
echo ""
echo "   a) Chunk 'optim-table-main' (~line 1558):"
echo "      Change:  #| eval=FALSE"
echo "      To:      #| eval=TRUE"
echo ""
echo "   b) Chunk 'optim-summ-table-robust' (~line 1889):"
echo "      Change:  #| eval=FALSE"
echo "      To:      #| eval=TRUE"
echo "      AND uncomment the table-building code inside the chunk."
echo ""
echo "2. Run both chunks (RStudio or knitr::knit). They write:"
echo "     ${TABLE_DIR}/optim-summ-table.tex"
echo "     ${TABLE_DIR}/optim-summ-robust-table.tex"
echo ""
echo "3. Copy the main table to the TeX-referenced filename:"
echo "   cp ${TABLE_DIR}/optim-summ-table.tex \\"
echo "      ${TABLE_DIR}/manual-optim-summ-table.tex"
echo ""
echo "4. Update 'ECM ReStud.tex' line 1402 — comp-dist-plot fit version:"
echo "   FROM: {misc-figures/comp-dist-plot3-fit86-util-identity-${MODEL}.pdf}"
echo "   TO:   {misc-figures/comp-dist-plot3-fit${VERSION}-util-identity-${MODEL}.pdf}"
echo ""
echo "5. Mark all optimization rows as 'done' in:"
echo "   recreate_structural_robustness_and_optimization.csv"
echo ""
echo "══════════════════════════════════════════════════════════════"
