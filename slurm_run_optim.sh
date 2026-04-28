#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=takeup-optim-105
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28G
#SBATCH --time=36:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/optim-105-%j.log
#SBATCH --error=temp/log/optim-105-%j.log
#SBATCH --export=IN_SLURM=1

set -euo pipefail

# ── Smoketest mode ────────────────────────────────────────────────────────────
# Set SMOKETEST=1 to run a single scenario with 5 draws at dist=3500 only.
# Use for interactive testing before submitting the full sbatch job:
#   sinteractive --partition=broadwl --cpus-per-task=4 --mem=20G --time=01:00:00
#   cd ~/projects/takeup && module load gurobi/9.2 R/4.2.0
#   SMOKETEST=1 bash slurm_run_optim.sh
SMOKETEST=${SMOKETEST:-0}

# ── Configuration ─────────────────────────────────────────────────────────────
VERSION=105
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
NUM_CORES=12          # used for sequential steps (gurobi, postprocess)
NUM_CORES_PARALLEL=1  # Step 2b already runs scenarios in parallel; avoid nested R workers
PREDICT_PARALLELISM="${PREDICT_PARALLELISM:-1}"
WELFARE="identity"
CTYPE="agg"
COUNTY="full"
SOLVER="gurobi"
STATIC_DIST=500

if [[ $SMOKETEST -eq 1 ]]; then
    NUM_POST_DRAWS=5
    DIST_CONSTRAINTS=(3500)
    SCENARIOS=("control:control:")
    MEDIAN_SCENARIOS=("control:control:")
    echo "[smoketest] 1 scenario, 5 draws, dist=3500 only"
else
    NUM_POST_DRAWS=200
    DIST_CONSTRAINTS=(2500 3500 4500 5500 10000)
    SCENARIOS=(
        "control:control:"
        "control:bracelet:"
        "control:control:static-"
        "control:bracelet:static-"
        "control:control:suppress-rep-"
    )
    MEDIAN_SCENARIOS=(
        "control:control:"
        "control:bracelet:"
        "control:bracelet:static-"
        "control:control:suppress-rep-"
    )
fi

STAN_INPUT="${STAN_INPUT:-/project/akaring/takeup-data/data/stan_analysis_data}"
DATA_INPUT="${COUNTY}-many-pots-experiment.rds"
DATA_DIR="optim/data/${MODEL}/${CTYPE}-${COUNTY}-many-pots"
PLOT_DIR="optim/plots/${MODEL}/${CTYPE}-${COUNTY}-many-pots"
PDF_DIR="presentations/misc-figures"
OPTIM_FIG_DIR="presentations/optim-figures"

# ── Environment ───────────────────────────────────────────────────────────────
cd ~/projects/takeup
module load -f gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0
export GUROBI_HOME="${HOME}/gurobi952/linux64"
export PATH="${GUROBI_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
# Use cluster token server license (works on compute nodes; not login node)
export GRB_LICENSE_FILE=/software/gurobi-9.2-el7-x86_64/gurobi.lic
# scratch-midway2 is also noexec at GPFS policy level; home is exec-allowed
export TMPDIR="${TMPDIR:-~/tmp-rcpp}"
mkdir -p "${TMPDIR}" "${DATA_DIR}" "${PLOT_DIR}" "${PDF_DIR}" "${OPTIM_FIG_DIR}" temp/log

# ── Helper: run an R script, log to file, echo progress ───────────────────────
rrun() {
    local label=$1
    local logfile="temp/log/optim-105-${label}.log"
    shift
    echo "[$(date +%H:%M:%S)] ${label}"
    Rscript --no-save --no-restore "$@" >> "${logfile}" 2>&1
}

# ── Step 1: Village/PoT distance data ─────────────────────────────────────────
if [[ ! -f "optim/data/${DATA_INPUT}" ]]; then
    rrun "create-distance-data" optim/create-distance-data.R \
        --output-name="${DATA_INPUT}" \
        --num-extra-pots=100 \
        --county-subset="${COUNTY}" \
        --distance-cutoff=Inf
else
    echo "[$(date +%H:%M:%S)] create-distance-data skipped (optim/data/${DATA_INPUT} exists)"
fi

# ── Step 2a: Extract posterior parameter draws once (1.7GB CSV → param CSV) ───
# Loading as_cmdstan_fit on a 1.7GB chain uses ~15GB RAM; do this once rather
# than once per scenario to avoid OOM and a write-race on the shared param CSV.
PARAM_CSV="${STAN_INPUT}/param_posterior_draws_dist_fit${VERSION}_${MODEL}.csv"
if [[ ! -f "${PARAM_CSV}" ]]; then
    echo "[$(date +%H:%M:%S)] Step 2a: extracting param draws to CSV..."
    rrun "extract-params" optim/predict-takeup-for-optim.R \
        "${VERSION}" "control" "control" \
        --output-name="cutoff-b-control-mu-control-${MODEL}" \
        --to-csv \
        --num-post-draws="${NUM_POST_DRAWS}" \
        --rep-cutoff=Inf \
        --dist-cutoff=3500 \
        --num-cores="${NUM_CORES}" \
        --type-lb=-Inf --type-ub=Inf \
        --data-input-name="${DATA_INPUT}" \
        --output-path="${DATA_DIR}" \
        --input-path="${STAN_INPUT}" \
        --model="${MODEL}" \
        --single-chain
else
    echo "[$(date +%H:%M:%S)] Step 2a: param CSV already exists, skipping extraction"
fi

# ── Step 2b: Demand prediction from pre-extracted CSV ─────────────────────────
echo "[$(date +%H:%M:%S)] Step 2b: demand prediction (${PREDICT_PARALLELISM} scenario process(es) at a time)..."

predict_demand() {
    local b_z=$1 mu_z=$2 prefix=$3 label=$4
    shift 4
    local out_csv="${DATA_DIR}/pred-demand-dist-fit${VERSION}-${prefix}cutoff-b-${b_z}-mu-${mu_z}-${MODEL}.csv"
    if [[ -f "${out_csv}" ]]; then
        echo "[$(date +%H:%M:%S)] predict-${label} skipped (output exists)"
        return
    fi
    Rscript --no-save --no-restore optim/predict-takeup-for-optim.R \
        "${VERSION}" "${b_z}" "${mu_z}" \
        --output-name="${prefix}cutoff-b-${b_z}-mu-${mu_z}-${MODEL}" \
        --from-csv \
        --num-post-draws="${NUM_POST_DRAWS}" \
        --rep-cutoff=Inf \
        --dist-cutoff=3500 \
        --num-cores="${NUM_CORES_PARALLEL}" \
        --type-lb=-Inf --type-ub=Inf \
        --data-input-name="${DATA_INPUT}" \
        --output-path="${DATA_DIR}" \
        --input-path="${STAN_INPUT}" \
        --model="${MODEL}" \
        --run-estimation \
        "$@" \
        >> "temp/log/optim-105-predict-${label}.log" 2>&1
}

predict_pids=()
wait_for_predict_slot() {
    if [[ ${#predict_pids[@]} -lt ${PREDICT_PARALLELISM} ]]; then
        return
    fi

    wait "${predict_pids[0]}"
    predict_pids=("${predict_pids[@]:1}")
}

wait_for_predict_slot
predict_demand "control" "control" "" "control-control" &
predict_pids+=("$!")
if [[ $SMOKETEST -eq 0 ]]; then
    wait_for_predict_slot
    predict_demand "control" "bracelet" "" "control-bracelet" &
    predict_pids+=("$!")
    wait_for_predict_slot
    predict_demand "control" "control"  "static-"     "static-control"  \
        --static-signal-pm --static-signal-distance="${STATIC_DIST}"    &
    predict_pids+=("$!")
    wait_for_predict_slot
    predict_demand "control" "bracelet" "static-"     "static-bracelet" \
        --static-signal-pm --static-signal-distance="${STATIC_DIST}"    &
    predict_pids+=("$!")
    wait_for_predict_slot
    predict_demand "control" "control"  "suppress-rep-" "suppress-rep"  \
        --suppress-reputation                                            &
    predict_pids+=("$!")
fi
for pid in "${predict_pids[@]}"; do
    wait "${pid}"
done
echo "[$(date +%H:%M:%S)] Demand prediction complete."

# ── Step 3: Experiment target constraint (posterior median, base scenario) ────
EXPERIMENT_TARGET_CSV="${DATA_DIR}/summ-${CTYPE}-${WELFARE}-experiment-target-constraint.csv"
if [[ ! -f "${EXPERIMENT_TARGET_CSV}" ]]; then
    rrun "experiment-target" optim/create-experiment-target.R \
        --constraint-type="${CTYPE}" \
        --welfare-function="${WELFARE}" \
        --min-cost \
        --output-path="${DATA_DIR}" \
        --output-basename="summ-${CTYPE}-${WELFARE}" \
        --cutoff-type=cutoff \
        --data-input-name="${DATA_INPUT}" \
        --posterior-median \
        --demand-input-path="${DATA_DIR}" \
        --demand-input-filename="pred-demand-dist-fit${VERSION}-cutoff-b-control-mu-control-${MODEL}.csv"
else
    echo "[$(date +%H:%M:%S)] experiment-target skipped (${EXPERIMENT_TARGET_CSV} exists)"
fi

# ── Step 4: Optimize + postprocess (full posterior draws) ─────────────────────
echo "[$(date +%H:%M:%S)] Step 4: Gurobi optimization (post-draws)..."

for SCENARIO in "${SCENARIOS[@]}"; do
    IFS=':' read -r b_z mu_z prefix <<< "${SCENARIO}"
    demand_file="pred-demand-dist-fit${VERSION}-${prefix}cutoff-b-${b_z}-mu-${mu_z}-${MODEL}.csv"
    label="${prefix}b-${b_z}-mu-${mu_z}"

    for dist in "${DIST_CONSTRAINTS[@]}"; do
        subdir="${DATA_DIR}/dist-constraint-${dist}"
        alloc_base="target-rep-distconstraint-${dist}-util-${WELFARE}-${prefix}cutoff-b-${b_z}-mu-${mu_z}-${MODEL}"
        mkdir -p "${subdir}"
        cp -n "${DATA_DIR}/${demand_file}" "${subdir}/" 2>/dev/null || true

        if [[ -f "${subdir}/${alloc_base}-post-draws-optimal-allocation.rds" ]]; then
            echo "[$(date +%H:%M:%S)]   optimize+postprocess ${label} @ ${dist}m skipped (output exists)"
            continue
        fi

        echo "[$(date +%H:%M:%S)]   optimize ${label} @ ${dist}m"
        Rscript --no-save --no-restore optim/optimal_allocation.R \
            --num-cores="${NUM_CORES}" \
            --min-cost \
            --constraint-type="${CTYPE}" \
            --target-constraint="summ-${CTYPE}-${WELFARE}-experiment-target-constraint.csv" \
            --output-path="${subdir}" \
            --output-filename="${alloc_base}" \
            --input-path="${DATA_DIR}" \
            --data-input-path="optim/data" \
            --data-input-name="${DATA_INPUT}" \
            --time-limit=10000 \
            --demand-input-filename="${demand_file}" \
            --welfare-function="${WELFARE}" \
            --solver="${SOLVER}" \
            --distance-constraint="${dist}" \
            >> "temp/log/optim-105-gurobi-${label}-dist${dist}.log" 2>&1

        echo "[$(date +%H:%M:%S)]   postprocess ${label} @ ${dist}m"
        Rscript --no-save --no-restore optim/postprocess_allocation.R \
            --min-cost \
            --constraint-type="${CTYPE}" \
            --welfare-function="${WELFARE}" \
            --optim-input-path="${subdir}" \
            --optim-input-a-filename="${alloc_base}-post-draws-optimal-allocation.rds" \
            --data-input-name="${DATA_INPUT}" \
            --output-path="${PLOT_DIR}" \
            --output-basename="${CTYPE}-${alloc_base}-post-draws" \
            --cutoff-type=cutoff \
            --pdf-output-path="${PDF_DIR}" \
            >> "temp/log/optim-105-postprocess-${label}-dist${dist}.log" 2>&1
    done
done

# ── Step 5: Optimize (posterior median) for figure generation ─────────────────
echo "[$(date +%H:%M:%S)] Step 5: Gurobi optimization (posterior median for figures)..."

for SCENARIO in "${MEDIAN_SCENARIOS[@]}"; do
    IFS=':' read -r b_z mu_z prefix <<< "${SCENARIO}"
    dist=3500
    subdir="${DATA_DIR}/dist-constraint-${dist}"
    demand_file="pred-demand-dist-fit${VERSION}-${prefix}cutoff-b-${b_z}-mu-${mu_z}-${MODEL}.csv"
    alloc_base="target-rep-distconstraint-${dist}-util-${WELFARE}-${prefix}cutoff-b-${b_z}-mu-${mu_z}-${MODEL}"
    label="${prefix}b-${b_z}-mu-${mu_z}"

    if [[ -f "${subdir}/${alloc_base}-median-optimal-allocation.rds" ]]; then
        echo "[$(date +%H:%M:%S)]   median optimize ${label} skipped (output exists)"
        continue
    fi

    echo "[$(date +%H:%M:%S)]   median optimize ${label}"
    Rscript --no-save --no-restore optim/optimal_allocation.R \
        --posterior-median \
        --num-cores="${NUM_CORES}" \
        --min-cost \
        --constraint-type="${CTYPE}" \
        --target-constraint="summ-${CTYPE}-${WELFARE}-experiment-target-constraint.csv" \
        --output-path="${subdir}" \
        --output-filename="${alloc_base}" \
        --input-path="${DATA_DIR}" \
        --data-input-path="optim/data" \
        --data-input-name="${DATA_INPUT}" \
        --time-limit=10000 \
        --demand-input-filename="${demand_file}" \
        --welfare-function="${WELFARE}" \
        --solver="${SOLVER}" \
        --distance-constraint="${dist}" \
        >> "temp/log/optim-105-median-gurobi-${label}.log" 2>&1
done

# ── Step 6: Aggregate all scenarios into posterior-clean-summ-optim.csv ───────
rrun "compare-optim" optim/compare-optim.R \
    --input-path="${DATA_DIR}" \
    --output-path="${DATA_DIR}" \
    --many-pots \
    --model="${MODEL}" \
    --welfare-function="${WELFARE}"

# ── Step 7: Figures ───────────────────────────────────────────────────────────
echo "[$(date +%H:%M:%S)] Step 7: generating figures..."

# Figure 1: Demand curves under alternative social image assumptions
# → presentations/figures/plot-scaled-MODEL-agg-identity-full-many-pots-pred-demand-vstar-comp-all.pdf
rrun "misc-optim-plots" optim/misc-optim-plots.R \
    --output-path="${DATA_DIR}" \
    --model="${MODEL}" \
    --fit-version="${VERSION}" \
    --welfare-function="${WELFARE}"

# Figure 2: Optimal placement panel
# → presentations/optim-figures/panel-scenarios-compare-optimal-allocation-plot-distconstraint-3500.pdf
rrun "optim-paper-panel" optim/create-optim-paper-panel.R \
    --output-path="${OPTIM_FIG_DIR}" \
    --model="${MODEL}" \
    --fit-version="${VERSION}" \
    --welfare-function="${WELFARE}" \
    --input-path="${DATA_DIR}" \
    --distance-constraint=3500

# Figure 3: Distribution of community-PoT distances under alternative allocations
# → presentations/misc-figures/comp-dist-plot3-fit105-util-identity-MODEL.pdf
# Writes experimental-control-allocation-data.rds (needed by compare-optim.R if re-run)
rrun "presentation-plots" optim/create-presentation-plots.R \
    --constraint-type="${CTYPE}" \
    --welfare-function="${WELFARE}" \
    --min-cost \
    --output-path="${DATA_DIR}/dist-constraint-3500" \
    --output-basename="target-rep-agg-${WELFARE}-cutoff-b-control-mu-control-${MODEL}-median" \
    --cutoff-type=cutoff \
    --data-input-name="${DATA_INPUT}" \
    --posterior-median \
    --pdf-output-path="${PDF_DIR}" \
    --demand-input-path="${DATA_DIR}" \
    --demand-input-filename="pred-demand-dist-fit${VERSION}-cutoff-b-control-mu-control-${MODEL}.csv" \
    --model="${MODEL}" \
    --fit-version="${VERSION}" \
    --distance-constraint=3500

echo "[$(date +%H:%M:%S)] Done."
echo ""
echo "Outputs on cluster:"
echo "  Summary CSV:   ${DATA_DIR}/posterior-clean-summ-optim.csv"
echo "  Demand figure: presentations/figures/plot-scaled-${MODEL}-agg-identity-full-many-pots-pred-demand-vstar-comp-all.pdf"
echo "  Panel figure:  ${OPTIM_FIG_DIR}/panel-scenarios-compare-optimal-allocation-plot-distconstraint-3500.pdf"
echo "  Dist figure:   ${PDF_DIR}/comp-dist-plot3-fit${VERSION}-util-identity-${MODEL}.pdf"
echo ""
echo "Sync to local, then run: bash run_optim_tables.sh"
