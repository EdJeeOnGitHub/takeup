#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=policy-model-robust
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=28G
#SBATCH --time=04:00:00
#SBATCH --output=temp/log/policy-model-robust-%A-%a.log
#SBATCH --error=temp/log/policy-model-robust-%A-%a.log

set -euo pipefail

STAGE=${STAGE:?Set STAGE to prepare, predict, optimize, or summarize}
MODEL_ID=${MODEL_ID:?Set MODEL_ID}
NUM_CORES=${NUM_CORES:-12}
optimizer_cpu_budget=${SLURM_CPUS_PER_TASK:-${NUM_CORES}}
OPTIMIZE_CORES=${OPTIMIZE_CORES:-$(( optimizer_cpu_budget < 8 ? optimizer_cpu_budget : 8 ))}
DRAW_BATCH_SIZE=${DRAW_BATCH_SIZE:-25}
SOLVER_THREADS=${SOLVER_THREADS:-1}
SOLVER_SEED=${SOLVER_SEED:-0}
POLICY_SOLVER=${POLICY_SOLVER:-gurobi}
POLICY_SCRATCH=${POLICY_SCRATCH:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

MAX_DRAWS=${MAX_DRAWS:-0}
ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
STREAMLINED_ROOT=${STREAMLINED_ROOT:-${ROOT}/streamlined-active-robustness}
OUTPUT_PATH=${OUTPUT_PATH:-/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-model-robustness/${MODEL_ID}}
COMPACT_CSV=${OUTPUT_PATH}/compact-policy-draws.csv
TARGET_CSV=${TARGET_CSV:-/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv}
DISTANCE_DATA=${DISTANCE_DATA:-/project/akaring/takeup-data/optim/data/full-many-pots-experiment.rds}
REPO_ROOT=${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}
DISTANCE_CAP=${DISTANCE_CAP:-3500}
POPULATION_WEIGHTING=${POPULATION_WEIGHTING:-equal-community}
TARGET_MODE=${TARGET_MODE:-legacy-fixed}
CENSUS_DATA=${CENSUS_DATA:-${POLICY_CENSUS:-${REPO_ROOT}/data/takeup_census.RData}}
export POLICY_CENSUS="$CENSUS_DATA"
if [[ "$TARGET_MODE" == draw-specific-experimental-control ]]; then
  TARGET_CSV="${OUTPUT_PATH}/policy-experimental-targets.csv"
fi

use_streamlined_fits() {
  local spec_id=$1
  FITS=()
  for chain in 1 2 3 4; do
    FITS+=("${STREAMLINED_ROOT}/${spec_id}/fits/${spec_id}-slim-chain${chain}-1.csv")
  done
}

case "${MODEL_ID}" in
  cluster-weighted)
    MODEL_LABEL="Exponential cluster-weighted modes"
    MODEL_FAMILY=gaussian
    LAMBDA_STRUCTURE=common
    FITS=()
    EXTRACT_OPTIONS=()
    ;;
  benchmark)
    MODEL_LABEL="Benchmark"
    MODEL_FAMILY=gaussian
    LAMBDA_STRUCTURE=common
    FITS=()
    EXTRACT_OPTIONS=()
    ;;
  private-distance-community-image)
    MODEL_LABEL="Individual travel costs"
    MODEL_FAMILY=private_distance_community_image
    LAMBDA_STRUCTURE=common
    use_streamlined_fits private-distance-community-image
    EXTRACT_OPTIONS=()
    ;;
  full-information)
    MODEL_LABEL="Individual distance observed by peers"
    MODEL_FAMILY=full_information
    LAMBDA_STRUCTURE=common
    use_streamlined_fits full-information
    EXTRACT_OPTIONS=()
    ;;
  exclude-dispersed)
    MODEL_LABEL="Excluding geographically dispersed communities"
    MODEL_FAMILY=gaussian
    LAMBDA_STRUCTURE=common
    use_streamlined_fits exclude-dispersed
    EXTRACT_OPTIONS=()
    ;;
  tight-multinomial)
    MODEL_LABEL="Correct classification of take-up"
    MODEL_FAMILY=asymmetric_conditional
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/main-core-report-distance-priors/tight/dist_fit106_MAIN_CORE_chain{1,2,3,4}-1.csv)
    EXTRACT_OPTIONS=(--include-asymmetric)
    ;;
  finite-mixture)
    MODEL_LABEL="Mixture v distribution"
    MODEL_FAMILY=finite_mixture
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/main-core-finite-mixture-robustness-800/fits/finite-mixture/finite-mixture-chain{1,2,3,4}-1.csv)
    EXTRACT_OPTIONS=(--include-finite-mixture)
    ;;
  correct-observability)
    MODEL_LABEL="Correct classification"
    MODEL_FAMILY=gaussian
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/dist_fit106_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS-{1,2,3,4}.csv)
    EXTRACT_OPTIONS=()
    ;;
  second-order-observability)
    MODEL_LABEL="Perceived observability (second order)"
    MODEL_FAMILY=gaussian
    LAMBDA_STRUCTURE=common
    use_streamlined_fits second-order-observability
    EXTRACT_OPTIONS=(--beliefs-order 2)
    ;;
  grouped-lambda)
    MODEL_LABEL="Any-signal/no-signal social image weight"
    MODEL_FAMILY=gaussian
    LAMBDA_STRUCTURE=grouped
    FITS=("${ROOT}"/main-core-lambda-identification/fits/grouped-sd0p25/grouped-sd0p25-chain{1,2,3,4}-1.csv)
    EXTRACT_OPTIONS=(--include-lambda grouped)
    ;;
  arm-lambda)
    MODEL_LABEL="Treatment-specific social image weight"
    MODEL_FAMILY=gaussian
    LAMBDA_STRUCTURE=arm
    FITS=("${ROOT}"/main-core-lambda-identification/fits/arm-sd0p25/arm-sd0p25-chain{1,2,3,4}-1.csv)
    EXTRACT_OPTIONS=(--include-lambda arm)
    ;;
  student-t5)
    MODEL_LABEL="Student-t(5) intrinsic-motivation types"
    MODEL_FAMILY=student_t5
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/main-core-student-t-robustness/fits/student-t5/student-t5-chain{1,2,3,4}-1.csv)
    EXTRACT_OPTIONS=()
    ;;
  cluster-shock)
    MODEL_LABEL="Cluster random shock"
    MODEL_FAMILY=gaussian
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/main-core-cluster-shock-production/dist_fit106_MAIN_CORE_SHOCK_SD0.1_chain{1,2,3,4}-1.csv)
    EXTRACT_OPTIONS=(--include-cluster-shock 144)
    ;;
  asymmetric-conditional)
    MODEL_LABEL="Asymmetric reports, conditional on recognition"
    MODEL_FAMILY=asymmetric_conditional
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/main-core-asym-conditional-production/dist_fit106_MAIN_CORE_chain{1,2,3,4}-1.csv)
    EXTRACT_OPTIONS=(--include-asymmetric)
    ;;
  asymmetric-unconditional)
    MODEL_LABEL="Asymmetric reports, unrecognized as null signal"
    MODEL_FAMILY=asymmetric_unconditional
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/main-core-asym-unconditional-production/dist_fit106_MAIN_CORE_chain{1,2,3,4}-1.csv)
    EXTRACT_OPTIONS=(--include-asymmetric)
    ;;
  asymmetric-f1|asymmetric-f2|asymmetric-f3|asymmetric-u3)
    LADDER_ID=${MODEL_ID#asymmetric-}
    MODEL_LABEL="Observability ladder ${LADDER_ID^^}"
    MODEL_FAMILY="asymmetric_${LADDER_ID}"
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/main-core-observability-ladder/"${LADDER_ID}"/dist_fit106_MAIN_CORE_chain{1,2,3,4}-1.csv)
    EXTRACT_OPTIONS=(--asymmetric-structure "${LADDER_ID}")
    ;;
  *)
    echo "Unknown MODEL_ID=${MODEL_ID}" >&2
    exit 2
    ;;
esac

cd "${REPO_ROOT}"
module load -f R/4.2.0
if [[ -n "${POLICY_GUROBI_ROOT:-}" ]]; then
  export GUROBI_HOME="$POLICY_GUROBI_ROOT"
  export GRB_LICENSE_FILE="$POLICY_GUROBI_ROOT/gurobi.lic"
else
  export GUROBI_HOME="${HOME}/gurobi952/linux64"
  if [[ "${SLURM_JOB_PARTITION:-}" == caslake ]]; then
    module load -f gurobi/11.0
    GUROBI_HOME=$(dirname "$(dirname "$(command -v gurobi_cl)")")
  else
    export GRB_LICENSE_FILE=/software/gurobi-9.2-el7-x86_64/gurobi.lic
  fi
fi
export PATH="${GUROBI_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
mkdir -p temp/log "${OUTPUT_PATH}"

case "${STAGE}" in
  prepare)
    if [[ "$MODEL_ID" == benchmark ]]; then
      Rscript --vanilla scripts/policy/prepare-baseline-posterior.R \
        "--fit-path=${BENCHMARK_FIT_PATH:?Set BENCHMARK_FIT_PATH to assigned slim chains}" \
        "--draws-per-chain=${BENCHMARK_DRAWS_PER_CHAIN:-400}" \
        "--distance-data=${DISTANCE_DATA}" "--output-path=${OUTPUT_PATH}"
      exit 0
    elif [[ "$MODEL_ID" == cluster-weighted ]]; then
      Rscript --vanilla scripts/policy/prepare-cluster-bootstrap.R \
        "--weighted-path=${WEIGHTED_PATH:?Set assigned-distance WEIGHTED_PATH}" \
        "--output-path=${OUTPUT_PATH}" "--num-replicates=${NUM_REPLICATES:-999}" \
        --method=exponential --distance-definition=assigned
      Rscript --vanilla scripts/policy/standardize-weighted-modes.R "--input-path=${OUTPUT_PATH}"
      exit 0
    fi
    python3 scripts/policy/extract-cmdstan-draws.py \
      --output "${COMPACT_CSV}" \
      ${EXTRACT_OPTIONS[@]+"${EXTRACT_OPTIONS[@]}"} "${FITS[@]}"
    Rscript --no-save --no-restore scripts/policy/prepare-model-robustness.R \
      "--input-csv=${COMPACT_CSV}" "--output-path=${OUTPUT_PATH}" \
      "--model-id=${MODEL_ID}" "--model-label=${MODEL_LABEL}" \
      "--model-family=${MODEL_FAMILY}" "--lambda-structure=${LAMBDA_STRUCTURE}" \
      "--workspace=${CLUSTER_WORKSPACE:-${REPO_ROOT}/data/stan_analysis_data/dist_fit106.RData}"
    ;;
  predict)
    PREDICT_OPTIONS=()
    if [[ "${MODEL_ID}" == "private-distance-community-image" || "${MODEL_ID}" == "full-information" ]]; then
      PREDICT_OPTIONS+=(--household-workspace="${ROOT}/dist_fit104.RData")
    fi
    Rscript --no-save --no-restore scripts/policy/predict-model-robustness.R \
      "--parameter-rds=${OUTPUT_PATH}/policy-model-parameters.rds" \
      "--distance-data=${DISTANCE_DATA}" "--output-path=${OUTPUT_PATH}" \
      "--distance-cap=${DISTANCE_CAP}" \
      "--population-weighting=${POPULATION_WEIGHTING}" "--census-data=${CENSUS_DATA}" \
      "--num-cores=${NUM_CORES}" "--max-draws=${MAX_DRAWS}" \
      ${PREDICT_OPTIONS[@]+"${PREDICT_OPTIONS[@]}"}
    ;;
  optimize)
    : "${SLURM_ARRAY_TASK_ID:?Optimize requires scenario array 1-5}"
    Rscript --no-save --no-restore scripts/policy/optimize-cluster-bootstrap.R \
      "--input-path=${OUTPUT_PATH}" "--target-csv=${TARGET_CSV}" \
      "--population-weighting=${POPULATION_WEIGHTING}" "--target-mode=${TARGET_MODE}" \
      "--distance-data=${DISTANCE_DATA}" "--census-data=${CENSUS_DATA}" \
      "--allocation-root=${ALLOCATION_ROOT:-${OUTPUT_PATH}}" \
      "--draw-start=${DRAW_START:-1}" "--draw-end=${DRAW_END:-100000}" \
      "--num-cores=${OPTIMIZE_CORES}" "--draw-batch-size=${DRAW_BATCH_SIZE}" \
      "--solver=${POLICY_SOLVER}" "--solver-threads=${SOLVER_THREADS}" \
      "--solver-seed=${SOLVER_SEED}" "--scratch-path=${POLICY_SCRATCH}" \
      "--time-limit=${SOLVER_TIME_LIMIT:-300}" \
      "--scenario-id=${SLURM_ARRAY_TASK_ID}" --num-replicates=100000
    if [[ -n "${ALLOCATION_ROOT:-}" ]]; then
      python3 scripts/policy/archive-policy-solver-logs.py \
        "--run-root=${ALLOCATION_ROOT}" \
        "--archive=${ALLOCATION_ROOT}/solver-logs-${SLURM_JOB_ID:-$$}.tar.gz"
      touch "${ALLOCATION_ROOT}/_SUCCESS"
    fi
    ;;
  collect)
    Rscript --vanilla scripts/policy/collect-policy-shards.R \
      "--input-path=${OUTPUT_PATH}" "--shard-root=${SHARD_ROOT:?Set SHARD_ROOT}"
    ;;
  summarize)
    Rscript --no-save --no-restore scripts/policy/summarize-model-results.R \
      "--input-path=${OUTPUT_PATH}"
    {
      echo "model_id=${MODEL_ID}"
      echo "model_label=${MODEL_LABEL}"
      echo "distance_definition=assigned"
      echo "git_commit=${POLICY_CODE_REVISION:-$(git rev-parse HEAD)}"
      echo "source_fit_directory=$(dirname "${FITS[0]:-prepared-balanced-assigned-distance-slim-chains}")"
      echo "source_fit_files=$(IFS=,; echo "${FITS[*]:-prepared-balanced-assigned-distance-slim-chains}")"
      echo "extract_options=${EXTRACT_OPTIONS[*]:-}"
      if [[ ${FITS[0]:-} == "${STREAMLINED_ROOT}"/* ]]; then
        echo "streamlined_active_robustness=true"
      else
        echo "streamlined_active_robustness=false"
      fi
      echo "structural_refit_performed=no"
      echo "candidate_sites=1451"
      echo "distance_cap_m=${DISTANCE_CAP}"
      echo "population_weighting=${POPULATION_WEIGHTING}"
      echo "target_mode=${TARGET_MODE}"
      echo "generated_utc=$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
    } > "${OUTPUT_PATH}/provenance.txt"
    ;;
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
