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
ROOT=/project/akaring/takeup-data/data/stan_analysis_data
OUTPUT_PATH=${OUTPUT_PATH:-optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-model-robustness/${MODEL_ID}}
COMPACT_CSV=${OUTPUT_PATH}/compact-policy-draws.csv
TARGET_CSV=${TARGET_CSV:-optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv}
REPO_ROOT=${REPO_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}

case "${MODEL_ID}" in
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
    FITS=("${ROOT}"/dist_fit106_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB-{1,2,3,4}.csv)
    EXTRACT_OPTIONS=()
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
  *)
    echo "Unknown MODEL_ID=${MODEL_ID}" >&2
    exit 2
    ;;
esac

cd "${REPO_ROOT}"
module load -f gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0
export GUROBI_HOME="${HOME}/gurobi952/linux64"
export PATH="${GUROBI_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export GRB_LICENSE_FILE=/software/gurobi-9.2-el7-x86_64/gurobi.lic
mkdir -p temp/log "${OUTPUT_PATH}"

case "${STAGE}" in
  prepare)
    python3 optim/extract-cmdstan-policy-draws.py \
      --output "${COMPACT_CSV}" \
      ${EXTRACT_OPTIONS[@]+"${EXTRACT_OPTIONS[@]}"} "${FITS[@]}"
    Rscript --no-save --no-restore optim/prepare-policy-model-robustness.R \
      "--input-csv=${COMPACT_CSV}" "--output-path=${OUTPUT_PATH}" \
      "--model-id=${MODEL_ID}" "--model-label=${MODEL_LABEL}" \
      "--model-family=${MODEL_FAMILY}" "--lambda-structure=${LAMBDA_STRUCTURE}"
    ;;
  predict)
    Rscript --no-save --no-restore optim/predict-policy-model-robustness.R \
      "--parameter-rds=${OUTPUT_PATH}/policy-model-parameters.rds" \
      "--output-path=${OUTPUT_PATH}" --distance-cap=3500 \
      "--num-cores=${NUM_CORES}"
    ;;
  optimize)
    : "${SLURM_ARRAY_TASK_ID:?Optimize requires scenario array 1-5}"
    Rscript --no-save --no-restore optim/optimize-policy-cluster-bootstrap.R \
      "--input-path=${OUTPUT_PATH}" "--target-csv=${TARGET_CSV}" \
      "--scenario-id=${SLURM_ARRAY_TASK_ID}" --num-replicates=100000
    ;;
  summarize)
    Rscript --no-save --no-restore optim/summarize-policy-model-results.R \
      "--input-path=${OUTPUT_PATH}"
    ;;
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
