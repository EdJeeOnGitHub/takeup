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
TARGET_CSV=${TARGET_CSV:-/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv}
DISTANCE_DATA=${DISTANCE_DATA:-/project/akaring/takeup-data/optim/data/full-many-pots-experiment.rds}
REPO_ROOT=${REPO_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}

case "${MODEL_ID}" in
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
    FITS=("${ROOT}"/dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS-{1,2,3,4}.csv)
    EXTRACT_OPTIONS=()
    ;;
  full-information)
    MODEL_LABEL="Individual distance observed by peers"
    MODEL_FAMILY=full_information
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP-{1,2,3,4}.csv)
    EXTRACT_OPTIONS=()
    ;;
  exclude-dispersed)
    MODEL_LABEL="Excluding geographically dispersed communities"
    MODEL_FAMILY=gaussian
    LAMBDA_STRUCTURE=common
    FITS=("${ROOT}"/dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS-{1,2,3,4}.csv)
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
module load -f gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export GUROBI_HOME="${HOME}/gurobi952/linux64"
export PATH="${GUROBI_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export GRB_LICENSE_FILE=/software/gurobi-9.2-el7-x86_64/gurobi.lic
mkdir -p temp/log "${OUTPUT_PATH}"

case "${STAGE}" in
  prepare)
    python3 scripts/policy/extract-cmdstan-draws.py \
      --output "${COMPACT_CSV}" \
      ${EXTRACT_OPTIONS[@]+"${EXTRACT_OPTIONS[@]}"} "${FITS[@]}"
    Rscript --no-save --no-restore scripts/policy/prepare-model-robustness.R \
      "--input-csv=${COMPACT_CSV}" "--output-path=${OUTPUT_PATH}" \
      "--model-id=${MODEL_ID}" "--model-label=${MODEL_LABEL}" \
      "--model-family=${MODEL_FAMILY}" "--lambda-structure=${LAMBDA_STRUCTURE}"
    ;;
  predict)
    PREDICT_OPTIONS=()
    if [[ "${MODEL_ID}" == "private-distance-community-image" || "${MODEL_ID}" == "full-information" ]]; then
      PREDICT_OPTIONS+=(--household-workspace="${ROOT}/dist_fit104.RData")
    fi
    Rscript --no-save --no-restore scripts/policy/predict-model-robustness.R \
      "--parameter-rds=${OUTPUT_PATH}/policy-model-parameters.rds" \
      "--distance-data=${DISTANCE_DATA}" "--output-path=${OUTPUT_PATH}" \
      --distance-cap=3500 \
      "--num-cores=${NUM_CORES}" "${PREDICT_OPTIONS[@]}"
    ;;
  optimize)
    : "${SLURM_ARRAY_TASK_ID:?Optimize requires scenario array 1-5}"
    Rscript --no-save --no-restore scripts/policy/optimize-cluster-bootstrap.R \
      "--input-path=${OUTPUT_PATH}" "--target-csv=${TARGET_CSV}" \
      "--scenario-id=${SLURM_ARRAY_TASK_ID}" --num-replicates=100000
    ;;
  summarize)
    Rscript --no-save --no-restore scripts/policy/summarize-model-results.R \
      "--input-path=${OUTPUT_PATH}"
    {
      echo "model_id=${MODEL_ID}"
      echo "model_label=${MODEL_LABEL}"
      echo "distance_definition=assigned"
      echo "git_commit=$(git rev-parse HEAD)"
      echo "source_fit_directory=$(dirname "${FITS[0]:-prepared-balanced-assigned-distance-slim-chains}")"
      echo "structural_refit_performed=no"
      echo "candidate_sites=1451"
      echo "distance_cap_m=3500"
      echo "generated_utc=$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
    } > "${OUTPUT_PATH}/provenance.txt"
    ;;
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
