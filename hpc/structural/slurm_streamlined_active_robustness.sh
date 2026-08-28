#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=slim-active-robust
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=1-12:00:00
#SBATCH --array=1-16
#SBATCH --output=temp/log/slim-active-robust-%A-%a.log
#SBATCH --error=temp/log/slim-active-robust-%A-%a.log

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
PREPARED_ROOT=${PREPARED_ROOT:-/project/akaring/takeup-data/candidate-multiplier-1250-20260825/workspaces}
OUTPUT_ROOT=${OUTPUT_ROOT:-${ANALYSIS_ROOT}/streamlined-active-robustness}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.35.0}
ITER_WARMUP=${ITER_WARMUP:-1000}
ITER_SAMPLING=${ITER_SAMPLING:-1000}
ADAPT_DELTA=${ADAPT_DELTA:-0.999}
MAX_TREEDEPTH=${MAX_TREEDEPTH:-12}
task=${SLURM_ARRAY_TASK_ID:?Submit as an array from 1 to 16}
spec_index=$(( (task - 1) / 4 + 1 ))
chain_id=$(( (task - 1) % 4 + 1 ))

case ${spec_index} in
  1)
    spec_id=private-distance-community-image
    model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
    workspace=${ANALYSIS_ROOT}/dist_fit104.RData
    stan_file=takeup_struct_private_info_slim.stan
    preprocessed=0
    ;;
  2)
    spec_id=full-information
    model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP
    workspace=${ANALYSIS_ROOT}/dist_fit105.RData
    stan_file=takeup_struct_indiv_fp_slim.stan
    preprocessed=0
    ;;
  3)
    spec_id=exclude-dispersed
    model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS
    workspace=${PREPARED_ROOT}/dist_fit1250.RData
    stan_file=takeup_struct_indiv_fp_slim.stan
    preprocessed=1
    ;;
  4)
    spec_id=second-order-observability
    model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB
    workspace=${PREPARED_ROOT}/dist_fit1252.RData
    stan_file=takeup_struct_indiv_fp_slim.stan
    preprocessed=1
    ;;
  *) exit 2 ;;
esac

output_path=${OUTPUT_ROOT}/${spec_id}/fits
output_basename=${spec_id}-slim-chain${chain_id}

cd "${PROJECT_ROOT}"
mkdir -p temp/log "${output_path}"
module load gcc/10.2.0
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

Rscript --no-save --no-restore --no-init-file \
  scripts/structural/sample-main-core.R \
  --model="${model}" \
  --workspace="${workspace}" \
  --workspace-already-preprocessed="${preprocessed}" \
  --output-path="${output_path}" \
  --output-basename="${output_basename}" \
  --stan-path=stan_models \
  --stan-file="${stan_file}" \
  --cmdstan-path="${CMDSTAN_PATH}" \
  --chains=1 --parallel-chains=1 --threads-per-chain=4 \
  --iter-warmup="${ITER_WARMUP}" --iter-sampling="${ITER_SAMPLING}" \
  --adapt-delta="${ADAPT_DELTA}" --max-treedepth="${MAX_TREEDEPTH}" --metric=diag_e \
  --seed="$((20260828 + chain_id + 100 * spec_index))" \
  --distance-definition=assigned \
  --refresh=25

cat > "${output_path}/${output_basename}-provenance.txt" <<EOF
spec_id=${spec_id}
model=${model}
stan_file=${stan_file}
workspace=${workspace}
workspace_already_preprocessed=${preprocessed}
distance_definition=assigned
git_commit=$(git rev-parse HEAD)
slurm_job_id=${SLURM_JOB_ID:-}
slurm_array_task_id=${task}
chain_id=${chain_id}
EOF
