#!/usr/bin/env bash
#SBATCH --partition=caslake
#SBATCH --account=pi-akaring
#SBATCH --job-name=cluster-shock-gq-1250
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=temp/log/cluster-shock-gq-1250-%j.log
#SBATCH --error=temp/log/cluster-shock-gq-1250-%j.log

# Re-run only the compact generated quantities for the existing community
# random-shock posterior. No Stan sampling or structural refitting occurs.
set -euo pipefail

REPO_ROOT=${REPO_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}
FIT_ROOT=${FIT_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-cluster-shock-production}
WORKSPACE=${WORKSPACE:-${PROJECT_ROOT}/build/structural-workspace/main-core-input.RData}
OUTPUT_PATH=${OUTPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-cluster-shock-population-gq-1250}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.35.0}
THREADS=${THREADS:-1}
PARALLEL_CHAINS=${PARALLEL_CHAINS:-4}

cd "${REPO_ROOT}"
mkdir -p temp/log "${OUTPUT_PATH}"

fits=(
  "${FIT_ROOT}/dist_fit106_MAIN_CORE_SHOCK_SD0.1_chain1-1.csv"
  "${FIT_ROOT}/dist_fit106_MAIN_CORE_SHOCK_SD0.1_chain2-1.csv"
  "${FIT_ROOT}/dist_fit106_MAIN_CORE_SHOCK_SD0.1_chain3-1.csv"
  "${FIT_ROOT}/dist_fit106_MAIN_CORE_SHOCK_SD0.1_chain4-1.csv"
)
for fit in "${fits[@]}"; do
  if [[ ! -s "${fit}" ]]; then
    echo "Missing random-shock fit: ${fit}" >&2
    exit 2
  fi
done

module load gcc/10.2.0
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/R/x86_64-pc-linux-gnu-library/4.2}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

fit_csvs=$(IFS=,; echo "${fits[*]}")
Rscript --no-save --no-restore --no-init-file \
  scripts/structural/generate-compact-gq.R \
  "--workspace=${WORKSPACE}" \
  "--fit-csvs=${fit_csvs}" \
  "--output-path=${OUTPUT_PATH}" \
  --output-basename=main-core-cluster-shock-population-1250 \
  --use-core-cluster-shock=1 \
  --core-cluster-shock-sd-prior=0.1 \
  --distance-definition=assigned \
  --sm-evaluation-distance-m=1250 \
  "--parallel-chains=${PARALLEL_CHAINS}" \
  "--threads=${THREADS}"

gq_files=("${OUTPUT_PATH}"/main-core-cluster-shock-population-1250-*.csv)
if [[ ${#gq_files[@]} -ne 4 ]]; then
  echo "Expected four compact-GQ outputs; found ${#gq_files[@]}." >&2
  exit 2
fi
gq_csvs=$(IFS=,; echo "${gq_files[*]}")
Rscript --no-save --no-restore --no-init-file \
  scripts/appendix/summarize-main-core-cluster-shock-population.R \
  "--gq-files=${gq_csvs}" \
  --distance-index=1 \
  --distance-m=1250 \
  "--output-path=${OUTPUT_PATH}/summary"

echo "Completed exact-1250m cluster-shock multiplier summaries: ${OUTPUT_PATH}/summary"
