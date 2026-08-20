#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=indiv-fp-prod-gq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=0-12:00:00
#SBATCH --output=temp/log/indiv-fp-prod-gq-%j.log
#SBATCH --error=temp/log/indiv-fp-prod-gq-%j.log

set -euo pipefail

CHAINS=${CHAINS:-2,3}
INPUT_PATH=${INPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data}
OUTPUT_PATH=${OUTPUT_PATH:-/project/akaring/takeup-data/temp-data/indiv-fp-reduce-sum-production-gq}
FIT_CSVS=${FIT_CSVS:-${INPUT_PATH}/indiv-fp-reduce-sum-production/dist_fit106_INDIV_FP_REDUCE_SUM_PROD_chain2-1.csv,${INPUT_PATH}/indiv-fp-reduce-sum-production/dist_fit106_INDIV_FP_REDUCE_SUM_PROD_chain3-1.csv}

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}

Rscript --no-save --no-restore --no-init-file \
  scripts/appendix/generate-compact-individual-sm-gq.R \
  --model=indiv_dist_indiv_fp \
  "--input-path=${INPUT_PATH}" \
  "--output-path=${OUTPUT_PATH}" \
  --stan-path=stan_models_fit105 \
  "--chains=${CHAINS}" \
  --parallel-chains=2 \
  --threads-per-chain=8 \
  "--fit-csvs=${FIT_CSVS}"
