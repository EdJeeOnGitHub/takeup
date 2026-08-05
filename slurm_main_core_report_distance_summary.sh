#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=report-dist-summary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=temp/log/report-dist-summary-%j.log
#SBATCH --error=temp/log/report-dist-summary-%j.log

set -euo pipefail
ROOT=${ROOT:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-report-distance-priors}
mkdir -p temp/log "${ROOT}/summary"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTANR_NO_VER_CHECK=TRUE
Rscript --no-save --no-restore --no-init-file \
  scratch/summarize-main-core-report-distance-priors.R \
  "--root=${ROOT}" "--output-path=${ROOT}/summary"
