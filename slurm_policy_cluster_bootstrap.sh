#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=policy-cluster-bs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=28G
#SBATCH --time=04:00:00
#SBATCH --output=temp/log/policy-cluster-bs-%A-%a.log
#SBATCH --error=temp/log/policy-cluster-bs-%A-%a.log

set -euo pipefail

STAGE=${STAGE:?Set STAGE to prepare, predict, optimize, or summarize}
NUM_REPLICATES=${NUM_REPLICATES:-210}
NUM_CORES=${NUM_CORES:-12}
MODEL=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
OUTPUT_PATH=${OUTPUT_PATH:-optim/data/${MODEL}/agg-full-many-pots-cluster-bootstrap}
WEIGHTED_PATH=${WEIGHTED_PATH:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-weighted-modes}
DISTANCE_DATA=${DISTANCE_DATA:-optim/data/full-many-pots-experiment.rds}
TARGET_CSV=${TARGET_CSV:-optim/data/${MODEL}/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv}
TABLE_PATH=${TABLE_PATH:-presentations/tables/fit105/optim-summ-cluster-bootstrap.tex}

cd ~/projects/takeup
module load -f gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0
export GUROBI_HOME="${HOME}/gurobi952/linux64"
export PATH="${GUROBI_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export GRB_LICENSE_FILE=/software/gurobi-9.2-el7-x86_64/gurobi.lic
mkdir -p temp/log "${OUTPUT_PATH}"

case "${STAGE}" in
  prepare)
    Rscript --no-save --no-restore optim/prepare-policy-cluster-bootstrap.R \
      "--weighted-path=${WEIGHTED_PATH}" "--output-path=${OUTPUT_PATH}" \
      "--num-replicates=${NUM_REPLICATES}"
    ;;
  predict)
    Rscript --no-save --no-restore optim/predict-policy-cluster-bootstrap.R \
      "--parameter-csv=${OUTPUT_PATH}/policy-bootstrap-parameters.csv" \
      "--distance-data=${DISTANCE_DATA}" "--output-path=${OUTPUT_PATH}" \
      --distance-cap=3500 "--num-cores=${NUM_CORES}" \
      "--num-replicates=${NUM_REPLICATES}"
    ;;
  optimize)
    : "${SLURM_ARRAY_TASK_ID:?Optimize requires scenario array 1-5}"
    Rscript --no-save --no-restore optim/optimize-policy-cluster-bootstrap.R \
      "--input-path=${OUTPUT_PATH}" "--target-csv=${TARGET_CSV}" \
      "--scenario-id=${SLURM_ARRAY_TASK_ID}" \
      "--num-replicates=${NUM_REPLICATES}"
    ;;
  summarize)
    Rscript --no-save --no-restore optim/summarize-policy-cluster-bootstrap.R \
      "--input-path=${OUTPUT_PATH}" "--table-path=${TABLE_PATH}" \
      "--num-replicates=${NUM_REPLICATES}"
    ;;
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
