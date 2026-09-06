#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --job-name=policy-speedup-verify
#SBATCH --output=/project/akaring/takeup-data/scratch/policy-speedup-20260906/verify-%j.log
set -euo pipefail
module load R/4.2.0
export GUROBI_HOME="$HOME/gurobi952/linux64"
export PATH="$GUROBI_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GUROBI_HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export GRB_LICENSE_FILE=/software/gurobi-9.2-el7-x86_64/gurobi.lic
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /project/akaring/takeup-data/scratch/policy-speedup-20260906/candidate-final
Rscript --vanilla scripts/checks/verify-final-policy.R
Rscript --vanilla scripts/checks/test-policy-speedup-execution.R --output-path=../execution-contract-final
