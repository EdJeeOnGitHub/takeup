#!/usr/bin/env bash
# Run the optimization smoketest locally.
# Overrides cluster-specific settings from archive/code/policy-v1/slurm_run_optim.sh.

set -euo pipefail

export SMOKETEST=1
export STAN_INPUT="data/stan_analysis_data"
export TMPDIR=/tmp/takeup-rcpp
export GUROBI_HOME=/opt/gurobi1000/linux64
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

mkdir -p "$TMPDIR"

# Stub out HPC module system
module() { echo "[local] skipping: module $*"; }
export -f module

bash archive/code/policy-v1/slurm_run_optim.sh
