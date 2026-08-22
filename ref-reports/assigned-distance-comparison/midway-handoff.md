# Midway handoff: assigned-distance policy and structural robustness

Date: 2026-08-22

## Objective

Produce the server-only assigned-distance artifacts needed to finish the
assigned-versus-realized manuscript comparison. Do not overwrite or promote
the current paper outputs. Keep all new results in isolated output directories
and return them with provenance.

The local container is handling the baseline structural generated quantities,
post-processing, tables, and figures for both assigned and realized Close/Far.
Midway is needed for the full optimal-policy workflow and any missing
structural robustness fits.

## Code and safety checks

Use the `ed-refine-todos` branch at commit
`cd5f2955ab5d0b31408e3e66b68518a5d54132c9` or a descendant containing that
commit. Before doing anything, report:

```bash
cd /home/edjee/projects/takeup-ed-refine-todos
git branch --show-current
git rev-parse HEAD
git status --short
squeue -u "$USER" -o '%.18i %.30j %.10T %.10M %.20R'
```

Do not reset, clean, or overwrite a dirty checkout. If the checkout is dirty or
the required commit is unavailable, stop and report the state.

All jobs below must use:

```bash
export DISTANCE_DEFINITION=assigned
export PROJECT_ROOT=/home/edjee/projects/takeup-ed-refine-todos
cd "$PROJECT_ROOT"
```

## 1. Full optimal-policy workflow

First check whether the complete assigned exponential cluster-weight policy
outputs already exist and are valid:

```bash
POLICY_ROOT=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots-exponential-cluster-weights
find "$POLICY_ROOT" -maxdepth 3 -type f -printf '%TY-%Tm-%Td %TH:%TM %s %p\n' 2>/dev/null | sort
for scenario in control bracelet static-control static-bracelet suppress-reputation; do
  test -f "$POLICY_ROOT/allocations/$scenario/status.csv" && \
    wc -l "$POLICY_ROOT/allocations/$scenario/status.csv"
done
```

If the output is absent or incomplete, submit the checkpointed 999-refit
workflow:

```bash
DISTANCE_DEFINITION=assigned NUM_REPLICATES=999 \
  bash hpc/policy/submit_policy_cluster_bootstrap.sh
```

Record the printed job IDs. Monitor the dependency chain through the final
`population` job. Existing per-replicate allocation files are checkpoints and
should be reused; do not delete them before rerunning.

Expected core outputs include:

```text
optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots-exponential-cluster-weights/policy-bootstrap-parameters.csv
optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots-exponential-cluster-weights/policy-cluster-bootstrap-replicates.csv
optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots-exponential-cluster-weights/policy-cluster-bootstrap-summary.csv
presentations/tables/fit105/optim-summ-exponential-cluster-weights.tex
ref-reports/policy-cost-sensitivity/
```

The comparison report also needs the paper-facing optimal-policy table and
figures. Determine whether the current sparse workflow renders these legacy
paper artifact names:

```text
tables/fit105/optim-summ-table.tex
optim-figures/fit105/plot-scaled-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-agg-identity-full-many-pots-pred-demand-vstar-comp-all.pdf
optim-figures/fit105/panel-scenarios-compare-optimal-allocation-plot-distconstraint-3500.pdf
```

If it does not, do not silently copy old realized outputs under new names.
Report which modern summaries were produced and which renderer/compatibility
step is still missing.

## 2. Policy model robustness

Check whether the six assigned model-robustness summaries already exist under:

```text
optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-model-robustness/
```

Required model IDs are `correct-observability`,
`second-order-observability`, `grouped-lambda`, `arm-lambda`, `student-t5`, and
`cluster-shock`. If any are incomplete, submit:

```bash
DISTANCE_DEFINITION=assigned REPO_ROOT="$PROJECT_ROOT" \
  bash hpc/policy/submit_policy_model_robustness.sh
```

Expected result from each model directory:

```text
policy-model-summary.csv
```

## 3. Missing structural robustness artifacts

The local comparison still lacks newly generated assigned versions of the
paper's structural robustness tables L1--L3 and L5. Before submitting new fits,
inventory the expected artifacts and their provenance:

```bash
for path in \
  presentations/tables/fit105/indiv-dist-community-fp-indiv-vis-robust-struct-overall-te-table.tex \
  presentations/tables/fit105/indiv-dist-indiv-fp-robust-struct-overall-te-table.tex \
  presentations/tables/fit105/struct-robustness-nooutliers-overall-te-table.tex \
  presentations/tables/fit105/struct-social-multiplier-robustness-table.tex; do
  if test -f "$path"; then sha256sum "$path"; stat "$path"; else echo "MISSING $path"; fi
done
```

If these are not already assigned-distance outputs, use the corresponding
commands in `build/candidate-hpc/assigned/export/job-manifest.csv`. In
particular, the requested server workflows are `individual-distance`,
`no-outliers`, and the model inputs used to assemble the social-multiplier
robustness table. Do not rerun completed fits merely to change table formatting;
postprocess existing valid chain CSVs where possible.

## 4. Return package

Do not commit generated results or replace current paper files. Create a
timestamped return directory outside the Git checkout containing:

1. the optimal-policy parameter, replicate, and summary CSVs;
2. the optimal-policy `.tex` table and every newly rendered policy PDF;
3. the six policy-model-robustness summary CSVs;
4. assigned structural robustness tables L1--L3 and L5;
5. all submission/job IDs and final `sacct` statuses;
6. a SHA-256 manifest with source paths; and
7. a text file recording branch, commit, distance definition, generation UTC,
   and any missing artifacts.

Suggested layout:

```text
assigned-distance-midway-YYYYMMDD-HHMMSS/
  provenance.txt
  jobs.txt
  artifact-manifest.sha256
  policy/
  policy-model-robustness/
  structural-robustness/
```

Return the absolute path to this directory. Do not delete the original Midway
outputs after packaging them.
