# Assigned-distance paper update runbook

This document describes how to generate and review a complete paper update using the randomized assigned Close/Far definition without replacing the currently approved paper artifacts.

The existing command `make paper DISTANCE_SPEC=assigned` is not, by itself, a complete assigned-distance regeneration. It regenerates ordinary reduced-form and balance outputs but validates and copies several approved or frozen structural and appendix artifacts that may still contain results based on realized Close/Far status.

## Current status

| Component | Status | Required action |
|---|---|---|
| Distance crosswalk and audit | Complete | Validate again during the final candidate build |
| Main reduced form, including Tables 1 and 2 | Complete | The 500-draw assigned and realized bootstrap runs are available |
| Balance tables | Complete locally | Six validated sections and candidate tables use 500 RI draws |
| Reduced-form appendix results | Complete locally | Campaign-day and multiplier use 1,000 bootstrap draws; RI uses 99,999 reallocations |
| Baseline structural ATEs | Complete locally | Four-chain slim fit, compact GQ, and candidate tables are available |
| Other structural specifications | Pending | Reaggregate saved generated quantities using assigned groups |
| Cluster structural inference | Pending | Reprocess the 999 weighted fits and cluster-random-shock model |
| Continuous-distance multiplier and policy results | Expected to be invariant | Verify rather than re-estimate when possible |
| Manuscript prose | Pending | Replace hard-coded numbers only in a staged manuscript copy |
| Candidate PDF | Pending | Compile separately under `build/` |

## 1. Use the non-promoting candidate-paper target

The implemented entrypoint is:

```bash
make paper-candidate DISTANCE_SPEC=assigned
```

It should write only below:

```text
build/paper-candidate/assigned/
```

The target:

- redirect newly generated tables and figures into `build/`;
- copy the manuscript into the candidate directory;
- overlay all assigned-distance candidate artifacts;
- compile the copied manuscript;
- produce a manifest comparing candidate and approved files;
- never modify `presentations/`, the stable structural-robustness appendix, artifact checksums, Overleaf, or the source manuscript.

This isolation is necessary because some existing appendix scripts copy their results directly into `appendix/structural-robustness/` or other stable paper-facing locations.

## 2. Assigned balance results

The full assigned balance build completed successfully. To reproduce it, run:

```bash
make balance \
  DISTANCE_SPEC=assigned \
  BALANCE_SECTIONS=all \
  RI_DRAWS=500 \
  TAKEUP_THREADS=8

make balance-tables DISTANCE_SPEC=assigned
```

The cached sections can instead be run independently:

```bash
make balance DISTANCE_SPEC=assigned BALANCE_SECTIONS=main RI_DRAWS=500 TAKEUP_THREADS=8
make balance DISTANCE_SPEC=assigned BALANCE_SECTIONS=orig RI_DRAWS=500 TAKEUP_THREADS=8
make balance DISTANCE_SPEC=assigned BALANCE_SECTIONS=attrition RI_DRAWS=500 TAKEUP_THREADS=8
make balance DISTANCE_SPEC=assigned BALANCE_SECTIONS=monitored-attrition RI_DRAWS=500 TAKEUP_THREADS=8
make balance DISTANCE_SPEC=assigned BALANCE_SECTIONS=sms RI_DRAWS=500 TAKEUP_THREADS=8
make balance DISTANCE_SPEC=assigned BALANCE_SECTIONS=fit-ri RI_DRAWS=500 TAKEUP_THREADS=8
make balance-tables DISTANCE_SPEC=assigned
```

Render the tables only after all required sections have completed. Completed
sections are cached and validated independently, so a repeat invocation skips
them.

## 3. Reduced-form results

The full 500-draw assigned reduced-form run is already complete, including the main take-up and observability results. It need not be rerun unless its code or inputs change.

For a clean end-to-end invocation, use:

```bash
make reduced-form \
  DISTANCE_SPEC=assigned \
  BOOTSTRAP_DRAWS=500 \
  TAKEUP_THREADS=8
```

The targets pipeline should recognize and skip completed work.

The local candidate workflow now covers the auxiliary reduced-form outputs
outside the main target, including:

- campaign-day results;
- randomization inference;
- reduced-form distance-multiplier results;
- continuous-distance specifications whose output also displays binary Close/Far estimates;
- `sms-TE-by-dist-incentive.pdf`;
- `dist-ri-plot.pdf`.

The production local settings are 1,000 bootstrap draws for campaign-day and
the reduced-form distance multiplier, and 99,999 conditional arm reallocations
for randomization inference.

## 4. Baseline and alternative structural results

The baseline structural model does not need to be refitted solely because of
the binary-distance change. For the candidate audit, however, a fresh
four-chain fit of the latest slim main-core model was completed and retained
under `build/structural-fit/assigned/`.

Compact generated quantities and modern renderer inputs are available under:

```text
build/assigned/structural/compact-gq/
build/candidate-components/assigned/structural-data/
build/candidate-components/assigned/tables/
```

Regenerate and render them with:

```bash
make structural-postprocess DISTANCE_SPEC=assigned
make structural-candidate-render DISTANCE_SPEC=assigned
```

The fit wrapper reuses the completed fit only after checking the manifest and
diagnostics. Set `FORCE_STRUCTURAL_FIT=1` when a deliberate refit is required.

Apply the same assigned-group reaggregation to every robustness output with Combined, Close, Far, or Far-minus-Close columns, including:

- individual-distance/community fixed-point robustness;
- individual-distance/individual fixed-point robustness;
- no-outliers robustness;
- relevant lambda or nested-model ATE tables;
- continuous-distance robustness tables that also report binary Close/Far groups.

When their saved draws or compact generated quantities are available, these should be postprocessing jobs rather than structural refits.

## 5. Cluster structural robustness

The current assigned structural comparison does not yet include all cluster robustness. Assigned Close/Far treatment-effect aggregation is still required for:

- the 999 exponential cluster-weighted fits;
- the cluster-random-shock specification;
- their appendix structural treatment-effect tables.

If the external compute directory retains the fitted parameter vectors or compact generated-quantity inputs for all 999 runs, rerun only compact generated quantities and aggregation. Do not repeat 999 optimizations unless the required fitted outputs are unavailable or incompatible.

The relevant external workflow corresponds to:

```text
main-core-exponential-cluster-weight-999/
```

The current local summary CSVs are not sufficient to reconstruct every assigned Close/Far estimate. Before launching the external workflow, ensure that the distance definition is passed explicitly and consistently through:

```text
hpc/structural/slurm_main_core_cluster_bootstrap_999.sh
scripts/appendix/summarize-main-core-cluster-robustness.R
scripts/appendix/generate-main-core-cluster-bootstrap-ate-table.R
```

Cluster multiplier contrasts evaluated at fixed continuous distances should remain unchanged. Regenerate only outputs that group observations into binary Close and Far categories, while validating continuous-distance outputs for identity.

## 6. Policy and design outputs

The optimal-policy workflow uses continuous distances, so changing the binary Close/Far label should not require rerunning the optimization.

For the complete candidate update:

- compare the policy inputs used by the assigned and realized builds;
- verify that policy tables do not group results by binary Close/Far;
- retain existing policy outputs if their inputs and results are identical;
- record this identity check in the candidate manifest.

The assignment-density design figure describes the experimental assignment mechanism and should normally remain unchanged. It provides evidence for using design-based assigned Close/Far status.

## 7. Stage, compile, and audit the candidate

After all numerical outputs are ready, run:

```bash
make paper-candidate DISTANCE_SPEC=assigned
make check DISTANCE_SPEC=assigned
```

Audit the staged manuscript for Close/Far estimates, standard errors, confidence intervals, and p-values that are hard-coded in prose rather than inserted from generated tables. Initially update these only in the staged TeX copy.

The candidate bundle should contain at least:

```text
build/paper-candidate/assigned/ECM ReStud.pdf
build/paper-candidate/assigned/artifact-comparison.csv
build/paper-candidate/assigned/numeric-text-audit.csv
build/paper-candidate/assigned/build-manifest.csv
```

For every changed artifact, the comparison should record:

- the current approved file;
- the assigned-distance candidate file;
- whether its contents changed;
- old and new key estimates;
- whether the result required estimation, postprocessing, or validation only.

## Recommended execution order

1. Export and run the Midway bundle for the remaining structural robustness jobs.
2. Reprocess the 999 cluster-weighted fits and cluster-random-shock output.
3. Generate and validate policy outputs from the canonical compact parameters.
4. Import the checksummed Midway result bundle.
5. Stage and compile the candidate paper.
6. Audit hard-coded manuscript numbers and create the old-versus-new comparison bundle.
7. Review the bundle before separately authorizing any promotion into the paper or Overleaf.

The critical path is therefore the isolated candidate target, full balance run, and complete structural Close/Far reaggregation—especially cluster inference. No approved result should be replaced until the candidate PDF and comparison manifest have been reviewed.
