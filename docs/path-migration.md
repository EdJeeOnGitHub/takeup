# Repository path migration

The Make target names and command-line options are unchanged. Direct script
callers must use the new canonical paths; no compatibility shims are kept.

| Old location | New location |
|---|---|
| `R/distance-spec.R` | `R/distance/spec.R` |
| `R/pipeline.R` | `R/workflow/pipeline.R` |
| `scratch/reduced-form-functions.R` | `R/reduced-form/functions.R` |
| `scratch/reduced-form-setup.R` | `R/reduced-form/context.R` |
| `scratch/reduced-form-bootstrap.R` | `scripts/reduced-form/bootstrap.R` |
| `scratch/main-core-data.R` | `R/structural/main-core-data.R` |
| `scratch/sample-slim-individual-fp.R` | `scripts/structural/sample-main-core.R` |
| `scratch/generate-main-core-compact-gq.R` | `scripts/structural/generate-compact-gq.R` |
| `structural_tables.R` | `scripts/structural/render-paper.R` |
| `balance.R` | `scripts/balance/run.R` |
| `optim/policy-bootstrap-functions.R` | `R/policy/bootstrap.R` |
| `optim/*policy-cluster-bootstrap.R` | `scripts/policy/*-cluster-bootstrap.R` |
| root `slurm_main_core*` and `submit_main_core*` | `hpc/structural/` |
| root policy Slurm and submit scripts | `hpc/policy/` |
| supported robustness launchers | `hpc/appendix/` |
| unsupported root and scratch code | `archive/code/` (see its manifest) |

Compact policy inputs moved from `optim/data/...cluster-weights/` to
`replication/inputs/policy/`.
