# TakeUp: analysis and replication code

This repository contains the R, Stan, and workflow code for the TakeUp study.
The active paper combines a randomized intervention varying social signals and
travel distance with a structural model of deworming take-up.

The top-level `Makefile` is the supported entry point. The primary paper
specification uses randomized assigned Close/Far status; the Close/Far group
reconstructed from realized travel distance is retained as a robustness
specification. Generated files are written below `build/`, so the two
definitions can be regenerated and compared without replacing approved paper
artifacts.

## Quick start

```sh
make setup
make paper DISTANCE_SPEC=assigned
make check DISTANCE_SPEC=assigned
make replication-package DISTANCE_SPEC=assigned
```

Use `make help` for the individual reduced-form, balance, structural, policy,
and audit targets. See `replication/README.md` for the complete reproduction
contract and restricted-input requirements.

## Repository layout

- `R/`: reusable analysis functions, with no command-line entry points.
- `scripts/`: supported local entry points, organized by analysis component.
- `hpc/`: Slurm launchers and workers for expensive structural, policy, and
  appendix workflows.
- `stan_models/`: current Stan models and model-level robustness variants.
- `replication/`: the data manifest, artifact contract, and compact depositable
  inputs.
- `appendix/` and `ref-reports/`: supported appendix and referee-response
  analyses and their stable paper-facing artifacts.
- `archive/code/`: unsupported historical code retained for provenance and
  indexed by `archive/code/manifest.csv`.
- `scratch/`: ignored local experimentation only; production code must not
  depend on it.
- `data/`: a private data submodule. The survey data are not currently public.
- `multilvlr/`: a public utility-code submodule.

The organization rules and the old-to-new path map are documented in
`docs/repository-layout.md` and `docs/path-migration.md`.

## Supported workflows

The ordinary paper build regenerates reduced-form and balance results, then
checksum-validates approved structural, policy, design, and other frozen
artifacts before staging the manuscript under `build/paper/<spec>/`.

The current structural workflow uses:

- `stan_models/takeup_struct_main_core.stan` for sampling;
- `stan_models/takeup_struct_main_core_compact_gq.stan` for compact generated
  quantities;
- `scripts/workflow/run-structural-fit.sh` for a local four-chain refit;
- `scripts/structural/render-paper.R` for focused paper tables and figures.

The current optimal-policy workflow is `make optimal-policy`. It uses compact
cluster-weighted parameters and sparse feasible-edge predictions, and requires
either GLPK (`glpsol`) or Gurobi. The older dense workflow is retained only as
the explicitly named `make optimal-policy-legacy` target under
`archive/code/policy-v1/`.

Restricted inputs and large saved fits are not committed. Their expected paths,
owners, and replication disposition are listed in
`replication/data-manifest.csv`. On a compute server, keep those inputs and
generated results outside the Git checkout and pass the relevant input/output
flags to the supported scripts.

## Output policy

New results are never promoted over approved paper outputs automatically.
Ordinary generated outputs live in `build/`; stable frozen artifacts are
declared in `replication/paper-artifact-contract.csv`. Run `make paper-audit`
to report default-build and full-regeneration coverage, and `make stan-inventory`
to inventory Stan sources and saved fits without moving or deleting them.
