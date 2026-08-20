# Repository layout and maintenance rules

The repository has a production spine and an explicit historical archive.
New code should follow these boundaries.

## Production code

- Put reusable R functions in `R/<component>/`. Files here should define
  functions and objects, not parse command-line arguments or launch analyses.
- Put runnable local entry points in `scripts/<component>/`.
- Put Slurm launchers and workers in `hpc/<component>/`; keep the computational
  implementation in `R/` or `scripts/` whenever practical.
- Put maintained Stan sources in `stan_models/` and document canonical models
  in `stan_models/README.md`.
- Write disposable and reviewable outputs under `build/`. Do not overwrite a
  frozen manuscript artifact as a side effect of an analysis target.

The supported scope is the main paper plus maintained appendix and
referee-response workflows. The top-level `Makefile` is the public interface;
`_targets.R` provides dependency tracking for the ordinary analysis pipeline.

## Scratch and archive

`scratch/` is ignored and is only for new local experiments. No supported
script, Make target, or paper dependency may source a file there.

Unsupported historical code belongs in `archive/code/`. Files are retained,
not deleted, and every archived file is recorded in
`archive/code/manifest.csv` with its former path and reason. Regenerate that
index with:

```sh
Rscript scripts/workflow/build-code-archive-manifest.R
```

Archived code may be invoked only through an explicitly legacy-named target;
it must not be a hidden dependency of the default paper build.

## Data and reproducibility

Private, licensed, or large external inputs are described in
`replication/data-manifest.csv`. Compact depositable inputs belong under
`replication/inputs/`. The replication staging script copies supported code,
declared inputs, Stan dependencies, and contracted paper artifacts into a
checksummed directory below `build/`.

Run `make check` after moving files. It checks the distance specification,
paper coverage, Stan inventory, replication dependencies, and these layout
invariants.
