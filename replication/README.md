# Take-up replication workflow

The public entrypoint is the top-level `Makefile`. The pipeline uses the R
`targets` package to track dependencies and avoid rerunning unchanged analyses.

## Standard reproduction

```sh
make setup
make paper DISTANCE_SPEC=assigned
make check DISTANCE_SPEC=assigned
```

`make paper` regenerates ordinary reduced-form and balance outputs, validates
all remaining approved outputs against
`replication/paper-artifact-contract.csv`, and stages a complete manuscript in
`build/paper/<spec>/`. It does not start structural sampling or production
policy optimization.

The canonical sampling model is `stan_models/takeup_struct_main_core.stan`,
its compact-GQ companion is
`stan_models/takeup_struct_main_core_compact_gq.stan`, and the default input
workspace is `data/stan_analysis_data/dist_fit104.RData`. These defaults are
shared by the Make pipeline and direct helper scripts.

To build and compare both Close/Far definitions:

```sh
make compare-distance
```

Balance is split into independently cached sections. For a quick iteration,
run only the section being changed, for example:

```sh
make balance DISTANCE_SPEC=assigned BALANCE_SECTIONS=main
make balance DISTANCE_SPEC=assigned BALANCE_SECTIONS=orig
make balance DISTANCE_SPEC=assigned BALANCE_SECTIONS=fit-ri RI_DRAWS=50
```

Valid sections are `main`, `orig`, `fit-ri`, `attrition`,
`monitored-attrition`, and `sms`. The default is all sections and 500 RI draws.
Independent balance sections and reduced-form bootstrap draws use up to
`TAKEUP_THREADS` workers (eight by default in the project container). Production
defaults remain 500 draws; use `BOOTSTRAP_DRAWS=50 RI_DRAWS=50` only for a
quick development smoke test.

Eligible linear reduced-form specifications prebuild their observed and four
treatment-counterfactual design matrices, then refit each Bayesian-bootstrap
draw with weighted cross-products. Before using this path, every specification
checks its first draw against the original Fixest fit and stops if collapsed
predictions differ by more than `1e-9`. IV, Lee-bound, and other unsupported
specifications retain the Fixest implementation. Set `TAKEUP_FAST_WLS=0` to
force the complete legacy Fixest path for comparison. The reproducible timing
and numerical-equivalence benchmark is `scratch/bench-rf-direct-wls.R`.

## Paper-output coverage

The default workflow is deliberately freeze-first. Each active dependency is
`generated` in the isolated build, a version-controlled `static` source asset,
or a checksum-validated `frozen` approved output in
`replication/paper-artifacts/`. New structural or policy results stay in their
build directories until they have been compared and deliberately promoted by
updating the frozen artifact and checksum.

Run `make paper-audit` to write
`build/manifest/paper-pipeline-coverage.csv`. It reports complete default-build
coverage separately from full-regeneration coverage. The latter continues to
identify the focused structural-render and isolated legacy-output work that
remains for a public replication package.

`make paper-full POLICY_SMOKETEST=0` refreshes the expensive compact structural
GQ and production policy workflows before staging. It still does not promote
new results over approved paper outputs automatically. `make design-paper-full`
runs the historical design simulation and PAP map source; both require the
restricted design inputs.

Production policy optimization is not silently run by `make paper`: it needs
Gurobi, posterior inputs, and substantially more compute. Its current entry
points are `make optimal-policy` for estimation/optimization/figures and
`make policy-tables` for validating synced outputs and rendering tables.
`make optimal-policy` is a five-draw smoke test by default; pass
`POLICY_SMOKETEST=0` only in a production Gurobi environment.

## Full structural refit

```sh
make structural-fit DISTANCE_SPEC=assigned
```

This uses the stripped main-core Stan model, runs four chains across at most
eight CPU threads, and writes outputs to `build/structural-fit/assigned/`.
It is intentionally separate from ordinary paper reproduction.

## Outputs and data

All generated outputs live below `build/`; existing results are not overwritten.
See `replication/data-manifest.csv` for the intended disposition of inputs.
Restricted inputs must be placed at their documented paths before running the
pipeline. Disposable `_targets/`, `build/`, `temp-data/`, logs, and compiled Stan
binaries are not replication inputs.

`make stan-inventory` writes a non-destructive inventory to
`build/manifest/stan-artifacts.csv`. It does not move or delete archived fits.
