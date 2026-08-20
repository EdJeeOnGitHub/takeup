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
and numerical-equivalence benchmark is `tests/benchmarks/rf-direct-wls.R`.

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

`make structural-render` runs the focused fit-105 table/figure renderer into
`build/structural-paper/fit105/` and writes `comparison.csv` against the frozen
paper artifacts. It reads the small postprocessed RDS summaries in
`temp-data/struct-postprocess`; it does not load the 12GB social-multiplier draw
object. `make structural-postprocess` remains the separate compact-GQ step from
the four saved Stan chains. New results are never promoted over approved paper
outputs automatically.

`make paper-full` refreshes compact structural GQ, focused structural renders,
and the fast sparse policy workflow before staging. The policy default is a
five-refit smoke test; use `POLICY_REPLICATES=999` for the complete exponential
cluster-weighted exercise. `make design-paper-full` runs the historical design
simulation and PAP map source; both require the restricted design inputs.

Production policy optimization is not silently run by `make paper`. The
current `make optimal-policy` entrypoint is the low-memory cluster-weighted
workflow: it reads the compact `policy-bootstrap-parameters.csv`, predicts only
the 1,252 feasible village--PoT edges, optimizes five scenarios, and writes
everything below `build/policy/cluster-weighted/`. It needs either `glpsol`
(Ubuntu package `glpk-utils`) or `gurobi_cl`. Five refits are used by default;
set `POLICY_REPLICATES=999` for production. The older fit-105 dense posterior
workflow is retained as `make optimal-policy-legacy` and still supports
`POLICY_SMOKETEST=1`.

The policy inputs have three layers. The external cluster-weighted mode fits
are needed only to rebuild the compact parameter CSV with
`optim/prepare-policy-cluster-bootstrap.R`. Thereafter the local workflow needs
that CSV, `optim/data/full-many-pots-experiment.rds`, and the fixed experimental
welfare-target CSV under `optim/data/.../agg-full-many-pots/`. The design
workflow is not downstream of the structural model: it separately needs the
restricted randomization and geospatial RDS inputs documented in
`replication/data-manifest.csv`.

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
