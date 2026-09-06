# Policy optimizer speedups: implementation and validation

Implemented 2026-09-06. The production population-weighting rerun has not been launched. This implements the execution-speed portion of [the plan](policy-optimizer-speedup-plan-2026-09-06.md); the adult-population/target corrections remain a separate task.

## Changes

- `predict-model-robustness.R` now writes compact edge-demand matrices for ordinary models as well as household-distance and cluster-shock models. Ordinary models still evaluate unique distances with the same equations, then map to ordered edges. Prediction diagnostics and a checksummed cache manifest are retained.
- `cost-sensitivity.R` caches invariant non-pooling assignment/open-site constraint text. The generated LP preserves the old constraint order and numerical formatting. Other callers retain the uncached interface, including pooling. Solver diagnostics now include termination, optimality, objective, bound, and gap; optional persistent logs survive temporary-file cleanup.
- `optimize-cluster-bootstrap.R` selects a scenario's matrix block once, caches village/edge grouping, processes bounded batches with `mclapply`, and collects status rows in lists. It supports explicit worker, batch, solver-thread, solver-seed, cache-format, and scratch settings.
- Each draw's assignment is published by atomic rename. The parent writes status checkpoints. Failed workers cause failed records and a nonzero job exit. A solver exit alone no longer counts as proof of optimality. Resume checks bind results to input, target, code, and solver settings; ambiguous or changed caches fail explicitly.
- Solver scratch files use the chosen temporary path; final assignments, diagnostic logs, and manifests remain in the output tree. Execution timing distinguishes cache/loading, lookup, model writing, solver execution, and output writing. Per-stage times summed over workers are work totals, not elapsed wall time.
- Both Slurm stage workers expose the execution controls. The model launcher no longer forces optimization into 4 GB; its default is 28 GB, overridable through `OPTIMIZE_MEMORY`. This is conservative headroom, not a claim that every run needs 28 GB.

## Default behavior and the thread-count finding

**Parallel draws are now the default, following user acceptance of the thread-count tie-breaking difference.** The optimizer uses up to eight workers, capped by allocated CPUs, with `Threads=1, Seed=0`. Choose serial execution with `--num-cores=1` or `OPTIMIZE_CORES=1` in the Slurm wrappers; this keeps the same solver settings as parallel mode. To reproduce the earlier automatic-thread reference specifically, also use `--solver-threads=auto --solver-seed=auto` (or `SOLVER_THREADS=auto SOLVER_SEED=auto`).

During the original validation, the default was serial with automatic solver threading. With those settings, the final code matched all 40 old-code test allocations exactly: two draws × five scenarios × four models (benchmark, cluster-shock, full-information, tight-multinomial). Predictions, scenario classifications, assignments, coverage, distance, and accounting inputs were checked, not just site counts. Optimal cases were certified from Gurobi logs/diagnostics; tight-multinomial cases included historical target infeasibility.

Parallel execution uses explicit solver thread and seed settings by default. At fixed `Threads=1, Seed=0`, the final parallel code matched all 40 single-thread serial reference allocations. The 32-draw performance pilot also matched serial versus four/eight workers.

Changing the solver thread count itself is **not allocation-equivalent in every case**. Benchmark draw 1, static-Control, changed six community assignments when switching from automatic Gurobi threading to one thread. Both solutions have a certified optimum of 99 sites and meet the same target of about 47.12069, but:

| Quantity | Old thread setting | One solver thread |
|---|---:|---:|
| Achieved welfare (legacy unit weights) | 47.20511 | 47.15554 |
| Mean community distance, metres | 1,461.214 | 1,468.464 |

The comparison correctly fails this case. No secondary objective or tolerance relaxation was introduced to hide the difference. The user accepted this tie-breaking consequence and authorized parallel execution as the default. Cached predictions and LP structure already provide a speedup without changing the solver-thread setting.

## Measured performance

Midway2, `broadwl`, eight allocated CPUs; Gurobi 9.5.2, R 4.2.0. Each row averages two fresh optimization runs on the same 32 benchmark Control draws, 3.5 km, existing unit-weight/common-target inputs. Prediction is performed before each variant's optimization measurements. These are end-to-end optimizer process times, including loading and output; queue time and prediction are excluded.

| Configuration | Mean seconds | Speed relative to old code |
|---|---:|---:|
| Frozen old code, automatic solver threads | 16.910 | 1.00× |
| New serial, automatic solver threads | 8.579 | 1.97× |
| New serial, one solver thread | 8.108 | 2.09× |
| Four draw workers, one solver thread each | 2.892 | 5.85× |
| Eight draw workers, one solver thread each | 2.034 | 8.32× |

The four/eight-worker pilots used batches of four draws; the configurable default batch size remains 25. With a tiny test set, choose a small batch size to actually exercise multiple workers. The default 25 is intended for longer draw lists and should be tuned with representative workloads.

GNU `time` recorded maximum resident size around 106–168 MiB for these small optimizer processes; this is not the sum of all forked workers' memory. Slurm reported approximately 1.24 GiB maximum RSS for the final multi-model validation batch. Retain conservative production memory requests pending a larger-draw/larger-cap measurement. These small pilots do not justify extrapolating an eightfold improvement to all models or the 10 km exercise. Scratch locality was not separately benchmarked, so no speedup is attributed specifically to local disk.

## Additional verification

- Byte-identical cached versus uncached LP text on bounded fixtures.
- After enabling parallel defaults, a 200-draw GLPK fixture comparison passed with no explicit worker/thread settings; the run respected a two-CPU allocation. The explicit one-worker serial override was also checked. These are synthetic execution checks, not additional empirical draws.
- Real interruption: sent SIGTERM via `timeout` after 84 of 200 synthetic draws; restart recovered all 200 assignments identically to an uninterrupted run.
- Synthetic infeasible and undefined-equilibrium cases retain their classifications.
- Changed settings, corrupted saved contracts, ambiguous caches, and worker errors produce failing runs rather than silently accepted results.
- R parsing and shell syntax checks pass for the changed implementation and launch scripts. Scoped whitespace checks pass. Unrelated pre-existing whitespace in `ref-reports/edited-feedback-refine.md` was left alone.

Corrected adult-weighted serial/parallel checks and the corrected 10 km pilot depend on the population/target implementation and have not been represented as completed here.

## Run log and artifacts

All server tests used an isolated directory:

`/project/akaring/takeup-data/scratch/policy-speedup-20260906`

Final tested code is in its `candidate-final/` directory. The regular Midway2 production checkout was not overwritten. Local implementation hashes are saved with the evidence.

| Job | Purpose | Result |
|---|---|---|
| 48979796 | Initial four-model comparison ladder | Failed only at the intentionally strict thread-change comparison; old/new same-settings and serial/parallel comparisons passed |
| 48979797 | Two-repeat 32-draw performance pilot and repeated comparison ladder | Performance comparisons passed; job failed at the same thread-change gate |
| 48979799 | Final-code default and parallel assignment verification; execution-contract checks | Completed, 1m53s |

Literal submitted batch scripts, the final verification driver, frozen pre-change source, CSV reports, cache/run manifests, and solver logs are in `temp-data/policy-speedup-validation-20260906/`. `server-evidence.tgz` preserves the detailed reports and logs. Full assignment and prediction RDS files remain in the isolated server directory. The frozen helper used by server comparisons was instrumented only to retain solver logs and set `OutputFlag=1`; the archive preserves the original unmodified source and hashes.

Review entry points:

- `final-verification.csv`: all 80 final-code assignment comparisons passed (40 then-default serial versus old; 40 parallel versus fixed-thread serial).
- `comparisons.csv`: the full settings-change ladder, including the deliberately visible thread-change failure.
- `performance.csv`: both timing repetitions and process RSS measurements.
- `interruption-comparison.csv`: the 200-draw interruption/restart comparison.

Submission commands executed from the local workspace were:

```bash
ssh -F /home/ed/.ssh/config midway 'sbatch --parsable' < /tmp/policy-speedup-pilot.sh
ssh -F /home/ed/.ssh/config midway 'sbatch --parsable' < /tmp/policy-speedup-final.sh
ssh -F /home/ed/.ssh/config midway 'sbatch --parsable' < /tmp/policy-speedup-verify.sh
```

The exact scripts are archived under the evidence directory with the same basenames. For future bounded experiments, explicitly select parallel settings, for example `OPTIMIZE_CORES=8 SOLVER_THREADS=1 SOLVER_SEED=0 DRAW_BATCH_SIZE=4`, within an eight-CPU allocation. The audit's production submission commands still require their separate population/target/benchmark-wrapper corrections.
