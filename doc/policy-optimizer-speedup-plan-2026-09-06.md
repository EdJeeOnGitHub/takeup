# Policy optimizer: low-risk speedup implementation plan

Date: 2026-09-06. Implementation and bounded validation completed; see [results and the thread-count compatibility finding](policy-optimizer-speedup-results-2026-09-06.md). The production rerun and population/target corrections remain separate.

## Objective and scope

Reduce elapsed time and memory for the policy rerun while preserving predictions, optimization problems, and reported results. Implement the changes below in small, independently checked steps. Start from the newer Midway2 checkout (`~/projects/takeup-ed-refine-todos`); its core policy files matched the local checkout during the server inventory.

This plan complements [the population-weighting audit](policy-population-weighting-rerun-audit-2026-09-06.md). Speedup equivalence must hold weights and targets fixed. The old incorrect-population/fixed-target results are not the expected answers for the corrected economic exercise. First verify execution changes against frozen old semantics, then verify serial versus parallel execution under the corrected adult-weighted, draw-specific target contract.

Defer an in-process Gurobi API rewrite, warm starts, geographic-floor shortcuts, C++ prediction kernels, and changes to solver tolerances or tie-breaking objectives. These require additional validation and are not necessary for this first pass.

## 1. Freeze reference inputs and instrument the serial baseline

- Preserve the exact pre-change scripts/helpers, commit and working-file hashes. Use fresh, separate reference/candidate output roots; do not overwrite or resume historical production outputs.
- Record selected model, chain/iteration/draw/replicate IDs, parameter and geometry hashes, cap, population vector, exact targets, prediction values, solver executable/version, parameters, R version, and relevant environment settings.
- Use the same solver version for comparisons. Explicitly select Gurobi rather than `auto`, which currently chooses GLPK if installed. First reproduce the old settings; assess the single-thread setting as a separate change.
- Measure total process wall time, prediction/cache loading, per-draw demand lookup, model writing, solver execution, and output writing. Existing `elapsed_seconds` excludes demand lookup and other overhead; do not use it alone to claim a speedup. Record peak RSS through Slurm accounting where available.
- Retain solver logs and parse termination status, incumbent objective, bound, and gap for the comparison. The current optimizer's hard-coded successful status and the CLI process exit code do not establish optimality.

## 2. Store compact predictions for every model

Files: `scripts/policy/predict-model-robustness.R`, `scripts/policy/optimize-cluster-bootstrap.R`.

- Extend the existing `policy-edge-demand-matrix.rds` and draw-map interface to the models currently producing long demand tables. Keep existing household and cluster-shock prediction adapters.
- For ordinary models, keep evaluating unique distances as today, then map predictions to the ordered feasible edges. Produce the numeric row inside each prediction worker rather than constructing a full long table for all draws and converting afterward.
- Preserve the established layout: rows follow the draw map; columns are five scenario blocks, each in exactly the saved feasible-edge order. Validate scenario IDs, duplicate distances, draw IDs, and edge alignment explicitly.
- Preserve undefined-equilibrium masks and fallback diagnostics in a small companion manifest. Leave experimental community-level predictions and their model-specific geography/mapping intact.
- Load the matrix once per optimizer process. Select the scenario block once and precompute draw-row indices; remove repeated whole-table filtering and repeated edge-distance matching.
- Retain reading of legacy long caches for the reference comparison. Select the cache format explicitly and validate its manifest so an old matrix cannot silently override newly generated predictions.

## 3. Cache invariant geometry and LP structure

Files: `R/policy/cost-sensitivity.R`, `scripts/policy/optimize-cluster-bootstrap.R`.

- Precompute village-to-edge and site-to-edge index lists, variable names, site IDs, population-to-edge mapping, and scenario column indices once per geometry/population combination.
- Let the solver helper accept a prepared structure, retaining its existing default behavior for other callers.
- Cache fixed assignment/open-site constraint text while keeping row/column order, coefficient formatting, and mathematical constraints unchanged. Update only draw-dependent welfare coefficients, target, and any applicable variable costs.
- Key caches by ordered geometry, population, and formulation, including pooling status. Pooling constraints must never reuse a non-pooling template by accident.
- Replace repeated `rbind()` of growing status tables with a preallocated list and combine at checkpoints/end.

## 4. Parallelize independent draws within a scenario

Files: `scripts/policy/optimize-cluster-bootstrap.R`, `R/policy/cost-sensitivity.R`, both policy Slurm stage workers and launchers.

- Extract a `solve_one_draw()` function used identically by serial and parallel execution. Add `--num-cores`, `--draw-batch-size`, `--solver-threads`, and `--solver-seed`; preserve serial operation as a diagnostic mode.
- Use a bounded process pool with dynamically dispatched small batches. An initial batch size of 25 draws is a tuning candidate, not a required constant. Avoid launching a new process for every subsecond draw; initialize each persistent worker once, or dispatch batch tasks rather than individual draws.
- Load read-only numeric inputs before forking on Midway/Linux where practical. Keep worker inputs immutable to preserve copy-on-write memory sharing. If a different backend copies data, measure its memory cost before production.
- Start with eight workers and one solver thread per worker on an eight-CPU allocation. Pass `Threads=1` to Gurobi explicitly, along with a reproducible seed. Keep BLAS/OpenMP at one thread. Enforce `workers * solver_threads <= allocated CPUs`.
- Workers write only their own draw-specific assignment files, using temporary files followed by atomic rename. The parent alone writes the scenario status manifest, sorted by original draw ID. Worker failures must produce an explicit failed record and a failing job, not disappear from summaries.
- Checkpoints must support interruption/restart without duplicate or missing draws. Validate saved results against the input/settings manifest before reusing them.
- Propagate controls through both model and bootstrap stage scripts. Keep model/scenario Slurm arrays; do not also launch unbounded draw arrays. Remove or revise launcher memory overrides (currently 4 GB for optimization) based on measured worker RSS.

## 5. Move disposable solver files off shared storage

- Add an explicit scratch-path control. Prefer the site's allocated node-local scratch when available; otherwise use a unique job/worker temporary directory. Verify that the chosen path is actually local before claiming a benefit.
- Keep LP, solution, and transient log files there. Preserve final assignments, status, provenance, and required solver diagnostics in the persistent output root. Copy failed-solve logs before cleanup.
- Retain the current assignment layout consumed by renderers. Avoid a separate output-format migration in this pass.

## Matched test cases and acceptance criteria

Create a small comparison driver and machine-readable report, proposed as `scripts/checks/check-policy-speedup-equivalence.R`. It should consume explicit reference/candidate output roots and fail when required comparisons do not pass.

### Minimum cases

1. **Benchmark, 3.5 km:** two explicitly recorded retained draws, all five scenarios. Exercises compact prediction conversion, ordinary allocations, and the no-social-image path.
2. **Cluster-shock, 3.5 km:** two recorded draws, all five scenarios. Exercises edge-specific demand and the community-shock mapping.

Also run two full-information draws to ensure household predictions survive the shared execution refactor. Select one known tight-multinomial target-infeasible case and, if available, an undefined-equilibrium case from historical statuses; retain their outcomes exactly. If no suitable case exists in the chosen inputs, state that explicitly and test the exceptional branch with a clearly labeled constructed fixture. Do not present a constructed fixture as an empirical result.

Run a two-draw benchmark comparison at 10 km once the correctness work has supplied cap-aware, matched targets; this is also an early performance check for the denser problem.

### Comparison ladder

1. Frozen old implementation versus new serial implementation, using identical old weights, targets, input draws, solver version, and solver settings.
2. New serial with old solver settings versus new serial with explicit single-thread settings.
3. New serial single-thread versus new parallel single-thread, using identical per-draw seeds and inputs.
4. Corrected adult-weighted serial reference versus corrected parallel implementation, after the population/target audit changes are implemented. Independently verify the 144-community/39,301-adult contract in this phase.

At each step compare:

- Exact draw/scenario accounting, parameter IDs, ordered edge IDs, population vector, targets, and undefined masks. Prediction conversion should be bit-for-bit equal where it merely rearranges existing values; allow at most `1e-12` absolute numerical difference if serialization/computation requires it, and explain any such difference.
- Identical classification of complete, target-infeasible, and equilibrium-undefined cases. Independently recompute maximum attainable coverage for the non-pooling infeasibility check.
- Exactly one legal assignment per community, active-site consistency, cap compliance, and coverage recomputed directly from the saved assignment. Use the same documented target tolerance on both sides; do not weaken it to pass a comparison.
- Exact integer site objective for optimal solutions, plus certified solver optimality. For these small tests, require a proven optimum; a time-limited incumbent is inconclusive rather than evidence of equivalence.
- Exact community-to-site assignments as the first comparison. Also compare expected takers, target slack, adult/community mean distances, and downstream cost inputs; use a stated `1e-8` absolute/relative tolerance for recomputed numeric summaries.

**Multiple optima are not a blanket pass.** If site count and feasibility match but assignments differ, record both assignments and all distance/cost consequences. Label this “objective equivalent, allocation/reporting different.” Investigate ordering, seed, thread count, and solver settings. Do not claim the same reported answer or enable the change for production until the difference is resolved or an explicitly documented reporting decision is made. Adding a new secondary objective is outside this speedup plan.

Include one interrupted/restarted run and confirm that its draw inventory, assignments, and numerical summary match the uninterrupted candidate run. Confirm a worker failure cannot yield a falsely complete manifest.

## Performance pilot and rollout

- After equivalence passes, run a fixed 32–64 draw subset on allocated compute nodes with 1, 4, and 8 workers and one solver thread each. Use identical inputs and fresh output roots; distinguish cold cache/loading costs from solve-loop throughput. Repeat enough to identify shared-filesystem or node-load noise.
- Report draws/minute, end-to-end elapsed time, peak memory, solver-time distribution, and output volume. Tune batch size only if scheduling overhead or stragglers are material. The primary target is lower elapsed time within the allocation, not an assumed eightfold gain.
- Keep the single-thread solver default for this pass. Test 4 workers × 2 solver threads only if the 10 km pilot shows solver-dominated draws, and rerun the assignment/reporting comparison for that setting.
- Choose concurrency and memory from measured results, allowing headroom for the largest cap. Model/scenario concurrency is additionally limited by scheduler/account capacity.
- Return the code diff, frozen reference manifest, comparison CSV, mismatch details if any, performance table, and literal pilot commands/job IDs. Update the production runtime estimate from measured throughput before launching the full rerun.
- Production requires both this equivalence gate and the separate adult-population correctness audit to pass. Existing historical results remain archived; use new output roots for corrected runs.
