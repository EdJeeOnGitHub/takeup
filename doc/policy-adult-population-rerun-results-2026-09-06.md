# Corrected adult-population policy rerun

All 75 model/cap/inference/scenario combinations are complete and independently audited: 152,995 draw-level evaluations, comprising 123,578 certified target-preserving optima, 29,385 target-infeasible outcomes, and 32 undefined counterfactual equilibria. No outcomes were silently discarded.

The exercise uses 39,301 unique census adults in 144 communities and all 1,451 candidate sites. Each model/draw preserves its own experimental-Control adult target across the five scenarios. The four benchmark caps have byte-identical parameter CSVs and identical adult targets. A maximum difference of approximately 1.07e-13 in the auxiliary community-weighted welfare column reflects floating-point arithmetic; that column is not the adult target used by these allocations.

## Review artifacts

Local return directory: `temp-data/policy-adult-rerun/return/`.

- `review.pdf` / `review.tex`: nine-page review packet, including the full alternative-model table, larger-cap table, short-cap diagnostic, allocation/cost/pooling tables, break-even figure, allocation maps, and distance/demand figures.
- `reports/all-model-cap-summary.csv`: all 75 scenario cells, with target-attainment and undefined/infeasible counts. Site-count summaries in this file condition on target attainment.
- `reports/all-model-cap-replicates.csv` and `reports/all-model-cap-paired-contrasts.csv`: full numerical rows, including experimental references and explicit fallback/target-preserving contrasts.
- `reports/policy-amplification-paired-contrasts.csv` and `reports/policy-amplification-summary.csv`: paired endogenous-minus-fixed-return site-saving contrasts.
- `reports/weighting-change-comparison.csv`: corrected adult weights versus equal community weights on matched draws, plus separately labeled historical incorrect-weight results.
- `reports/independent-audit-summary.json` and `reports/independent-allocation-audit.csv`: complete independent validation evidence.
- `reports/input-source-manifest.csv`, `reports/selected-draw-manifest.csv`, and `run-records/`: source hashes, draw selection, exact submitted scripts, job IDs, settings, and recovery history.
- `policy-numerical-outputs.tar.gz`: all model/cap numerical CSV/RDS outputs and shard/attempt provenance. Its checksum is in `policy-numerical-outputs.sha256`.
- `executed-code-snapshots.tar.gz`: exact main, reporting, and precision-retry code snapshots.
- `legacy-incorrect-population-reference.tar.gz`: preserved historical outputs, clearly separated from corrected results.

The full server run remains at `/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-adult-population-20260906/production-v3`.

## Main results

| Inference / weighting | Median Control sites | Median Bracelet sites | Median paired sites saved | Central 95% interval for savings |
|---|---:|---:|---:|---:|
| Benchmark, corrected adult weights | 102 | 77 | 25 | 15.975–31 |
| Benchmark, matched equal community weights | 106 | 79 | 27 | 16–35 |
| 999 assigned-distance bootstrap modes, corrected adult weights | 102 | 78 | 24 | 7–32 |
| Bootstrap, matched equal community weights | 106 | 80 | 26 | 7–36 |

Paired savings are computed within each draw, not by subtracting separately summarized medians. Benchmark endogenous-minus-fixed-return savings have median 18 sites (95% interval 10–26); the corresponding bootstrap median is 16 (3–31).

Median corrected break-even fixed site costs:

| Inference | No travel cost | $0.10 per participant round-trip kilometre |
|---|---:|---:|
| Benchmark posterior | $103.96 | $192.63 |
| Assigned-distance bootstrap | $108.80 | $208.22 |

The historical benchmark $260.58/$482.52 thresholds used incorrect population weights and older fit/draw provenance. The old and new median saving both happen to be 25 sites, but the cost figures do not survive correction. This historical before/after comparison is not a pure weighting effect; use the matched equal-community rows to isolate the planner-weight comparison under a common target convention.

Bootstrap median site savings under 0/50/100% attenuation toward Control observability at reassigned sites are 24, 10, and 3, with intervals 7–32, 3–20.05, and -1–9 respectively. Zero attenuation reuses canonical allocations.

At the representative posterior-median parameter vector, the 3.5 km Control/Bracelet allocations use 102/76 sites. These differ from posterior distribution summaries because a solution evaluated at median parameters is a different estimand. The 2.5 km diagnostic still reaches the geographic floor of 114 sites in both regimes; at 3.5 km the geographic floor is 71, so equality with the floor is not assumed.

The distance-density figure weights communities equally. The underlying allocations preserve adult-weighted coverage; allocation tables report distances using the row's stated weights. Cost accounting uses participant-weighted round-trip travel.

## Draw provenance and verification

- Benchmark: 400 draws per assigned-distance slim chain, 1,600 total. The exact previous 50-per-chain (200 total) subset is retained and independently checked against its original parameter CSV. The older cost input was not spliced into the new draw set.
- Alternatives: all retained draws in the handoff, using their existing model-specific extraction/prediction adapters, household geography, and community-shock mappings.
- Bootstrap: all 999 complete exponential assigned-distance modes from `candidate-hpc-cd5f295-assigned/work/cluster-weight/modes`.
- Every main saved assignment was checked against original source edges, census counts, draw-specific targets, demand joins, adult-weighted summaries, infeasibility bounds, and integer-optimality certificates.
- Both full cost runs were independently checked at assignment level, including source edges, weights, targets, site counts, participant travel, and every cost-grid value. Main Control/Bracelet allocations are reused for accounting.
- Code regression checks preserve complete allocations and summaries in two old/new cost-accounting cases. Three successful finite-mixture cases also retain identical assignments, targets, and optimal site counts under tighter solver tolerances. Controller tests cover missing observations, failure exits, cancellation races, submission reconciliation, and migration destination handling.

## Execution and recoveries

The corrected production DAG took approximately 43.9 minutes, including recovery delays: 445 successful task outcomes over 467 attempts (240 on Midway2, 227 on Midway3). Recorded allocated CPU time for this DAG, including retries, was approximately 48.2 core-hours; this excludes smoke tests, reporting, auditing, and packaging. Concurrency was raised from 16 to 32 eight-core jobs per cluster after checking available capacity. Each optimizer job runs independent draws in parallel with one Gurobi thread per draw and seed zero.

Three issues were resolved and retained in the execution record:

1. The shared default distance file contained 1,092 sites. Its initial runs were invalidated and stopped. The verified 1,451-site file was staged separately and a full-geography guard added.
2. Finite-mixture draw 1327 exposed near-binary numerical residue that lost about 0.001 expected adults on extraction. Tighter solver feasibility tolerances recovered the same 103-site optimum while satisfying the unchanged target check.
3. The project reached its file-count quota. Completed solver logs were archived and byte-verified before individual copies were removed. Quota-affected tasks were retried in fresh directories. Future shard workers archive logs automatically.

Production primarily used `code-v3`, with tighter-feasibility retries recorded as `code-precision-v7`; reporting used `code-report-v6` plus documented final review/audit scripts. Saved assignments are authoritative. Equal-objective solver tie choices can affect distances and costs even when site counts agree, so the returned costs are tied to those actual saved assignments.

No structural estimation, reduced-form reruns, or Overleaf/manuscript updates were performed. See `policy-adult-rerun-execution-2026-09-06.md` for detailed hashes and recovery evidence.
