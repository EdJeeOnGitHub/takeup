# Policy cost and pooling sensitivity (review draft)

This directory is intentionally separate from the paper and the assembled
structural-robustness appendix. Nothing here is input by `ECM ReStud.tex`.

The local pilot extends the sparse geographic PoT optimization in four ways:

1. it weights the 144 experimental communities by adult population from
   `data/takeup_census.RData`;
2. it retains the paper's 3.5 km cap and separately documents why a 2.5 km cap
   leaves essentially no scope for site consolidation in this geography;
3. it adds an illustrative fixed site cost, the documented $0.20 Bracelet
   cost per expected taker, and a grid of household travel resource costs; and
4. it lets Bracelet observability for assignments away from a community's
   experimental site move 0%, 50%, or 100% toward the Control observability
   schedule. Each experimental site is matched to its nearest candidate school;
   all 144 matches are unique. This is a conservative, pre-specified proxy for
   the information change caused by site consolidation.

The optimization preserves predicted take-up under the experimental
allocation. Population-weighted specifications preserve the predicted share
of census adults treated. The cost scenarios minimize site, signal, and
travel resource costs subject to that target. Travel cost is applied to
expected takers and round-trip distance. It is not added to the demand index,
where the estimated distance cost already affects participation.

The baseline results propagate all retained fit-105 posterior draws. The
preferred cluster analysis propagates 999 county-stratified exponential
cluster-weighted modes through the same allocation code. The conventional
multinomial cluster bootstrap is not promoted.

The cost scenarios retain the 3.5 km cap so that site consolidation remains
feasible; the 2.5 km row separately shows the consequence of staying within
experimental distance support. The consolidation rows minimize site count without
adding signal or travel costs, cleanly isolating the observability assumption.

The $100 site cost and $0.10/$0.25 per round-trip kilometre travel values are
illustrative scenarios, not estimated welfare parameters. The fixed site cost
roughly spans the paper's stipend-only implementation calculation but excludes
training, supervision, and transport. The final presentation should emphasize
cost grids or break-even thresholds rather than a single calibration.

Run locally from the repository root with:

```bash
python scripts/policy/extract-cmdstan-draws.py \
  --output=temp-data/policy-cost-sensitivity/fit105-compact.csv \
  data/stan_analysis_data/dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-{1,2,3,4}.csv
Rscript scripts/policy/run-cost-sensitivity.R
```

The local runner writes a standard LP and uses the installed open-source GLPK
solver. The same model can use `gurobi_cl` on a licensed machine.

Outputs:

- `policy-allocation-draws-<analysis>.csv`: draw-level equal- and
  population-weighted allocations;
- `policy-break-even-draws-<analysis>.csv`: draw-level resource-cost
  thresholds;
- `policy-pooling-draws-<analysis>.csv`: observability-under-consolidation
  sensitivity;
- `policy-run-audit-<analysis>.csv` and
  `policy-input-audit-<analysis>.csv`: success and geography audits; and
- the assembled tables and figure under
  `appendix/structural-robustness/{tables,figures}/`.
