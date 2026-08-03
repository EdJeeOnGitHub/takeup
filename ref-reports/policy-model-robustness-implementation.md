# Policy robustness across structural models

The fast sparse policy workflow now applies the same 3.5 km feasible-site set,
experimental-allocation welfare target, and five information counterfactuals to
six fitted alternative structural models:

1. correct-classification observability;
2. perceived second-order observability;
3. Any-Signal/No-Signal social-image weights (log-ratio prior SD 0.25);
4. treatment-specific social-image weights (pairwise log-ratio prior SD 0.25);
5. variance-normalized Student-t(5) intrinsic-motivation types; and
6. the cluster random-shock model (shock prior SD 0.1).

All retained HMC draws are used: 800 draws for each legacy observability fit,
1,600 for each lambda and Student-t fit, and 4,000 for the cluster-shock fit.
The baseline comparison remains the 210 valid county-stratified cluster-
bootstrap mode refits. Thus, the baseline interval measures bootstrap
variation while alternative-model intervals measure posterior uncertainty;
the table labels this distinction explicitly.

The Gaussian/common-lambda models use the existing vectorized fixed-point
solver. Signal-specific lambda draws are transformed with the exact grouped or
Helmert parameterization used in Stan. The Student-t predictor ports the same
12-component generalized Gauss-Laguerre normal-variance mixture and eight-step
Newton solve used by the fitted model. Cluster-shock predictions retain every
village's posterior shock and therefore evaluate the 1,252 sparse village-site
edges directly rather than treating demand as only a function of distance.

If a posterior draw cannot attain the fixed welfare target, the workflow saves
the best feasible allocation and marks the draw `target_infeasible`; it does
not discard the draw. The assembled table reports the share of such draws.
If an alternative-model draw has no fixed point under a new counterfactual,
the workflow records `equilibrium_undefined` rather than inserting a
pseudo-root. It retains that draw in the diagnostic denominator and reports
the undefined share separately.

Optimization uses a `1e-5` tolerance on the aggregate welfare target, aligned
with the MIP solver's numerical feasibility scale. In the production audit the
triggering edge case missed the aggregate target by `8.3e-7`, or less than
`6e-9` predicted takers per village; it is therefore treated as numerically
feasible rather than as target-infeasible.

In the production run, Gaussian alternative-model prediction took 14--27
seconds, the 1,600-draw Student-t prediction took 29 seconds, and the full
4,000-draw cluster-shock prediction took 40 seconds. Compact numeric storage
reduced the cluster-shock prediction object to 191 MB and peak job memory to
about 1.9 GB. Any policy-edge counterfactual that crosses a saddle-node
boundary is reported as undefined rather than being replaced by the nearby
minimum-residual point. The correctly mapped production matrix contained one
such value among 25,040,000 predictions, affecting 0.025 percent of the
static-Bracelet draws and none of the other policy scenarios.

The cluster-shock adapter constructs an explicit crosswalk from Stan's
fit-106 `cluster_id` order to the policy distance file's external `cluster.id`.
This is required because the two row orders are not identical. Preparation
fails if the crosswalk is incomplete, duplicated, or does not cover all 144
policy villages. The audit confirms that all external IDs match but 109 of 144
positions differ, so positional indexing would be invalid.

Run all models with:

```bash
bash submit_policy_model_robustness.sh
```

After the model summaries have been copied back locally, assemble the combined
CSV and portable LaTeX table with:

```bash
Rscript optim/assemble-policy-model-robustness.R
```

The principal outputs are:

- `temp-data/policy-model-robustness-summary.csv`;
- `presentations/tables/fit105/optim-policy-model-robustness.tex`; and
- the model-level draw, prediction, allocation, and summary artifacts below
  `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-model-robustness/`
  on Midway. Compact model summaries are copied to
  `temp-data/policy-model-robustness/` for table assembly and version control.

Validation is in `scratch/test-policy-model-robustness.R`, alongside the
existing fast-policy equivalence and benchmark tests.

## Production results

The paired median reduction in assigned PoTs when the policymaker uses
Bracelet rather than Control observability is 25 (95 percent interval 10--42)
under the cluster-bootstrap baseline. It remains positive and economically
large under correct classification (18, 6--29), second-order observability
(21, 13--30), grouped signal-specific lambda (28, 18--36), Student-t types
(27, 17--35), and the cluster random shock (12, 1--28). The fully
treatment-specific lambda model has a similar median (25) but a much wider
interval (-27--37); 1.31 percent of its Bracelet draws cannot attain the fixed
welfare target, exposing the expected weak-identification cost of that very
flexible specification.

When social image returns are held fixed at their 0.5 km level, paired median
savings shrink to 1--11 PoTs across models and most intervals include zero.
This contrast supports the interpretation that the larger dynamic-policy
savings operate through the endogenous distance-social-image feedback rather
than a mechanical level difference between observability regimes.
