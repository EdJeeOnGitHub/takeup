# Remaining Refine TODOs

The editorial items below remain unassigned. The exponential cluster-weight
production follow-up is an explicit computational handoff for Ed/the analysis
workflow.

- [ ] Reframe the 12 table-opening paragraphs so they lead with the main takeaway.
- [ ] Simplify the Mermaid diagram and clarify that the SMS treatment is excluded specifically from the take-up analysis sample.
- [ ] Clarify whether baseline correspondents were intended to be sampled at endline.
- [ ] Explain why the structural model is necessary without a zero-observability arm.
- [ ] Add discussion of the observability robustness checks to the text.
- [ ] Add the site-consolidation footnote explaining that observability could move in either direction.

## Exponential cluster-weight production follow-up

- [x] Run
  `hpc/structural/submit_main_core_cluster_bootstrap_999.sh`. It now defaults to 999
  county-stratified exponential/Dirichlet cluster-weighted refits. Confirm the
  attempted, successful, included, and unique-weight counts before replacing
  the current 200-refit appendix artifacts. All 999 attempted refits completed,
  all were included, and all 999 weight vectors were unique.
- [x] Regenerate the exponential-weight compact generated quantities and add
  the finite 500--2500m No-Signal minus Any-Signal multiplier contrast. The
  older 200-refit artifacts contain the three point multipliers but did not
  save this finite-change quantity.
- [x] Run `hpc/policy/submit_policy_cluster_bootstrap.sh` against the completed
  exponential-weight modes. The pipeline now defaults to exponential weights
  and all 999 refits. The production allocation and summary completed for all
  five policy scenarios.
- [x] Run the population/cost/pooling policy analysis on the 999 canonical
  exponential-weight modes (included as the final stage of the same submitted
  production chain; all 999 draws completed):
  `Rscript scripts/policy/run-population-cost.R
  --parameter-csv=replication/inputs/policy/policy-bootstrap-parameters.csv
  --parameter-type=canonical --analysis-id=exponential-cluster-weights
  --cores=12`. Copy the allocation, break-even, pooling, and audit CSVs to
  `ref-reports/policy-cost-sensitivity/`, rerun
  `scripts/policy/assemble-population-cost.R`, and validate with
  `scripts/policy/validate-population-cost.R`. The assembler will add the
  cluster-weighted rows and uncertainty band automatically.
- [x] Put the regenerated tables in
  `appendix/structural-robustness/tables/`, and the stated replication counts
  have been updated. The standalone appendix compiles and the new policy pages
  have been visually checked.
- [ ] After final substantive review, copy the verified bundle into the
  paper/Overleaf tree.
- [ ] Do not relabel the older multinomial cluster-bootstrap policy output as
  exponential weighting. Keep conventional multinomial results archived as
  an internal check, not as the preferred paper-facing specification.
