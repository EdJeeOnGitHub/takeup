# Cluster-bootstrap policy exercise

## Design

The policy robustness exercise uses all 210 valid county-stratified multinomial
cluster-bootstrap mode refits. The experimental welfare target is held fixed at
the existing baseline target. Results are medians and 2.5th--97.5th percentile
intervals across modes, not posterior-within-bootstrap credible intervals.

Only village--PoT edges within the paper's 3.5km policy cap enter prediction
and optimization. This leaves 1,252 feasible edges and 1,246 unique distances,
compared with 157,248 pairs in the legacy dense file. The exact distance-level
demand output remains sufficient to regenerate the demand-over-distance plot.

Dynamic Control and Bracelet fixed points use vectorized bracketed bisection.
Non-bracketed points fall back to the legacy `nleqslv` solve initialized at
minus private benefit, preserving its root-selection convention. Fixed-signal
and no-social-image predictions use their exact closed forms.

## Validation and performance

- Five-mode validation covered both dynamic scenarios at all 1,246 feasible
  distances. Maximum demand disagreement with legacy `nleqslv` was
  `3.75e-9`; median solver-only speedup was 14.7x.
- The 210-mode, five-scenario demand stage took 5.98 seconds, peaked at 402MB,
  avoided 163,802,100 dense rows, and used 229 legacy fallbacks.
- The five original 200-mode production optimizer tasks took 22--37 seconds
  each; adding the final ten modes used the checkpointed outputs rather than
  repeating completed solves. All 840 feasible-target problems were optimal.
  All 210 no-social-image cases were
  explicitly classified as target-infeasible and report the best/closest-site
  benchmark rather than a silent optimizer fallback.
- Bootstrap replicate 205 had a status-marked mode CSV containing a header but
  no optimizer row. It is recorded in `policy-bootstrap-excluded-modes.csv`;
  all remaining 210 valid modes are retained.

## Reproduction

Submit the complete dependency chain from the Midway repository root:

```bash
bash hpc/policy/submit_policy_cluster_bootstrap.sh
```

For a five-mode recovery test:

```bash
NUM_REPLICATES=5 OUTPUT_PATH=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots-cluster-bootstrap-smoke \
  bash hpc/policy/submit_policy_cluster_bootstrap.sh
```

The production table is
`presentations/tables/fit105/optim-summ-cluster-bootstrap.tex`. Compact local
summaries and manifests are stored in
`temp-data/policy-cluster-bootstrap/`; full distance-demand and checkpointed
allocation files remain on Midway under the matching `optim/data` directory.
