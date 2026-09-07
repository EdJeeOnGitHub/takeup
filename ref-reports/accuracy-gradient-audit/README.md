# Accuracy-gradient audit — 2026-09-07

## Question and findings

Do conditional reporting-accuracy gradients vary by treatment, and is the
constant-accuracy restriction in the historical two-stage model plausible?

The full-sample descriptive estimates show improving correct No reports among
definite answers about nonparticipants as distance increases. A joint test of
zero accuracy slopes across the four arms rejects (p=0.000495). The corresponding
test for actual participants does not reject (p=0.544).

The evidence does not establish different accuracy slopes across treatment arms:
community-clustered joint equality tests give p=0.188 for nonparticipants and
p=0.605 for participants. These are not equivalence tests and do not demonstrate
that common slopes are correct. In particular, uncertainty permits economically
meaningful differences.

Fitted changes in conditional accuracy from 0.5 to 2.5 km, percentage points
(95% community-clustered delta-method intervals):

| True status | Control | Calendar | Ink | Bracelet |
|---|---:|---:|---:|---:|
| Nonparticipant | 33.5 [8.8,58.2] | 21.6 [-2.7,45.9] | 9.6 [-4.2,23.4] | 30.6 [12.3,48.9] |
| Participant | -8.9 [-28.2,10.5] | -5.1 [-18.6,8.5] | 2.9 [-10.4,16.1] | -8.3 [-22.3,5.7] |

Excluding the five dispersed communities does not change the broad conclusion:
zero-slope p=0.00139 for nonparticipants and 0.585 for participants; equal-slope
p=0.281 and 0.830 respectively.

## Sample and estimation

Input: `data/clean-data/clean-endline-know-table-data-long.rds`, linked to
`build/structural-workspace/main-core-input.RData` for active communities and
community-centroid distance to the treatment site. The existing structural data
builder independently verifies recognition/definite counts against the saved
full input before estimation.

Named Table-A peers, SMS-control respondents, at least one recognized peer;
measurement regressions restrict to recognized peers with administrative truth.
The full sample has 144 communities, 1,141 respondents, 11,410 peer records,
10,456 linked records, 4,962 recognized linked records, and 3,926 definite
linked records. The historical exclusion recovers 139 communities, 1,098
respondents and 10,079 linked records, matching the earlier fit's documented audit.

Separate logistic regressions for each true-status group and each stage:
definite answer among recognized linked peers; correct answer among definite
linked responses. Each regression has unrestricted arm intercepts and arm-specific
distance slopes. Distance is continuous centroid-to-site distance in kilometres.
No county or demographic controls are added: this is a descriptive diagnostic
of the reporting schedule used by the structural model, not a causal distance
estimate or a replication of Table B16's adjusted reduced form.

Uncertainty uses `sandwich::vcovCL`, community clustering, HC1 adjustment;
intervals use t(G-1), joint Wald tests use F(q,G-1). This permits dependence
among respondent and peer reports within a community. P-values are not adjusted
for multiple comparisons. Accuracy conditions on giving a definite answer;
changes can reflect selection into answering, not improvements for fixed dyads.

## Implication for structural work

The empirical pattern motivates relaxing constant conditional accuracy. It does
not establish that Bracelet accuracy improves uniquely with distance, nor that
pooling will recover benchmark multipliers. Compare a full-sample two-stage fit
with arm-specific accuracy intercepts and truth-specific common distance slopes
against a version allowing arm-specific slopes (partial pooling is a sensitivity).
Keep the definite-answer specification, other priors, sample and evaluation
weights fixed. Evaluate reporting and take-up fit before multiplier/policy claims.

This audit does not estimate structural parameters or multipliers, and does not
recover the historical two-stage ATEs. Full structural refits remain server work.
The existing fitted report matrices are joint structural predictions; these
unrestricted measurement-only estimates need not reproduce them.

## Files and reproduction

Run `Rscript ref-reports/accuracy-gradient-audit/audit.R` from the repository root.
Outputs: `sample-counts.csv`, `slopes.csv`, `predictions.csv`, `joint-tests.csv`,
`accuracy-gradients.pdf`, and `session-info.txt`. No manuscript files are changed.
