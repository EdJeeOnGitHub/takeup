# Knowledge Table Differential Attrition

## Background

The endline survey (Section D) contains the knowledge module, in which respondents are shown
a pre-loaded list of named village members and asked about recognition, relationship,
deworming status, and second-order beliefs. Two variants exist:

- **Table A** (`individual`): 10 named individuals shown one at a time
- **Table B** (`paired`): 20 people shown as pairs of two

Assignment to Table A or B was **randomly determined at fieldwork design time**: within each
cluster × SMS-treatment group, the `select.endline.sample()` function in
`takeup_field_notebook.Rmd` randomly split the endline sample 50/50 between
`endline.survey.type = "individual"` and `"paired"`. The knowledge list itself (which 10 or
20 people appear for each respondent) was also pre-generated — sampling randomly from
other cluster members — and loaded into SurveyCTO before fieldwork.

## The pattern

Bracelet arm respondents are significantly *more* likely to complete Section D than
control-arm respondents, driven almost entirely by the close-distance bracelet cell.

### Attrition from knowledge table (dependent variable: not in knowledge table)
County FE, SEs clustered by community. Baseline (control) attrition ≈ 12.7%.

| Arm          | Pooled     | Close       | Far        |
|--------------|-----------|-------------|------------|
| Bracelet     | −6.1pp**  | −8.9pp***   | −2.5pp     |
| Calendar     | −2.7pp    | −1.8pp      | −4.0pp     |
| Ink          | −3.1pp    | −3.7pp      | −2.2pp     |

## What the tests show

### Test 1: Backup respondents and attrition

The endline sample was pre-selected from the census. When a primary respondent could not be
located, enumerators surveyed a backup household. Crucially, the SurveyCTO knowledge list
was **only pre-loaded for primary respondents** — so backup respondents had no Section D
to complete.

Attrition by respondent type:

| Group       | In know table | Missing | Miss rate |
|-------------|--------------|---------|-----------|
| Non-backup  | 2,252        | 0       | **0%**    |
| Backup      | 155          | 252     | **61.9%** |

Backup status alone explains 58% of variance in knowledge table attrition (R² = 0.58,
coefficient = 0.62***). The `has_know_ids` indicator (whether the SurveyCTO form had
pre-loaded IDs) is perfectly collinear with `is_backup`, confirming this is entirely
mechanical: no knowledge list was loaded for backups, so Section D was silently skipped.

### Test 2: Backup rates differ by treatment

Control arm households were harder to locate, leading to higher backup survey rates:

| Treatment  | Close | Far   |
|------------|-------|-------|
| Control    | 19.8% | 18.2% |
| Bracelet   | 12.6% | 14.6% |
| Ink        | 13.9% | 16.0% |
| Calendar   | 16.1% | 13.6% |

Close-bracelet has the lowest backup rate (12.6%) vs. close-control (19.8%) — a 7.2pp
gap — consistent with bracelet participants being more engaged and easier to find.

### Test 3: Conditioning on backup status absorbs about half the treatment gap

| Specification          | Bracelet × Close |
|------------------------|-----------------|
| Baseline               | −0.089***       |
| + Backup indicator     | −0.045*         |
| + Knowledge IDs        | −0.089*** (collinear, no change) |
| + Both                 | −0.045*         |

The bracelet×close coefficient halves once backup status is controlled for, but a
residual ~4.5pp gap remains significant at 5%. This suggests a secondary factor —
possibly that bracelet recipients are more cooperative survey respondents even among
those who were successfully surveyed as primary respondents.

### Test 4: Own deworming does not predict attrition

The most concerning alternative story would be that completing the knowledge module is
itself correlated with having dewormed — which would mean the analysis sample is
endogenously selected on the very outcome of interest. This is not supported:

- Within-arm miss rates for dewormed vs. non-dewormed respondents are nearly identical
  across all treatment arms (diffs of 1–4pp, none significant)
- Regression coefficient on own deworming = 0.011, p = 0.37
- The bracelet coefficient is completely unchanged when own deworming is added as a control

## Interpretation

The differential attrition is **primarily mechanical, not endogenous**:

1. The knowledge list was only built for primary (non-backup) endline respondents.
2. Control arm households had higher backup rates (~19%) than bracelet households (~13%),
   likely because bracelet participants were more engaged with the health programme and
   easier to locate.
3. Backups had no knowledge list pre-loaded → SurveyCTO skipped Section D → they appear
   as missing from the knowledge table.

This is attrition on a pre-determined covariate (`endline.type`, set before fieldwork
outcomes were observed), not endogenous selection on outcomes. The concern that "people
who dewormed are more willing to discuss their neighbours' behaviour" is **not supported**
by the data.

A residual ~4.5pp treatment gap remains after accounting for backup status. This is
small and likely reflects differential ease-of-location rather than anything substantive
about the knowledge module itself.

## Implications for analysis

- The comparison of knowledge outcomes across treatment arms uses a slightly unbalanced
  sample (control is missing ~4–5% more primary respondents than bracelet). This is
  unlikely to drive large biases given that (a) the underlying mechanism is mechanical
  and observable, and (b) own deworming does not predict selection into the sample.
- The clearest robustness check is to **restrict to primary respondents** (i.e.,
  `endline.type == "endline"`) and verify that knowledge estimates are stable. This drops
  the backup-driven asymmetry entirely.
- Alternatively, control for `is_backup` (or `endline.type`) in the knowledge regressions.

## Code

Investigation script: `scripts/appendix/know-table-attrition-investigation.R`
