# Monitoring Status Check: `monitored` vs `true.monitored`

## Background

The census data distinguishes between two monitoring variables:

- **`monitored`**: Whether an individual was *eligible* for monitoring (SMS control, or SMS-treated who consented, or in baseline sample pool)
- **`true.monitored`**: Whether we *actually have PoT survey data* for them

The distinction was introduced by Karim in commit `e41212e` (Apr 2017):

> "Identify true monitoring; for the wave 1 some individuals did not make it to the PoT monitoring list (these are name matched)"

For **wave 2**, all monitored individuals are automatically `true.monitored`. For **wave 1**, an individual is only `true.monitored` if their `KEY.individ` appears in `takeup_survey_dict_wave1.csv`, a lookup table mapping numeric `person_key` values from the field survey to census identifiers.

## Findings

### 1. 991 individuals are `monitored` but not `true.monitored`

All are in wave 1:

| wave | monitored | true.monitored | n      |
|------|-----------|----------------|--------|
| 1    | TRUE      | FALSE          | 991    |
| 1    | TRUE      | TRUE           | 5,494  |
| 2    | TRUE      | TRUE           | 7,478  |

### 2. 147 valid survey entries don't match the dict

Of 38,462 rows in the wave 1 PoT survey CSV, most are empty/placeholder rows (`person_key = 0` or `NA`). After filtering to valid entries (`cluster_key` non-NA, `person_key > 0`), there are **2,687 valid rows**, of which **147 across 19 clusters** have no match in `takeup_survey_dict_wave1.csv`. 74 of these fall in clusters that also contain `monitored & !true.monitored` individuals.

The wave 1 survey uses numeric `person_key` values that require the dict to link back to census `KEY.individ`. The Kakamega survey uses string keys that are already `KEY.individ` — zero unmatched entries there.

### 3. Name-matched individuals are in the analysis data

The `prepare.analysis.data()` function in `R/common/analysis.R` handles the gap via fuzzy name matching. For anyone with unknown or negative deworming status, it runs `name.match.monitored()` which computes string distances between census names and PoT survey names within each cluster.

In the final `analysis.data` (SMS control sample):

| name_matched | monitored | n      | dewormed.any = TRUE |
|--------------|-----------|--------|---------------------|
| FALSE        | TRUE      | 9,805  | 3,767               |
| TRUE         | FALSE     | 25,190 | 4,150               |

Name-matched individuals have `dewormed = NA` (since `dewormed` is only set for `true.monitored`) but `dewormed.any` and `dewormed.matched` are populated from the fuzzy match.

## Conclusions

1. **Not a coding error.** The `monitored`/`true.monitored` distinction is intentional and documented. Wave 1 had a known gap where some eligible individuals weren't on the PoT monitoring roster.

2. **Name matching compensates.** The 991 gap individuals (and others) are recovered via fuzzy name matching. 4,150 name-matched individuals were identified as dewormed through this process.

3. **Name-matched individuals are excluded from the structural model.** `run_takeup.R` filters to `mon_status == "monitored"` (i.e., `true.monitored == TRUE`) and passes `dewormed` (not `dewormed.any`) to Stan. So the 25,190 name-matched individuals are not used in estimation. Their `dewormed.any` status is available in the analysis data but is not currently used.

4. **147 orphan survey entries remain.** These are valid PoT survey rows with no dict match. They may represent walk-ins, roster errors, or data entry mistakes in `person_key`. The name matching may have recovered some of these, but we cannot confirm a 1-to-1 correspondence without further investigation.
