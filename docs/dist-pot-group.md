# Randomization/Distance Group Generation

## Step 1: Drawing RCT Cluster Locations

**File:** `rct-design-fieldwork/takeup_rct_prep.R`<br>
**Inputs:** `rct.schools.data` from `~/Dropbox/Data/TakeUp/Kenya_Primary_Schools.csv`, `~/Dropbox/data/TakeUp/Kakamega schools geolocation.csv`<br>
**Outputs:** `rct_cluster_selection_2.0.rds`, `rct_targetable_schools_2.0.rds`

Uses a greedy spatial packing algorithm to generate 158 RCT clusters using school/PoT location data. These are guaranteed to be 50-50 split between bracelet/calendar and ink/control with different sized buffers and a validity checking algorithm for overlapping clusters that takes chunks out of the buffer depending on nearby schools.

This produces `data/rct_cluster_selection_2.0.rds` (saved interactively, not by a committed `write_rds` call) which has the 158 clusters and their bracelet.airtime / control.ink group assignment.

Next, the file creates `rct.targetable.schools`. This takes the 158 clusters and categorises them into `close`, `far`, or `mixed` based on how many other schools (i.e. potential PoTs) there are nearby:

- **Close:** other schools within the cluster buffer and <= 1.25km away
- **Far:** other schools within the buffer and > 1.25km and <= 2.5km
- **Mixed:** both close and far schools present

The algorithm assigns close clusters to close and far to far. Then, it randomizes the mixed clusters to ensure total close/far are balanced (N.B. an off-by-one bug means actual assignment is 80-78, not 79-79). This is saved as `rct_targetable_schools_2.0.rds`.

Key variables in `rct_targetable_schools_2.0.rds`:
- `assigned.dist.cat` -- cluster-level randomized distance assignment (close/far)
- `pot.cluster.id` -- cluster ID
- `cluster.dist.cat` -- close/far/mixed classification based on schools in the area

## Step 2: RCT Cluster Survey

**File:** `rct-design-fieldwork/takeup_cluster_survey.R`<br>
**Inputs:** `rct_cluster_selection_2.0.rds`, `data/Cluster Survey V1.dta`, `data/Cluster Survey 27.07.2016.dta`<br>
**Outputs:** `takeup_cluster_survey.rds`

Enumerators are sent to each cluster and survey villages in the surrounding area. This file takes these village surveys and joins them onto the RCT cluster (PoT), saving the output as `takeup_cluster_survey.rds`.

## Step 3: Assigning RCT Villages

**File:** `rct-design-fieldwork/takeup_rct_target_villages.R`<br>
**Inputs:** `rct_cluster_selection_2.0.rds`, `rct_targetable_schools_2.0.rds`, `takeup_cluster_survey.rds`<br>
**Outputs:** `rct_target_villages_2.0.rds`

Loads the actual village data (`takeup_cluster_survey.rds`) and calculates the actual distance between the cluster PoT and the villages. Then classifies each village-PoT relationship as close/far depending on distance.

The script has swapping logic (lines 100-102) that flips `village.dist.cat` for clusters where `dist.switchable=TRUE` (i.e. the only valid target villages are in the "wrong" distance category). However, this swap modifies `rct.cluster.selection@data` in memory only -- it is never saved to disk. The `rct.villages` object that gets saved to `rct_target_villages_2.0.rds` is built from `cluster.survey.data`, which carries the original `village.dist.cat` from `rct_targetable_schools_2.0.rds` (joined at line 68). This is why `village.dist.cat` in `rct_target_villages_2.0.rds` always matches `assigned.dist.cat` -- the swap was computed but never persisted.

`rct_target_villages_2.0-3.rds` also exists in `data/` but is not referenced anywhere in the codebase; it appears to be an orphaned backup file (possibly from a manual re-run).

## Step 4: Treatment Arm Randomization

**File:** `rct-design-fieldwork/takeup_field_notebook.Rmd`<br>
**Inputs:** `hh.census.data` (census household GPS), `vill.pot.dist` (computed earlier in notebook), `takeup_cluster_randomization_1.0.rds`<br>
**Outputs:** `takeup_processed_cluster_strat.rds`

Performs stratified randomization of the four treatment arms (control, calendar, ink, bracelet) by county x `dist.pot.group`.

Critically, the `dist.pot.group` used here is **not** the same as `assigned.dist.cat` from Step 1. It is recomputed from scratch via `vill.pot.dist` (section "Village-PoT Distances", line 1210):

1. **`pot.info`** (line 102): PoT coordinates from the Cluster Survey V3 (`Cluster Survey V3 July 04.csv`) and V1 (`Cluster Survey V1.dta`). The `alt.pot.lon`/`alt.pot.lat` fields come from enumerator-collected GPS (`gps2-Longitude`/`gps2-Latitude`), with fallbacks to `read_lon_ant`/`read_lat_ant` and then `manual_long2`/`manual_lat2`. These are joined with a separate PoT verification survey (`POT verification.csv`) which provides `lon.verify`/`lat.verify`.
2. **`village.centers`** (line 1215): For each cluster, takes **all census household GPS coordinates** from `hh.census.data`, converts to spatial points, and computes the `gCentroid` — i.e., a single (lon, lat) per cluster representing the centroid of surveyed households. Then joins `pot.info` to get the PoT coordinates (`alt.pot.lon`/`alt.pot.lat` and `lon.verify`/`lat.verify`).
3. **`get.vill.pot.dist`** (line 1228): Computes a full pairwise distance matrix (`gDistance` with Kenya projection) between all village centers and all PoTs. Extracts `dist.to.own.pot` as the diagonal (distance from each village center to its own PoT), plus `closest.other` (nearest other PoT) and `num.within.*` counts (how many PoTs within various radii).
4. **`vill.pot.dist`** (line 1246): Applies `get.vill.pot.dist` to `village.centers` filtered to clusters with known PoT locations. The default uses `alt.pot.lon`/`alt.pot.lat`; a second version (`verify.vill.pot.dist`, line 1250) uses `lon.verify`/`lat.verify` from the verification survey.
5. **`dist.pot.group`** (line 1237): Assigned by `factor.dist.pot` (line 1213): close if `dist.to.own.pot` <= 1250m, far otherwise.
6. Both `vill.pot.dist` and `verify.vill.pot.dist` are saved to `data/takeup_village_pot_dist.RData` (line 1254).

The randomization chunk (line 1422, `eval=FALSE`) joins `vill.pot.dist` onto `hh.census.data` and stratifies by `county x dist.pot.group`. This was run interactively and saved as `takeup_cluster_randomization_1.0.rds`. The live chunk (line 1444) reads this file back, adds wave/county info, and writes `takeup_processed_cluster_strat.rds`.

This is the source of the 26 mismatches documented below: `assigned.dist.cat` (Step 1) is based on distances between schools/potential PoTs, while `dist.pot.group` (Step 4) is based on census household centroid-to-PoT distance. These two distance measures disagree for 26 clusters.

## Feasibility Drops (158 → 144)

Source: `scratch/randomization-checks.R`

14 clusters were dropped before the analysis, reducing the sample from 158 to 144. Reasons are documented in `scratch/distance/distance-recreation-attempt.Rmd` (lines 49–53) and `rct-design-fieldwork/takeup_rct_target_villages.R` (lines 6–7).

| Cluster | Original cat | Reason |
|---|---|---|
| 1 | far | Village dispute about PoT |
| 99 | close | Problems at cluster survey stage |
| 201 | close | No longer valid: too close to another PoT |
| 239 | far | Problems at cluster survey stage |
| 277 | far | Too close to another cluster (503) |
| 491 | close | Problematic urban cluster |
| 492 | close | Problematic urban cluster |
| 678 | far | Hostile community member |
| 691 | far | Problems at cluster survey stage |
| 737 | close | Data fabrication and medication theft |
| 853 | far | No longer valid: too close to another PoT |
| 892 | close | Problems at cluster survey stage |
| 1293 | far | Problems at cluster survey stage |
| 1402 | far | No longer valid: too close to another PoT |

This removes 6 close and 8 far clusters. Combined with the 26 category switches (see below), the full count progression is:

| Stage | Close | Far |
|---|---|---|
| Design (original randomization) | 80 | 78 |
| Dropped (feasibility) | −6 | −8 |
| After drops | 74 | 70 |
| Category changes: close → far | −10 | +10 |
| Category changes: far → close | +16 | −16 |
| **Final analysis** | **80** | **64** |

Note the coincidence: 6 close clusters dropped but 16 far→close switches, leaving close at exactly 80 again.

## Randomization Mismatches

Source: `scratch/randomization-checks.R` (output saved to `docs/randomization-mismatches.csv`)

26 clusters where `processed_cat` (from `takeup_processed_cluster_strat.rds`) differs from `clust_randomization_cat` (the original `assigned.dist.cat` from `rct_targetable_schools_2.0.rds`): 10 close -> far, 16 far -> close. In every case, `swap_cat` (from `rct_target_villages_2.0.rds`) matches the original randomization, confirming the swap in Step 3 was never persisted.

The root cause is that `dist.pot.group` in the processed data is computed from a **different distance measure** than the original `assigned.dist.cat` (see Step 4): census household centroid-to-PoT distance vs. school-to-school distance. These two methods disagree for 26 clusters near the 1250m threshold.

| Cluster | Classification | Randomization | Swap | Processed | Analysis |
|---|---|---|---|---|---|
| 108 | mixed | far | far | close | close |
| 141 | mixed | far | far | close | close |
| 149 | mixed | far | far | close | close |
| 164 | far | far | far | close | close |
| 213 | close | close | close | far | far |
| 215 | mixed | close | close | far | far |
| 235 | far | far | far | close | close |
| 238 | mixed | far | far | close | close |
| 408 | close | close | close | far | far |
| 587 | far | far | far | close | close |
| 631 | mixed | close | close | far | far |
| 648 | mixed | close | close | far | far |
| 735 | close | close | close | far | far |
| 738 | far | far | far | close | close |
| 796 | close | close | close | far | far |
| 830 | mixed | far | far | close | close |
| 835 | far | far | far | close | close |
| 910 | mixed | close | close | far | far |
| 954 | far | far | far | close | close |
| 1014 | close | close | close | far | far |
| 1059 | mixed | far | far | close | close |
| 1120 | far | far | far | close | close |
| 1219 | mixed | far | far | close | close |
| 1272 | close | close | close | far | far |
| 1347 | mixed | far | far | close | close |
| 1370 | far | far | far | close | close |

## Dataset Reference

| Dataset | Created by | Key variables | Description |
|---|---|---|---|
| `rct_cluster_selection_2.0.rds` | `takeup_rct_prep.R` (interactive) | `cluster.id`, `cluster.group`, `county` | 158 RCT clusters with spatial buffers and bracelet.airtime/control.ink group assignment |
| `rct_targetable_schools_2.0.rds` | `takeup_rct_prep.R` | `pot.cluster.id`, `assigned.dist.cat`, `cluster.dist.cat`, `village.dist.cat`, `selected.targeted` | Clusters classified into close/far/mixed distance groups with randomized assignment and selected target villages |
| `takeup_cluster_survey.rds` | `takeup_cluster_survey.R` | `cluster.id`, `target.village.id`, `target.lat`, `target.lon` | Village-level survey data joined to cluster PoTs |
| `rct_target_villages_2.0.rds` | `takeup_rct_target_villages.R` | `cluster.id`, `target.village.id`, `dist.to.pot`, `valid.target.village` | Final selected target villages with actual distances to PoT |
| `takeup_village_pot_dist.RData` | `takeup_field_notebook.Rmd` | `dist.to.own.pot`, `dist.pot.group`, `closest.other` | Village center-to-PoT distances computed from census household centroid GPS |
| `takeup_cluster_randomization_1.0.rds` | `takeup_field_notebook.Rmd` (interactive) | `cluster.id`, `dist.pot.group`, `assigned.treatment` | Randomization output; `dist.pot.group` comes from `vill.pot.dist` |
| `takeup_processed_cluster_strat.rds` | `takeup_field_notebook.Rmd` | `assigned.treatment`, `dist.pot.group` | Processed version of randomization with wave/county joined |
