# RCT Design and Fieldwork

This directory contains the scripts that implement the experimental design for the deworming take-up RCT: cluster selection, village targeting, and treatment randomization.

## Pipeline Overview

```
rct_targetable_schools.rds        (spatial school/village data)
        │
        ▼
takeup_rct_assign_clusters.R      Step 1: Select RCT clusters spatially
        │
        ▼
takeup_rct_prep.R                 Step 2: Assign villages to distance groups, select target villages
        │
        ▼
takeup_field_notebook.Rmd         Step 3: Clean census, randomize treatment arms, prepare fieldwork data
        │
        ▼
takeup_processed_cluster_strat.rds   (final cluster stratification with treatment assignments)
```

## Step 1: Cluster Selection (`takeup_rct_assign_clusters.R`)

`get.rct.clusters()` (line 180) selects non-overlapping clusters of schools/villages from the spatial data using an iterative buffer-based algorithm.

- Schools are grouped into two spatial groups based on required buffer radii:
  - `control.ink` — 3000m outer radius
  - `bracelet.airtime` — 4000m outer radius
- Inner radius for all clusters: 2500m
- Target: 79 clusters per group (158 total)
- Selection ensures clusters don't overlap spatially and meet minimum population/area thresholds

Output: `rct_cluster_selection_2.0.rds`

`get.cluster.villages.data()` (line 387) then identifies candidate villages within each selected cluster and classifies them by distance to their Point of Treatment (PoT).

## Step 2: Distance Group Assignment (`takeup_rct_prep.R`)

Lines 449–484 take the selected clusters and assign each to a distance condition.

### Key variables

| Variable | Level | Values | Description |
|---|---|---|---|
| `village.dist.cat` | Village | `"close"` / `"far"` | Objective distance classification: ≤1250m from PoT = close, >1250m = far |
| `cluster.dist.cat` | PoT cluster | `"close"` / `"far"` / `"mixed"` | Summary of all villages in the cluster: all close, all far, or a mix |
| `assigned.dist.cat` | PoT cluster | `"close"` / `"far"` | **Randomized assignment.** Clusters that are unambiguously close/far keep their label. Mixed clusters are randomly assigned close or far to balance counts within each `cluster.group` (roughly 50/50). |

After `assigned.dist.cat` is determined, one village is randomly selected from those matching the assigned distance category (`village.dist.cat == assigned.dist.cat`), marked as `selected.targeted = TRUE`.

Output: `rct_targetable_schools_2.0.rds`

## Step 3: Treatment Arm Randomization (`takeup_field_notebook.Rmd`)

Lines ~1422–1463 perform stratified randomization of the four treatment arms.

### Treatment arms

| Arm | Signal? |
|---|---|
| `control` | No |
| `calendar` | No |
| `ink` | Yes |
| `bracelet` | Yes |

### Stratification

Randomization is stratified by **county** × **distance-to-PoT group** (`dist.pot.group`: close/far, derived from village-centroid distance with a 1250m threshold). Within each stratum, clusters are assigned to one of the four arms with balanced allocation (remainders handled by random sampling).

### Process

1. Build `cluster.strat.data` by joining household counts, village-PoT distances, cluster density, and phone ownership
2. Stratify by `county` × `dist.pot.group`
3. Within each stratum, assign `assigned.treatment` via balanced random allocation
4. Save to `takeup_cluster_randomization_1.0.rds`, then process and re-save as `takeup_processed_cluster_strat.rds`
5. Merge `assigned.treatment` back into `hh.census.data` and `census.data`

## Other Files

- **`clean-baseline-data.R`** — Cleans baseline survey data (household and individual level)
- **`takeup_assigned_distance.qmd`** — Analysis of assigned vs actual distances to PoT
