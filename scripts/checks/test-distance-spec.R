#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(readr)
})
source("R/distance/spec.R")

crosswalk <- takeup_distance_crosswalk()
analysis_crosswalk <- filter(crosswalk, in_main_analysis)
stopifnot(
  nrow(analysis_crosswalk) == 144L,
  sum(analysis_crosswalk$switched) == 26L,
  identical(
    as.integer(table(analysis_crosswalk$assigned_dist_group,
                     analysis_crosswalk$realized_dist_group)),
    c(64L, 16L, 10L, 54L)
  )
)

data <- read_csv("temp-data/analysis-cluster-covariate-data.csv",
                 show_col_types = FALSE)
cluster_expected <- read_csv("data/cluster_expected_dist.csv",
                             show_col_types = FALSE) |>
  transmute(cluster_id = as.integer(cluster.id), clust_expected_dist = dist)
data <- data |>
  left_join(cluster_expected, by = "cluster_id") |>
  mutate(mu_d = clust_expected_dist / sd(cluster.dist.to.pot))
assigned <- takeup_apply_distance_spec(data, crosswalk, "assigned") |>
  mutate(
    assigned_treatment = factor(assigned.treatment),
    cluster.id = factor(cluster_id),
    wt = 1
  )
realized <- takeup_apply_distance_spec(data, crosswalk, "realized") |>
  mutate(
    assigned_treatment = factor(assigned.treatment),
    cluster.id = factor(cluster_id),
    wt = 1
  )

continuous_fit <- function(x) {
  feols(
    dewormed ~ 0 + assigned_treatment + standard_cluster.dist.to.pot +
      i(assigned_treatment, cluster.dist.to.pot, "control") + mu_d +
      female + age.census | county,
    data = x, weights = ~wt, warn = FALSE, notes = FALSE
  )
}
assigned_coef <- coef(continuous_fit(assigned))
realized_coef <- coef(continuous_fit(realized))
stopifnot(
  identical(names(assigned_coef), names(realized_coef)),
  isTRUE(all.equal(unname(assigned_coef), unname(realized_coef),
                   tolerance = 1e-12)),
  sum(
    assigned |>
      distinct(cluster.id.x, analysis_dist_group) |>
      arrange(cluster.id.x) |>
      pull(analysis_dist_group) !=
    realized |>
      distinct(cluster.id.x, analysis_dist_group) |>
      arrange(cluster.id.x) |>
      pull(analysis_dist_group)
  ) == 26L
)

cat("Distance specification tests passed.\n")
