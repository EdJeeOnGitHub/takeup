#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
})
options(fixest.notes = FALSE, dplyr.summarise.inform = FALSE)
fixest::setFixest_notes(FALSE)

args <- commandArgs(trailingOnly = TRUE)
draws <- if (length(args)) as.integer(args[[1L]]) else 30L
stopifnot(!is.na(draws), draws > 1L)

data <- read.csv("temp-data/analysis-cluster-recentered-covariate-data.csv") %>%
  select(dewormed, assigned_treatment, assigned_dist_group, female, age.census,
         mu_d, county, cluster.id, cluster_id_rank) %>%
  filter(if_all(everything(), ~ !is.na(.x))) %>%
  mutate(
    assigned_treatment = factor(
      assigned_treatment, c("control", "ink", "calendar", "bracelet")
    ),
    assigned_dist_group = factor(assigned_dist_group, c("close", "far")),
    county = factor(county)
  )

direct_formula <- ~ 0 + assigned_treatment * assigned_dist_group +
  female + age.census + mu_d + county
X <- model.matrix(direct_formula, data)
y <- data$dewormed
stopifnot(qr(X)$rank == ncol(X))

counterfactual <- bind_rows(lapply(levels(data$assigned_treatment), function(arm) {
  mutate(data, assigned_treatment = factor(
    arm, levels = levels(data$assigned_treatment)
  ))
}))
X_counterfactual <- model.matrix(direct_formula, counterfactual)
counterfactual$arm <- rep(levels(data$assigned_treatment), each = nrow(data))

draw_weights <- function(seed) {
  set.seed(seed)
  cluster_weights <- rgamma(length(unique(data$cluster.id)), shape = 1)
  cluster_weights <- cluster_weights / sum(cluster_weights)
  cluster_weights[data$cluster_id_rank]
}

collapse_predictions <- function(prediction) {
  tibble(
    arm = counterfactual$arm,
    distance = rep(data$assigned_dist_group,
                   times = length(levels(data$assigned_treatment))),
    prediction = prediction
  ) %>%
    group_by(arm, distance) %>%
    summarise(value = mean(prediction), .groups = "drop") %>%
    arrange(arm, distance) %>%
    pull(value)
}

fixest_draw <- function(seed) {
  data$wt <- draw_weights(seed)
  fit <- feols(
    dewormed ~ 0 + assigned_treatment * assigned_dist_group +
      female + age.census + mu_d | county,
    data = data, weights = ~wt, nthreads = 1
  )
  collapse_predictions(predict(fit, newdata = counterfactual))
}

crossproduct_draw <- function(seed) {
  weights <- draw_weights(seed)
  xtwx <- crossprod(X, weights * X)
  xtwy <- crossprod(X, weights * y)
  coefficients <- solve(xtwx, xtwy)
  collapse_predictions(drop(X_counterfactual %*% coefficients))
}

# Warm libraries, formula parsing, and BLAS before timing.
invisible(fixest_draw(1L))
invisible(crossproduct_draw(1L))

fixest_time <- system.time({
  fixest_results <- lapply(seq_len(draws), fixest_draw)
})[["elapsed"]]
direct_time <- system.time({
  direct_results <- lapply(seq_len(draws), crossproduct_draw)
})[["elapsed"]]

fixest_matrix <- do.call(rbind, fixest_results)
direct_matrix <- do.call(rbind, direct_results)
max_difference <- max(abs(fixest_matrix - direct_matrix))

cat(sprintf("Observations: %d; coefficients: %d; draws: %d\n",
            nrow(data), ncol(X), draws))
cat(sprintf("fixest refit + predict: %.3f seconds\n", fixest_time))
cat(sprintf("cross-products + predict: %.3f seconds\n", direct_time))
cat(sprintf("speedup: %.2fx\n", fixest_time / direct_time))
cat(sprintf("maximum collapsed-prediction difference: %.3g\n", max_difference))
if (max_difference > 1e-8) {
  stop("Direct WLS did not reproduce fixest predictions closely enough.")
}
