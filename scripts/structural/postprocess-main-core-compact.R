#!/usr/bin/env Rscript

# Convert the lean main-core generated quantities into the RDS contracts used
# by the paper renderer. This replaces the legacy cluster-by-draw arrays with
# the sufficient aggregated counterfactual cells saved by the compact GQ.

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (!length(hit)) default else substring(hit, nchar(prefix) + 1L)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
})

gq_path <- option_value(
  "--compact-gq-path", "build/assigned/structural/compact-gq"
)
fit_path <- option_value("--fit-path", "build/structural-fit/assigned")
output_path <- option_value(
  "--output-path", "build/candidate-components/assigned/structural-data"
)
fit_version <- option_value("--fit-version", "105")
model <- option_value(
  "--model", "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
)

gq_files <- list.files(
  gq_path, full.names = TRUE,
  pattern = "^main-core-baseline-compact-[1-4][.]csv$"
)
fit_files <- list.files(
  fit_path, full.names = TRUE,
  pattern = paste0("^dist_fit", fit_version, "_", model, "-[1-4][.]csv$")
)
if (length(gq_files) != 4L || length(fit_files) != 4L) {
  stop("Expected four compact-GQ files and four slim fit files.", call. = FALSE)
}

gq <- as_draws_rvars(read_cmdstan_csv(gq_files)$generated_quantities)
required <- c(
  "core_compact_takeup_cf_level",
  "core_compact_belief_prob_1ord",
  "core_compact_belief_prob_2ord",
  "core_compact_prob_prefer_calendar"
)
missing <- setdiff(required, names(gq))
if (length(missing)) {
  stop("Compact GQ lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}
fit <- as_draws_rvars(
  read_cmdstan_csv(fit_files, variables = "hyper_wtp_mu")$post_warmup_draws
)

as_rvar_vector <- function(values) {
  draw_matrix <- do.call(cbind, lapply(values, draws_of))
  rvar(draw_matrix)
}
write_contract <- function(data, values, name) {
  data$value <- as_rvar_vector(values)
  path <- file.path(output_path, name)
  saveRDS(data, path)
  path
}

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
treatments <- c("Control", "Ink", "Calendar", "Bracelet")
distances <- c("close", "far", "combined")
distance_index <- c(close = 2L, far = 3L, combined = 1L)
cf <- gq$core_compact_takeup_cf_level

ate_rows <- list()
ate_values <- list()
add_ate <- function(treatment = NA_character_, dist_group, estimand,
                    mu_treatment = NA_character_, value) {
  index <- length(ate_rows) + 1L
  ate_rows[[index]] <<- data.frame(
    fit_type = "fit", model = model, fit_version = fit_version,
    treatment = treatment, dist_group = dist_group, estimand = estimand,
    mu_treatment = mu_treatment, variable = "pr_takeup",
    stringsAsFactors = FALSE
  )
  ate_values[[index]] <<- value
}

for (dist_group in distances) {
  group <- distance_index[[dist_group]]
  control <- cf[group, 1L, 1L]
  for (treatment in seq_along(treatments)) {
    level <- cf[group, treatment, treatment]
    add_ate(
      treatment = treatments[treatment], dist_group = dist_group,
      estimand = "overall",
      value = if (treatment == 1L) level else level - control
    )
  }
  for (signal in seq_along(treatments)) {
    level <- cf[group, 1L, signal]
    add_ate(
      dist_group = dist_group, estimand = "signal",
      mu_treatment = treatments[signal],
      value = if (signal == 1L) level else level - control
    )
  }
}
for (treatment in seq_along(treatments)) {
  level <- cf[1L, treatment, 1L]
  control <- cf[1L, 1L, 1L]
  add_ate(
    treatment = treatments[treatment], dist_group = "combined",
    estimand = "private",
    value = if (treatment == 1L) level else level - control
  )
}
ates <- bind_rows(ate_rows) |>
  mutate(
    treatment = factor(treatment, levels = treatments),
    mu_treatment = factor(mu_treatment, levels = treatments),
    dist_group = factor(dist_group, levels = c("far", "close", "combined"))
  )
ate_path <- write_contract(
  ates, ate_values,
  paste0(
    "rvar_processed_dist_fit", fit_version, "_ates_", model, "_1-4.rds"
  )
)

level_rows <- list()
level_values <- list()
for (dist_group in distances) {
  group <- distance_index[[dist_group]]
  for (treatment in seq_along(treatments)) {
    index <- length(level_rows) + 1L
    level_rows[[index]] <- data.frame(
      fit_type = "fit", model = model, fit_version = fit_version,
      treatment = treatments[treatment], dist_group = dist_group,
      estimand = "overall", variable = "pr_takeup",
      stringsAsFactors = FALSE
    )
    level_values[[index]] <- cf[group, treatment, treatment]
  }
}
levels <- bind_rows(level_rows) |>
  mutate(
    treatment = factor(treatment, levels = treatments),
    dist_group = factor(dist_group, levels = c("far", "close", "combined"))
  )
level_path <- write_contract(
  levels, level_values,
  paste0(
    "rvar_processed_dist_fit", fit_version, "_levels_", model, "_1-4.rds"
  )
)

belief_rows <- list()
belief_values <- list()
belief_matrices <- list(
  prob_1ord = gq$core_compact_belief_prob_1ord,
  prob_2ord = gq$core_compact_belief_prob_2ord
)
for (dist_group in c("close", "far")) {
  group <- if (dist_group == "close") 1L else 2L
  for (treatment in seq_along(treatments)) {
    for (variable in names(belief_matrices)) {
      index <- length(belief_rows) + 1L
      belief_rows[[index]] <- data.frame(
        dist_treat_idx = (group - 1L) * 4L + treatment,
        model = model, fit_version = fit_version, fit_type = "fit",
        treatment = treatments[treatment], dist_group = dist_group,
        variable = variable, stringsAsFactors = FALSE
      )
      belief_values[[index]] <- belief_matrices[[variable]][group, treatment]
    }
  }
}
for (treatment in seq_along(treatments)) {
  for (variable in names(belief_matrices)) {
    index <- length(belief_rows) + 1L
    belief_rows[[index]] <- data.frame(
      dist_treat_idx = NA_integer_, model = model, fit_version = fit_version,
      fit_type = "fit", treatment = treatments[treatment],
      dist_group = "combined", variable = variable,
      stringsAsFactors = FALSE
    )
    belief_values[[index]] <-
      (belief_matrices[[variable]][1L, treatment] +
       belief_matrices[[variable]][2L, treatment]) / 2
  }
}
belief_probs <- bind_rows(belief_rows) |>
  mutate(treatment = factor(treatment, levels = treatments))
belief_prob_path <- write_contract(
  belief_probs, belief_values,
  paste0(
    "rvar_processed_dist_fit", fit_version, "_belief_probs_", model,
    "_1-4.rds"
  )
)

belief_ate_rows <- list()
belief_ate_values <- list()
for (dist_group in distances) {
  for (treatment in 2:4) {
    for (order in c("1ord", "2ord")) {
      variable <- paste0("prob_", order)
      matrix <- belief_matrices[[variable]]
      value <- if (dist_group == "combined") {
        ((matrix[1L, treatment] - matrix[1L, 1L]) +
         (matrix[2L, treatment] - matrix[2L, 1L])) / 2
      } else {
        group <- if (dist_group == "close") 1L else 2L
        matrix[group, treatment] - matrix[group, 1L]
      }
      index <- length(belief_ate_rows) + 1L
      belief_ate_rows[[index]] <- data.frame(
        model = model, fit_version = fit_version, fit_type = "fit",
        treatment = treatments[treatment], dist_group = dist_group,
        variable = paste0("ate_", order), stringsAsFactors = FALSE
      )
      belief_ate_values[[index]] <- value
    }
  }
}
belief_ates <- bind_rows(belief_ate_rows) |>
  mutate(treatment = factor(treatment, levels = treatments))
belief_ate_path <- write_contract(
  belief_ates, belief_ate_values,
  paste0(
    "rvar_processed_dist_fit", fit_version, "_belief_ates_", model,
    "_1-4.rds"
  )
)

wtp_rows <- list()
wtp_values <- list()
for (index in seq_len(dim(gq$core_compact_prob_prefer_calendar)[1L])) {
  for (variable in c("prob_prefer_calendar", "hyper_wtp_mu")) {
    row <- length(wtp_rows) + 1L
    wtp_rows[[row]] <- data.frame(
      pref_value_diff_idx = index, model = model, fit_version = fit_version,
      fit_type = "fit", variable = variable, stringsAsFactors = FALSE
    )
    wtp_values[[row]] <- if (variable == "prob_prefer_calendar") {
      gq$core_compact_prob_prefer_calendar[index]
    } else {
      fit$hyper_wtp_mu
    }
  }
}
wtp_path <- write_contract(
  bind_rows(wtp_rows), wtp_values,
  paste0(
    "rvar_processed_dist_fit", fit_version, "_wtp_params_", model,
    "_1-4.rds"
  )
)

paths <- c(
  ate_path, level_path, belief_prob_path, belief_ate_path, wtp_path
)
cat(paths, sep = "\n")
