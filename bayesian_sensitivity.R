#!/usr/bin/env Rscript

args <- docopt::docopt(
  "
Calculate Bayesian data-to-parameter sensitivities.

Usage:
  bayesian_sensitivity.R <fit-version> --model=<model> [options]

Options:
  --fit-path=<path>       Directory containing fitted CmdStan CSV files [default: data/stan_analysis_data]
  --gq-path=<path>        Directory containing sensitivity GQ files [default: temp-data/bayesian-sensitivity]
  --output-path=<path>    Directory for tables and data [default: temp-data/bayesian-sensitivity]
  --chains=<chains>       Comma-separated chain numbers [default: 1,2,3,4]
  --digits=<n>            Digits in formatted tables [default: 2]
  --derived               Include social multipliers and social-image returns
  --num-clusters=<n>      Clusters in the fitted model [default: 144]
  ",
  args = commandArgs(trailingOnly = TRUE)
)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(readr)
  library(tidyr)
})

fit_version <- args$fit_version
model_name <- args$model
chains <- as.integer(strsplit(args$chains, ",", fixed = TRUE)[[1]])
digits <- as.integer(args$digits)
num_clusters <- as.integer(args$num_clusters)

fit_files <- file.path(
  args$fit_path,
  sprintf("dist_fit%s_%s-%d.csv", fit_version, model_name, chains)
)
gq_files <- file.path(
  args$gq_path,
  sprintf("dist_fit%s_%s_sensitivity-%d.csv", fit_version, model_name, chains)
)

missing_files <- c(fit_files[!file.exists(fit_files)], gq_files[!file.exists(gq_files)])
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

parameter_map <- c(
  lambda = "base_mu_rep",
  beta_control = "beta[1]",
  beta_ink = "beta[2]",
  beta_calendar = "beta[3]",
  beta_bracelet = "beta[4]",
  delta = "dist_beta_v[1]",
  sigma_u = "u_sd[1]",
  visibility_level_control = "centered_cluster_beta_beliefs[1,1]",
  visibility_level_ink = "centered_cluster_beta_beliefs[1,2]",
  visibility_level_calendar = "centered_cluster_beta_beliefs[1,3]",
  visibility_level_bracelet = "centered_cluster_beta_beliefs[1,4]",
  visibility_distance_control = "centered_cluster_dist_beta_beliefs[1,1]",
  visibility_distance_ink = "centered_cluster_dist_beta_beliefs[1,2]",
  visibility_distance_calendar = "centered_cluster_dist_beta_beliefs[1,3]",
  visibility_distance_bracelet = "centered_cluster_dist_beta_beliefs[1,4]",
  psi = "hyper_wtp_mu"
)

likelihood_map <- c(
  takeup = "sensitivity_log_lik_takeup",
  beliefs = "sensitivity_log_lik_beliefs",
  wtp = "sensitivity_log_lik_wtp",
  distance = "sensitivity_log_lik_distance"
)

parameter_labels <- c(
  lambda = "\\lambda",
  beta_control = "\\beta_{Control}",
  beta_ink = "\\beta_{Ink}",
  beta_calendar = "\\beta_{Calendar}",
  beta_bracelet = "\\beta_{Bracelet}",
  delta = "\\delta",
  sigma_u = "\\sigma_u",
  visibility_level_control = "\\beta^O_{Control}",
  visibility_level_ink = "\\beta^O_{Ink}",
  visibility_level_calendar = "\\beta^O_{Calendar}",
  visibility_level_bracelet = "\\beta^O_{Bracelet}",
  visibility_distance_control = "\\gamma^O_{Control}",
  visibility_distance_ink = "\\gamma^O_{Ink}",
  visibility_distance_calendar = "\\gamma^O_{Calendar}",
  visibility_distance_bracelet = "\\gamma^O_{Bracelet}",
  psi = "\\psi"
)

parameter_groups <- c(
  rep("Takeup Parameters", 7),
  rep("Visibility Parameters", 8),
  "WTP Parameter"
)

derived_spec <- tidyr::crossing(
  measure = c("social_multiplier", "social_image_return"),
  distance_m = c(500L, 1500L, 2500L),
  treatment_id = seq_len(4)
) %>%
  mutate(
    distance_index = distance_m %/% 100L + 1L,
    treatment = c("Control", "Ink", "Calendar", "Bracelet")[treatment_id],
    stan_base = if_else(
      measure == "social_multiplier",
      "cluster_social_multiplier",
      "cluster_rep_return"
    ),
    quantity = paste(measure, tolower(treatment), distance_m, sep = "_"),
    quantity_label = if_else(
      measure == "social_multiplier",
      sprintf("M_{%s}(%0.1f\\,\\mathrm{km})", treatment, distance_m / 1000),
      sprintf("SI_{%s}(%0.1f\\,\\mathrm{km})", treatment, distance_m / 1000)
    ),
    quantity_group = if_else(
      measure == "social_multiplier",
      "Social Multiplier",
      "Social-Image Return"
    )
  )

derived_variables <- if (args$derived) {
  unlist(lapply(seq_len(nrow(derived_spec)), function(row_index) {
    row <- derived_spec[row_index, ]
    sprintf(
      "%s[%d,%d,%d]",
      row$stan_base,
      row$distance_index,
      seq_len(num_clusters),
      row$treatment_id
    )
  }))
} else {
  character()
}

parameter_draws_by_chain <- vector("list", length(fit_files))
likelihood_draws_by_chain <- vector("list", length(gq_files))
derived_draws_by_chain <- if (args$derived) vector("list", length(fit_files)) else NULL

# The structural CSVs are large. Process chains sequentially and retain only
# the requested columns and cluster averages.
for (chain_index in seq_along(fit_files)) {
  requested_fit_variables <- c(unname(parameter_map), derived_variables)
  fit_draws <- read_cmdstan_csv(
    fit_files[chain_index],
    variables = requested_fit_variables
  )$post_warmup_draws %>%
    as_draws_matrix()

  parameter_draws_by_chain[[chain_index]] <-
    fit_draws[, unname(parameter_map), drop = FALSE]

  if (args$derived) {
    chain_derived <- matrix(
      NA_real_,
      nrow = nrow(fit_draws),
      ncol = nrow(derived_spec),
      dimnames = list(NULL, derived_spec$quantity)
    )
    for (row_index in seq_len(nrow(derived_spec))) {
      row <- derived_spec[row_index, ]
      cluster_variables <- sprintf(
        "%s[%d,%d,%d]",
        row$stan_base,
        row$distance_index,
        seq_len(num_clusters),
        row$treatment_id
      )
      cluster_average <-
        rowMeans(fit_draws[, cluster_variables, drop = FALSE])
      # Match the paper convention in structural_tables.R:
      # reported multiplier = -cluster_social_multiplier / dist_beta_v.
      chain_derived[, row_index] <- if (row$measure == "social_multiplier") {
        -cluster_average / fit_draws[, "dist_beta_v[1]"]
      } else {
        cluster_average
      }
    }
    derived_draws_by_chain[[chain_index]] <- chain_derived
  }

  rm(fit_draws)
  likelihood_draws_by_chain[[chain_index]] <- read_cmdstan_csv(
    gq_files[chain_index],
    variables = unname(likelihood_map)
  )$generated_quantities %>%
    as_draws_matrix()
  gc(verbose = FALSE)
}

parameter_draws <- do.call(rbind, parameter_draws_by_chain)
likelihood_draws <- do.call(rbind, likelihood_draws_by_chain)
derived_draws <- if (args$derived) do.call(rbind, derived_draws_by_chain) else NULL

if (nrow(parameter_draws) != nrow(likelihood_draws) ||
    (args$derived && nrow(derived_draws) != nrow(likelihood_draws))) {
  stop("Parameter, derived-quantity, and likelihood draws are not aligned.")
}

colnames(parameter_draws) <- names(parameter_map)
colnames(likelihood_draws) <- names(likelihood_map)

calculate_sensitivity <- function(draws, labels, groups) {
  posterior_mean <- colMeans(draws)
  posterior_sd <- apply(draws, 2, sd)
  posterior_variance <- posterior_sd^2

  invalid <- !is.finite(posterior_sd) | posterior_sd <= 0
  if (any(invalid)) {
    stop(
      "Nonpositive or nonfinite posterior SD for: ",
      paste(names(posterior_sd)[invalid], collapse = ", ")
    )
  }

  mean_sensitivity <- cov(draws, likelihood_draws)
  mean_sensitivity <- sweep(mean_sensitivity, 1, posterior_sd, "/")

  centered_squared <- sweep(draws, 2, posterior_mean, "-")^2
  log_sd_sensitivity <- cov(centered_squared, likelihood_draws)
  log_sd_sensitivity <- sweep(
    log_sd_sensitivity,
    1,
    2 * posterior_variance,
    "/"
  )

  mean_sensitivity <- as.data.frame(mean_sensitivity)
  names(mean_sensitivity) <- paste0("mean_", names(likelihood_map))
  log_sd_sensitivity <- as.data.frame(log_sd_sensitivity)
  names(log_sd_sensitivity) <- paste0("log_sd_", names(likelihood_map))

  bind_cols(
    tibble(
      quantity = colnames(draws),
      quantity_label = unname(labels[colnames(draws)]),
      quantity_group = unname(groups)
    ),
    tibble(
      posterior_mean = unname(posterior_mean),
      posterior_sd = unname(posterior_sd)
    ),
    mean_sensitivity,
    log_sd_sensitivity
  )
}

parameter_sensitivity <- calculate_sensitivity(
  parameter_draws,
  parameter_labels,
  parameter_groups
)

derived_sensitivity <- if (args$derived) {
  derived_labels <- setNames(
    derived_spec$quantity_label,
    derived_spec$quantity
  )
  calculate_sensitivity(
    derived_draws,
    derived_labels,
    derived_spec$quantity_group
  )
} else {
  NULL
}

dir.create(args$output_path, recursive = TRUE, showWarnings = FALSE)
output_stem <- sprintf(
  "bayesian-parameter-sensitivity-fit%s-%s",
  fit_version,
  model_name
)

format_number <- function(x) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

write_sensitivity_table <- function(data, stem, type = c("mean", "log_sd")) {
  type <- match.arg(type)
  columns <- paste0(type, "_", names(likelihood_map))
  table_data <- data %>%
    transmute(
      quantity_label,
      quantity_group,
      takeup = .data[[columns[1]]],
      beliefs = .data[[columns[2]]],
      wtp = .data[[columns[3]]],
      distance = .data[[columns[4]]]
    )

  table_rows <- table_data %>%
    mutate(
      across(c(takeup, beliefs, wtp, distance), format_number),
      row = paste0(
        "$", quantity_label, "$ & ",
        takeup, " & ", beliefs, " & ", wtp, " & ", distance, " \\\\"
      )
    )

  latex_lines <- c(
    "\\begin{tabular}{lrrrr}",
    "\\toprule",
    "Quantity & Take-up & Beliefs & WTP & Distance \\\\",
    "\\midrule"
  )
  for (group_name in unique(table_rows$quantity_group)) {
    latex_lines <- c(
      latex_lines,
      "\\addlinespace[0.3em]",
      sprintf("\\multicolumn{5}{l}{\\textbf{%s}} \\\\", group_name),
      table_rows$row[table_rows$quantity_group == group_name]
    )
  }
  latex_lines <- c(latex_lines, "\\bottomrule", "\\end{tabular}")
  writeLines(latex_lines, paste0(stem, ".tex"))

  markdown <- table_data %>%
    transmute(
      Quantity = paste0("$", quantity_label, "$"),
      `Take-up` = format_number(takeup),
      Beliefs = format_number(beliefs),
      WTP = format_number(wtp),
      Distance = format_number(distance)
    )
  markdown_lines <- c(
    "| Quantity | Take-up | Beliefs | WTP | Distance |",
    "|:--|--:|--:|--:|--:|",
    apply(markdown, 1, function(row) {
      paste0("| ", paste(row, collapse = " | "), " |")
    })
  )
  writeLines(markdown_lines, paste0(stem, ".md"))
}

write_csv(
  parameter_sensitivity,
  file.path(args$output_path, paste0(output_stem, ".csv"))
)
write_sensitivity_table(
  parameter_sensitivity,
  file.path(args$output_path, output_stem),
  "mean"
)
write_sensitivity_table(
  parameter_sensitivity,
  file.path(args$output_path, paste0(output_stem, "-uncertainty")),
  "log_sd"
)

if (args$derived) {
  derived_stem <- sprintf(
    "bayesian-derived-sensitivity-fit%s-%s",
    fit_version,
    model_name
  )
  write_csv(
    derived_sensitivity,
    file.path(args$output_path, paste0(derived_stem, ".csv"))
  )
  write_sensitivity_table(
    derived_sensitivity,
    file.path(args$output_path, derived_stem),
    "mean"
  )
  write_sensitivity_table(
    derived_sensitivity,
    file.path(args$output_path, paste0(derived_stem, "-uncertainty")),
    "log_sd"
  )
}

message(
  "Wrote posterior-mean and posterior-uncertainty sensitivity tables to ",
  args$output_path
)
