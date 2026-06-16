#!/usr/bin/env Rscript

library(tidyverse)
library(cmdstanr)
library(posterior)

script_options <- docopt::docopt(
  "Usage:
  stan-hmc-diagnostics.R <fit-version> --model=<model> [options]

Options:
  --input-path=<path>       Stan CSV input path [default: data/stan_analysis_data]
  --output-path=<path>      Diagnostics output path [default: temp-data/stan-diagnostics]
  --chains=<chains>         Comma-separated chains [default: 1,2,3,4]
  --fit-prefix=<prefix>     Fit prefix [default: dist_fit]
  --variables=<variables>   Comma-separated Stan variables to retain
",
  args = if (interactive()) "
    105
    --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS
  " else commandArgs(trailingOnly = TRUE)
)

default_variables <- c(
  "beta",
  "u_sd",
  "dist_beta_v",
  "base_mu_rep",
  "centered_cluster_beta_beliefs",
  "centered_cluster_dist_beta_beliefs",
  "hyper_wtp_mu",
  "wtp_value_utility"
)

chains <- str_split(script_options$chains, ",", simplify = TRUE) %>%
  as.character() %>%
  str_trim()

variables <- if (is.null(script_options$variables)) {
  default_variables
} else {
  str_split(script_options$variables, ",", simplify = TRUE) %>%
    as.character() %>%
    str_trim()
}

fit_files <- file.path(
  script_options$input_path,
  str_glue(
    "{script_options$fit_prefix}{script_options$fit_version}_{script_options$model}-{chains}.csv"
  )
)

missing_files <- fit_files[!file.exists(fit_files)]
if (length(missing_files) > 0) {
  stop("Missing Stan CSV(s): ", paste(missing_files, collapse = ", "))
}

fit <- cmdstanr::read_cmdstan_csv(fit_files, variables = variables)

draw_summary <- posterior::summarise_draws(fit$post_warmup_draws) %>%
  filter(!is.na(rhat), !is.na(ess_bulk), !is.na(ess_tail))

sampler_draws <- posterior::as_draws_df(fit$post_warmup_sampler_diagnostics)

diagnostic_summary <- tibble(
  fit_version = script_options$fit_version,
  model = script_options$model,
  chains = paste(chains, collapse = ","),
  variables = paste(variables, collapse = ","),
  n_parameters = nrow(draw_summary),
  max_rhat = max(draw_summary$rhat),
  min_bulk_ess = min(draw_summary$ess_bulk),
  min_tail_ess = min(draw_summary$ess_tail),
  divergences = sum(sampler_draws$divergent__, na.rm = TRUE)
)

dir.create(script_options$output_path, recursive = TRUE, showWarnings = FALSE)

safe_model <- str_replace_all(script_options$model, "[^A-Za-z0-9_-]", "_")
summary_path <- file.path(
  script_options$output_path,
  str_glue("hmc-diagnostics-fit{script_options$fit_version}-{safe_model}.csv")
)
parameter_path <- file.path(
  script_options$output_path,
  str_glue("hmc-diagnostics-fit{script_options$fit_version}-{safe_model}-parameters.csv")
)

readr::write_csv(diagnostic_summary, summary_path)
readr::write_csv(draw_summary, parameter_path)

cat(str_glue(
  "fit {script_options$fit_version} {script_options$model}: ",
  "max split Rhat={round(diagnostic_summary$max_rhat, 3)}, ",
  "min bulk ESS={round(diagnostic_summary$min_bulk_ess)}, ",
  "min tail ESS={round(diagnostic_summary$min_tail_ess)}, ",
  "divergences={diagnostic_summary$divergences}\n"
))
cat("Wrote ", summary_path, "\n", sep = "")
cat("Wrote ", parameter_path, "\n", sep = "")
