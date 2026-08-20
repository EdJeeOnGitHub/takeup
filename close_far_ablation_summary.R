#!/usr/bin/env Rscript

args <- docopt::docopt(
  "
Summarize matched continuous-distance and Close/Far-only structural fits.

Usage:
  close_far_ablation_summary.R --continuous=<prefix> --close-far=<prefix> [options]

Options:
  --model-continuous=<model>  Continuous-distance model name [default: STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP]
  --model-close-far=<model>   Close/Far-only model name [default: STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CLOSE_FAR_ONLY]
  --input-path=<path>         Fit directory [default: data/stan_analysis_data]
  --output-path=<path>        Output directory [default: temp-data/bayesian-sensitivity]
  --chains=<chains>           Comma-separated chains [default: 1,2,3,4]
  ",
  args = commandArgs(trailingOnly = TRUE)
)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(readr)
})

chains <- as.integer(strsplit(args$chains, ",", fixed = TRUE)[[1]])
continuous_files <- file.path(
  args$input_path,
  sprintf(
    "%s_%s-%d.csv",
    args$continuous,
    args$model_continuous,
    chains
  )
)
close_far_files <- file.path(
  args$input_path,
  sprintf(
    "%s_%s-%d.csv",
    args$close_far,
    args$model_close_far,
    chains
  )
)
missing_files <- c(
  continuous_files[!file.exists(continuous_files)],
  close_far_files[!file.exists(close_far_files)]
)
if (length(missing_files) > 0) {
  stop("Missing fit files:\n", paste(missing_files, collapse = "\n"))
}

read_fit_data <- function(files) {
  parsed_by_chain <- lapply(files, function(path) {
    read_cmdstan_csv(
      path,
      variables = c("u_sd[1]", "base_mu_rep", "dist_beta_v[1]")
    )
  })
  draws_by_chain <- lapply(parsed_by_chain, `[[`, "post_warmup_draws")
  chain_draw_counts <- vapply(
    draws_by_chain,
    dim,
    FUN.VALUE = integer(3)
  )[1, ]
  list(
    draws = posterior::bind_draws(draws_by_chain, along = "chain"),
    chain_draw_counts = chain_draw_counts,
    divergences = sum(vapply(
      parsed_by_chain,
      function(parsed) {
        sum(
          as_draws_matrix(parsed$post_warmup_sampler_diagnostics)[, "divergent__"]
        )
      },
      numeric(1)
    ))
  )
}

continuous_fit <- read_fit_data(continuous_files)
close_far_fit <- read_fit_data(close_far_files)
continuous_draws <- continuous_fit$draws
close_far_draws <- close_far_fit$draws

if (any(continuous_fit$chain_draw_counts != 400L) ||
    any(close_far_fit$chain_draw_counts != 400L)) {
  stop(
    "Expected 400 retained draws per chain; continuous counts are ",
    paste(continuous_fit$chain_draw_counts, collapse = ", "),
    " and Close/Far counts are ",
    paste(close_far_fit$chain_draw_counts, collapse = ", ")
  )
}

summarize_specification <- function(draws, specification) {
  q025 <- function(x) unname(quantile(x, 0.025))
  q975 <- function(x) unname(quantile(x, 0.975))
  summarise_draws(
    draws,
    mean,
    median,
    sd,
    q025,
    q975,
    rhat,
    ess_bulk,
    ess_tail
  ) %>%
    mutate(specification = specification, .before = 1)
}

summary <- bind_rows(
  summarize_specification(continuous_draws, "Continuous distance"),
  summarize_specification(close_far_draws, "Close/Far only")
)

sigma_summary <- summary %>% filter(variable == "u_sd[1]")
continuous_sigma <- sigma_summary %>%
  filter(specification == "Continuous distance")
close_far_sigma <- sigma_summary %>%
  filter(specification == "Close/Far only")

comparison <- tibble(
  continuous_mean = continuous_sigma$mean,
  close_far_mean = close_far_sigma$mean,
  mean_difference = close_far_sigma$mean - continuous_sigma$mean,
  continuous_sd = continuous_sigma$sd,
  close_far_sd = close_far_sigma$sd,
  sd_ratio_close_far_to_continuous = close_far_sigma$sd / continuous_sigma$sd,
  variance_ratio_close_far_to_continuous =
    (close_far_sigma$sd / continuous_sigma$sd)^2
)

diagnostics <- bind_rows(
  tibble(
    specification = "Continuous distance",
    divergences = continuous_fit$divergences,
    draws_per_chain = paste(
      continuous_fit$chain_draw_counts,
      collapse = ";"
    )
  ),
  tibble(
    specification = "Close/Far only",
    divergences = close_far_fit$divergences,
    draws_per_chain = paste(
      close_far_fit$chain_draw_counts,
      collapse = ";"
    )
  )
)

dir.create(args$output_path, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(args$output_path, "close-far-distance-ablation")
write_csv(summary, paste0(stem, "-posterior-summary.csv"))
write_csv(comparison, paste0(stem, "-sigma-comparison.csv"))
write_csv(diagnostics, paste0(stem, "-diagnostics.csv"))

sigma_markdown <- c(
  "| Specification | Mean | SD | 95% interval | Rhat | Bulk ESS |",
  "|:--|--:|--:|:--|--:|--:|",
  apply(sigma_summary, 1, function(row) {
    sprintf(
      "| %s | %.3f | %.3f | [%.3f, %.3f] | %.3f | %.0f |",
      row[["specification"]],
      as.numeric(row[["mean"]]),
      as.numeric(row[["sd"]]),
      as.numeric(row[["q025"]]),
      as.numeric(row[["q975"]]),
      as.numeric(row[["rhat"]]),
      as.numeric(row[["ess_bulk"]])
    )
  }),
  "",
  sprintf(
    "Close/Far-to-continuous posterior-SD ratio: **%.3f**.",
    comparison$sd_ratio_close_far_to_continuous
  )
)
writeLines(sigma_markdown, paste0(stem, "-sigma-comparison.md"))

message("Wrote Close/Far ablation summaries to ", args$output_path)
