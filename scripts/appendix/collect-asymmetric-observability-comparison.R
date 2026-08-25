#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (length(hit) > 1L) stop("Duplicate option: ", name)
  if (length(hit)) substring(hit, nchar(name) + 2L) else default
}

ladder_root <- option_value(
  "--ladder-root",
  "/project/akaring/takeup-data/data/stan_analysis_data/main-core-observability-ladder"
)
prior_root <- option_value(
  "--prior-root",
  "/project/akaring/takeup-data/data/stan_analysis_data/main-core-report-distance-priors"
)
output_path <- option_value(
  "--output-path", "temp-data/asymmetric-observability-comparison"
)
git_commit <- option_value("--git-commit", "unknown")
gq_1250_root <- option_value("--gq-1250-root", "not-yet-submitted")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(purrr)
  library(tidyr)
})

specifications <- tibble::tribble(
  ~id, ~specification, ~fit_dir, ~gq_dir, ~options,
  "tight", "Tight multinomial",
  file.path(prior_root, "tight"), file.path(prior_root, "tight", "gq"),
  "observation_model=1; recognition_structure=0; report_structure=0; report_arm_dist_prior_scale=0.10",
  "f3", "F3 two-stage conditional",
  file.path(ladder_root, "f3"), file.path(ladder_root, "f3", "gq"),
  "observation_model=1; recognition_structure=2; report_structure=2; report_arm_dist_prior_scale=0.25",
  "u3", "U3 two-stage unconditional",
  file.path(ladder_root, "u3"), file.path(ladder_root, "u3", "gq"),
  "observation_model=2; recognition_structure=1; report_structure=2; report_arm_dist_prior_scale=0.25"
)
treatments <- c("Control", "Ink", "Calendar", "Bracelet")
distances <- c(500L, 1500L, 2500L)
truth_states <- c("Nontaker", "Taker")
categories <- c("Yes", "No", "Don't know", "Unrecognized")

csvs <- function(path, fit = FALSE) {
  files <- sort(list.files(path, pattern = "[.]csv$", full.names = TRUE))
  files <- files[!grepl("peer-link-audit|profile|status", basename(files))]
  if (fit) files <- files[grepl("chain[1-4]-1[.]csv$", basename(files))]
  files[file.info(files)$size > 10000]
}
summarize_finite <- function(x) {
  finite <- x[is.finite(x)]
  c(
    median = median(finite), q025 = unname(quantile(finite, 0.025)),
    q975 = unname(quantile(finite, 0.975)),
    probability_below_one = mean(finite < 1),
    probability_above_one = mean(finite > 1),
    draws = length(x), finite_draw_share = length(finite) / length(x)
  )
}

multiplier_draw_rows <- multiplier_curve_rows <- response_rows <- diagnostic_rows <- list()
for (i in seq_len(nrow(specifications))) {
  spec <- specifications[i, ]
  fit_files <- csvs(spec$fit_dir, fit = TRUE)
  gq_files <- csvs(spec$gq_dir)
  if (length(fit_files) != 4L || length(gq_files) != 4L) {
    stop("Expected four fit and four GQ files for ", spec$id)
  }
  fit_csv <- read_cmdstan_csv(fit_files)
  fit_summary <- summarise_draws(fit_csv$post_warmup_draws)
  sampler <- as_draws_df(fit_csv$post_warmup_sampler_diagnostics)
  fit_draws <- as_draws_df(fit_csv$post_warmup_draws)
  diagnostic_rows[[i]] <- data.frame(
    specification = spec$specification,
    max_rhat = max(fit_summary$rhat, na.rm = TRUE),
    min_bulk_ess = min(fit_summary$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(fit_summary$ess_tail, na.rm = TRUE),
    divergences = sum(sampler$divergent__, na.rm = TRUE),
    maximum_treedepth_transitions = sum(sampler$treedepth__ >= 12, na.rm = TRUE),
    retained_draws = nrow(fit_draws)
  )

  gq <- as_draws_df(read_cmdstan_csv(gq_files)$generated_quantities)
  for (arm in seq_along(treatments)) {
    value_1500 <- -gq[[sprintf("core_compact_sm_rescaled[2,%d]", arm)]]
    multiplier_draw_rows[[length(multiplier_draw_rows) + 1L]] <- data.frame(
      specification = spec$specification, treatment = treatments[arm],
      draw = gq$.draw, multiplier = value_1500
    )
    for (distance_index in seq_along(distances)) {
      value <- -gq[[sprintf(
        "core_compact_sm_rescaled[%d,%d]", distance_index, arm
      )]]
      stats <- summarize_finite(value)
      multiplier_curve_rows[[length(multiplier_curve_rows) + 1L]] <- data.frame(
        specification = spec$specification, treatment = treatments[arm],
        distance = distances[distance_index],
        median = stats[["median"]], q025 = stats[["q025"]],
        q975 = stats[["q975"]]
      )
      response_row <- (distance_index - 1L) * length(treatments) + arm
      sensitivity <- gq[[sprintf(
        "core_compact_sensitivity[%d,%d]", distance_index, arm
      )]]
      specificity <- gq[[sprintf(
        "core_compact_specificity[%d,%d]", distance_index, arm
      )]]
      information <- gq[[sprintf(
        "core_compact_information_factor[%d,%d]", distance_index, arm
      )]]
      for (truth in seq_along(truth_states)) {
        probabilities <- vapply(seq_along(categories), function(category) {
          median(gq[[sprintf(
            "core_compact_response_matrix[%d,%d,%d]",
            truth, response_row, category
          )]], na.rm = TRUE)
        }, numeric(1))
        recognition_name <- if (truth == 1L) {
          "core_compact_recognition_nontaker"
        } else {
          "core_compact_recognition_taker"
        }
        recognition <- gq[[sprintf(
          "%s[%d,%d]", recognition_name, distance_index, arm
        )]]
        response_rows[[length(response_rows) + 1L]] <- data.frame(
          specification = spec$specification, treatment = treatments[arm],
          distance = distances[distance_index], true_takeup_state = truth_states[truth],
          p_yes = probabilities[1], p_no = probabilities[2],
          p_dont_know = probabilities[3], p_unrecognized = probabilities[4],
          recognition = median(recognition, na.rm = TRUE),
          sensitivity = median(sensitivity, na.rm = TRUE),
          specificity = median(specificity, na.rm = TRUE),
          information_factor = median(information, na.rm = TRUE),
          check.names = FALSE
        )
      }
    }
  }
}

multiplier_draws <- bind_rows(multiplier_draw_rows)
multiplier_summary <- multiplier_draws |>
  group_by(specification, treatment) |>
  summarise(stats = list(summarize_finite(multiplier)), .groups = "drop") |>
  unnest_wider(stats) |>
  rename(`Pr(M<1)` = probability_below_one, `Pr(M>1)` = probability_above_one)
response_summary <- bind_rows(response_rows)
diagnostics <- bind_rows(diagnostic_rows)
curve <- bind_rows(multiplier_curve_rows)

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(multiplier_summary, file.path(output_path, "multiplier-summary-1500.csv"))
readr::write_csv(multiplier_draws, file.path(output_path, "multiplier-draws-1500.csv"))
readr::write_csv(response_summary, file.path(output_path, "response-channel-summary.csv"))
readr::write_csv(diagnostics, file.path(output_path, "fit-diagnostics.csv"))

pdf(file.path(output_path, "multiplier-comparison.pdf"), width = 11, height = 4.5)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1), las = 1)
colors <- c("#333333", "#0072B2", "#E69F00", "#CC79A7")
offsets <- c(-45, -15, 15, 45)
for (spec_name in specifications$specification) {
  block <- curve[curve$specification == spec_name, ]
  yrange <- range(c(block$q025, block$q975, 1), finite = TRUE)
  plot(NA, xlim = c(350, 2650), ylim = yrange, xlab = "Distance (metres)",
       ylab = "Social multiplier", main = spec_name, xaxt = "n")
  axis(1, at = distances, labels = distances)
  abline(h = 1, lty = 2, col = "grey55")
  for (arm in seq_along(treatments)) {
    arm_data <- block[block$treatment == treatments[arm], ]
    x <- arm_data$distance + offsets[arm]
    segments(x, arm_data$q025, x, arm_data$q975, col = colors[arm], lwd = 2)
    lines(x, arm_data$median, col = colors[arm], lwd = 1.5)
    points(x, arm_data$median, col = colors[arm], pch = 16)
  }
  legend("topright", treatments, col = colors, pch = 16, lty = 1,
         bty = "n", cex = 0.72)
}
dev.off()

manifest <- c(
  "# Asymmetric observability comparison",
  "",
  paste0("- Git commit: `", git_commit, "`"),
  "- Distance definition: `assigned`",
  "- Posterior sampling rerun: no",
  "- 1.5 km GQ rerun: no; existing compact GQs were used",
  paste0("- Pre-emptive exact-1.25 km GQ root: `", gq_1250_root, "`"),
  "",
  "## Specifications and source directories",
  ""
)
for (i in seq_len(nrow(specifications))) {
  manifest <- c(
    manifest,
    paste0("### ", specifications$specification[i]), "",
    paste0("- Options: `", specifications$options[i], "`"),
    paste0("- Fit directory: `", specifications$fit_dir[i], "`"),
    paste0("- Saved 500/1500/2500m GQ directory: `", specifications$gq_dir[i], "`"),
    ""
  )
}
writeLines(manifest, file.path(output_path, "manifest.md"))
