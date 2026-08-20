#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")

output_path <- main_core_option_value(
  args, "--output-path", "temp-data/main-core-prior-grid"
)
baseline_u_mean <- 1.1 / (3.3 - 1)

baseline <- list(
  mu_rep_sd = 1,
  dist_beta_v_sd = 0.25,
  raw_u_sd_alpha = 3.3,
  raw_u_sd_beta = 1.1,
  visibility_prior_multiplier = 1,
  beta_intercept_sd = 1,
  beta_ink_effect_sd = 0.25,
  beta_calendar_effect_sd = 0.25,
  beta_bracelet_effect_sd = 0.25,
  wtp_value_utility_mean = 0,
  wtp_value_utility_sd = 0.0001,
  lnorm_wtp_value_utility_prior = 0L,
  wtp_mu_prior_sd = 2
)
tight <- modifyList(baseline, list(
  mu_rep_sd = 0.5,
  dist_beta_v_sd = 0.125,
  raw_u_sd_alpha = 8,
  raw_u_sd_beta = baseline_u_mean * (8 - 1),
  visibility_prior_multiplier = 0.5,
  beta_intercept_sd = 0.5,
  beta_ink_effect_sd = 0.125,
  beta_calendar_effect_sd = 0.125,
  beta_bracelet_effect_sd = 0.125,
  wtp_value_utility_sd = 0.00005,
  wtp_mu_prior_sd = 1
))
diffuse <- modifyList(baseline, list(
  mu_rep_sd = 2,
  dist_beta_v_sd = 0.5,
  raw_u_sd_alpha = 2.5,
  raw_u_sd_beta = baseline_u_mean * (2.5 - 1),
  visibility_prior_multiplier = 2,
  beta_intercept_sd = 2,
  beta_ink_effect_sd = 0.5,
  beta_calendar_effect_sd = 0.5,
  beta_bracelet_effect_sd = 0.5,
  wtp_value_utility_mean = -10,
  wtp_value_utility_sd = 4,
  lnorm_wtp_value_utility_prior = 1L,
  wtp_mu_prior_sd = 4
))

make_spec <- function(spec_id, label, panel, setting, values) {
  data.frame(
    spec_id = spec_id,
    label = label,
    panel = panel,
    setting = setting,
    as.data.frame(values, stringsAsFactors = FALSE),
    seed = 20260820L + 100L * spec_id,
    sample_status = "pending",
    gq_status = "pending",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
vary <- function(base, alternative, fields) {
  base[fields] <- alternative[fields]
  base
}
private_fields <- c(
  "beta_intercept_sd", "beta_ink_effect_sd", "beta_calendar_effect_sd",
  "beta_bracelet_effect_sd", "wtp_value_utility_mean",
  "wtp_value_utility_sd", "lnorm_wtp_value_utility_prior",
  "wtp_mu_prior_sd"
)

specifications <- list(
  make_spec(1, "baseline", "Baseline", "Baseline", baseline),
  make_spec(2, "image-tight", "Social-image weight", "Tight",
            vary(baseline, tight, "mu_rep_sd")),
  make_spec(3, "image-diffuse", "Social-image weight", "Diffuse",
            vary(baseline, diffuse, "mu_rep_sd")),
  make_spec(4, "distance-tight", "Distance cost", "Tight",
            vary(baseline, tight, "dist_beta_v_sd")),
  make_spec(5, "distance-diffuse", "Distance cost", "Diffuse",
            vary(baseline, diffuse, "dist_beta_v_sd")),
  make_spec(6, "heterogeneity-tight", "Idiosyncratic heterogeneity", "Tight",
            vary(baseline, tight, c("raw_u_sd_alpha", "raw_u_sd_beta"))),
  make_spec(7, "heterogeneity-diffuse", "Idiosyncratic heterogeneity", "Diffuse",
            vary(baseline, diffuse, c("raw_u_sd_alpha", "raw_u_sd_beta"))),
  make_spec(8, "visibility-tight", "Visibility schedule", "Tight",
            vary(baseline, tight, "visibility_prior_multiplier")),
  make_spec(9, "visibility-diffuse", "Visibility schedule", "Diffuse",
            vary(baseline, diffuse, "visibility_prior_multiplier")),
  make_spec(10, "private-tight", "Private utility", "Tight",
            vary(baseline, tight, private_fields)),
  make_spec(11, "private-diffuse", "Private utility", "Diffuse",
            vary(baseline, diffuse, private_fields)),
  make_spec(12, "joint-tight", "Joint stress test", "Tight", tight),
  make_spec(13, "joint-diffuse", "Joint stress test", "Diffuse", diffuse)
)
manifest <- do.call(rbind, specifications)

if (nrow(manifest) != 13L || anyDuplicated(manifest$label) ||
    any(!is.finite(as.matrix(manifest[, names(baseline)])))) {
  stop("Invalid prior-grid manifest.", call. = FALSE)
}
u_mean <- manifest$raw_u_sd_beta / (manifest$raw_u_sd_alpha - 1)
if (max(abs(u_mean - baseline_u_mean)) > 1e-12) {
  stop("Heterogeneity priors do not preserve the baseline mean.", call. = FALSE)
}

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(output_path, "prior-grid-manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
message("Wrote 13-specification prior grid: ", manifest_path)
print(manifest[, c("spec_id", "label", "panel", "setting")], row.names = FALSE)
