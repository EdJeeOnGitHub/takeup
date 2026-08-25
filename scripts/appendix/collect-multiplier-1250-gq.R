#!/usr/bin/env Rscript

# Collect exact-1250m compact generated quantities. This script never samples.
args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (length(hit) > 1L) stop("Duplicate option: ", name)
  if (length(hit)) substring(hit, nchar(name) + 2L) else default
}

input_root <- value("--input-root")
output_csv <- value("--output-csv", file.path(input_root, "multiplier-draws-1250.csv"))
if (is.null(input_root) || !dir.exists(input_root)) stop("--input-root must exist.")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(purrr)
  library(tidyr)
})

catalog <- tibble::tribble(
  ~table_id, ~spec_id, ~directory, ~schema,
  "main", "baseline", "baseline", "core",
  "prior", "baseline", "baseline", "core",
  "main", "private-distance-community-image", "private-distance-community-image", "legacy",
  "main", "full-information", "full-information", "legacy",
  "main", "exclude-dispersed", "exclude-dispersed", "legacy",
  "main", "correct-classification", "correct-classification", "legacy",
  "main", "second-order-observability", "second-order-observability", "legacy",
  "main", "grouped-lambda", "grouped-lambda", "core",
  "main", "arm-specific-lambda", "arm-specific-lambda", "core",
  "main", "student-t5", "student-t5", "core",
  "main", "finite-mixture", "finite-mixture", "core",
  "prior", "image-tight", "image-tight", "core",
  "prior", "image-diffuse", "image-diffuse", "core",
  "prior", "distance-tight", "distance-tight", "core",
  "prior", "distance-diffuse", "distance-diffuse", "core",
  "prior", "heterogeneity-tight", "heterogeneity-tight", "core",
  "prior", "heterogeneity-diffuse", "heterogeneity-diffuse", "core",
  "prior", "visibility-tight", "visibility-tight", "core",
  "prior", "visibility-diffuse", "visibility-diffuse", "core",
  "prior", "private-tight", "private-tight", "core",
  "prior", "private-diffuse", "private-diffuse", "core",
  "prior", "joint-tight", "joint-tight", "core",
  "prior", "joint-diffuse", "joint-diffuse", "core"
)

treatment <- c("Control", "Ink", "Calendar", "Bracelet")
read_spec <- function(table_id, spec_id, directory, schema) {
  files <- list.files(file.path(input_root, directory), pattern = "\\.csv$", full.names = TRUE)
  if (length(files) != 4L) stop("Expected four GQ CSVs in ", directory, "; found ", length(files))
  prefix <- if (schema == "core") "core_compact_sm_rescaled" else "compact_sm_rescaled"
  variables <- sprintf("%s[1,%d]", prefix, 1:4)
  draws <- cmdstanr::as_cmdstan_fit(sort(files))$draws(variables, format = "draws_df")
  as.data.frame(draws) |>
    transmute(draw = .draw, across(all_of(variables))) |>
    pivot_longer(-draw, names_to = "variable", values_to = "raw_multiplier") |>
    mutate(
      table_id = table_id,
      spec_id = spec_id,
      treatment = treatment[match(variable, variables)],
      distance_m = 1250,
      multiplier = -raw_multiplier
    ) |>
    select(table_id, spec_id, treatment, draw, distance_m, multiplier)
}

result <- pmap_dfr(catalog, read_spec) |>
  mutate(treatment = factor(treatment, levels = c("Control", "Calendar", "Ink", "Bracelet"))) |>
  arrange(table_id, spec_id, treatment, draw) |>
  mutate(treatment = as.character(treatment))

expected_rows <- 23L * 4L * length(unique(result$draw[result$table_id == "main" & result$spec_id == "baseline"]))
if (nrow(result) != expected_rows) {
  warning("Draw counts vary across specifications; retained all available draws.")
}
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(result, output_csv)
message("Wrote ", nrow(result), " rows to ", output_csv)
