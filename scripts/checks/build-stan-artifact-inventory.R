#!/usr/bin/env Rscript

# Non-destructive inventory of Stan sources and saved fit artifacts. This script
# deliberately does not move, compress, or delete anything.

roots <- c(
  "stan_models",
  "data/stan_analysis_data",
  "temp-data/archive",
  "optim/data/archive"
)
roots <- roots[dir.exists(roots)]

files <- unlist(lapply(
  roots,
  list.files,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
), use.names = FALSE)
files <- files[file.info(files)$isdir %in% FALSE]

main_tracked <- system2("git", "ls-files", stdout = TRUE)
data_tracked <- if (dir.exists("data/.git") || file.exists("data/.git")) {
  system2("git", c("-C", "data", "ls-files"), stdout = TRUE)
} else {
  character()
}

tracked <- vapply(files, function(path) {
  if (startsWith(path, "data/")) {
    sub("^data/", "", path) %in% data_tracked
  } else {
    path %in% main_tracked
  }
}, logical(1))

code_roots <- c("R", "scripts", "scratch", "optim", "presentations")
code_roots <- code_roots[dir.exists(code_roots)]
root_sources <- c("README.md", "Makefile", "_targets.R")
root_sources <- root_sources[file.exists(root_sources)]
reference_tokens <- system2(
  "rg",
  c(
    "-o", "--no-filename",
    "--glob", shQuote("*.R"), "--glob", shQuote("*.Rmd"),
    "--glob", shQuote("*.qmd"), "--glob", shQuote("*.sh"),
    "--glob", shQuote("*.md"), "--glob", shQuote("*.yml"),
    "--glob", shQuote("*.yaml"),
    shQuote("[A-Za-z0-9_.-]+\\.(stan|RData|rds|zip|csv)"),
    code_roots, root_sources
  ),
  stdout = TRUE,
  stderr = FALSE
)
reference_tokens <- unique(reference_tokens)

stan_dir <- "stan_models"
collect_stan_dependencies <- function(entrypoints) {
  queue <- entrypoints
  found <- character()
  while (length(queue)) {
    path <- queue[[1L]]
    queue <- queue[-1L]
    if (path %in% found || !file.exists(path)) next
    found <- c(found, path)
    lines <- readLines(path, warn = FALSE)
    includes <- sub(
      ".*#include[[:space:]]+([^[:space:]]+).*", "\\1",
      grep("#include[[:space:]]+", lines, value = TRUE)
    )
    includes <- file.path(stan_dir, basename(includes))
    queue <- c(queue, includes)
  }
  unique(found)
}

canonical_stan <- collect_stan_dependencies(c(
  file.path(stan_dir, "takeup_struct_main_core.stan"),
  file.path(stan_dir, "takeup_struct_main_core_compact_gq.stan")
))
canonical_chains <- paste0(
  "data/stan_analysis_data/",
  "dist_fit104_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-", 1:4, ".csv"
)
canonical_workspace <- "build/structural-workspace/main-core-input.RData"

extension <- tools::file_ext(files)
kind <- ifelse(
  extension == "stan", "stan_source",
  ifelse(extension == "hpp", "generated_cpp",
  ifelse(extension == "csv" & grepl("-[1-4]\\.csv$", files), "stan_chain_csv",
  ifelse(extension == "RData", "fit_workspace",
  ifelse(extension == "rds", "processed_fit",
  ifelse(extension == "zip", "compressed_archive",
  ifelse(extension == "" & startsWith(files, "stan_models/"),
         "compiled_executable", "other")))))))

referenced <- basename(files) %in% reference_tokens

role <- rep("historical_or_unclassified", length(files))
role[files %in% canonical_stan] <- "canonical_stan_source"
role[files == canonical_workspace] <- "canonical_workspace"
role[files %in% canonical_chains] <- "canonical_saved_chain"
role[kind %in% c("generated_cpp", "compiled_executable")] <-
  "generated_build_product"
role[startsWith(files, "temp-data/archive/") |
     startsWith(files, "optim/data/archive/")] <- "explicit_archive"
role[role == "historical_or_unclassified" & referenced] <-
  "referenced_support_or_robustness"

disposition <- rep("retain pending workflow review", length(files))
disposition[role %in% c(
  "canonical_stan_source", "canonical_workspace", "canonical_saved_chain",
  "referenced_support_or_robustness"
)] <- "retain"
disposition[role == "generated_build_product"] <-
  "regenerable; safe cleanup candidate after review"
disposition[role == "explicit_archive"] <-
  "retain in place; no cleanup authorized"

inventory <- data.frame(
  path = files,
  bytes = unname(file.info(files)$size),
  size_mib = round(unname(file.info(files)$size) / 1024^2, 3),
  kind = kind,
  git_tracked = tracked,
  referenced_by_text_source = referenced,
  workflow_role = role,
  recommended_disposition = disposition,
  stringsAsFactors = FALSE
)
inventory <- inventory[order(-inventory$bytes, inventory$path), ]

dir.create("build/manifest", recursive = TRUE, showWarnings = FALSE)
output <- "build/manifest/stan-artifacts.csv"
write.csv(inventory, output, row.names = FALSE)

summary <- aggregate(
  inventory$bytes,
  list(workflow_role = inventory$workflow_role),
  sum
)
summary$size_gib <- round(summary$x / 1024^3, 3)
summary$x <- NULL
write.csv(summary, "build/manifest/stan-artifacts-summary.csv", row.names = FALSE)

cat("Wrote", nrow(inventory), "artifacts to", output, "\n")
print(summary, row.names = FALSE)
