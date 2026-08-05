#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("scratch/main-core-data.R")

workspace_path <- main_core_option_value(
  args, "--workspace", "data/stan_analysis_data/dist_fit105.RData"
)
output_path <- main_core_option_value(
  args, "--output-path", "temp-data/main-core-cluster-weights"
)
replicates <- as.integer(main_core_option_value(args, "--replicates", "200"))
start_replicate <- as.integer(
  main_core_option_value(args, "--start-replicate", "1")
)
manifest_name <- main_core_option_value(
  args, "--manifest-name", "weight-manifest.csv"
)
seed <- as.integer(main_core_option_value(args, "--seed", "20260802"))
methods <- strsplit(
  main_core_option_value(args, "--methods", "multinomial"),
  ",", fixed = TRUE
)[[1L]]
if (length(methods) == 0L ||
    any(!methods %in% c("multinomial", "exponential"))) {
  stop("--methods must contain multinomial and/or exponential.", call. = FALSE)
}
if (replicates < 1L) stop("--replicates must be positive.", call. = FALSE)
if (start_replicate < 1L) {
  stop("--start-replicate must be positive.", call. = FALSE)
}

fit_env <- new.env(parent = emptyenv())
load(workspace_path, envir = fit_env)
cluster_county_id <- as.integer(fit_env$stan_data$cluster_county_id)
num_clusters <- length(cluster_county_id)
county_clusters <- split(seq_len(num_clusters), cluster_county_id)

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)
manifest <- vector("list", length(methods) * replicates)
manifest_index <- 1L

for (method in methods) {
  method_path <- file.path(output_path, method)
  dir.create(method_path, recursive = TRUE, showWarnings = FALSE)
  for (replicate_id in seq.int(
      start_replicate, length.out = replicates
  )) {
    weight <- numeric(num_clusters)
    for (cluster_ids in county_clusters) {
      if (method == "multinomial") {
        sampled <- sample(cluster_ids, length(cluster_ids), replace = TRUE)
        weight[cluster_ids] <- tabulate(
          match(sampled, cluster_ids), nbins = length(cluster_ids)
        )
      } else {
        raw_weight <- rexp(length(cluster_ids))
        weight[cluster_ids] <-
          length(cluster_ids) * raw_weight / sum(raw_weight)
      }
    }
    weight_file <- file.path(
      method_path, sprintf("weights-%04d.csv", replicate_id)
    )
    write.csv(
      data.frame(cluster_id = seq_len(num_clusters), weight = weight),
      weight_file,
      row.names = FALSE
    )
    manifest[[manifest_index]] <- data.frame(
      method = method,
      replicate = replicate_id,
      seed = seed,
      weight_file = normalizePath(weight_file),
      min_weight = min(weight),
      max_weight = max(weight),
      zero_clusters = sum(weight == 0),
      stringsAsFactors = FALSE
    )
    manifest_index <- manifest_index + 1L
  }
}
manifest <- do.call(rbind, manifest)
manifest$weight_hash <- unname(tools::md5sum(manifest$weight_file))
duplicate_hash <- duplicated(manifest[c("method", "weight_hash")]) |
  duplicated(manifest[c("method", "weight_hash")], fromLast = TRUE)
if (any(duplicate_hash)) {
  stop("Generated duplicate cluster-bootstrap weight vectors.", call. = FALSE)
}
write.csv(manifest, file.path(output_path, manifest_name), row.names = FALSE)
message("Wrote ", nrow(manifest), " stratified cluster-weight files to ", output_path)
