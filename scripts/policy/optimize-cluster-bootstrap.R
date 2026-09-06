#!/usr/bin/env Rscript

process_started <- proc.time()[[3L]]
args <- commandArgs(trailingOnly = TRUE)
source("R/policy/bootstrap.R")
source("R/policy/cost-sensitivity.R")

input_path <- policy_option_value(args, "--input-path")
target_csv <- policy_option_value(args, "--target-csv")
scenario_id <- as.integer(policy_option_value(args, "--scenario-id"))
num_replicates <- as.integer(policy_option_value(args, "--num-replicates", "999"))
time_limit <- as.numeric(policy_option_value(args, "--time-limit", "10000"))
target_tolerance <- as.numeric(policy_option_value(args, "--target-tolerance", "1e-5"))
solver <- policy_option_value(args, "--solver", "auto")
population_weighting <- policy_option_value(args, "--population-weighting", "equal-community")
target_mode <- policy_option_value(args, "--target-mode", "legacy-fixed")
distance_data <- policy_option_value(args, "--distance-data")
census_data <- policy_option_value(args, "--census-data", Sys.getenv("POLICY_CENSUS", "data/takeup_census.RData"))
allocation_root <- policy_option_value(args, "--allocation-root", input_path)
draw_start <- as.integer(policy_option_value(args, "--draw-start", "1"))
draw_end <- as.integer(policy_option_value(args, "--draw-end", as.character(num_replicates)))
if (!population_weighting %in% c("adult-census", "equal-community") ||
    !target_mode %in% c("legacy-fixed", "draw-specific-experimental-control")) stop("Unknown population/target convention.")
if (population_weighting == "adult-census" && target_mode != "draw-specific-experimental-control") stop("Adult weighting requires explicit draw-specific experimental targets.")
if (anyNA(c(draw_start, draw_end)) || draw_start < 1L || draw_end < draw_start) stop("Invalid draw shard.")
allocated <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", as.character(parallel::detectCores())))
if (is.na(allocated) || allocated < 1L) allocated <- 1L
num_cores <- as.integer(policy_option_value(args, "--num-cores", as.character(min(8L, allocated))))
batch_size <- as.integer(policy_option_value(args, "--draw-batch-size", "25"))
thread_option <- policy_option_value(args, "--solver-threads", "1")
solver_threads <- if (thread_option == "auto") NULL else as.integer(thread_option)
seed_option <- policy_option_value(args, "--solver-seed", "0")
solver_seed <- if (seed_option == "auto") NULL else as.integer(seed_option)
cache_format <- policy_option_value(args, "--cache-format", "auto")
scratch_path <- policy_option_value(args, "--scratch-path", Sys.getenv("SLURM_TMPDIR", tempdir()))
if (anyNA(c(num_cores, batch_size)) || num_cores < 1L || batch_size < 1L ||
    (!is.null(solver_threads) && (is.na(solver_threads) || solver_threads < 1L)) ||
    (!is.null(solver_seed) && (is.na(solver_seed) || solver_seed < 0L)) ||
    !cache_format %in% c("auto", "matrix", "legacy")) stop("Invalid execution controls.")
if (num_cores > 1L && (is.null(solver_threads) || is.null(solver_seed))) stop("Parallel draws require explicit solver threads and seed.")
if (!is.null(solver_threads) && num_cores * solver_threads > allocated) stop("Workers exceed allocated CPUs.")
if (num_cores > 1L) Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (solver == "auto") solver <- if (nzchar(Sys.which("glpsol"))) "glpk" else "gurobi"
if (!solver %in% c("glpk", "gurobi")) stop("Unknown solver.")
solver_executable <- Sys.which(if (solver == "gurobi") "gurobi_cl" else "glpsol")
if (!nzchar(solver_executable)) stop("Selected solver is not on PATH.")
if (anyNA(c(time_limit, target_tolerance)) || time_limit <= 0 || target_tolerance < 0) stop("Invalid solver limits.")
if (is.na(num_replicates) || num_replicates < 1L) stop("Invalid replicate count.")
atomic_rds <- function(value, path) {
  tmp <- tempfile(pattern = ".pending-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  saveRDS(value, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}
atomic_csv <- function(value, path) {
  tmp <- tempfile(pattern = ".pending-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  write.csv(value, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}
if (is.null(input_path) || is.null(target_csv) || !scenario_id %in% policy_scenarios$scenario_id) {
  stop("--input-path, --target-csv, and a valid --scenario-id are required.", call. = FALSE)
}

scenario <- policy_scenarios[policy_scenarios$scenario_id == scenario_id, ]
edges <- readRDS(file.path(input_path, "policy-feasible-edges.rds"))
edge_demand_path <- file.path(input_path, "policy-edge-demand-matrix.rds")
legacy_path <- file.path(input_path, "policy-demand-curves.rds")
manifest_path <- file.path(input_path, "policy-cache-manifest.rds")
if (cache_format == "auto") {
  if (file.exists(manifest_path)) cache_format <- readRDS(manifest_path)$format else {
    if (file.exists(edge_demand_path) && file.exists(legacy_path)) stop("Ambiguous caches: select --cache-format explicitly.")
    cache_format <- if (file.exists(edge_demand_path)) "matrix" else "legacy"
  }
}
if (cache_format == "matrix" && file.exists(manifest_path)) {
  m <- readRDS(manifest_path)
  actual <- tools::md5sum(file.path(input_path, names(m$hashes)))
  if (!identical(unname(actual), unname(m$hashes))) stop("Prediction cache checksum mismatch.")
}
if (cache_format == "matrix") {
  edge_demand <- readRDS(edge_demand_path)
  draw_map <- read.csv(
    file.path(input_path, "policy-edge-demand-draw-map.csv"),
    stringsAsFactors = FALSE
  )
  if (nrow(edge_demand) != nrow(draw_map) || anyDuplicated(draw_map$draw) ||
      anyNA(draw_map$draw) || ncol(edge_demand) != nrow(edges) * nrow(policy_scenarios)) {
    stop("Compact edge-demand matrix and draw map disagree.", call. = FALSE)
  }
  curves <- NULL
  draws <- draw_map$draw
  scenario_columns <- (scenario_id - 1L) * nrow(edges) + seq_len(nrow(edges))
  edge_demand <- edge_demand[, scenario_columns, drop = FALSE]
} else {
  edge_demand <- NULL
  curves <- readRDS(file.path(input_path, "policy-demand-curves.rds"))
  curves <- curves[curves$scenario_id == scenario_id, ]
  draws <- sort(unique(curves$draw))
}
if (!length(draws) || anyNA(draws) || anyDuplicated(draws)) stop("Invalid draw inventory.")
draws <- draws[seq_len(min(num_replicates, length(draws)))]
if (draw_start > length(draws)) stop("Draw shard starts beyond input draws.")
draws <- draws[seq.int(draw_start, min(draw_end, length(draws)))]

target_data <- read.csv(target_csv, stringsAsFactors = FALSE)
population_table <- if (population_weighting == "adult-census") {
  if (is.null(distance_data)) stop("Adult weighting requires --distance-data.")
  policy_adult_population(readRDS(distance_data), census_data)
} else NULL
population <- if (is.null(population_table)) rep(1, length(unique(edges$village_i))) else population_table$population
if (target_mode == "draw-specific-experimental-control") {
  required <- c("draw", "replicate", "target_expected_adults", "target_community_welfare", "population_total", "target_mode", "population_checksum")
  if (!all(required %in% names(target_data)) || anyDuplicated(target_data$draw) ||
      any(!target_data$target_mode %in% target_mode)) stop("Invalid draw-specific target manifest.")
  target_index <- match(draws, target_data$draw)
  if (anyNA(target_index)) stop("Missing targets for selected draws.")
  if (population_weighting == "adult-census" &&
      (any(target_data$population_total != sum(population)) ||
       any(target_data$population_checksum != policy_object_hash(population_table)))) stop("Target and census population disagree.")
  targets <- target_data[[if (population_weighting == "adult-census") "target_expected_adults" else "target_community_welfare"]][target_index]
  if (!is.null(edge_demand) && any(draw_map$replicate[match(draws, draw_map$draw)] != target_data$replicate[target_index])) stop("Target replicate IDs disagree with predictions.")
} else {
  if (!"social_welfare" %in% names(target_data)) stop("Target CSV lacks social_welfare.")
  target_by_draw <- if ("draw" %in% names(target_data)) vapply(split(target_data$social_welfare, target_data$draw), function(value) {
    finite <- unique(value[is.finite(value)])
    if (length(finite) != 1L) stop("Ambiguous target within a posterior draw.")
    finite
  }, numeric(1)) else target_data$social_welfare
  target <- mean(target_by_draw, na.rm = TRUE)
  if (!is.finite(target)) stop("Non-finite legacy target.")
  targets <- rep(target, length(draws))
}
names(targets) <- as.character(draws)

village_ids <- sort(unique(edges$village_i))
num_villages <- length(village_ids)
num_edges <- nrow(edges)
village_index <- match(edges$village_i, village_ids)
prepared <- policy_cost_prepare(edges)
if (!identical(as.integer(village_ids), seq_len(num_villages))) stop("Village IDs must be contiguous.")
input_files <- c(file.path(input_path, "policy-feasible-edges.rds"), target_csv,
                 if (cache_format == "matrix") c(edge_demand_path, file.path(input_path, "policy-edge-demand-draw-map.csv")) else legacy_path)
solver_version <- system2(solver_executable, "--version", stdout = TRUE, stderr = TRUE,
  env = if (solver == "gurobi") paste0("LD_LIBRARY_PATH=", file.path(dirname(dirname(solver_executable)), "lib")) else character())
contract <- list(version = 1L, hashes = unname(tools::md5sum(input_files)),
                 scenario_id = scenario_id, targets = targets, population = population, population_weighting = population_weighting, target_mode = target_mode,
                 solver = solver, solver_version = solver_version, solver_executable_hash = unname(tools::md5sum(solver_executable)),
                 threads = solver_threads, seed = solver_seed, time_limit = time_limit,
                 target_tolerance = target_tolerance,
                 code_hashes = unname(tools::md5sum(c("scripts/policy/optimize-cluster-bootstrap.R", "R/policy/cost-sensitivity.R", "R/policy/population.R"))))
contract_file <- tempfile(); saveRDS(contract, contract_file)
contract_hash <- unname(tools::md5sum(contract_file)); unlink(contract_file)
load_seconds <- proc.time()[[3L]] - process_started

scenario_dir <- file.path(allocation_root, "allocations", scenario$scenario)
dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(scenario_dir, "status.csv")
run_manifest <- file.path(scenario_dir, "run-manifest.rds")
if (file.exists(run_manifest) && !identical(readRDS(run_manifest)$contract, contract)) stop("Existing run has different inputs/settings; use a fresh output root.")
atomic_rds(list(contract = contract, contract_hash = contract_hash, R = R.version.string,
                draws = draws, workers = num_cores, batch_size = batch_size,
                scratch_path = scratch_path, solver_version = solver_version,
                input_files = normalizePath(input_files),
                environment = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "SLURM_JOB_ID", "SLURM_CPUS_PER_TASK"))), run_manifest)

solve_one_draw <- function(draw) {
  target <- targets[as.character(draw)]
  output_file <- file.path(scenario_dir, sprintf("replicate-%04d.rds", draw))
  if (file.exists(output_file)) {
    saved <- readRDS(output_file)
    if (!identical(saved$contract_hash, contract_hash) || saved$status$draw != draw ||
        saved$status$status %in% c("failed", "solver_incomplete")) stop("Invalid resumable draw: ", draw)
    return(saved$status)
  }
  lookup_started <- proc.time()[[3L]]
  if (!is.null(edge_demand)) {
    draw_row <- match(draw, draw_map$draw)
    demand <- edge_demand[draw_row, ]
    replicate_value <- draw_map$replicate[draw_row]
  } else {
    prediction <- curves[curves$draw == draw, ]
    edge_specific <- "village_i" %in% names(prediction) &&
      all(!is.na(prediction$village_i))
    if (edge_specific) {
      if (nrow(prediction) != nrow(edges) ||
          !identical(as.integer(prediction$village_i), as.integer(edges$village_i)) ||
          !isTRUE(all.equal(prediction$distance, edges$distance, tolerance = 0))) {
        stop("Edge-specific predictions do not align with feasible edges.", call. = FALSE)
      }
      demand <- prediction$demand
    } else {
      demand <- prediction$demand[match(edges$distance, prediction$distance)]
    }
    replicate_value <- prediction$replicate[1L]
  }
  lookup_seconds <- proc.time()[[3L]] - lookup_started
  if (any(!is.finite(demand)) || !is.finite(target)) {
    status <- data.frame(
      draw = draw, replicate = replicate_value, scenario_id = scenario_id,
      scenario = scenario$scenario, scenario_label = scenario$label,
      status = "equilibrium_undefined", solver_status = NA_integer_,
      elapsed_seconds = 0, n_pot = NA_real_, mean_demand = NA_real_,
      mean_distance = NA_real_, achieved_welfare = NA_real_,
      target_welfare = target, stringsAsFactors = FALSE
    )
    output_started <- proc.time()[[3L]]
    atomic_rds(list(status = status, allocation = edges[FALSE, ], contract_hash = contract_hash), output_file)
    attr(status, "execution_timing") <- c(lookup_seconds = lookup_seconds,
      output_write_seconds = proc.time()[[3L]] - output_started)
    return(status)
  }

  started <- Sys.time()
  best_edge <- unlist(lapply(prepared$village_edges, function(index) {
    index[which.max(demand[index])]
  }), use.names = FALSE)
  maximum_achievable <- sum(population[edges$village_i[best_edge]] * demand[best_edge])
  fit <- NULL
  if (maximum_achievable + target_tolerance < target) {
    # Retain the best feasible allocation and record target infeasibility. This
    # is expected for the no-social-image benchmark and can also occur in tail
    # posterior draws of flexible alternative structural models.
    selected <- best_edge
    status_code <- NA_integer_
    run_status <- "target_infeasible"
  } else {
    fit <- policy_cost_solve(
      edges = edges, demand = demand,
      population = population,
      target_rate = target / sum(population),
      site_cost = 1, solver = solver, time_limit = time_limit,
      work_path = file.path(scratch_path, paste0("policy-worker-", Sys.getpid())),
      prepared = prepared, solver_threads = solver_threads, solver_seed = solver_seed,
      log_file = file.path(scenario_dir, "solver-logs", sprintf("replicate-%04d.log", draw))
    )
    status_code <- 0L
    allocation <- fit$allocation
    run_status <- if (isTRUE(fit$diagnostics$optimal)) "complete" else "solver_incomplete"
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (run_status == "target_infeasible") {
    allocation <- edges[selected, ]
    allocation$demand <- demand[selected]
  }
  if (nrow(allocation) != num_villages || anyDuplicated(allocation$village_i)) {
    stop("Invalid sparse allocation for draw ", draw, call. = FALSE)
  }
  allocation$population <- population[allocation$village_i]
  allocation$expected_takers <- allocation$population * allocation$demand
  achieved <- sum(allocation$expected_takers)
  if (run_status == "complete" && achieved + target_tolerance < target) {
    stop(
      "Allocation missed fixed target by ", target - achieved, ".",
      call. = FALSE
    )
  }
  status <- data.frame(
    draw = draw,
    replicate = replicate_value,
    scenario_id = scenario_id,
    scenario = scenario$scenario,
    scenario_label = scenario$label,
    status = run_status,
    solver_status = status_code,
    elapsed_seconds = elapsed,
    n_pot = length(unique(allocation$pot_j)),
    mean_demand = weighted.mean(allocation$demand, allocation$population),
    mean_distance = weighted.mean(allocation$distance, allocation$population),
    achieved_welfare = achieved,
    target_welfare = target,
    stringsAsFactors = FALSE
  )
  output_started <- proc.time()[[3L]]
  atomic_rds(list(status = status, allocation = allocation, contract_hash = contract_hash,
                 solver = if (!is.null(fit)) fit$diagnostics else NULL,
                 timing = c(lookup_seconds = lookup_seconds, if (!is.null(fit)) fit$timing)), output_file)
  attr(status, "execution_timing") <- c(lookup_seconds = lookup_seconds,
    if (!is.null(fit)) fit$timing,
    output_write_seconds = proc.time()[[3L]] - output_started)
  if (run_status == "solver_incomplete") stop("Solver did not prove optimality for draw ", draw)
  status
}

failed_draw <- function(draw, message) {
  target <- targets[as.character(draw)]
  data.frame(draw = draw, replicate = NA_integer_, scenario_id = scenario_id,
             scenario = scenario$scenario, scenario_label = scenario$label,
             status = "failed", solver_status = NA_integer_, elapsed_seconds = NA_real_,
             n_pot = NA_real_, mean_demand = NA_real_, mean_distance = NA_real_,
             achieved_welfare = NA_real_, target_welfare = target,
             error = message, stringsAsFactors = FALSE)
}
safe_draw <- function(draw) {
  value <- tryCatch(solve_one_draw(draw), error = function(e) failed_draw(draw, conditionMessage(e)))
  value$population_weighting <- population_weighting
  value$population_total <- sum(population)
  value$target_mode <- target_mode
  value
}
# Batch tasks bound process creation and keep large inputs read-only after fork.
batches <- split(draws, ceiling(seq_along(draws) / batch_size))
results <- list()
for (first in seq.int(1L, length(batches), by = max(1L, num_cores * 4L))) {
  wave <- batches[seq.int(first, min(length(batches), first + num_cores * 4L - 1L))]
  computed <- parallel::mclapply(wave, function(batch) lapply(batch, safe_draw),
                                mc.cores = num_cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
  for (i in seq_along(computed)) {
    if (inherits(computed[[i]], "try-error") || is.null(computed[[i]])) {
      computed[[i]] <- lapply(wave[[i]], function(draw) {
        value <- failed_draw(draw, "Worker terminated unexpectedly; inspect saved assignments before resuming.")
        value$population_weighting <- population_weighting
        value$population_total <- sum(population)
        value$target_mode <- target_mode
        value
      })
    }
  }
  results <- c(results, unlist(computed, recursive = FALSE))
  statuses <- do.call(rbind, lapply(results, function(x) {
    if (!"error" %in% names(x)) x$error <- NA_character_
    x
  }))
  statuses <- statuses[order(statuses$draw), ]
  atomic_csv(statuses, status_path)
  message(scenario$scenario, ": processed ", nrow(statuses), "/", length(draws), " draws")
}
if (nrow(statuses) != length(draws) || anyDuplicated(statuses$draw) || !setequal(statuses$draw, draws)) stop("Incomplete scenario status manifest.")
stage_totals <- vapply(c("lookup_seconds", "model_write_seconds", "solver_seconds", "output_write_seconds"), function(stage) {
  sum(vapply(results, function(x) { value <- attr(x, "execution_timing")[stage];
    if (length(value) && is.finite(value)) value else 0 }, numeric(1)))
}, numeric(1))
atomic_csv(data.frame(workers = num_cores, batch_size = batch_size,
                      as.list(stage_totals),
                      load_seconds = load_seconds,
                      elapsed_seconds = proc.time()[[3L]] - process_started),
           file.path(scenario_dir, "execution-timing.csv"))
if (any(statuses$status == "failed")) stop("Failed draws: ", paste(statuses$draw[statuses$status == "failed"], collapse = ", "))
message("Completed ", scenario$scenario, " for ", nrow(statuses), " draws.")
