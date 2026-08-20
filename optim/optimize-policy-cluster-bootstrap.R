#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("optim/policy-bootstrap-functions.R")
source("optim/policy-cost-sensitivity-functions.R")

input_path <- policy_option_value(args, "--input-path")
target_csv <- policy_option_value(args, "--target-csv")
scenario_id <- as.integer(policy_option_value(args, "--scenario-id"))
num_replicates <- as.integer(policy_option_value(args, "--num-replicates", "999"))
time_limit <- as.numeric(policy_option_value(args, "--time-limit", "10000"))
target_tolerance <- as.numeric(policy_option_value(args, "--target-tolerance", "1e-5"))
solver <- policy_option_value(args, "--solver", "auto")
if (is.null(input_path) || is.null(target_csv) || !scenario_id %in% policy_scenarios$scenario_id) {
  stop("--input-path, --target-csv, and a valid --scenario-id are required.", call. = FALSE)
}

scenario <- policy_scenarios[policy_scenarios$scenario_id == scenario_id, ]
edges <- readRDS(file.path(input_path, "policy-feasible-edges.rds"))
edge_demand_path <- file.path(input_path, "policy-edge-demand-matrix.rds")
if (file.exists(edge_demand_path)) {
  edge_demand <- readRDS(edge_demand_path)
  draw_map <- read.csv(
    file.path(input_path, "policy-edge-demand-draw-map.csv"),
    stringsAsFactors = FALSE
  )
  if (nrow(edge_demand) != nrow(draw_map)) {
    stop("Compact edge-demand matrix and draw map disagree.", call. = FALSE)
  }
  curves <- NULL
  draws <- draw_map$draw
} else {
  edge_demand <- NULL
  curves <- readRDS(file.path(input_path, "policy-demand-curves.rds"))
  curves <- curves[curves$scenario_id == scenario_id, ]
  draws <- sort(unique(curves$draw))
}
draws <- draws[seq_len(min(num_replicates, length(draws)))]

target_data <- read.csv(target_csv, stringsAsFactors = FALSE)
if (!"social_welfare" %in% names(target_data)) {
  stop("Target CSV lacks social_welfare.", call. = FALSE)
}
# The legacy file repeats the draw-level target on every village row and has
# one isolated missing cell. Recover one finite target per draw, then take the
# same across-draw mean intended by the existing optimizer.
if ("draw" %in% names(target_data)) {
  target_by_draw <- vapply(split(target_data$social_welfare, target_data$draw), function(value) {
    finite <- unique(value[is.finite(value)])
    if (length(finite) != 1L) stop("Ambiguous target within a posterior draw.", call. = FALSE)
    finite
  }, numeric(1))
  target <- mean(target_by_draw)
} else {
  target <- mean(target_data$social_welfare, na.rm = TRUE)
}
if (!is.finite(target)) stop("Non-finite policy target.", call. = FALSE)

village_ids <- sort(unique(edges$village_i))
num_villages <- length(village_ids)
num_edges <- nrow(edges)
village_index <- match(edges$village_i, village_ids)

scenario_dir <- file.path(input_path, "allocations", scenario$scenario)
dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(scenario_dir, "status.csv")
statuses <- data.frame()

for (draw in draws) {
  output_file <- file.path(scenario_dir, sprintf("replicate-%04d.rds", draw))
  if (file.exists(output_file)) {
    saved <- readRDS(output_file)
    statuses <- rbind(statuses, saved$status)
    next
  }
  if (!is.null(edge_demand)) {
    draw_row <- match(draw, draw_map$draw)
    scenario_columns <- (scenario_id - 1L) * nrow(edges) + seq_len(nrow(edges))
    demand <- edge_demand[draw_row, scenario_columns]
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
  if (any(!is.finite(demand))) {
    status <- data.frame(
      draw = draw, replicate = replicate_value, scenario_id = scenario_id,
      scenario = scenario$scenario, scenario_label = scenario$label,
      status = "equilibrium_undefined", solver_status = NA_integer_,
      elapsed_seconds = 0, n_pot = NA_real_, mean_demand = NA_real_,
      mean_distance = NA_real_, achieved_welfare = NA_real_,
      target_welfare = target, stringsAsFactors = FALSE
    )
    saveRDS(
      list(status = status, allocation = edges[FALSE, ]), output_file,
      compress = FALSE
    )
    statuses <- rbind(statuses, status)
    if (draw %% 25L == 0L || draw == max(draws)) {
      write.csv(statuses, status_path, row.names = FALSE)
      message(scenario$scenario, ": processed draw ", draw, "/", max(draws))
    }
    next
  }
  started <- Sys.time()
  best_edge <- unlist(lapply(split(seq_len(num_edges), village_index), function(index) {
    index[which.max(demand[index])]
  }), use.names = FALSE)
  maximum_achievable <- sum(demand[best_edge])
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
      population = rep(1, num_villages),
      target_rate = target / num_villages,
      site_cost = 1, solver = solver, time_limit = time_limit,
      work_path = scenario_dir
    )
    status_code <- 0L
    allocation <- fit$allocation
    run_status <- "complete"
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (run_status == "target_infeasible") {
    allocation <- edges[selected, ]
    allocation$demand <- demand[selected]
  }
  if (nrow(allocation) != num_villages || anyDuplicated(allocation$village_i)) {
    stop("Invalid sparse allocation for draw ", draw, call. = FALSE)
  }
  achieved <- sum(allocation$demand)
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
    mean_demand = mean(allocation$demand),
    mean_distance = mean(allocation$distance),
    achieved_welfare = achieved,
    target_welfare = target,
    stringsAsFactors = FALSE
  )
  saveRDS(list(status = status, allocation = allocation), output_file, compress = FALSE)
  statuses <- rbind(statuses, status)
  if (draw %% 25L == 0L || draw == max(draws)) {
    write.csv(statuses, status_path, row.names = FALSE)
    message(scenario$scenario, ": processed draw ", draw, "/", max(draws))
  }
}

statuses <- statuses[order(statuses$draw), ]
if (nrow(statuses) != length(draws) || anyDuplicated(statuses$draw)) {
  stop("Incomplete scenario status manifest.", call. = FALSE)
}
write.csv(statuses, status_path, row.names = FALSE)
message("Completed ", scenario$scenario, " for ", nrow(statuses), " draws.")
