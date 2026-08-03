#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("optim/policy-bootstrap-functions.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(ROI)
  library(ROI.plugin.gurobi)
  library(slam)
})

input_path <- policy_option_value(args, "--input-path")
target_csv <- policy_option_value(args, "--target-csv")
scenario_id <- as.integer(policy_option_value(args, "--scenario-id"))
num_replicates <- as.integer(policy_option_value(args, "--num-replicates", "210"))
time_limit <- as.numeric(policy_option_value(args, "--time-limit", "10000"))
if (is.null(input_path) || is.null(target_csv) || !scenario_id %in% policy_scenarios$scenario_id) {
  stop("--input-path, --target-csv, and a valid --scenario-id are required.", call. = FALSE)
}

scenario <- policy_scenarios[policy_scenarios$scenario_id == scenario_id, ]
edges <- readRDS(file.path(input_path, "policy-feasible-edges.rds"))
curves <- readRDS(file.path(input_path, "policy-demand-curves.rds"))
curves <- curves[curves$scenario_id == scenario_id, ]
draws <- sort(unique(curves$draw))
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
pot_ids <- sort(unique(edges$pot_j))
num_villages <- length(village_ids)
num_pots <- length(pot_ids)
num_edges <- nrow(edges)
village_index <- match(edges$village_i, village_ids)
pot_index <- match(edges$pot_j, pot_ids)

build_sparse_problem <- function(demand) {
  assignment_rows <- village_index
  availability_rows <- num_villages + seq_len(num_edges)
  target_row <- num_villages + num_edges + 1L
  x_columns <- num_pots + seq_len(num_edges)
  constraint_matrix <- sparseMatrix(
    i = c(
      assignment_rows,
      availability_rows, availability_rows,
      rep.int(target_row, num_edges)
    ),
    j = c(
      x_columns,
      x_columns, pot_index,
      x_columns
    ),
    x = c(
      rep(1, num_edges),
      rep(1, num_edges), rep(-1, num_edges),
      demand
    ),
    dims = c(target_row, num_pots + num_edges)
  )
  OP(
    objective = L_objective(c(rep(1, num_pots), rep(0, num_edges))),
    constraints = L_constraint(
      L = as.simple_triplet_matrix(constraint_matrix),
      dir = c(rep("==", num_villages), rep("<=", num_edges), ">="),
      rhs = c(rep(1, num_villages), rep(0, num_edges), target)
    ),
    types = rep("B", num_pots + num_edges),
    maximum = FALSE
  )
}

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
  prediction <- curves[curves$draw == draw, ]
  demand <- prediction$demand[match(edges$distance, prediction$distance)]
  if (any(!is.finite(demand))) stop("Missing edge demand for draw ", draw, call. = FALSE)
  started <- Sys.time()
  best_edge <- unlist(lapply(split(seq_len(num_edges), village_index), function(index) {
    index[which.max(demand[index])]
  }), use.names = FALSE)
  maximum_achievable <- sum(demand[best_edge])
  if (maximum_achievable + 1e-7 < target) {
    if (!scenario$suppress_reputation) {
      stop("Policy target is infeasible for scenario ", scenario$scenario,
           ", draw ", draw, call. = FALSE)
    }
    # This is the paper's explicit no-social-image benchmark: the target cannot
    # be attained even with every village using its best/closest feasible site.
    selected <- best_edge
    status_code <- NA_integer_
    run_status <- "target_infeasible"
  } else {
    fit <- ROI_solve(
      build_sparse_problem(demand), solver = "gurobi",
      control = list(TimeLimit = time_limit, OutputFlag = 0)
    )
    status_code <- fit$status$code
    if (!identical(as.integer(status_code), 0L)) {
      stop("Gurobi failed for scenario ", scenario$scenario, ", draw ", draw,
           " (status ", status_code, ").", call. = FALSE)
    }
    solution <- fit$solution
    selected <- which(solution[num_pots + seq_len(num_edges)] > 0.5)
    run_status <- "complete"
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  allocation <- edges[selected, ]
  allocation$demand <- demand[selected]
  if (nrow(allocation) != num_villages || anyDuplicated(allocation$village_i)) {
    stop("Invalid sparse allocation for draw ", draw, call. = FALSE)
  }
  achieved <- sum(allocation$demand)
  if (run_status == "complete" && achieved + 1e-7 < target) {
    stop("Allocation missed fixed target.", call. = FALSE)
  }
  status <- data.frame(
    draw = draw,
    replicate = prediction$replicate[1L],
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
  write.csv(statuses, status_path, row.names = FALSE)
  message(scenario$scenario, ": ", run_status, " draw ", draw, "/", max(draws))
}

statuses <- statuses[order(statuses$draw), ]
if (nrow(statuses) != length(draws) || anyDuplicated(statuses$draw)) {
  stop("Incomplete scenario status manifest.", call. = FALSE)
}
write.csv(statuses, status_path, row.names = FALSE)
message("Completed ", scenario$scenario, " for ", nrow(statuses), " bootstrap modes.")
