# Helpers for the local, review-only policy cost and pooling sensitivity.
#
# This module deliberately does not depend on ROI or an R Gurobi package.  It
# writes a standard LP file and calls gurobi_cl, which makes the exercise
# runnable on the local workstation while Midway is unavailable.

policy_cost_weighted_quantile <- function(x, weights, probability) {
  stopifnot(length(x) == length(weights), length(probability) == 1L)
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  if (!any(keep)) return(NA_real_)
  x <- x[keep]
  weights <- weights[keep]
  ordering <- order(x)
  x <- x[ordering]
  weights <- weights[ordering]
  x[which(cumsum(weights) / sum(weights) >= probability)[1L]]
}

policy_cost_lp_terms <- function(coefficients, variables, tolerance = 1e-12) {
  keep <- is.finite(coefficients) & abs(coefficients) > tolerance
  coefficients <- coefficients[keep]
  variables <- variables[keep]
  if (!length(coefficients)) return("0")
  pieces <- paste0(
    ifelse(coefficients >= 0, "+ ", "- "),
    format(abs(coefficients), scientific = TRUE, digits = 15, trim = TRUE),
    " ", variables
  )
  paste(pieces, collapse = " ")
}

policy_cost_write_lp <- function(
    path, edges, objective_x, objective_y, welfare_x, target,
    pooled_welfare_change = NULL, pooled_objective_change = NULL) {
  village_ids <- sort(unique(edges$village_i))
  pot_ids <- sort(unique(edges$pot_j))
  x_names <- paste0("x_", seq_len(nrow(edges)))
  y_names <- paste0("y_", pot_ids)
  pooling <- !is.null(pooled_welfare_change)
  if (pooling) {
    p_names <- paste0("p_", seq_len(nrow(edges)))
    z_names <- paste0("z_", pot_ids)
  } else {
    p_names <- z_names <- character()
  }

  connection <- file(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write_lp_line <- function(value) {
    writeLines(strwrap(value, width = 180, exdent = 2), connection)
  }
  writeLines("Minimize", connection)
  objective <- policy_cost_lp_terms(objective_y, y_names)
  objective_x_terms <- policy_cost_lp_terms(objective_x, x_names)
  if (objective_x_terms != "0") objective <- paste(objective, objective_x_terms)
  if (pooling) {
    pooled_objective_terms <- policy_cost_lp_terms(pooled_objective_change, p_names)
    if (pooled_objective_terms != "0") {
      objective <- paste(objective, pooled_objective_terms)
    }
  }
  write_lp_line(paste(" obj:", objective))
  writeLines("Subject To", connection)

  for (village in village_ids) {
    index <- which(edges$village_i == village)
    write_lp_line(paste0(
      " assign_", village, ": ",
      policy_cost_lp_terms(rep(1, length(index)), x_names[index]), " = 1"
    ))
  }
  for (edge in seq_len(nrow(edges))) {
    write_lp_line(paste0(
      " open_", edge, ": + 1 ", x_names[edge], " - 1 y_",
      edges$pot_j[edge], " <= 0"
    ))
  }

  welfare_terms <- policy_cost_lp_terms(welfare_x, x_names)
  if (pooling) {
    pooled_welfare_terms <- policy_cost_lp_terms(pooled_welfare_change, p_names)
    if (pooled_welfare_terms != "0") {
      welfare_terms <- paste(welfare_terms, pooled_welfare_terms)
    }
  }
  write_lp_line(paste0(" target: ", welfare_terms, " >= ",
                       format(target, scientific = TRUE, digits = 15)))

  if (pooling) {
    big_m <- length(village_ids) - 1L
    for (pot in pot_ids) {
      index <- which(edges$pot_j == pot)
      # z_j equals one exactly when two or more villages use site j.
      write_lp_line(paste0(
        " pool_upper_", pot, ": ",
        policy_cost_lp_terms(rep(1, length(index)), x_names[index]),
        " - ", big_m, " z_", pot, " <= 1"
      ))
      write_lp_line(paste0(
        " pool_lower_", pot, ": ",
        policy_cost_lp_terms(rep(1, length(index)), x_names[index]),
        " - 2 z_", pot, " >= 0"
      ))
    }
    for (edge in seq_len(nrow(edges))) {
      pot <- edges$pot_j[edge]
      write_lp_line(paste0(
        " product_x_", edge, ": + 1 p_", edge,
        " - 1 x_", edge, " <= 0"
      ))
      write_lp_line(paste0(
        " product_z_", edge, ": + 1 p_", edge,
        " - 1 z_", pot, " <= 0"
      ))
      write_lp_line(paste0(
        " product_lower_", edge, ": + 1 p_", edge,
        " - 1 x_", edge, " - 1 z_", pot, " >= -1"
      ))
    }
  }

  writeLines("Binary", connection)
  binary_names <- c(x_names, y_names, p_names, z_names)
  for (start in seq(1L, length(binary_names), by = 10L)) {
    end <- min(start + 9L, length(binary_names))
    writeLines(paste(" ", paste(binary_names[start:end], collapse = " ")), connection)
  }
  writeLines("End", connection)
  invisible(list(x = x_names, y = y_names, p = p_names, z = z_names))
}

policy_cost_read_solution <- function(path) {
  if (!file.exists(path)) stop("Gurobi did not write a solution: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  pieces <- strsplit(trimws(lines), "[[:space:]]+")
  values <- vapply(pieces, function(value) as.numeric(value[2L]), numeric(1))
  names(values) <- vapply(pieces, `[`, character(1), 1L)
  values
}

policy_cost_read_glpk_solution <- function(solution_path, problem_path) {
  problem_lines <- strsplit(trimws(readLines(problem_path, warn = FALSE)), "[[:space:]]+")
  column_lines <- problem_lines[vapply(
    problem_lines,
    function(value) length(value) >= 4L && value[1L] == "n" && value[2L] == "j",
    logical(1)
  )]
  column_index <- vapply(column_lines, function(value) as.integer(value[3L]), integer(1))
  column_name <- vapply(column_lines, function(value) value[4L], character(1))
  names(column_name) <- as.character(column_index)

  solution_lines <- strsplit(trimws(readLines(solution_path, warn = FALSE)), "[[:space:]]+")
  status_line <- solution_lines[vapply(
    solution_lines,
    function(value) length(value) >= 6L && value[1L] == "s" && value[2L] == "mip",
    logical(1)
  )]
  if (length(status_line) != 1L || !status_line[[1L]][5L] %in% c("o", "f")) {
    stop("GLPK did not return a feasible integer solution.", call. = FALSE)
  }
  value_lines <- solution_lines[vapply(
    solution_lines,
    function(value) length(value) >= 3L && value[1L] == "j",
    logical(1)
  )]
  index <- vapply(value_lines, function(value) as.integer(value[2L]), integer(1))
  values <- vapply(value_lines, function(value) as.numeric(value[3L]), numeric(1))
  names(values) <- unname(column_name[as.character(index)])
  values
}

policy_cost_solve <- function(
    edges, demand, population, target_rate, site_cost,
    signal_cost_per_taker = 0, travel_cost_per_roundtrip_km = 0,
    pooled_demand = NULL,
    gurobi = if (nzchar(Sys.which("gurobi_cl"))) {
      Sys.which("gurobi_cl")
    } else "/opt/gurobi1000/linux64/bin/gurobi_cl",
    gurobi_library = file.path(dirname(dirname(gurobi)), "lib"),
    time_limit = 300, work_path = tempdir(), solver = "auto",
    glpsol = if (nzchar(Sys.which("glpsol"))) {
      Sys.which("glpsol")
    } else "/usr/bin/glpsol") {
  stopifnot(nrow(edges) == length(demand), all(edges$village_i %in% seq_along(population)))
  pooling <- !is.null(pooled_demand)
  if (pooling && length(pooled_demand) != nrow(edges)) {
    stop("Pooled demand does not align with edges.", call. = FALSE)
  }
  variable_cost_rate <- signal_cost_per_taker +
    travel_cost_per_roundtrip_km * 2 * edges$distance_km
  population_edge <- population[edges$village_i]
  welfare_x <- population_edge * demand
  objective_x <- welfare_x * variable_cost_rate
  target <- target_rate * sum(population)
  pooled_welfare_change <- pooled_objective_change <- NULL
  if (pooling) {
    pooled_welfare_change <- population_edge * (pooled_demand - demand)
    pooled_objective_change <- pooled_welfare_change * variable_cost_rate
  }
  pot_ids <- sort(unique(edges$pot_j))
  objective_y <- rep(site_cost, length(pot_ids))

  dir.create(work_path, recursive = TRUE, showWarnings = FALSE)
  identifier <- paste0(Sys.getpid(), "-", format(Sys.time(), "%H%M%OS6"))
  identifier <- gsub("[^0-9-]", "", identifier)
  lp_path <- file.path(work_path, paste0("policy-", identifier, ".lp"))
  solution_path <- file.path(work_path, paste0("policy-", identifier, ".sol"))
  log_path <- file.path(work_path, paste0("policy-", identifier, ".log"))
  problem_path <- file.path(work_path, paste0("policy-", identifier, ".glp"))
  names_map <- policy_cost_write_lp(
    lp_path, edges, objective_x, objective_y, welfare_x, target,
    pooled_welfare_change, pooled_objective_change
  )
  on.exit(unlink(c(lp_path, solution_path, log_path, problem_path)), add = TRUE)
  if (solver == "auto") solver <- if (file.exists(glpsol)) "glpk" else "gurobi"
  if (!solver %in% c("glpk", "gurobi")) stop("Unknown MILP solver: ", solver)
  if (solver == "glpk") {
    output <- system2(
      glpsol,
      args = c(
        "--lp", lp_path, "--tmlim", as.integer(time_limit),
        "--write", solution_path, "--wglp", problem_path, "--log", log_path
      ),
      stdout = TRUE, stderr = TRUE
    )
  } else {
    output <- system2(
      gurobi,
      args = c(
        paste0("ResultFile=", solution_path), paste0("LogFile=", log_path),
        "OutputFlag=0", paste0("TimeLimit=", time_limit), lp_path
      ),
      stdout = TRUE, stderr = TRUE,
      env = paste0("LD_LIBRARY_PATH=", gurobi_library)
    )
  }
  status <- attr(output, "status") %||% 0L
  if (status != 0L) {
    stop(solver, " failed (", status, "): ", paste(output, collapse = "\n"), call. = FALSE)
  }
  solution <- if (solver == "glpk") {
    policy_cost_read_glpk_solution(solution_path, problem_path)
  } else {
    policy_cost_read_solution(solution_path)
  }
  x <- unname(solution[names_map$x])
  x[is.na(x)] <- 0
  selected <- which(x > 0.5)
  if (length(selected) != length(unique(edges$village_i))) {
    stop("Invalid assignment returned by Gurobi.", call. = FALSE)
  }
  allocation <- edges[selected, , drop = FALSE]
  allocation$base_demand <- demand[selected]
  allocation$pooled <- FALSE
  if (pooling) {
    z <- unname(solution[names_map$z])
    z[is.na(z)] <- 0
    names(z) <- sort(unique(edges$pot_j))
    allocation$pooled <- z[as.character(allocation$pot_j)] > 0.5
    allocation$demand <- ifelse(
      allocation$pooled, pooled_demand[selected], demand[selected]
    )
  } else {
    allocation$demand <- allocation$base_demand
    allocation$pooled <- ave(
      allocation$pot_j, allocation$pot_j, FUN = length
    ) > 1L
  }
  allocation$population <- population[allocation$village_i]
  allocation$expected_takers <- allocation$population * allocation$demand
  allocation$signal_cost <- allocation$expected_takers * signal_cost_per_taker
  allocation$travel_cost <- allocation$expected_takers *
    travel_cost_per_roundtrip_km * 2 * allocation$distance_km
  total_population <- sum(population)
  site_count <- length(unique(allocation$pot_j))
  site_cost_total <- site_count * site_cost
  result <- data.frame(
    sites = site_count,
    pooled_sites = length(unique(allocation$pot_j[allocation$pooled])),
    pooled_population_share = sum(allocation$population[allocation$pooled]) / total_population,
    unweighted_takeup = mean(allocation$demand),
    population_weighted_takeup = sum(allocation$expected_takers) / total_population,
    expected_takers = sum(allocation$expected_takers),
    population_mean_distance_km = weighted.mean(
      allocation$distance_km, allocation$population
    ),
    population_p90_distance_km = policy_cost_weighted_quantile(
      allocation$distance_km, allocation$population, 0.9
    ),
    site_cost = site_cost_total,
    signal_cost = sum(allocation$signal_cost),
    travel_cost = sum(allocation$travel_cost),
    total_cost = site_cost_total + sum(allocation$signal_cost) +
      sum(allocation$travel_cost),
    target_rate = target_rate,
    target_slack_takers = sum(allocation$expected_takers) - target,
    stringsAsFactors = FALSE
  )
  list(summary = result, allocation = allocation, solver_output = output)
}

`%||%` <- function(left, right) if (is.null(left)) right else left
