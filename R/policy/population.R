# One population/experimental-target contract for all policy runners.
policy_object_hash <- function(object) {
  path <- tempfile(); on.exit(unlink(path))
  saveRDS(object, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

policy_adult_population <- function(distance_object, census_path = Sys.getenv(
    "POLICY_CENSUS", "data/takeup_census.RData")) {
  villages <- distance_object$village_df
  sites <- distance_object$pot_df
  geography <- distance_object$long_distance_mat
  if (!identical(distance_object$candidate_site_mode, "all") || nrow(sites) != 1451L ||
      anyDuplicated(sites$cluster.id) || nrow(geography) != 144L * 1451L ||
      anyDuplicated(geography[c("index_i", "index_j")])) {
    stop("Adult policy requires the complete canonical 144-by-1451 geography.")
  }
  if (nrow(villages) != 144L || anyNA(villages$cluster.id) ||
      anyDuplicated(villages$cluster.id) ||
      !identical(as.integer(villages$id), seq_len(144L))) {
    stop("Population contract requires ordered canonical 144 communities.")
  }
  env <- new.env(parent = emptyenv()); load(census_path, envir = env)
  census <- env$census.data
  required <- c("KEY.individ", "cluster.id", "age.census")
  if (!all(required %in% names(census))) stop("Missing census identity/age fields.")
  selected <- census[census$cluster.id %in% villages$cluster.id, required]
  if (anyNA(selected) || any(!is.finite(selected$age.census))) stop("Missing census identity/age.")
  if (anyDuplicated(selected$KEY.individ)) stop("Duplicate census adults or inconsistent community assignments.")
  selected <- selected[selected$age.census >= 18, ]
  counts <- tabulate(match(selected$cluster.id, villages$cluster.id), nbins = 144L)
  if (nrow(selected) != 39301L || sum(counts) != 39301L || any(counts <= 0)) {
    stop("Expected 39,301 unique census adults in 144 policy communities.")
  }
  data.frame(village_i = villages$id, cluster.id = villages$cluster.id,
             population = counts, stringsAsFactors = FALSE)
}

policy_write_experimental_targets <- function(experimental, population, output_path,
                                               model_id, source_hash) {
  if (anyDuplicated(experimental[c("draw", "village_i")])) stop("Repeated experimental community/draw.")
  rows <- lapply(split(experimental, experimental$draw), function(x) {
    index <- match(population$village_i, x$village_i)
    if (nrow(x) != nrow(population) || anyNA(index) || length(unique(x$replicate)) != 1L) {
      stop("Incomplete experimental predictions within draw.")
    }
    demand <- x$demand[index]
    expected <- sum(population$population * demand)
    data.frame(draw = x$draw[1L], replicate = x$replicate[1L],
      target_expected_adults = expected, target_rate = expected / sum(population$population),
      target_community_welfare = sum(demand), target_community_rate = mean(demand),
      population_total = sum(population$population), model_id = model_id,
      parameter_hash = source_hash, population_checksum = policy_object_hash(population),
      target_mode = "draw-specific-experimental-control",
      target_defined = all(is.finite(demand)), stringsAsFactors = FALSE)
  })
  targets <- do.call(rbind, rows); targets <- targets[order(targets$draw), ]
  write.csv(population, file.path(output_path, "policy-population.csv"), row.names = FALSE)
  write.csv(targets, file.path(output_path, "policy-experimental-targets.csv"), row.names = FALSE)
  invisible(targets)
}

policy_experimental_summary <- function(input_path, scenario_results) {
  x <- readRDS(file.path(input_path, "policy-experimental-demand.rds"))
  path <- file.path(input_path, "policy-population.csv")
  population <- if (file.exists(path)) read.csv(path) else NULL
  weighting <- if ("population_weighting" %in% names(scenario_results)) unique(scenario_results$population_weighting) else "equal-community"
  if (length(weighting) != 1L) stop("Mixed population conventions in scenario summaries.")
  rows <- lapply(split(x, x$draw), function(d) {
    weights <- if (weighting == "adult-census") population$population[match(d$village_i, population$village_i)] else rep(1, nrow(d))
    if (anyNA(weights) || anyDuplicated(d$village_i)) stop("Invalid experimental summary join.")
    draw_target <- unique(scenario_results$target_welfare[scenario_results$draw == d$draw[1L]])
    if (length(draw_target) != 1L) stop("Scenarios disagree on experimental target.")
    data.frame(draw = d$draw[1L], replicate = d$replicate[1L],
      mean_demand = weighted.mean(d$demand, weights), mean_distance = weighted.mean(d$distance, weights),
      scenario_id = 0L, scenario = "experimental", scenario_label = "Experimental allocation",
      status = if (all(is.finite(d$demand))) "observed_allocation" else "equilibrium_undefined",
      solver_status = NA_integer_, elapsed_seconds = NA_real_, n_pot = nrow(d),
      achieved_welfare = sum(weights * d$demand), target_welfare = draw_target,
      population_weighting = weighting, population_total = sum(weights), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
