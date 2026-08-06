# Shared data preparation for the minimal main structural model and its compact
# generated quantities. Callers must load dplyr, purrr, and rlang first.

source("scratch/main-core-asymmetric-observability-data.R")

main_core_option_value <- function(args, name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
}

# Change scalar controls in an existing Stan JSON file without parsing and
# reserializing the full document. jsonlite simplifies length-one arrays to
# scalars, which changes Stan data dimensions; a textual scalar patch preserves
# every original array and matrix shape.
main_core_patch_stan_json_scalars <- function(input_path, output_path, values) {
  if (!file.exists(input_path)) {
    stop("Stan JSON not found: ", input_path, call. = FALSE)
  }
  if (is.null(names(values)) || any(!nzchar(names(values)))) {
    stop("JSON scalar replacements must be named.", call. = FALSE)
  }
  json_text <- paste(readLines(input_path, warn = FALSE), collapse = "\n")
  json_number <- paste0(
    "[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)",
    "(?:[eE][+-]?[0-9]+)?"
  )
  for (field in names(values)) {
    if (length(values[[field]]) != 1L || !is.finite(values[[field]])) {
      stop("Replacement for ", field, " must be one finite scalar.", call. = FALSE)
    }
    pattern <- paste0(
      "(\\\"", field, "\\\"[[:space:]]*:[[:space:]]*)(", json_number, ")"
    )
    hits <- gregexpr(pattern, json_text, perl = TRUE)[[1L]]
    if (identical(hits, -1L) || length(hits) != 1L) {
      stop("Expected exactly one scalar JSON field named ", field, ".", call. = FALSE)
    }
    encoded <- as.character(jsonlite::toJSON(
      values[[field]], auto_unbox = TRUE, digits = NA
    ))
    json_text <- sub(pattern, paste0("\\1", encoded), json_text, perl = TRUE)
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(json_text, output_path, useBytes = TRUE)
  output_path
}

main_core_student_t_mixture <- function(df = 5, components = 12L) {
  if (!is.finite(df) || df <= 2) stop("Student-t df must exceed 2.", call. = FALSE)
  components <- as.integer(components)
  if (components < 2L) stop("At least two mixture components are required.", call. = FALSE)
  shape <- df / 2
  alpha <- shape - 1
  index <- seq_len(components)
  jacobi <- diag(2 * index - 1 + alpha)
  off_diagonal <- sqrt(index[-components] * (index[-components] + alpha))
  jacobi[cbind(index[-components], index[-1L])] <- off_diagonal
  jacobi[cbind(index[-1L], index[-components])] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(decomposition$values)
  list(
    df = df,
    components = components,
    scale_sq = (df - 2) / df,
    precision = decomposition$values[ordering] / shape,
    weight = decomposition$vectors[1L, ordering]^2
  )
}

read_main_core_weights <- function(path, num_clusters) {
  if (is.null(path)) return(rep(1, num_clusters))
  if (!file.exists(path)) stop("Weight file not found: ", path, call. = FALSE)
  weight_data <- read.csv(path, stringsAsFactors = FALSE)
  if (!identical(names(weight_data), c("cluster_id", "weight")) ||
      nrow(weight_data) != num_clusters ||
      anyDuplicated(weight_data$cluster_id) ||
      !setequal(weight_data$cluster_id, seq_len(num_clusters)) ||
      any(!is.finite(weight_data$weight)) ||
      any(weight_data$weight < 0)) {
    stop(
      "Weight CSV must contain one nonnegative cluster_id,weight row per cluster.",
      call. = FALSE
    )
  }
  weight_data$weight[order(weight_data$cluster_id)]
}

prepare_main_core_data <- function(
    workspace_path,
    model_name = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP",
    weight_file = NULL,
    use_cluster_shock = 0L,
    cluster_shock_sd_prior = 0.1,
    lambda_structure = 0L,
    lambda_log_ratio_sd_prior = 0.25,
    profile_group_lambda = 0L,
    profile_group_log_ratio = 0,
    gq_override_lambda = 0L,
    gq_lambda_override = NULL,
    type_distribution = 0L,
    student_t_df = 5,
    student_t_components = 12L,
    observation_model = 0L,
    recognition_structure = 0L,
    report_structure = 0L,
    report_arm_dist_hierarchical = 0L,
    report_arm_dist_prior_scale = 0.25,
    project_root = ".",
    peer_audit_path = NULL) {
  if (!file.exists(workspace_path)) {
    stop("Workspace not found: ", workspace_path, call. = FALSE)
  }
  fit_env <- new.env(parent = globalenv())
  load(workspace_path, envir = fit_env)
  if (!all(c("models", "stan_data") %in% ls(fit_env)) ||
      !model_name %in% names(fit_env$models)) {
    stop("Workspace lacks the requested model: ", model_name, call. = FALSE)
  }
  model_info <- fit_env$models[[model_name]]
  stan_data_preprocess <- model_info$stan_data_preprocess %||% identity
  model_info$stan_data_preprocess <- NULL
  model_info$model_file <- NULL

  sample_data <- stan_data_preprocess(fit_env$stan_data) |>
    map_at(
      c("cluster_treatment_map", "beliefs_ate_pairs"),
      \(x) mutate(x, across(everything(), as.integer)) |> as.matrix()
    ) |>
    list_modify(!!!model_info) |>
    map_if(is.factor, as.integer)
  sample_data$num_dist_group_mix <- 1L
  sample_data$use_belief_row_cluster_mu_rep <-
    sample_data$use_belief_row_cluster_mu_rep %||% 0L
  if (length(sample_data$num_belief_rows_by_cluster) != sample_data$num_clusters) {
    sample_data$num_belief_rows_by_cluster <- tabulate(
      sample_data$obs_cluster_id[sample_data$beliefs_obs_index],
      nbins = sample_data$num_clusters
    )
  }
  if (length(sample_data$belief_observed) != sample_data$num_beliefs_obs) {
    sample_data$belief_observed <- rep.int(1L, sample_data$num_beliefs_obs)
  }
  # Preserve a one-element Stan vector when older CmdStanR/jsonlite versions
  # would otherwise auto-unbox it to a scalar JSON number.
  if (sample_data$num_optim_distances == 1L) {
    sample_data$optim_distances <- array(sample_data$optim_distances, dim = 1L)
  }
  sample_data$core_cluster_weight <- read_main_core_weights(
    weight_file,
    sample_data$num_clusters
  )
  sample_data$use_core_cluster_shock <- as.integer(use_cluster_shock)
  sample_data$core_cluster_shock_sd_prior <- cluster_shock_sd_prior
  sample_data$core_lambda_structure <- as.integer(lambda_structure)
  sample_data$core_lambda_log_ratio_sd_prior <- lambda_log_ratio_sd_prior
  sample_data$core_profile_group_lambda <- as.integer(profile_group_lambda)
  sample_data$core_profile_group_log_ratio <- profile_group_log_ratio
  sample_data$core_gq_override_lambda <- as.integer(gq_override_lambda)
  sample_data$core_gq_lambda_override <- if (is.null(gq_lambda_override)) {
    rep(0, sample_data$num_treatments)
  } else {
    gq_lambda_override
  }
  type_distribution <- as.integer(type_distribution)
  if (!type_distribution %in% 0:2) {
    stop("type_distribution must be 0 (Gaussian), 1 (Student-t), or 2 (finite mixture).",
         call. = FALSE)
  }
  type_mixture <- main_core_student_t_mixture(
    student_t_df, student_t_components
  )
  sample_data$core_type_distribution <- type_distribution
  sample_data$core_student_t_df <- type_mixture$df
  sample_data$core_type_scale_sq <- type_mixture$scale_sq
  sample_data$core_type_mixture_components <- type_mixture$components
  sample_data$core_type_mixture_precision <- type_mixture$precision
  sample_data$core_type_mixture_weight <- type_mixture$weight
  observation_model <- as.integer(observation_model)
  if (!observation_model %in% 0:2) {
    stop("observation_model must be 0, 1, or 2.", call. = FALSE)
  }
  sample_data$core_observation_model <- observation_model
  recognition_structure <- as.integer(recognition_structure)
  report_structure <- as.integer(report_structure)
  if (!recognition_structure %in% 0:2 || !report_structure %in% 0:2) {
    stop("recognition_structure and report_structure must be 0, 1, or 2.",
         call. = FALSE)
  }
  if (observation_model == 0L &&
      (recognition_structure != 0L || report_structure != 0L)) {
    stop("Observation restrictions require an asymmetric model.", call. = FALSE)
  }
  if (observation_model == 2L && recognition_structure == 2L) {
    stop("The unconditional model cannot condition recognition out.",
         call. = FALSE)
  }
  sample_data$core_recognition_structure <- recognition_structure
  sample_data$core_report_structure <- report_structure
  report_arm_dist_hierarchical <- as.integer(report_arm_dist_hierarchical)
  if (!report_arm_dist_hierarchical %in% 0:2 ||
      !is.finite(report_arm_dist_prior_scale) ||
      report_arm_dist_prior_scale <= 0) {
    stop("Invalid report arm-distance hierarchy controls.", call. = FALSE)
  }
  if (report_arm_dist_hierarchical > 0L &&
      (observation_model == 0L || report_structure != 0L)) {
    stop("Hierarchical slopes require the full multinomial channel.",
         call. = FALSE)
  }
  sample_data$core_report_arm_dist_hierarchical <-
    report_arm_dist_hierarchical
  sample_data$core_report_arm_dist_prior_scale <-
    report_arm_dist_prior_scale
  peer_data <- if (observation_model == 0L) {
    main_core_empty_peer_response_data()
  } else {
    main_core_prepare_peer_response_data(
      sample_data,
      project_root = project_root,
      write_audit = peer_audit_path
    )
  }
  peer_data$core_peer_link_audit <- NULL
  sample_data <- modifyList(sample_data, peer_data)

  if (is.null(sample_data$wtp_cluster_id)) {
    if (!is.null(weight_file)) {
      stop(
        "Weighted refits require a regenerated workspace with wtp_cluster_id.",
        call. = FALSE
      )
    }
    sample_data$wtp_cluster_id <- rep.int(1L, sample_data$num_wtp_obs)
  }
  if (length(sample_data$wtp_cluster_id) != sample_data$num_wtp_obs ||
      any(!sample_data$wtp_cluster_id %in% seq_len(sample_data$num_clusters))) {
    stop("Invalid wtp_cluster_id in workspace.", call. = FALSE)
  }

  sample_data$control <- NULL
  sample_data$analysis_data <- NULL
  discard(
    sample_data,
    \(x) is.function(x) || is.character(x) || is.null(x)
  )
}

main_core_nudge_init_boundaries <- function(init, epsilon = 1e-4) {
  # Optimizers may return an exact boundary value (most commonly zero for
  # wtp_value_utility).  Such a value is a valid constrained-space optimum but
  # CmdStan cannot transform it back to the unconstrained space for HMC.
  positive_parameters <- c(
    "group_dist_sd", "county_dist_effect_sd", "cluster_dist_effect_sd",
    "strata_wtp_mu_tau", "wtp_sigma",
    "stratum_beta_1ord_sd", "cluster_beta_1ord_sd", "obs_beta_1ord_sd",
    "stratum_dist_beta_1ord_sd", "cluster_dist_beta_1ord_sd",
    "obs_dist_beta_1ord_sd", "stratum_beta_2ord_sd",
    "cluster_beta_2ord_sd", "obs_beta_2ord_sd",
    "stratum_dist_beta_2ord_sd", "cluster_dist_beta_2ord_sd",
    "obs_dist_beta_2ord_sd", "obs_beta_common_sd",
    "structural_beta_cluster_sd", "structural_beta_county_sd",
    "base_mu_rep", "raw_u_sd", "raw_cluster_sd_tilde",
    "dist_beta_cluster_sd", "dist_beta_county_sd",
    "dist_quadratic_beta_v", "wtp_value_utility",
    "core_cluster_shock_sd", "core_report_within_dist_sd"
  )
  nudge_value <- function(value) {
    if (is.list(value)) return(lapply(value, nudge_value))
    pmax(value, epsilon)
  }
  for (parameter in intersect(names(init), positive_parameters)) {
    init[[parameter]] <- nudge_value(init[[parameter]])
  }
  for (parameter in intersect(
    names(init),
    c("core_finite_mixture_weight", "core_finite_mixture_between_share")
  )) {
    init[[parameter]] <- pmin(1 - epsilon, pmax(epsilon, init[[parameter]]))
  }
  # These effects are lower-bounded only when the private-incentive
  # restrictions are active. Preserve negative unrestricted values, but move
  # an optimizer solution at exactly zero into the constrained space.
  for (parameter in intersect(
    names(init), c("beta_calendar_effect", "beta_bracelet_effect")
  )) {
    value <- init[[parameter]]
    if (is.numeric(value) && length(value) == 1L &&
        value >= 0 && value < epsilon) {
      init[[parameter]] <- epsilon
    }
  }
  init
}

write_mode_init_json <- function(fit, model, path) {
  available <- unique(sub("\\[.*$", "", fit$metadata()$model_params))
  parameters <- intersect(names(model$variables()$parameters), available)
  draws <- posterior::as_draws_rvars(
    fit$draws(variables = parameters)
  )
  init <- lapply(names(draws), function(parameter) {
    rv <- draws[[parameter]]
    values <- posterior::draws_of(rv)
    parameter_rank <- model$variables()$parameters[[parameter]]$dimensions
    if (parameter_rank == 0L) return(as.numeric(values[1]))
    value_dim <- dim(rv)
    index <- c(list(values, 1L), rep(list(TRUE), length(value_dim)))
    selected <- do.call(`[`, c(index, list(drop = FALSE)))
    array(as.numeric(selected), dim = value_dim)
  })
  names(init) <- names(draws)
  init <- main_core_nudge_init_boundaries(init)
  # posterior::rvar represents a length-one vector like a scalar, and scalar
  # pmin/pmax boundary nudges can also drop its dimension. Restore the declared
  # Stan rank last; otherwise CmdStan rejects vector[1] initial values as
  # scalars (notably the finite-mixture parameters).
  for (parameter in names(init)) {
    parameter_rank <- model$variables()$parameters[[parameter]]$dimensions
    if (parameter_rank > 0L && length(init[[parameter]]) == 1L &&
        is.null(dim(init[[parameter]]))) {
      init[[parameter]] <- array(
        as.numeric(init[[parameter]]),
        dim = rep.int(1L, parameter_rank)
      )
    }
  }
  cmdstanr::write_stan_json(init, path)
  invisible(path)
}
