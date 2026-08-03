# Shared data preparation for the minimal main structural model and its compact
# generated quantities. Callers must load dplyr, purrr, and rlang first.

main_core_option_value <- function(args, name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
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
    cluster_shock_sd_prior = 0.1) {
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
  sample_data$core_cluster_weight <- read_main_core_weights(
    weight_file,
    sample_data$num_clusters
  )
  sample_data$use_core_cluster_shock <- as.integer(use_cluster_shock)
  sample_data$core_cluster_shock_sd_prior <- cluster_shock_sd_prior

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
    "core_cluster_shock_sd"
  )
  nudge_value <- function(value) {
    if (is.list(value)) return(lapply(value, nudge_value))
    pmax(value, epsilon)
  }
  for (parameter in intersect(names(init), positive_parameters)) {
    init[[parameter]] <- nudge_value(init[[parameter]])
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
  cmdstanr::write_stan_json(init, path)
  invisible(path)
}
