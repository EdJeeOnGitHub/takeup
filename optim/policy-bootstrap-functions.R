# Shared helpers for the fast, sparse cluster-bootstrap policy exercise.

policy_option_value <- function(args, name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Duplicate option: ", name, call. = FALSE)
  if (length(hit) == 0L) default else substring(hit, nchar(prefix) + 1L)
}

policy_scenarios <- data.frame(
  scenario_id = 1:5,
  scenario = c(
    "control", "bracelet", "static-control", "static-bracelet",
    "suppress-reputation"
  ),
  label = c(
    "Control", "Bracelet", "Control social image returns at 0.5km",
    "Bracelet social image returns at 0.5km", "No social image returns"
  ),
  visibility = c("control", "bracelet", "control", "bracelet", "control"),
  static_signal = c(FALSE, FALSE, TRUE, TRUE, FALSE),
  suppress_reputation = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

read_cmdstan_mode_row <- function(path) {
  if (!file.exists(path)) stop("Mode CSV not found: ", path, call. = FALSE)
  value <- read.csv(path, comment.char = "#", check.names = FALSE)
  if (nrow(value) != 1L) {
    stop("Expected one optimizer row in: ", path, call. = FALSE)
  }
  value
}

mode_scalar <- function(mode, name) {
  candidates <- c(name, paste0(name, ".1"), paste0(name, "[1]"))
  hit <- intersect(candidates, names(mode))
  if (length(hit) != 1L) stop("Missing/ambiguous mode parameter: ", name, call. = FALSE)
  value <- as.numeric(mode[[hit]])
  if (length(value) != 1L || !is.finite(value)) {
    stop("Non-finite mode parameter: ", name, call. = FALSE)
  }
  value
}

mode_vector <- function(mode, name, length_out) {
  bracket <- paste0(name, "[", seq_len(length_out), "]")
  dotted <- paste0(name, ".", seq_len(length_out))
  columns <- if (all(bracket %in% names(mode))) bracket else dotted
  if (!all(columns %in% names(mode))) stop("Missing mode vector: ", name, call. = FALSE)
  value <- as.numeric(mode[1L, columns, drop = TRUE])
  if (length(value) != length_out || any(!is.finite(value))) {
    stop("Invalid mode vector: ", name, call. = FALSE)
  }
  value
}

canonical_policy_parameters <- function(mode, replicate, mode_csv) {
  hyper_beta <- mode_vector(mode, "hyper_beta_1ord", 4L)
  hyper_dist_beta <- mode_vector(mode, "hyper_dist_beta_1ord", 4L)
  bracelet_beta <- mode_scalar(mode, "beta_bracelet_effect")
  raw_u_sd <- mode_scalar(mode, "raw_u_sd")
  data.frame(
    draw = as.integer(replicate),
    replicate = as.integer(replicate),
    mode_csv = mode_csv,
    beta_control = mode_scalar(mode, "beta_intercept"),
    beta_ink = mode_scalar(mode, "beta_ink_effect"),
    beta_calendar = bracelet_beta +
      mode_scalar(mode, "wtp_value_utility") * mode_scalar(mode, "hyper_wtp_mu"),
    beta_bracelet = bracelet_beta,
    dist_beta = mode_scalar(mode, "dist_beta_v"),
    mu_control = hyper_beta[1L],
    mu_ink = hyper_beta[2L],
    mu_calendar = hyper_beta[3L],
    mu_bracelet = hyper_beta[4L],
    mu_dist_control = hyper_dist_beta[1L],
    mu_dist_ink = hyper_dist_beta[2L],
    mu_dist_calendar = hyper_dist_beta[3L],
    mu_dist_bracelet = hyper_dist_beta[4L],
    base_mu_rep = mode_scalar(mode, "base_mu_rep"),
    u_sd = raw_u_sd,
    total_error_sd = sqrt(1 + raw_u_sd^2),
    stringsAsFactors = FALSE
  )
}

policy_delta <- function(cutoff, u_sd) {
  total_sd <- sqrt(1 + u_sd^2)
  probability <- pnorm(cutoff, sd = total_sd)
  denominator <- pmax(probability * (1 - probability), .Machine$double.xmin)
  numerator <- dnorm(cutoff, sd = total_sd)
  numerator / denominator
}

policy_fixedpoint_residual <- function(cutoff, benefit, mu_rep, u_sd) {
  cutoff + benefit + mu_rep * policy_delta(cutoff, u_sd)
}

# Exact vectorized scalar solve. The stable-equilibrium restriction gives one
# crossing for retained draws. Brackets expand defensively before failing.
solve_policy_fixedpoint <- function(
    benefit, mu_rep, u_sd, tolerance = 1e-11, iterations = 60L) {
  stopifnot(length(benefit) == length(mu_rep), length(u_sd) == 1L)
  lower <- rep(-5.5, length(benefit))
  upper <- rep(5.5, length(benefit))
  lower_value <- policy_fixedpoint_residual(lower, benefit, mu_rep, u_sd)
  upper_value <- policy_fixedpoint_residual(upper, benefit, mu_rep, u_sd)
  unresolved <- lower_value * upper_value > 0
  for (bound in c(8, 12, 20)) {
    if (!any(unresolved)) break
    lower[unresolved] <- -bound
    upper[unresolved] <- bound
    lower_value[unresolved] <- policy_fixedpoint_residual(
      lower[unresolved], benefit[unresolved], mu_rep[unresolved], u_sd
    )
    upper_value[unresolved] <- policy_fixedpoint_residual(
      upper[unresolved], benefit[unresolved], mu_rep[unresolved], u_sd
    )
    unresolved <- lower_value * upper_value > 0
  }
  fallback <- unresolved
  bracketed <- !fallback
  for (step in seq_len(iterations)) {
    midpoint <- 0.5 * (lower + upper)
    midpoint_value <- policy_fixedpoint_residual(midpoint, benefit, mu_rep, u_sd)
    take_lower_half <- bracketed & lower_value * midpoint_value <= 0
    upper[take_lower_half] <- midpoint[take_lower_half]
    take_upper_half <- bracketed & !take_lower_half
    lower[take_upper_half] <- midpoint[take_upper_half]
    lower_value[take_upper_half] <- midpoint_value[take_upper_half]
  }
  solution <- 0.5 * (lower + upper)
  if (any(fallback)) {
    if (!requireNamespace("nleqslv", quietly = TRUE)) {
      stop(
        "Legacy nleqslv fallback is required for ", sum(fallback),
        " non-bracketed fixed points.", call. = FALSE
      )
    }
    fallback_index <- which(fallback)
    fallback_fit <- lapply(fallback_index, function(index) {
      nleqslv::nleqslv(
        x = -benefit[index],
        fn = function(cutoff) policy_fixedpoint_residual(
          cutoff, benefit[index], mu_rep[index], u_sd
        )
      )
    })
    fallback_code <- vapply(fallback_fit, `[[`, numeric(1), "termcd")
    if (any(fallback_code > 2)) {
      stop("Legacy fallback failed for ", sum(fallback_code > 2),
           " policy fixed points.", call. = FALSE)
    }
    solution[fallback_index] <- vapply(fallback_fit, function(fit) fit$x, numeric(1))
  }
  residual <- abs(policy_fixedpoint_residual(solution, benefit, mu_rep, u_sd))
  residual_tolerance <- if (any(fallback)) max(tolerance, 1e-7) else tolerance
  if (any(!is.finite(solution)) || max(residual) > residual_tolerance) {
    stop("Policy fixed-point residual exceeds tolerance: ", max(residual), call. = FALSE)
  }
  attr(solution, "fallback_count") <- sum(fallback)
  solution
}

policy_mu_rep <- function(distance_sd, parameter, visibility) {
  intercept <- parameter[[paste0("mu_", visibility)]]
  slope <- parameter[[paste0("mu_dist_", visibility)]]
  if (visibility == "control") {
    latent <- intercept + slope * distance_sd
  } else {
    latent <- parameter$mu_control + intercept +
      (parameter$mu_dist_control + slope) * distance_sd
  }
  parameter$base_mu_rep * plogis(latent)
}

predict_policy_draw <- function(parameter, distances, scenarios = policy_scenarios) {
  distance_sd <- distances / parameter$sd_of_dist
  benefit <- parameter$beta_control - parameter$dist_beta * distance_sd
  dynamic_cache <- list()
  for (visibility in unique(scenarios$visibility[!scenarios$suppress_reputation])) {
    mu <- policy_mu_rep(distance_sd, parameter, visibility)
    cutoff <- solve_policy_fixedpoint(benefit, mu, parameter$u_sd)
    dynamic_cache[[visibility]] <- list(mu = mu, cutoff = cutoff)
  }
  static_signal <- lapply(unique(scenarios$visibility[scenarios$static_signal]), function(visibility) {
    distance_500_sd <- 500 / parameter$sd_of_dist
    benefit_500 <- parameter$beta_control - parameter$dist_beta * distance_500_sd
    mu_500 <- policy_mu_rep(distance_500_sd, parameter, visibility)
    cutoff_500 <- solve_policy_fixedpoint(benefit_500, mu_500, parameter$u_sd)
    signal <- mu_500 * policy_delta(cutoff_500, parameter$u_sd)
    attr(signal, "fallback_count") <- attr(cutoff_500, "fallback_count")
    signal
  })
  names(static_signal) <- unique(scenarios$visibility[scenarios$static_signal])

  pieces <- lapply(seq_len(nrow(scenarios)), function(index) {
    scenario <- scenarios[index, ]
    if (scenario$suppress_reputation) {
      cutoff <- -benefit
      fallback_count <- 0L
    } else if (scenario$static_signal) {
      cutoff <- -(benefit + static_signal[[scenario$visibility]])
      fallback_count <- attr(static_signal[[scenario$visibility]], "fallback_count") %||% 0L
    } else {
      cutoff <- dynamic_cache[[scenario$visibility]]$cutoff
      fallback_count <- attr(cutoff, "fallback_count") %||% 0L
    }
    data.frame(
      draw = parameter$draw,
      replicate = parameter$replicate,
      scenario_id = scenario$scenario_id,
      scenario = scenario$scenario,
      scenario_label = scenario$label,
      visibility = scenario$visibility,
      distance = distances,
      distance_km = distances / 1000,
      demand = 1 - pnorm(cutoff / parameter$total_error_sd),
      fixedpoint_fallbacks = fallback_count,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

`%||%` <- function(left, right) if (is.null(left)) right else left

summarize_policy_values <- function(values) {
  c(
    estimate = median(values, na.rm = TRUE),
    conf_low = unname(quantile(values, 0.025, na.rm = TRUE)),
    conf_high = unname(quantile(values, 0.975, na.rm = TRUE))
  )
}
