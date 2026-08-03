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

# Convert one retained HMC draw from any of the compact main-core models to the
# policy parameterization.  Optional model-specific quantities are represented
# explicitly so the prediction/optimization code does not need to know Stan
# parameter names.
canonical_policy_draw <- function(
    draw, draw_id, chain, model_id, model_label, model_family = "gaussian",
    lambda_structure = "common", lambda_log_ratio_sd_prior = 0.25,
    source_csv = NA_character_) {
  value <- canonical_policy_parameters(draw, draw_id, source_csv)
  value$chain <- as.integer(chain)
  value$model_id <- model_id
  value$model_label <- model_label
  value$model_family <- model_family
  value$lambda_control <- value$base_mu_rep
  value$lambda_ink <- value$base_mu_rep
  value$lambda_calendar <- value$base_mu_rep
  value$lambda_bracelet <- value$base_mu_rep
  if (lambda_structure == "grouped") {
    ratio <- lambda_log_ratio_sd_prior *
      mode_scalar(draw, "core_lambda_group_log_ratio_raw")
    value[c("lambda_control", "lambda_calendar")] <-
      value$base_mu_rep * exp(-0.5 * ratio)
    value[c("lambda_ink", "lambda_bracelet")] <-
      value$base_mu_rep * exp(0.5 * ratio)
  } else if (lambda_structure == "arm") {
    helmert_basis <- matrix(
      c(
        1 / sqrt(2), 1 / sqrt(6), 1 / sqrt(12),
        -1 / sqrt(2), 1 / sqrt(6), 1 / sqrt(12),
        0, -2 / sqrt(6), 1 / sqrt(12),
        0, 0, -3 / sqrt(12)
      ),
      nrow = 4, byrow = TRUE
    )
    raw <- mode_vector(draw, "core_lambda_arm_log_ratio_raw", 3L)
    lambda <- value$base_mu_rep * exp(
      lambda_log_ratio_sd_prior / sqrt(2) * as.vector(helmert_basis %*% raw)
    )
    value[c("lambda_control", "lambda_ink", "lambda_calendar", "lambda_bracelet")] <-
      as.list(lambda)
  } else if (lambda_structure != "common") {
    stop("Unknown lambda structure: ", lambda_structure, call. = FALSE)
  }
  value
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
    fallback_index <- which(fallback)
    if (requireNamespace("nleqslv", quietly = TRUE)) {
      fallback_fit <- lapply(fallback_index, function(index) {
        nleqslv::nleqslv(
          x = -benefit[index],
          fn = function(cutoff) policy_fixedpoint_residual(
            cutoff, benefit[index], mu_rep[index], u_sd
          )
        )
      })
      # nleqslv termination codes are advisory here: validate the equation
      # directly below instead of discarding accurate small-step solutions.
      solution[fallback_index] <- vapply(
        fallback_fit, function(fit) fit$x, numeric(1)
      )
    } else {
      # Dependency-free diagnostic fallback. This also distinguishes a
      # near-tangent root from a counterfactual with no equilibrium.
      solution[fallback_index] <- vapply(fallback_index, function(index) {
        optimize(
          function(cutoff) policy_fixedpoint_residual(
            cutoff, benefit[index], mu_rep[index], u_sd
          )^2,
          interval = c(-20, 20), tol = tolerance
        )$minimum
      }, numeric(1))
    }
  }
  residual <- abs(policy_fixedpoint_residual(solution, benefit, mu_rep, u_sd))
  residual_tolerance <- if (any(fallback)) max(tolerance, 1e-7) else tolerance
  if (any(!is.finite(solution))) stop("Non-finite policy fixed point.", call. = FALSE)
  undefined <- residual > residual_tolerance
  # A flexible posterior draw can cross the saddle-node boundary under a new
  # policy even when the fitted, observed allocation has a valid equilibrium.
  # Preserve the draw and mark that counterfactual undefined rather than
  # substituting a numerical pseudo-root.
  solution[undefined] <- NA_real_
  attr(solution, "fallback_count") <- sum(fallback)
  attr(solution, "undefined_count") <- sum(undefined)
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
  lambda_name <- paste0("lambda_", visibility)
  lambda <- if (lambda_name %in% names(parameter)) {
    parameter[[lambda_name]]
  } else {
    parameter$base_mu_rep
  }
  lambda * plogis(latent)
}

policy_student_t_mixture <- function(df = 5, components = 12L) {
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
    scale_sq = (df - 2) / df,
    precision = decomposition$values[ordering] / shape,
    weight = decomposition$vectors[1L, ordering]^2
  )
}

policy_student_t_moments <- function(cutoff, u_sd, mixture) {
  type_variance <- mixture$scale_sq / mixture$precision
  index_variance <- type_variance + u_sd^2
  standardized <- outer(cutoff, sqrt(index_variance), "/")
  density_kernel <- dnorm(standardized)
  below <- as.vector(pnorm(standardized) %*% mixture$weight)
  below <- pmin(1 - 1e-12, pmax(1e-12, below))
  numerator <- as.vector(
    density_kernel %*% (mixture$weight * type_variance / sqrt(index_variance))
  )
  index_density <- as.vector(
    density_kernel %*% (mixture$weight / sqrt(index_variance))
  )
  numerator_derivative <- as.vector(
    density_kernel %*% (mixture$weight * type_variance /
      (sqrt(index_variance) * index_variance))
  ) * -cutoff
  denominator <- below * (1 - below)
  delta <- numerator / denominator
  delta_derivative <- numerator_derivative / denominator -
    numerator * index_density * (1 - 2 * below) / denominator^2
  list(probability_below = below, delta = delta, delta_derivative = delta_derivative)
}

solve_policy_student_t_fixedpoint <- function(benefit, mu_rep, u_sd, mixture) {
  cutoff <- pmin(8, pmax(-8, -benefit))
  for (step in seq_len(8L)) {
    moments <- policy_student_t_moments(cutoff, u_sd, mixture)
    residual <- cutoff + benefit + mu_rep * moments$delta
    derivative <- 1 + mu_rep * moments$delta_derivative
    safe_derivative <- ifelse(
      abs(derivative) < 0.1, ifelse(derivative < 0, -0.1, 0.1), derivative
    )
    cutoff <- pmin(8, pmax(-8, cutoff - pmin(1, pmax(-1, residual / safe_derivative))))
  }
  cutoff
}

predict_policy_draw <- function(
    parameter, distances, scenarios = policy_scenarios, village_ids = NULL) {
  distance_sd <- distances / parameter$sd_of_dist
  benefit <- parameter$beta_control - parameter$dist_beta * distance_sd
  if (!is.null(parameter$cluster_shock)) {
    if (is.null(village_ids) || any(!village_ids %in% seq_along(parameter$cluster_shock))) {
      stop("Cluster-shock prediction requires valid village indices.", call. = FALSE)
    }
    benefit <- benefit + parameter$cluster_shock[village_ids]
  }
  dynamic_cache <- list()
  for (visibility in unique(scenarios$visibility[!scenarios$suppress_reputation])) {
    mu <- policy_mu_rep(distance_sd, parameter, visibility)
    cutoff <- if (identical(parameter$model_family, "student_t5")) {
      solve_policy_student_t_fixedpoint(
        benefit, mu, parameter$u_sd, policy_student_t_mixture()
      )
    } else {
      solve_policy_fixedpoint(benefit, mu, parameter$u_sd)
    }
    dynamic_cache[[visibility]] <- list(mu = mu, cutoff = cutoff)
  }
  static_signal <- lapply(unique(scenarios$visibility[scenarios$static_signal]), function(visibility) {
    distance_500_sd <- 500 / parameter$sd_of_dist
    benefit_500 <- parameter$beta_control - parameter$dist_beta * distance_500_sd
    if (!is.null(parameter$cluster_shock)) {
      benefit_500 <- benefit_500 + parameter$cluster_shock[village_ids]
    }
    mu_500 <- policy_mu_rep(distance_500_sd, parameter, visibility)
    if (!is.null(parameter$cluster_shock)) {
      mu_500 <- rep(mu_500, length(benefit_500))
    }
    if (identical(parameter$model_family, "student_t5")) {
      mixture <- policy_student_t_mixture()
      cutoff_500 <- solve_policy_student_t_fixedpoint(
        benefit_500, mu_500, parameter$u_sd, mixture
      )
      signal <- mu_500 * policy_student_t_moments(
        cutoff_500, parameter$u_sd, mixture
      )$delta
    } else {
      cutoff_500 <- solve_policy_fixedpoint(benefit_500, mu_500, parameter$u_sd)
      signal <- mu_500 * policy_delta(cutoff_500, parameter$u_sd)
    }
    attr(signal, "fallback_count") <- attr(cutoff_500, "fallback_count")
    attr(signal, "undefined_count") <- attr(cutoff_500, "undefined_count")
    signal
  })
  names(static_signal) <- unique(scenarios$visibility[scenarios$static_signal])

  pieces <- lapply(seq_len(nrow(scenarios)), function(index) {
    scenario <- scenarios[index, ]
    if (scenario$suppress_reputation) {
      cutoff <- -benefit
      fallback_count <- 0L
      undefined_count <- 0L
    } else if (scenario$static_signal) {
      cutoff <- -(benefit + static_signal[[scenario$visibility]])
      fallback_count <- attr(static_signal[[scenario$visibility]], "fallback_count") %||% 0L
      undefined_count <- attr(static_signal[[scenario$visibility]], "undefined_count") %||% 0L
    } else {
      cutoff <- dynamic_cache[[scenario$visibility]]$cutoff
      fallback_count <- attr(cutoff, "fallback_count") %||% 0L
      undefined_count <- attr(cutoff, "undefined_count") %||% 0L
    }
    data.frame(
      draw = parameter$draw,
      replicate = parameter$replicate,
      scenario_id = scenario$scenario_id,
      scenario = scenario$scenario,
      scenario_label = scenario$label,
      visibility = scenario$visibility,
      village_i = if (is.null(village_ids)) NA_integer_ else village_ids,
      distance = distances,
      distance_km = distances / 1000,
      demand = if (identical(parameter$model_family, "student_t5")) {
        1 - policy_student_t_moments(
          cutoff, parameter$u_sd, policy_student_t_mixture()
        )$probability_below
      } else {
        1 - pnorm(cutoff / parameter$total_error_sd)
      },
      fixedpoint_fallbacks = fallback_count,
      fixedpoint_undefined = undefined_count,
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
