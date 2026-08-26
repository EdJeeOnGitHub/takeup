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
policy_asymmetric_families <- c(
  "asymmetric_conditional", "asymmetric_unconditional",
  "asymmetric_f1", "asymmetric_f2", "asymmetric_f3", "asymmetric_u3"
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

mode_matrix <- function(mode, name, nrow, ncol) {
  bracket <- as.vector(outer(
    seq_len(nrow), seq_len(ncol),
    function(i, j) paste0(name, "[", i, ",", j, "]")
  ))
  dotted <- as.vector(outer(
    seq_len(nrow), seq_len(ncol),
    function(i, j) paste(name, i, j, sep = ".")
  ))
  columns <- if (all(bracket %in% names(mode))) bracket else dotted
  if (!all(columns %in% names(mode))) stop("Missing mode matrix: ", name, call. = FALSE)
  matrix(as.numeric(mode[1L, columns, drop = TRUE]), nrow = nrow, ncol = ncol)
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
    # Optimizer modes used by the cluster bootstrap are draws from the
    # Gaussian-type baseline model.  Record that family explicitly so the
    # shared prediction code follows the same branch as retained HMC draws.
    model_family = "gaussian",
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
  if (identical(model_family, "finite_mixture")) {
    value$mixture_weight <- mode_scalar(draw, "core_finite_mixture_weight")
    value$mixture_between_share <-
      mode_scalar(draw, "core_finite_mixture_between_share")
  }
  if (model_family %in% policy_asymmetric_families) {
    if (model_family %in% c("asymmetric_f3", "asymmetric_u3")) {
      measurement <- c(
        if (model_family == "asymmetric_u3") {
          mode_vector(draw, "core_recognition_intercept", 2L)
        } else numeric(),
        mode_vector(draw, "core_definite_intercept", 2L),
        mode_vector(draw, "core_definite_dist_slope", 2L),
        as.vector(mode_matrix(
          draw, "core_definite_arm_intercept_raw", 2L, 3L
        )),
        mode_vector(draw, "core_definite_public_signal_dist_slope", 1L),
        mode_vector(draw, "core_accuracy_intercept", 2L),
        as.vector(mode_matrix(
          draw, "core_accuracy_arm_intercept_raw", 2L, 3L
        ))
      )
    } else {
      full_recognition <- model_family != "asymmetric_f2"
      full_report <- model_family %in% c(
        "asymmetric_conditional", "asymmetric_unconditional"
      )
      measurement <- c(
        if (full_recognition) {
          mode_vector(draw, "core_recognition_intercept", 2L)
        } else numeric(2L),
        if (full_recognition) {
          mode_vector(draw, "core_recognition_dist_slope", 2L)
        } else numeric(2L),
        if (full_recognition) as.vector(mode_matrix(
          draw, "core_recognition_arm_intercept_raw", 2L, 3L
        )) else numeric(6L),
        if (full_recognition) as.vector(mode_matrix(
          draw, "core_recognition_arm_dist_raw", 2L, 3L
        )) else numeric(6L),
        as.vector(mode_matrix(draw, "core_report_intercept", 2L, 2L)),
        as.vector(mode_matrix(draw, "core_report_dist_slope", 2L, 2L)),
        as.vector(mode_matrix(draw, "core_report_arm_intercept_raw", 2L, 6L)),
        if (full_report) as.vector(mode_matrix(
          draw, "core_report_arm_dist_raw", 2L, 6L
        )) else numeric(12L)
      )
    }
    names(measurement) <- paste0("observation_parameter_", seq_along(measurement))
    for (parameter_name in names(measurement)) {
      value[[parameter_name]] <- measurement[[parameter_name]]
    }
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

policy_noisy_channel <- function(distance_sd, parameter, visibility) {
  helmert_basis <- matrix(
    c(1 / sqrt(2), 1 / sqrt(6), 1 / sqrt(12),
      -1 / sqrt(2), 1 / sqrt(6), 1 / sqrt(12),
      0, -2 / sqrt(6), 1 / sqrt(12),
      0, 0, -3 / sqrt(12)),
    nrow = 4, byrow = TRUE
  )
  treatment <- match(visibility, c("control", "ink", "calendar", "bracelet"))
  if (is.na(treatment)) stop("Unknown policy visibility: ", visibility, call. = FALSE)
  if (parameter$model_family %in% c("asymmetric_f3", "asymmetric_u3")) {
    parameter_count <- if (parameter$model_family == "asymmetric_u3") 21L else 19L
    value <- unlist(parameter[paste0(
      "observation_parameter_", seq_len(parameter_count)
    )], use.names = FALSE)
    if (length(value) != parameter_count || any(!is.finite(value))) {
      stop("Missing two-stage observability policy parameters.", call. = FALSE)
    }
    cursor <- 0L
    take <- function(n) {
      index <- cursor + seq_len(n)
      cursor <<- cursor + n
      value[index]
    }
    recognition_intercept <- if (parameter$model_family == "asymmetric_u3") {
      take(2L)
    } else c(Inf, Inf)
    definite_intercept <- take(2L)
    definite_dist_slope <- take(2L)
    definite_arm_intercept <- matrix(take(6L), 2L, 3L)
    definite_public_dist_slope <- take(1L)
    accuracy_intercept <- take(2L)
    accuracy_arm_intercept <- matrix(take(6L), 2L, 3L)
    public_signal <- as.integer(visibility %in% c("ink", "bracelet"))
    q <- recognition <- vector("list", 2L)
    for (truth in 1:2) {
      recognition[[truth]] <- plogis(recognition_intercept[truth])
      definite <- plogis(
        definite_intercept[truth] +
          sum(definite_arm_intercept[truth, ] * helmert_basis[treatment, ]) +
          (definite_dist_slope[truth] +
           definite_public_dist_slope * public_signal) * distance_sd
      )
      accuracy <- plogis(
        accuracy_intercept[truth] +
          sum(accuracy_arm_intercept[truth, ] * helmert_basis[treatment, ])
      )
      conditional <- if (truth == 2L) {
        cbind(definite * accuracy, definite * (1 - accuracy), 1 - definite)
      } else {
        cbind(definite * (1 - accuracy), definite * accuracy, 1 - definite)
      }
      q[[truth]] <- if (parameter$model_family == "asymmetric_f3") {
        conditional
      } else {
        cbind(recognition[[truth]] * conditional, 1 - recognition[[truth]])
      }
    }
    return(list(nontaker = q[[1L]], taker = q[[2L]], recognition = recognition))
  }
  value <- unlist(parameter[paste0("observation_parameter_", 1:48)], use.names = FALSE)
  if (length(value) != 48L || any(!is.finite(value))) {
    stop("Missing asymmetric-observability policy parameters.", call. = FALSE)
  }
  cursor <- 0L
  take <- function(n) {
    index <- cursor + seq_len(n)
    cursor <<- cursor + n
    value[index]
  }
  recognition_intercept <- take(2L)
  recognition_dist_slope <- take(2L)
  recognition_arm_intercept <- matrix(take(6L), 2L, 3L)
  recognition_arm_dist <- matrix(take(6L), 2L, 3L)
  report_intercept <- matrix(take(4L), 2L, 2L)
  report_dist_slope <- matrix(take(4L), 2L, 2L)
  report_arm_intercept <- matrix(take(12L), 2L, 6L)
  report_arm_dist <- matrix(take(12L), 2L, 6L)
  q <- vector("list", 2L)
  recognition <- vector("list", 2L)
  for (truth in 1:2) {
    recognition_eta <- recognition_intercept[truth] +
      sum(recognition_arm_intercept[truth, ] * helmert_basis[treatment, ]) +
      (recognition_dist_slope[truth] +
       sum(recognition_arm_dist[truth, ] * helmert_basis[treatment, ])) *
      distance_sd
    recognition[[truth]] <- plogis(recognition_eta)
    report_eta <- vapply(1:2, function(category) {
      columns <- (category - 1L) * 3L + 1:3
      report_intercept[truth, category] +
        sum(report_arm_intercept[truth, columns] * helmert_basis[treatment, ]) +
        (report_dist_slope[truth, category] +
         sum(report_arm_dist[truth, columns] * helmert_basis[treatment, ])) *
        distance_sd
    }, numeric(length(distance_sd)))
    if (length(distance_sd) == 1L) {
      report_eta <- matrix(report_eta, nrow = 1L, ncol = 2L)
    }
    full_eta <- cbind(report_eta, 0)
    exp_eta <- exp(full_eta - apply(full_eta, 1L, max))
    conditional <- exp_eta / rowSums(exp_eta)
    q[[truth]] <- if (parameter$model_family %in% c(
      "asymmetric_conditional", "asymmetric_f1", "asymmetric_f2"
    )) {
      conditional
    } else {
      cbind(recognition[[truth]] * conditional, 1 - recognition[[truth]])
    }
  }
  list(nontaker = q[[1L]], taker = q[[2L]], recognition = recognition)
}

policy_noisy_information <- function(cutoff, total_error_sd, channel) {
  if (nrow(channel$taker) == 1L && length(cutoff) > 1L) {
    channel$taker <- channel$taker[rep.int(1L, length(cutoff)), , drop = FALSE]
    channel$nontaker <- channel$nontaker[
      rep.int(1L, length(cutoff)), , drop = FALSE
    ]
  }
  if (nrow(channel$taker) != length(cutoff)) {
    stop("Noisy policy channel and cutoff lengths differ.", call. = FALSE)
  }
  probability <- pnorm(-cutoff / total_error_sd)
  difference <- channel$taker - channel$nontaker
  signal_probability <- probability * channel$taker +
    (1 - probability) * channel$nontaker
  probability * (1 - probability) * rowSums(difference^2 / signal_probability)
}

policy_noisy_residual <- function(cutoff, benefit, lambda, u_sd, channel) {
  cutoff + benefit + lambda *
    policy_noisy_information(cutoff, sqrt(1 + u_sd^2), channel) *
    policy_delta(cutoff, u_sd)
}

solve_policy_noisy_fixedpoint <- function(
    benefit, lambda, u_sd, channel, tolerance = 1e-9, iterations = 70L) {
  lower <- rep(-8, length(benefit))
  upper <- rep(8, length(benefit))
  lower_value <- policy_noisy_residual(lower, benefit, lambda, u_sd, channel)
  upper_value <- policy_noisy_residual(upper, benefit, lambda, u_sd, channel)
  bracketed <- lower_value * upper_value <= 0
  for (step in seq_len(iterations)) {
    midpoint <- 0.5 * (lower + upper)
    midpoint_value <- policy_noisy_residual(midpoint, benefit, lambda, u_sd, channel)
    use_lower <- bracketed & lower_value * midpoint_value <= 0
    upper[use_lower] <- midpoint[use_lower]
    use_upper <- bracketed & !use_lower
    lower[use_upper] <- midpoint[use_upper]
    lower_value[use_upper] <- midpoint_value[use_upper]
  }
  solution <- 0.5 * (lower + upper)
  residual <- abs(policy_noisy_residual(solution, benefit, lambda, u_sd, channel))
  undefined <- !bracketed | !is.finite(residual) | residual > tolerance
  solution[undefined] <- NA_real_
  attr(solution, "fallback_count") <- sum(!bracketed)
  attr(solution, "undefined_count") <- sum(undefined)
  solution
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

policy_finite_mixture <- function(parameter) {
  weight <- parameter$mixture_weight
  between <- parameter$mixture_between_share
  if (!is.finite(weight) || !is.finite(between) || weight <= 0 || weight >= 1 ||
      between < 0 || between >= 1) {
    stop("Invalid finite-mixture shape parameters.", call. = FALSE)
  }
  separation <- sqrt(between / (weight * (1 - weight)))
  list(
    mean = c(-(1 - weight) * separation, weight * separation),
    variance = rep(1 - between, 2L),
    weight = c(weight, 1 - weight)
  )
}

policy_finite_mixture_moments <- function(cutoff, u_sd, mixture) {
  index_variance <- mixture$variance + u_sd^2
  index_sd <- sqrt(index_variance)
  z <- outer(cutoff, mixture$mean, "-")
  z <- sweep(z, 2L, index_sd, "/")
  density <- dnorm(z)
  component_below <- pnorm(z)
  below <- as.vector(component_below %*% mixture$weight)
  below <- pmin(1 - 1e-12, pmax(1e-12, below))
  numerator <- as.vector((
    sweep(density, 2L, mixture$variance / index_sd, "*") +
      sweep(component_below, 2L, -mixture$mean, "*")
  ) %*% mixture$weight)
  delta <- numerator / (below * (1 - below))
  list(probability_below = below, delta = delta)
}

solve_policy_finite_mixture_fixedpoint <- function(
    benefit, mu_rep, u_sd, mixture, tolerance = 1e-9, iterations = 70L) {
  residual <- function(cutoff) cutoff + benefit + mu_rep *
    policy_finite_mixture_moments(cutoff, u_sd, mixture)$delta
  lower <- rep(-12, length(benefit))
  upper <- rep(12, length(benefit))
  lower_value <- residual(lower)
  upper_value <- residual(upper)
  bracketed <- is.finite(lower_value) & is.finite(upper_value) &
    lower_value * upper_value <= 0
  for (step in seq_len(iterations)) {
    midpoint <- 0.5 * (lower + upper)
    midpoint_value <- residual(midpoint)
    left <- bracketed & lower_value * midpoint_value <= 0
    upper[left] <- midpoint[left]
    right <- bracketed & !left
    lower[right] <- midpoint[right]
    lower_value[right] <- midpoint_value[right]
  }
  solution <- 0.5 * (lower + upper)
  bad <- !bracketed | !is.finite(residual(solution)) |
    abs(residual(solution)) > tolerance
  solution[bad] <- NA_real_
  attr(solution, "fallback_count") <- sum(!bracketed)
  attr(solution, "undefined_count") <- sum(bad)
  solution
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
  # Older canonical parameter files predate the explicit family column and
  # are draws from the Gaussian baseline model.
  model_family <- parameter$model_family %||% "gaussian"
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
    asymmetric <- model_family %in% policy_asymmetric_families
    mu <- if (asymmetric) {
      rep(parameter$base_mu_rep, length(distance_sd))
    } else {
      policy_mu_rep(distance_sd, parameter, visibility)
    }
    cutoff <- if (asymmetric) {
      solve_policy_noisy_fixedpoint(
        benefit, parameter$base_mu_rep, parameter$u_sd,
        policy_noisy_channel(distance_sd, parameter, visibility)
      )
    } else if (identical(model_family, "student_t5")) {
      solve_policy_student_t_fixedpoint(
        benefit, mu, parameter$u_sd, policy_student_t_mixture()
      )
    } else if (identical(model_family, "finite_mixture")) {
      solve_policy_finite_mixture_fixedpoint(
        benefit, mu, parameter$u_sd, policy_finite_mixture(parameter)
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
    asymmetric <- model_family %in% policy_asymmetric_families
    mu_500 <- if (asymmetric) parameter$base_mu_rep else
      policy_mu_rep(distance_500_sd, parameter, visibility)
    if (!is.null(parameter$cluster_shock)) {
      mu_500 <- rep(mu_500, length(benefit_500))
    }
    if (asymmetric) {
      channel_500 <- policy_noisy_channel(
        distance_500_sd, parameter, visibility
      )
      cutoff_500 <- solve_policy_noisy_fixedpoint(
        benefit_500, parameter$base_mu_rep, parameter$u_sd, channel_500
      )
      signal <- parameter$base_mu_rep * policy_noisy_information(
        cutoff_500, parameter$total_error_sd, channel_500
      ) * policy_delta(cutoff_500, parameter$u_sd)
    } else if (identical(model_family, "student_t5")) {
      mixture <- policy_student_t_mixture()
      cutoff_500 <- solve_policy_student_t_fixedpoint(
        benefit_500, mu_500, parameter$u_sd, mixture
      )
      signal <- mu_500 * policy_student_t_moments(
        cutoff_500, parameter$u_sd, mixture
      )$delta
    } else if (identical(model_family, "finite_mixture")) {
      mixture <- policy_finite_mixture(parameter)
      cutoff_500 <- solve_policy_finite_mixture_fixedpoint(
        benefit_500, mu_500, parameter$u_sd, mixture
      )
      signal <- mu_500 * policy_finite_mixture_moments(
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
      demand = if (identical(model_family, "student_t5")) {
        1 - policy_student_t_moments(
          cutoff, parameter$u_sd, policy_student_t_mixture()
        )$probability_below
      } else if (identical(model_family, "finite_mixture")) {
        1 - policy_finite_mixture_moments(
          cutoff, parameter$u_sd, policy_finite_mixture(parameter)
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

policy_haversine_m <- function(lon1, lat1, lon2, lat2) {
  radians <- pi / 180
  lon1 <- lon1 * radians
  lat1 <- lat1 * radians
  lon2 <- lon2 * radians
  lat2 <- lat2 * radians
  delta_lon <- lon2 - lon1
  delta_lat <- lat2 - lat1
  value <- sin(delta_lat / 2)^2 +
    cos(lat1) * cos(lat2) * sin(delta_lon / 2)^2
  2 * 6371008.8 * asin(pmin(1, sqrt(value)))
}

prepare_policy_household_edges <- function(distance_object, workspace, edges) {
  if (!file.exists(workspace)) stop("Household workspace not found: ", workspace)
  fit_environment <- new.env(parent = emptyenv())
  suppressWarnings(load(workspace, envir = fit_environment))
  analysis <- fit_environment$stan_data$analysis_data
  required <- c("cluster.id", "lon", "lat")
  if (!all(required %in% names(analysis))) {
    stop("Household workspace lacks cluster.id/lon/lat.", call. = FALSE)
  }
  household <- data.frame(
    household_i = seq_len(nrow(analysis)),
    cluster.id = analysis$cluster.id,
    lon = as.numeric(analysis$lon), lat = as.numeric(analysis$lat),
    observed_distance = as.numeric(analysis$dist.to.pot)
  )
  if (any(!is.finite(household$lon) | !is.finite(household$lat))) {
    stop("Household coordinate frame contains missing coordinates.", call. = FALSE)
  }
  household$village_i <- match(
    household$cluster.id, distance_object$village_df$cluster.id
  )
  if (anyNA(household$village_i)) {
    stop("Household-to-policy-community mapping is incomplete.", call. = FALSE)
  }
  counts <- tabulate(household$village_i, nbins = nrow(distance_object$village_df))
  if (any(counts == 0L)) stop("At least one policy community has no households.")
  edge_households <- merge(
    transform(edges, edge_i = seq_len(nrow(edges))),
    household, by = "village_i", sort = FALSE
  )
  site_row <- match(edge_households$pot_j, distance_object$pot_df$id)
  if (anyNA(site_row)) stop("Candidate-site mapping is incomplete.", call. = FALSE)
  edge_households$household_distance <- policy_haversine_m(
    edge_households$lon, edge_households$lat,
    distance_object$pot_df$lon[site_row], distance_object$pot_df$lat[site_row]
  )
  edge_households <- edge_households[order(edge_households$edge_i), ]
  list(
    data = edge_households,
    status = data.frame(
      coordinate_source = normalizePath(workspace),
      observations = nrow(household), communities = length(unique(household$village_i)),
      candidate_sites = nrow(distance_object$pot_df), feasible_edges = nrow(edges),
      household_edge_rows = nrow(edge_households),
      missing_coordinates = 0L,
      distance_method = "great-circle haversine; radius 6371008.8m",
      cap_interpretation = "3.5km cap on community-centroid assignment; all mapped households retained",
      memory_mb = as.numeric(object.size(edge_households)) / 1024^2,
      stringsAsFactors = FALSE
    )
  )
}

predict_policy_household_draw <- function(parameter, edges, household_edges) {
  family <- parameter$model_family
  if (!family %in% c("private_distance_community_image", "full_information")) {
    stop("Unsupported household-distance policy family: ", family)
  }
  hh <- household_edges
  if (identical(family, "full_information")) {
    raw <- predict_policy_draw(parameter, hh$household_distance)
    pieces <- lapply(seq_len(nrow(policy_scenarios)), function(index) {
      rows <- ((index - 1L) * nrow(hh) + 1L):(index * nrow(hh))
      sums <- rowsum(raw$demand[rows], hh$edge_i, reorder = FALSE)
      counts <- rowsum(rep(1, nrow(hh)), hh$edge_i, reorder = FALSE)
      as.numeric(sums / counts)
    })
    return(matrix(unlist(pieces, use.names = FALSE), nrow = 1L))
  }

  centroid_distance_sd <- edges$distance / parameter$sd_of_dist
  household_distance_sd <- hh$household_distance / parameter$sd_of_dist
  centroid_benefit <- parameter$beta_control -
    parameter$dist_beta * centroid_distance_sd
  household_benefit <- parameter$beta_control -
    parameter$dist_beta * household_distance_sd
  values <- vector("list", nrow(policy_scenarios))
  for (scenario_index in seq_len(nrow(policy_scenarios))) {
    scenario <- policy_scenarios[scenario_index, ]
    if (scenario$suppress_reputation) {
      household_cutoff <- -household_benefit
    } else if (scenario$static_signal) {
      distance_500_sd <- 500 / parameter$sd_of_dist
      benefit_500 <- parameter$beta_control - parameter$dist_beta * distance_500_sd
      mu_500 <- policy_mu_rep(distance_500_sd, parameter, scenario$visibility)
      cutoff_500 <- solve_policy_fixedpoint(benefit_500, mu_500, parameter$u_sd)
      signal_500 <- mu_500 * policy_delta(cutoff_500, parameter$u_sd)
      household_cutoff <- -(household_benefit + signal_500)
    } else {
      mu <- policy_mu_rep(centroid_distance_sd, parameter, scenario$visibility)
      community_cutoff <- solve_policy_fixedpoint(
        centroid_benefit, mu, parameter$u_sd
      )
      # This is algebraically identical to the fitted private-info model:
      # -community cutoff - household cost + community-centroid cost.
      household_cutoff <- community_cutoff[hh$edge_i] +
        parameter$dist_beta *
          (household_distance_sd - centroid_distance_sd[hh$edge_i])
    }
    demand <- 1 - pnorm(household_cutoff / parameter$total_error_sd)
    sums <- rowsum(demand, hh$edge_i, reorder = FALSE)
    counts <- rowsum(rep(1, length(demand)), hh$edge_i, reorder = FALSE)
    values[[scenario_index]] <- as.numeric(sums / counts)
  }
  matrix(unlist(values, use.names = FALSE), nrow = 1L)
}

`%||%` <- function(left, right) if (is.null(left)) right else left

summarize_policy_values <- function(values) {
  c(
    estimate = median(values, na.rm = TRUE),
    conf_low = unname(quantile(values, 0.025, na.rm = TRUE)),
    conf_high = unname(quantile(values, 0.975, na.rm = TRUE))
  )
}
