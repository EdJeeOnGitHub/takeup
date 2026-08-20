#!/usr/bin/env Rscript

# Unit checks for the identified two-normal type mixture used by
# takeup_struct_main_core.stan. These test the normalization and the analytic
# inference-moment derivative independently of Stan.

mixture_components <- function(weight, between_share) {
  stopifnot(weight > 0, weight < 1, between_share >= 0, between_share < 1)
  separation <- sqrt(between_share / (weight * (1 - weight)))
  list(
    weight = c(weight, 1 - weight),
    mean = c(-(1 - weight) * separation, weight * separation),
    variance = rep(1 - between_share, 2)
  )
}

mixture_moments <- function(cutoff, u_sd, component, use_u = TRUE) {
  shock_variance <- as.numeric(use_u) * u_sd^2
  index_variance <- component$variance + shock_variance
  index_sd <- sqrt(index_variance)
  z <- (cutoff - component$mean) / index_sd
  below <- sum(component$weight * pnorm(z))
  density <- sum(component$weight * dnorm(z) / index_sd)
  numerator <- sum(component$weight * (
    component$variance / index_sd * dnorm(z) -
      component$mean * pnorm(z)
  ))
  numerator_deriv <- -sum(component$weight *
    (component$variance * cutoff + component$mean * shock_variance) /
      (index_sd * index_variance) * dnorm(z))
  denominator <- below * (1 - below)
  c(
    below = below,
    density = density,
    delta = numerator / denominator,
    delta_deriv = numerator_deriv / denominator -
      numerator * density * (1 - 2 * below) / denominator^2
  )
}

for (weight in c(0.2, 0.5, 0.8)) {
  for (between_share in c(0, 0.15, 0.6)) {
    component <- mixture_components(weight, between_share)
    weighted_mean <- sum(component$weight * component$mean)
    weighted_variance <- sum(component$weight * (
      component$variance + component$mean^2
    )) - weighted_mean^2
    stopifnot(abs(weighted_mean) < 1e-12, abs(weighted_variance - 1) < 1e-12)
  }
}

component <- mixture_components(0.37, 0.42)
for (cutoff in c(-1.5, -0.2, 0.7, 1.8)) {
  analytic <- mixture_moments(cutoff, 0.6, component)
  step <- 1e-5
  numeric_delta_deriv <- (
    mixture_moments(cutoff + step, 0.6, component)[["delta"]] -
      mixture_moments(cutoff - step, 0.6, component)[["delta"]]
  ) / (2 * step)
  stopifnot(abs(analytic[["delta_deriv"]] - numeric_delta_deriv) < 1e-7)
}

# At zero between-component variance, the mixture exactly equals N(0, 1),
# regardless of its otherwise unidentified weight.
gaussian_component <- mixture_components(0.23, 0)
for (cutoff in c(-1, 0, 1)) {
  mixture <- mixture_moments(cutoff, 0.4, gaussian_component)
  sd_index <- sqrt(1 + 0.4^2)
  stopifnot(
    abs(mixture[["below"]] - pnorm(cutoff / sd_index)) < 1e-12,
    abs(mixture[["density"]] - dnorm(cutoff / sd_index) / sd_index) < 1e-12
  )
}

message("Finite-mixture normalization and moment checks passed.")
