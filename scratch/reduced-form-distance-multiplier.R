#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(prefix, default) {
  hit <- args[str_starts(args, paste0(prefix, "="))]
  if (length(hit) == 0) return(default)
  str_sub(hit[[1]], str_length(prefix) + 2)
}

B <- as.integer(get_arg("--bootstrap-draws", "1000"))
bootstrap_seed <- as.integer(get_arg("--seed", "20260808"))
output_dir <- get_arg("--output-dir", "ref-reports/reduced-form-distance-multiplier")
if (is.na(B) || B < 100) stop("--bootstrap-draws must be at least 100.")
if (is.na(bootstrap_seed)) stop("--seed must be an integer.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

message("Reading the canonical 9,805-person no-SMS take-up sample")
d <- read_csv(
  "temp-data/analysis-cluster-recentered-covariate-data.csv",
  show_col_types = FALSE,
  progress = FALSE
) %>%
  transmute(
    person_id = KEY.individ,
    cluster_id = as.integer(cluster.id.x),
    county = factor(county),
    arm = factor(
      assigned.treatment,
      levels = c("control", "calendar", "ink", "bracelet")
    ),
    signal_group = factor(
      if_else(assigned.treatment %in% c("ink", "bracelet"),
              "Any Signal", "No Signal"),
      levels = c("No Signal", "Any Signal")
    ),
    distance_km = cluster.dist.to.pot / 1000,
    dewormed = as.numeric(as.logical(dewormed)),
    female,
    age = age.census,
    expected_distance = mu_d
  )

original_assignment <- read_rds("data/rct_targetable_schools_2.0.rds") %>%
  as_tibble() %>%
  transmute(
    cluster_id = as.integer(as.character(pot.cluster.id)),
    assigned_far = as.integer(assigned.dist.cat == "far")
  ) %>%
  distinct()
assert_true(
  nrow(original_assignment) == n_distinct(original_assignment$cluster_id),
  "Original distance assignment is not unique by community."
)
d <- d %>%
  left_join(original_assignment, by = "cluster_id", relationship = "many-to-one")

assert_true(nrow(d) == 9805, "Expected 9,805 observations.")
assert_true(n_distinct(d$cluster_id) == 144, "Expected 144 communities.")
assert_true(!anyDuplicated(d$person_id), "Person identifiers are not unique.")
assert_true(!anyNA(d), "Analysis variables contain missing values.")

cluster_data <- d %>%
  distinct(cluster_id, county, arm, signal_group, distance_km)
assert_true(nrow(cluster_data) == 144, "Community-level fields are not constant.")

support <- cluster_data %>%
  group_by(arm) %>%
  summarise(
    communities = n(),
    minimum_km = min(distance_km),
    p05_km = quantile(distance_km, 0.05),
    p95_km = quantile(distance_km, 0.95),
    maximum_km = max(distance_km),
    .groups = "drop"
  )
write_csv(support, file.path(output_dir, "data", "distance-support.csv"))

# The flexible analysis is deliberately restricted to a conservative interval
# lying inside the empirical support of all four arms.
grid <- seq(0.4, 2.0, by = 0.05)
finite_points <- c(0.5, 2.0)
assert_true(
  min(grid) >= max(support$p05_km) && max(grid) <= min(support$p95_km),
  "Chosen grid is outside the common 5th--95th percentile support."
)

# Fixed natural-spline basis. Knots are based on the 144 community distances,
# rather than the person-weighted empirical distribution.
interior_knots <- unname(quantile(cluster_data$distance_km, c(1 / 3, 2 / 3)))
boundary_knots <- range(cluster_data$distance_km)
basis <- ns(
  d$distance_km,
  knots = interior_knots,
  Boundary.knots = boundary_knots,
  intercept = FALSE
)
colnames(basis) <- paste0("s", seq_len(ncol(basis)))
d <- bind_cols(d, as_tibble(basis))

linear_formula <- dewormed ~ 0 + arm + signal_group:distance_km +
  female + age + expected_distance + county
arm_linear_formula <- dewormed ~ 0 + arm + arm:distance_km +
  female + age + expected_distance + county
spline_formula <- dewormed ~ 0 + arm + signal_group:(s1 + s2 + s3) +
  female + age + expected_distance + county
arm_spline_formula <- dewormed ~ 0 + arm + arm:(s1 + s2 + s3) +
  female + age + expected_distance + county
itt_formula <- dewormed ~ 0 + arm + signal_group:assigned_far +
  female + age + expected_distance + county
arm_itt_formula <- dewormed ~ 0 + arm + arm:assigned_far +
  female + age + expected_distance + county

X_linear <- model.matrix(linear_formula, d)
X_arm_linear <- model.matrix(arm_linear_formula, d)
X_spline <- model.matrix(spline_formula, d)
X_arm_spline <- model.matrix(arm_spline_formula, d)
X_itt <- model.matrix(itt_formula, d)
X_arm_itt <- model.matrix(arm_itt_formula, d)
y <- d$dewormed

cluster_levels <- cluster_data$cluster_id
cluster_index <- match(d$cluster_id, cluster_levels)
county_by_cluster <- cluster_data$county

weighted_fit <- function(X, weights) {
  qr.solve(crossprod(X, X * weights), crossprod(X, y * weights), tol = 1e-10)
}

make_newdata <- function(distance, arm_name) {
  nd <- d
  nd$distance_km <- distance
  nd$arm <- factor(arm_name, levels = levels(d$arm))
  nd$signal_group <- factor(
    if (arm_name %in% c("ink", "bracelet")) "Any Signal" else "No Signal",
    levels = levels(d$signal_group)
  )
  nb <- ns(
    nd$distance_km,
    knots = interior_knots,
    Boundary.knots = boundary_knots,
    intercept = FALSE
  )
  colnames(nb) <- c("s1", "s2", "s3")
  nd$s1 <- nb[, 1]
  nd$s2 <- nb[, 2]
  nd$s3 <- nb[, 3]
  nd
}

prediction_row <- function(formula, distance, arm_name, columns) {
  nd <- make_newdata(distance, arm_name)
  mm <- model.matrix(delete.response(terms(formula)), nd)
  colMeans(mm[, columns, drop = FALSE])
}

group_prediction_row <- function(formula, distance, group_name, columns) {
  group_arms <- if (group_name == "No Signal") {
    c("control", "calendar")
  } else {
    c("ink", "bracelet")
  }
  Reduce(`+`, map(group_arms, ~prediction_row(formula, distance, .x, columns))) /
    length(group_arms)
}

derivative_row <- function(formula, distance, group_name, columns, h = 1e-4) {
  (group_prediction_row(formula, distance + h, group_name, columns) -
     group_prediction_row(formula, distance - h, group_name, columns)) / (2 * h)
}

arm_derivative_row <- function(formula, distance, arm_name, columns, h = 1e-4) {
  (prediction_row(formula, distance + h, arm_name, columns) -
     prediction_row(formula, distance - h, arm_name, columns)) / (2 * h)
}

itt_prediction_row <- function(far, arm_name, columns) {
  nd <- d
  nd$assigned_far <- far
  nd$arm <- factor(arm_name, levels = levels(d$arm))
  nd$signal_group <- factor(
    if (arm_name %in% c("ink", "bracelet")) "Any Signal" else "No Signal",
    levels = levels(d$signal_group)
  )
  mm <- model.matrix(delete.response(terms(itt_formula)), nd)
  colMeans(mm[, columns, drop = FALSE])
}

itt_group_row <- function(group_name, columns) {
  group_arms <- if (group_name == "No Signal") {
    c("control", "calendar")
  } else {
    c("ink", "bracelet")
  }
  Reduce(`+`, map(group_arms, function(a) {
    itt_prediction_row(1, a, columns) - itt_prediction_row(0, a, columns)
  })) / length(group_arms)
}

arm_itt_row <- function(arm_name, columns) {
  prediction <- function(far) {
    nd <- d
    nd$assigned_far <- far
    nd$arm <- factor(arm_name, levels = levels(d$arm))
    mm <- model.matrix(delete.response(terms(arm_itt_formula)), nd)
    colMeans(mm[, columns, drop = FALSE])
  }
  prediction(1) - prediction(0)
}

groups <- levels(d$signal_group)
arms <- levels(d$arm)

P_curve <- map_dfr(groups, function(g) {
  map_dfr(grid, function(x) {
    as_tibble_row(group_prediction_row(spline_formula, x, g, colnames(X_spline))) %>%
      mutate(signal_group = g, distance_km = x, .before = 1)
  })
})
P_curve <- as.matrix(P_curve %>% select(all_of(colnames(X_spline))))

D_curve <- map_dfr(groups, function(g) {
  map_dfr(grid, function(x) {
    as_tibble_row(derivative_row(spline_formula, x, g, colnames(X_spline))) %>%
      mutate(signal_group = g, distance_km = x, .before = 1)
  })
})
D_curve <- as.matrix(D_curve %>% select(all_of(colnames(X_spline))))

D_arm_curve <- map_dfr(arms, function(a) {
  map_dfr(grid, function(x) {
    as_tibble_row(arm_derivative_row(
      arm_spline_formula, x, a, colnames(X_arm_spline)
    )) %>% mutate(arm = a, distance_km = x, .before = 1)
  })
})
D_arm_curve <- as.matrix(D_arm_curve %>% select(all_of(colnames(X_arm_spline))))

curve_index <- crossing(signal_group = groups, distance_km = grid) %>%
  arrange(match(signal_group, groups), distance_km)
arm_curve_index <- crossing(arm = arms, distance_km = grid) %>%
  arrange(match(arm, arms), distance_km)

finite_rows <- map_dfr(groups, function(g) {
  lo <- group_prediction_row(spline_formula, finite_points[1], g, colnames(X_spline))
  hi <- group_prediction_row(spline_formula, finite_points[2], g, colnames(X_spline))
  as_tibble_row((hi - lo) / diff(finite_points)) %>%
    mutate(signal_group = g, .before = 1)
})
finite_rows <- as.matrix(finite_rows %>% select(all_of(colnames(X_spline))))

arm_finite_rows <- map_dfr(arms, function(a) {
  lo <- prediction_row(arm_spline_formula, finite_points[1], a, colnames(X_arm_spline))
  hi <- prediction_row(arm_spline_formula, finite_points[2], a, colnames(X_arm_spline))
  as_tibble_row((hi - lo) / diff(finite_points)) %>%
    mutate(arm = a, .before = 1)
})
arm_finite_rows <- as.matrix(arm_finite_rows %>% select(all_of(colnames(X_arm_spline))))

linear_slope_rows <- map_dfr(groups, function(g) {
  as_tibble_row(derivative_row(linear_formula, 1.2, g, colnames(X_linear))) %>%
    mutate(signal_group = g, .before = 1)
})
linear_slope_rows <- as.matrix(linear_slope_rows %>% select(all_of(colnames(X_linear))))

arm_linear_slope_rows <- map_dfr(arms, function(a) {
  as_tibble_row(arm_derivative_row(arm_linear_formula, 1.2, a, colnames(X_arm_linear))) %>%
    mutate(arm = a, .before = 1)
})
arm_linear_slope_rows <- as.matrix(
  arm_linear_slope_rows %>% select(all_of(colnames(X_arm_linear)))
)

itt_rows <- map_dfr(groups, function(g) {
  as_tibble_row(itt_group_row(g, colnames(X_itt))) %>%
    mutate(signal_group = g, .before = 1)
})
itt_rows <- as.matrix(itt_rows %>% select(all_of(colnames(X_itt))))

arm_itt_rows <- map_dfr(arms, function(a) {
  as_tibble_row(arm_itt_row(a, colnames(X_arm_itt))) %>% mutate(arm = a, .before = 1)
})
arm_itt_rows <- as.matrix(arm_itt_rows %>% select(all_of(colnames(X_arm_itt))))

message("Estimating point curves")
beta_linear <- weighted_fit(X_linear, rep(1, nrow(d)))
beta_arm_linear <- weighted_fit(X_arm_linear, rep(1, nrow(d)))
beta_spline <- weighted_fit(X_spline, rep(1, nrow(d)))
beta_arm_spline <- weighted_fit(X_arm_spline, rep(1, nrow(d)))
beta_itt <- weighted_fit(X_itt, rep(1, nrow(d)))
beta_arm_itt <- weighted_fit(X_arm_itt, rep(1, nrow(d)))

point_curve <- as.numeric(P_curve %*% beta_spline)
point_derivative <- as.numeric(D_curve %*% beta_spline)
point_arm_derivative <- as.numeric(D_arm_curve %*% beta_arm_spline)
point_finite <- as.numeric(finite_rows %*% beta_spline)
point_arm_finite <- as.numeric(arm_finite_rows %*% beta_arm_spline)
point_linear <- as.numeric(linear_slope_rows %*% beta_linear)
point_arm_linear <- as.numeric(arm_linear_slope_rows %*% beta_arm_linear)
point_itt <- as.numeric(itt_rows %*% beta_itt)
point_arm_itt <- as.numeric(arm_itt_rows %*% beta_arm_itt)

message("Running ", B, " county-stratified exponential community-weight draws")
set.seed(bootstrap_seed)
draw_curve <- matrix(NA_real_, B, nrow(P_curve))
draw_derivative <- matrix(NA_real_, B, nrow(D_curve))
draw_arm_derivative <- matrix(NA_real_, B, nrow(D_arm_curve))
draw_finite <- matrix(NA_real_, B, length(groups))
draw_arm_finite <- matrix(NA_real_, B, length(arms))
draw_linear <- matrix(NA_real_, B, length(groups))
draw_arm_linear <- matrix(NA_real_, B, length(arms))
draw_itt <- matrix(NA_real_, B, length(groups))
draw_arm_itt <- matrix(NA_real_, B, length(arms))
failed <- logical(B)

for (b in seq_len(B)) {
  cw <- rexp(length(cluster_levels))
  for (county_i in levels(county_by_cluster)) {
    idx <- which(county_by_cluster == county_i)
    cw[idx] <- cw[idx] / mean(cw[idx])
  }
  ow <- cw[cluster_index]
  result <- tryCatch({
    bl <- weighted_fit(X_linear, ow)
    bal <- weighted_fit(X_arm_linear, ow)
    bs <- weighted_fit(X_spline, ow)
    ba <- weighted_fit(X_arm_spline, ow)
    bi <- weighted_fit(X_itt, ow)
    bai <- weighted_fit(X_arm_itt, ow)
    list(
      curve = as.numeric(P_curve %*% bs),
      derivative = as.numeric(D_curve %*% bs),
      arm_derivative = as.numeric(D_arm_curve %*% ba),
      finite = as.numeric(finite_rows %*% bs),
      arm_finite = as.numeric(arm_finite_rows %*% ba),
      linear = as.numeric(linear_slope_rows %*% bl),
      arm_linear = as.numeric(arm_linear_slope_rows %*% bal),
      itt = as.numeric(itt_rows %*% bi),
      arm_itt = as.numeric(arm_itt_rows %*% bai)
    )
  }, error = function(e) NULL)
  if (is.null(result)) {
    failed[b] <- TRUE
  } else {
    draw_curve[b, ] <- result$curve
    draw_derivative[b, ] <- result$derivative
    draw_arm_derivative[b, ] <- result$arm_derivative
    draw_finite[b, ] <- result$finite
    draw_arm_finite[b, ] <- result$arm_finite
    draw_linear[b, ] <- result$linear
    draw_arm_linear[b, ] <- result$arm_linear
    draw_itt[b, ] <- result$itt
    draw_arm_itt[b, ] <- result$arm_itt
  }
  if (b %% 100 == 0) message("  completed ", b, " of ", B)
}

keep <- which(!failed)
assert_true(length(keep) >= 0.98 * B, "More than two percent of draws failed.")
draw_curve <- draw_curve[keep, , drop = FALSE]
draw_derivative <- draw_derivative[keep, , drop = FALSE]
draw_arm_derivative <- draw_arm_derivative[keep, , drop = FALSE]
draw_finite <- draw_finite[keep, , drop = FALSE]
draw_arm_finite <- draw_arm_finite[keep, , drop = FALSE]
draw_linear <- draw_linear[keep, , drop = FALSE]
draw_arm_linear <- draw_arm_linear[keep, , drop = FALSE]
draw_itt <- draw_itt[keep, , drop = FALSE]
draw_arm_itt <- draw_arm_itt[keep, , drop = FALSE]

summarize_vector <- function(point, draws) {
  n_draws <- length(draws)
  lower_tail <- (sum(draws <= 0) + 1) / (n_draws + 1)
  upper_tail <- (sum(draws >= 0) + 1) / (n_draws + 1)
  tibble(
    estimate = point,
    std_error = sd(draws),
    conf_low = unname(quantile(draws, 0.025)),
    conf_high = unname(quantile(draws, 0.975)),
    p_value = min(1, 2 * min(lower_tail, upper_tail))
  )
}

curve_summary <- map_dfr(seq_len(nrow(curve_index)), function(j) {
  bind_cols(curve_index[j, ], summarize_vector(point_curve[j], draw_curve[, j]))
})
derivative_summary <- map_dfr(seq_len(nrow(curve_index)), function(j) {
  bind_cols(curve_index[j, ], summarize_vector(point_derivative[j], draw_derivative[, j]))
})
arm_derivative_summary <- map_dfr(seq_len(nrow(arm_curve_index)), function(j) {
  bind_cols(
    arm_curve_index[j, ],
    summarize_vector(point_arm_derivative[j], draw_arm_derivative[, j])
  )
})

no_idx <- which(curve_index$signal_group == "No Signal")
any_idx <- which(curve_index$signal_group == "Any Signal")
point_diff <- point_derivative[any_idx] - point_derivative[no_idx]
draw_diff <- draw_derivative[, any_idx, drop = FALSE] -
  draw_derivative[, no_idx, drop = FALSE]
diff_summary <- map_dfr(seq_along(grid), function(j) {
  bind_cols(tibble(distance_km = grid[j]), summarize_vector(point_diff[j], draw_diff[, j]))
})

# A simultaneous band for the derivative difference over the prespecified grid.
grid_sd <- apply(draw_diff, 2, sd)
studentized <- sweep(draw_diff, 2, point_diff, "-")
studentized <- sweep(studentized, 2, grid_sd, "/")
critical_value <- unname(quantile(apply(abs(studentized), 1, max), 0.95))
diff_summary <- diff_summary %>%
  mutate(
    simultaneous_low = estimate - critical_value * std_error,
    simultaneous_high = estimate + critical_value * std_error
  )

finite_diff_point <- point_finite[2] - point_finite[1]
finite_diff_draw <- draw_finite[, 2] - draw_finite[, 1]
finite_ratio_point <- point_finite[2] / point_finite[1]
finite_ratio_draw <- draw_finite[, 2] / draw_finite[, 1]
finite_attenuation_point <- 1 - finite_ratio_point
finite_attenuation_draw <- 1 - finite_ratio_draw

linear_diff_point <- point_linear[2] - point_linear[1]
linear_diff_draw <- draw_linear[, 2] - draw_linear[, 1]
linear_ratio_point <- point_linear[2] / point_linear[1]
linear_ratio_draw <- draw_linear[, 2] / draw_linear[, 1]

itt_diff_point <- point_itt[2] - point_itt[1]
itt_diff_draw <- draw_itt[, 2] - draw_itt[, 1]
itt_ratio_point <- point_itt[2] / point_itt[1]
itt_ratio_draw <- draw_itt[, 2] / draw_itt[, 1]

main_summary <- bind_rows(
  summarize_vector(point_itt[1], draw_itt[, 1]) %>%
    mutate(specification = "Original assignment", estimand = "No Signal distance effect"),
  summarize_vector(point_itt[2], draw_itt[, 2]) %>%
    mutate(specification = "Original assignment", estimand = "Any Signal distance effect"),
  summarize_vector(itt_diff_point, itt_diff_draw) %>%
    mutate(specification = "Original assignment", estimand = "Any Signal minus No Signal"),
  summarize_vector(itt_ratio_point, itt_ratio_draw) %>%
    mutate(specification = "Original assignment", estimand = "Relative distance response"),
  summarize_vector(1 - itt_ratio_point, 1 - itt_ratio_draw) %>%
    mutate(specification = "Original assignment", estimand = "Distance-response attenuation"),
  summarize_vector(point_linear[1], draw_linear[, 1]) %>%
    mutate(specification = "Linear", estimand = "No Signal slope"),
  summarize_vector(point_linear[2], draw_linear[, 2]) %>%
    mutate(specification = "Linear", estimand = "Any Signal slope"),
  summarize_vector(linear_diff_point, linear_diff_draw) %>%
    mutate(specification = "Linear", estimand = "Any Signal minus No Signal slope"),
  summarize_vector(linear_ratio_point, linear_ratio_draw) %>%
    mutate(specification = "Linear", estimand = "Relative distance response"),
  summarize_vector(1 - linear_ratio_point, 1 - linear_ratio_draw) %>%
    mutate(specification = "Linear", estimand = "Distance-response attenuation"),
  summarize_vector(point_finite[1], draw_finite[, 1]) %>%
    mutate(specification = "Spline, 0.5--2.0 km", estimand = "No Signal average slope"),
  summarize_vector(point_finite[2], draw_finite[, 2]) %>%
    mutate(specification = "Spline, 0.5--2.0 km", estimand = "Any Signal average slope"),
  summarize_vector(finite_diff_point, finite_diff_draw) %>%
    mutate(specification = "Spline, 0.5--2.0 km", estimand = "Any Signal minus No Signal slope"),
  summarize_vector(finite_ratio_point, finite_ratio_draw) %>%
    mutate(specification = "Spline, 0.5--2.0 km", estimand = "Relative distance response"),
  summarize_vector(finite_attenuation_point, finite_attenuation_draw) %>%
    mutate(specification = "Spline, 0.5--2.0 km", estimand = "Distance-response attenuation")
) %>%
  select(specification, estimand, everything())

write_csv(main_summary, file.path(output_dir, "data", "main-summary.csv"))

calendar_idx <- match("calendar", arms)
bracelet_idx <- match("bracelet", arms)
bracelet_calendar_summary <- function(specification, point, draws) {
  difference_point <- point[bracelet_idx] - point[calendar_idx]
  difference_draw <- draws[, bracelet_idx] - draws[, calendar_idx]
  ratio_point <- point[bracelet_idx] / point[calendar_idx]
  ratio_draw <- draws[, bracelet_idx] / draws[, calendar_idx]
  bind_rows(
    summarize_vector(point[calendar_idx], draws[, calendar_idx]) %>%
      mutate(estimand = "Calendar distance response"),
    summarize_vector(point[bracelet_idx], draws[, bracelet_idx]) %>%
      mutate(estimand = "Bracelet distance response"),
    summarize_vector(difference_point, difference_draw) %>%
      mutate(estimand = "Bracelet minus Calendar"),
    summarize_vector(ratio_point, ratio_draw) %>%
      mutate(estimand = "Bracelet relative distance response"),
    summarize_vector(1 - ratio_point, 1 - ratio_draw) %>%
      mutate(estimand = "Bracelet attenuation relative to Calendar")
  ) %>% mutate(specification = specification, .before = 1)
}

bracelet_calendar_summary_table <- bind_rows(
  bracelet_calendar_summary("Original assignment", point_arm_itt, draw_arm_itt),
  bracelet_calendar_summary("Linear", point_arm_linear, draw_arm_linear),
  bracelet_calendar_summary("Spline, 0.5--2.0 km", point_arm_finite, draw_arm_finite)
) %>% select(specification, estimand, everything())
write_csv(
  bracelet_calendar_summary_table,
  file.path(output_dir, "data", "bracelet-calendar-summary.csv")
)

bracelet_calendar_ratio_audit <- bind_rows(
  tibble(specification = "Original assignment", denominator = draw_arm_itt[, calendar_idx]),
  tibble(specification = "Linear", denominator = draw_arm_linear[, calendar_idx]),
  tibble(specification = "Spline, 0.5--2.0 km", denominator = draw_arm_finite[, calendar_idx])
) %>%
  group_by(specification) %>%
  summarise(
    draws = n(),
    calendar_positive_share = mean(denominator >= 0),
    calendar_near_zero_share = mean(abs(denominator) < 0.02),
    calendar_p025 = quantile(denominator, 0.025),
    calendar_p975 = quantile(denominator, 0.975),
    .groups = "drop"
  )
write_csv(
  bracelet_calendar_ratio_audit,
  file.path(output_dir, "data", "bracelet-calendar-ratio-audit.csv")
)
write_csv(curve_summary, file.path(output_dir, "data", "pooled-curves.csv"))
write_csv(derivative_summary, file.path(output_dir, "data", "pooled-derivatives.csv"))
write_csv(diff_summary, file.path(output_dir, "data", "derivative-contrast.csv"))
write_csv(arm_derivative_summary, file.path(output_dir, "data", "arm-derivatives.csv"))

bootstrap_audit <- tibble(
  requested_draws = B,
  successful_draws = length(keep),
  failed_draws = sum(failed),
  seed = bootstrap_seed,
  simultaneous_critical_value = critical_value,
  grid_minimum_km = min(grid),
  grid_maximum_km = max(grid)
)
write_csv(bootstrap_audit, file.path(output_dir, "data", "bootstrap-audit.csv"))

group_palette <- c("No Signal" = "#C44E52", "Any Signal" = "#4C72B0")
p_curve <- ggplot(curve_summary, aes(distance_km, estimate, color = signal_group, fill = signal_group)) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = group_palette) +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Community-centroid distance to treatment site (km)",
       y = "Adjusted take-up", color = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

p_derivative <- ggplot(diff_summary, aes(distance_km, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_ribbon(aes(ymin = simultaneous_low, ymax = simultaneous_high),
              fill = "#8172B2", alpha = 0.13) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high),
              fill = "#8172B2", alpha = 0.25) +
  geom_line(color = "#8172B2", linewidth = 1) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Community-centroid distance to treatment site (km)",
       y = "Any Signal minus No Signal derivative\n(per km)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(output_dir, "figures", "pooled-distance-curves.pdf"),
       p_curve, width = 6.4, height = 4.2)
ggsave(file.path(output_dir, "figures", "derivative-contrast.pdf"),
       p_derivative, width = 6.4, height = 4.2)

fmt <- function(x) sprintf("%.3f", x)
tex_rows <- main_summary %>%
  mutate(
    label = case_when(
      estimand %in% c("No Signal distance effect", "No Signal slope",
                      "No Signal average slope") ~ "No Signal distance response",
      estimand %in% c("Any Signal distance effect", "Any Signal slope",
                      "Any Signal average slope") ~ "Any Signal distance response",
      estimand %in% c("Any Signal minus No Signal",
                      "Any Signal minus No Signal slope") ~ "Any Signal $-$ No Signal",
      estimand == "Relative distance response" ~ "Relative distance response",
      estimand == "Distance-response attenuation" ~ "Distance-response attenuation",
      TRUE ~ estimand
    ),
    cell = paste0(fmt(estimate), " \\ (", fmt(conf_low), ", ", fmt(conf_high), ")")
  ) %>%
  select(specification, label, cell) %>%
  pivot_wider(names_from = specification, values_from = cell)

tex <- c(
  "\\begin{tabular}{lccc}",
  "\\toprule",
  " & Original assignment & Linear distance & Spline, 0.5--2.0 km \\\\",
  "\\midrule",
  paste0(tex_rows$label, " & ", tex_rows[["Original assignment"]],
         " & ", tex_rows$Linear, " & ",
         tex_rows[["Spline, 0.5--2.0 km"]], " \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(tex, file.path(output_dir, "tables", "reduced-form-distance-multiplier.tex"))

note <- c(
  "# Reduced-form distance response and relative multiplier",
  "",
  sprintf("This analysis uses %s adults in %s communities and %s successful county-stratified exponential community-weight draws.",
          format(nrow(d), big.mark = ","), n_distinct(d$cluster_id), length(keep)),
  "The primary comparison pools Control and Calendar as No Signal and Ink and Bracelet as Any Signal.",
  "",
  "## Main estimates",
  "",
  sprintf("- Linear No Signal slope: %.3f per km (95%% interval %.3f, %.3f).",
          point_linear[1], quantile(draw_linear[, 1], .025), quantile(draw_linear[, 1], .975)),
  sprintf("- Linear Any Signal slope: %.3f per km (95%% interval %.3f, %.3f).",
          point_linear[2], quantile(draw_linear[, 2], .025), quantile(draw_linear[, 2], .975)),
  sprintf("- Linear slope contrast, Any Signal minus No Signal: %.3f (95%% interval %.3f, %.3f).",
          linear_diff_point, quantile(linear_diff_draw, .025), quantile(linear_diff_draw, .975)),
  sprintf("- Linear relative response: %.3f; attenuation: %.3f.",
          linear_ratio_point, 1 - linear_ratio_point),
  sprintf("- Spline finite-change relative response over 0.5--2.0 km: %.3f; attenuation: %.3f.",
          finite_ratio_point, finite_attenuation_point),
  sprintf("- Original-assignment No Signal Far--Close effect: %.3f (95%% interval %.3f, %.3f).",
          point_itt[1], quantile(draw_itt[, 1], .025), quantile(draw_itt[, 1], .975)),
  sprintf("- Original-assignment Any Signal Far--Close effect: %.3f (95%% interval %.3f, %.3f).",
          point_itt[2], quantile(draw_itt[, 2], .025), quantile(draw_itt[, 2], .975)),
  sprintf("- Original-assignment relative response: %.3f; attenuation: %.3f.",
          itt_ratio_point, 1 - itt_ratio_point),
  "",
  "## Interpretation",
  "",
  "The arm-specific and pooled derivatives are reduced-form total equilibrium distance responses.",
  "Their pooled ratio identifies a relative multiplier if the direct access-cost response is common across randomized signal environments.",
  "The exercise does not identify an absolute multiplier relative to zero observability, because every experimental arm retains positive observability.",
  "The spline is restricted to 0.4--2.0 km, within the common 5th--95th percentile distance support of all four arms.",
  "Pointwise bands describe each distance; the derivative-contrast figure also shows a simultaneous 95 percent band over the prespecified grid."
)
writeLines(note, file.path(output_dir, "README.md"))

message("Wrote reduced-form distance multiplier review bundle to ", output_dir)
