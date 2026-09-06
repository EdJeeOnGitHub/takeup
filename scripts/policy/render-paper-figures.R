#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/policy/bootstrap.R")
value <- function(name, default = NULL) policy_option_value(args, name, default)
parameter_csv <- value("--parameter-csv", value("--median-parameter-csv"))
policy_path <- value("--policy-path")
distance_data <- value("--distance-data", "optim/data/full-many-pots-experiment.rds")
output_path <- value("--output-path", "build/policy/posterior/figures")
legacy_draw <- value("--draw", "1")
parameter_draw <- value("--parameter-draw", legacy_draw)
allocation_draw <- as.integer(value("--allocation-draw", legacy_draw))
median_allocation_path <- value("--median-allocation-path")
if (any(vapply(list(parameter_csv, policy_path), is.null, logical(1)))) {
  stop("--parameter-csv (or --median-parameter-csv) and --policy-path are required.")
}
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(readr) })
parameters <- read.csv(parameter_csv, stringsAsFactors = FALSE)
if (!"draw" %in% names(parameters)) stop("Parameter CSV lacks draw.")
if (parameter_draw == "median") {
  parameter <- parameters[1L, , drop = FALSE]
  numeric_columns <- names(parameters)[vapply(parameters, is.numeric, logical(1))]
  numeric_columns <- setdiff(numeric_columns, c("draw", "replicate"))
  for (column in numeric_columns) parameter[[column]] <- median(parameters[[column]])
  parameter$draw <- NA_integer_
  parameter$replicate <- NA_integer_
} else {
  selected_draw <- as.integer(parameter_draw)
  parameter <- parameters[parameters$draw == selected_draw, , drop = FALSE]
  if (nrow(parameter) != 1L) {
    stop("Expected exactly one parameter row for draw ", parameter_draw, ".")
  }
}
distance_object <- readRDS(distance_data)
# These maps use stored longitude/latitude columns, with no spatial operations.
for (name in c("village_df", "pot_df")) {
  x <- distance_object[[name]]
  columns <- setdiff(names(x), c("geometry", attr(x, "sf_column")))
  distance_object[[name]] <- as.data.frame(setNames(lapply(columns, function(column) x[[column]]), columns))
}
if (!"sd_of_dist" %in% names(parameter)) {
  parameter$sd_of_dist <- distance_object$sd_of_dist
}
distance <- seq(0, 3500, by = 25)
curves <- predict_policy_draw(
  parameter, distance,
  policy_scenarios[policy_scenarios$scenario_id %in% 1:4, ]
)
curves$visibility <- ifelse(grepl("bracelet", curves$scenario), "Bracelet", "Control")
curves$return_type <- ifelse(grepl("static", curves$scenario),
                             "Fixed at 0.5 km", "Distance dependent")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
p_demand <- ggplot(curves, aes(distance / 1000, demand, colour = visibility,
                               linetype = return_type)) +
  geom_line(linewidth = 1) + theme_minimal(base_size = 11) +
  labs(x = "Distance (km)", y = "Predicted take-up", colour = "Observability",
       linetype = "Social image returns")
demand_name <- paste0(
  "plot-scaled-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-agg-identity-",
  "full-many-pots-pred-demand-vstar-comp-all.pdf"
)
ggsave(file.path(output_path, demand_name), p_demand, width = 8, height = 6)

read_allocation <- function(scenario) {
  if (!is.null(median_allocation_path)) {
    if (parameter_draw != "median") stop("Median allocations require --parameter-draw=median.")
    fit <- readRDS(file.path(median_allocation_path, paste0("cap-3500-", scenario, ".rds")))
    if (!isTRUE(fit$diagnostics$optimal)) stop("Median allocation is not certified.")
    return(fit$allocation)
  }
  filename <- sprintf("replicate-%04d.rds", allocation_draw)
  candidates <- c(
    file.path(policy_path, "median-allocations", scenario, filename),
    file.path(policy_path, "allocations", scenario, filename)
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Missing median allocation for ", scenario)
  readRDS(hit[[1]])$allocation
}
distances <- bind_rows(
  data.frame(distance_km = distance_object$village_df$dist.to.pot / 1000,
             allocation = "Experimental"),
  transform(
    read_allocation("control"), allocation = "Optimal: Control"
  )[, c("distance_km", "allocation")],
  transform(
    read_allocation("bracelet"), allocation = "Optimal: Bracelet"
  )[, c("distance_km", "allocation")]
) |>
  mutate(allocation = factor(
    allocation,
    levels = c("Experimental", "Optimal: Control", "Optimal: Bracelet")
  ))
max_distance_km <- max(distances$distance_km, na.rm = TRUE)
p_distance <- ggplot(distances, aes(distance_km, fill = allocation)) +
  geom_density(colour = "black", alpha = 0.5) +
  annotate(
    "text", x = 0.5, y = 0.6, label = "Amplification",
    size = 5, alpha = 0.7
  ) +
  annotate(
    "text", x = 3, y = 0.25, label = "Mitigation",
    size = 5, alpha = 0.7
  ) +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, max_distance_km + 0.1)) +
  labs(x = "Distance Walked (km)", y = "Density", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
distance_name <- paste0(
  "comp-dist-plot3-fit105-util-identity-",
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.pdf"
)
write.csv(distances, file.path(output_path, "community-distance-figure-inputs.csv"), row.names = FALSE)
write.csv(curves, file.path(output_path, "median-demand-figure-inputs.csv"), row.names = FALSE)
ggsave(file.path(output_path, distance_name), p_distance, width = 8, height = 6)

resolve_allocation_geometry <- function(allocation) {
  # Candidate-site row numbers are local to the distance object used for an
  # optimization run.  Returned allocations can therefore not safely join
  # pot_j to a separately regenerated distance object's pot_df$id.  Recover
  # the site from the village and the stored edge distance, which jointly
  # identify the actual optimized edge.  Coincident duplicate sites are
  # harmless, but must have identical coordinates.
  edges <- distance_object$long_distance_mat[, c(
    "index_i", "index_j", "dist"
  )]
  allocation$.allocation_row <- seq_len(nrow(allocation))
  candidates <- merge(
    allocation, edges, by.x = "village_i", by.y = "index_i",
    suffixes = c("_allocation", "_candidate")
  )
  candidates$distance_gap <- abs(
    candidates$distance - candidates$dist
  )
  best_gap <- ave(
    candidates$distance_gap, candidates$.allocation_row, FUN = min
  )
  matched <- candidates[candidates$distance_gap == best_gap, ]

  pot_xy <- distance_object$pot_df[, c(
    "id", "lon", "lat"
  )]
  matched <- merge(
    matched, pot_xy, by.x = "index_j", by.y = "id", all.x = TRUE
  )
  coordinate_counts <- aggregate(
    cbind(lon, lat) ~ .allocation_row, matched,
    function(value) length(unique(value))
  )
  if (any(coordinate_counts$lon != 1L | coordinate_counts$lat != 1L)) {
    stop("Stored allocation distance maps to multiple candidate locations.")
  }
  matched <- matched[!duplicated(matched$.allocation_row), ]
  matched <- matched[order(matched$.allocation_row), ]
  if (nrow(matched) != nrow(allocation) ||
      max(matched$distance_gap) > 1e-4 ||
      any(!is.finite(matched$lon) | !is.finite(matched$lat))) {
    stop("Could not recover candidate-site geometry from allocation edges.")
  }

  village_xy <- distance_object$village_df[, c(
    "id", "lon", "lat"
  )]
  names(village_xy) <- c("village_i", "village_lon", "village_lat")
  matched <- merge(matched, village_xy, by = "village_i", all.x = TRUE)
  matched[order(matched$.allocation_row), ]
}

allocation_panel <- function(allocation, title) {
  village <- distance_object$village_df
  pot <- distance_object$pot_df
  line <- resolve_allocation_geometry(allocation)
  used <- unique(line[, c("index_j", "lon", "lat")])
  ggplot() +
    geom_point(data = pot, aes(lon, lat), colour = "#9DB7D7",
               shape = 17, size = .65, alpha = .65) +
    geom_segment(data = line, aes(village_lon, village_lat,
                                  xend = lon, yend = lat),
                 linewidth = .25, colour = "black") +
    geom_point(data = used, aes(lon, lat), colour = "#F04B36",
               shape = 17, size = 1.05) +
    geom_point(data = village, aes(lon, lat), colour = "black", size = .5) +
    coord_equal() + theme_void() + ggtitle(title) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("Assigned PoTs: ", nrow(used)),
      hjust = 1.1, vjust = 1.4, size = 3
    )
}
experimental_allocation <- data.frame(
  village_i = distance_object$village_df$id,
  pot_j = distance_object$village_df$id,
  distance = distance_object$village_df$dist.to.pot,
  distance_km = distance_object$village_df$dist.to.pot / 1000
)
# Experimental PoT identifiers need not equal candidate-site identifiers. Use
# the recorded experimental coordinates directly for this panel.
experimental_panel <- ggplot(distance_object$village_df) +
  geom_point(data = distance_object$pot_df, aes(lon, lat),
             colour = "#9DB7D7", shape = 17, size = .65, alpha = .65) +
  geom_segment(aes(lon, lat, xend = pot.lon, yend = pot.lat),
               linewidth = .25, colour = "black") +
  geom_point(aes(pot.lon, pot.lat), colour = "#F04B36",
             shape = 17, size = 1.05) +
  geom_point(aes(lon, lat), colour = "black", size = .5) +
  coord_equal() + theme_void() + ggtitle("Experimental") +
  annotate(
    "text", x = Inf, y = Inf, label = "Assigned PoTs: 144",
    hjust = 1.1, vjust = 1.4, size = 3
  )
control_panel <- allocation_panel(read_allocation("control"), "Control")
bracelet_panel <- allocation_panel(read_allocation("bracelet"), "Bracelet")
if (!requireNamespace("cowplot", quietly = TRUE)) stop("cowplot is required for policy panels.")
panel <- cowplot::plot_grid(experimental_panel, control_panel, bracelet_panel, nrow = 1)
panel_name <- "panel-scenarios-compare-optimal-allocation-plot-distconstraint-3500.pdf"
ggsave(file.path(output_path, panel_name), panel, width = 11, height = 4)

parameter_basis <- if (parameter_draw == "median" && nrow(parameters) > 1L) {
  paste0("componentwise median of ", nrow(parameters), " parameter draws")
} else if (nrow(parameters) == 1L && parameter_draw == "1") {
  "componentwise posterior median"
} else {
  paste0("parameter draw ", parameter_draw)
}
write.csv(data.frame(artifact = c(demand_name, distance_name, panel_name),
                     parameter_basis = parameter_basis,
                     allocation_basis = if (!is.null(median_allocation_path)) "componentwise posterior median, adult-weighted optimum" else paste0("representative draw ", allocation_draw),
                     distance_density_weighting = "equal-community",
                     source_parameter_csv = normalizePath(parameter_csv),
                     source_policy_path = normalizePath(policy_path)),
          file.path(output_path, "policy-figure-manifest.csv"), row.names = FALSE)
message("Wrote policy paper figures to ", output_path)
