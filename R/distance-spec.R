# Canonical Close/Far definitions shared by the replication pipeline.

takeup_distance_levels <- c("close", "far")

takeup_distance_spec <- function(value = Sys.getenv("TAKEUP_DISTANCE_SPEC", "assigned")) {
  value <- tolower(trimws(value))
  if (!value %in% c("assigned", "realized")) {
    stop(
      "Distance definition must be 'assigned' or 'realized'; got: ", value,
      call. = FALSE
    )
  }
  value
}

takeup_clean_cluster_id <- function(x, label = "cluster ID") {
  value <- trimws(as.character(x))
  bad <- is.na(value) | value == "" | !grepl("^[0-9]+$", value)
  if (any(bad)) {
    stop(label, " contains missing or non-integer values.", call. = FALSE)
  }
  as.integer(value)
}

takeup_distance_crosswalk <- function(
    assigned_path = file.path("data", "rct_targetable_schools_2.0.rds"),
    realized_path = file.path("data", "takeup_processed_cluster_strat.rds"),
    analysis_path = file.path("data", "clean-data",
                              "monitored-nosms-takeup-data.rds"),
    validate_counts = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  assigned_source <- readRDS(assigned_path)
  realized_source <- readRDS(realized_path)
  analysis_source <- readRDS(analysis_path)
  analysis_clusters <- dplyr::as_tibble(analysis_source) |>
    dplyr::transmute(
      cluster.id = takeup_clean_cluster_id(.data$cluster.id, "analysis cluster.id")
    ) |>
    dplyr::distinct()

  assigned <- dplyr::as_tibble(assigned_source) |>
    dplyr::transmute(
      cluster.id = takeup_clean_cluster_id(.data$pot.cluster.id, "pot.cluster.id"),
      assigned_dist_group = factor(.data$assigned.dist.cat,
                                   levels = takeup_distance_levels)
    ) |>
    dplyr::distinct()

  realized <- dplyr::as_tibble(realized_source) |>
    dplyr::transmute(
      cluster.id = takeup_clean_cluster_id(.data$cluster.id, "cluster.id"),
      realized_dist_group = factor(.data$dist.pot.group,
                                   levels = takeup_distance_levels),
      assigned_treatment = as.character(.data$assigned.treatment),
      realized_cluster_distance_m = .data$dist.to.own.pot
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(in_main_analysis = .data$cluster.id %in% analysis_clusters$cluster.id)

  if (anyDuplicated(assigned$cluster.id) || anyDuplicated(realized$cluster.id)) {
    stop("Distance sources are not unique by cluster ID.", call. = FALSE)
  }

  crosswalk <- dplyr::left_join(realized, assigned, by = "cluster.id") |>
    dplyr::mutate(
      switched = .data$assigned_dist_group != .data$realized_dist_group,
      switch_direction = ifelse(
        .data$switched,
        paste(.data$assigned_dist_group, "to", .data$realized_dist_group),
        "unchanged"
      )
    ) |>
    dplyr::arrange(.data$cluster.id)

  if (anyNA(crosswalk$assigned_dist_group) ||
      anyNA(crosswalk$realized_dist_group)) {
    stop("The distance crosswalk has unmatched or invalid groups.", call. = FALSE)
  }

  if (validate_counts) {
    analysis_crosswalk <- dplyr::filter(crosswalk, .data$in_main_analysis)
    observed <- table(
      factor(analysis_crosswalk$assigned_dist_group, levels = takeup_distance_levels),
      factor(analysis_crosswalk$realized_dist_group, levels = takeup_distance_levels)
    )
    expected <- matrix(c(64L, 10L, 16L, 54L), nrow = 2L, byrow = TRUE)
    if (nrow(analysis_crosswalk) != 144L ||
        !identical(as.integer(observed), as.integer(expected))) {
      stop(
        "Distance crosswalk does not reproduce the documented 144-cluster ",
        "64/10/16/54 assigned-by-realized split.",
        call. = FALSE
      )
    }
  }
  crosswalk
}

takeup_apply_distance_spec <- function(data, crosswalk = NULL,
                                       specification = takeup_distance_spec(),
                                       compatibility_column = TRUE) {
  specification <- takeup_distance_spec(specification)
  cluster_column <- intersect(c("cluster.id", "cluster.id.x"), names(data))[1L]
  if (is.na(cluster_column)) return(data)
  if (is.null(crosswalk)) crosswalk <- takeup_distance_crosswalk()

  data <- dplyr::as_tibble(data)
  if ("assigned_dist_group" %in% names(data)) {
    data$assigned_dist_group <- NULL
  }
  if ("realized_dist_group" %in% names(data)) {
    data$realized_dist_group <- NULL
  }
  raw_cluster_id <- trimws(as.character(data[[cluster_column]]))
  missing_cluster_id <- is.na(raw_cluster_id) | raw_cluster_id == ""
  if (any(!missing_cluster_id & !grepl("^[0-9]+$", raw_cluster_id))) {
    stop("cluster ID contains non-integer values.", call. = FALSE)
  }
  data$.takeup_cluster_id <- NA_integer_
  data$.takeup_cluster_id[!missing_cluster_id] <-
    as.integer(raw_cluster_id[!missing_cluster_id])
  joined <- dplyr::left_join(
    data,
    dplyr::select(
      crosswalk,
      "cluster.id", "assigned_dist_group", "realized_dist_group"
    ),
    by = c(".takeup_cluster_id" = "cluster.id"),
    relationship = "many-to-one"
  )
  matched_rows <- !is.na(joined$.takeup_cluster_id)
  if (any(is.na(joined$assigned_dist_group[matched_rows])) ||
      any(is.na(joined$realized_dist_group[matched_rows]))) {
    stop("Analysis data contain clusters absent from the distance crosswalk.",
         call. = FALSE)
  }
  joined$analysis_dist_group <- factor(
    joined[[paste0(specification, "_dist_group")]],
    levels = takeup_distance_levels
  )
  if (compatibility_column) {
    joined$dist.pot.group <- joined$analysis_dist_group
  }
  joined$.takeup_cluster_id <- NULL
  joined
}

takeup_write_distance_audit <- function(crosswalk, output_dir,
                                        specification = takeup_distance_spec()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    crosswalk,
    file.path(output_dir, "distance-crosswalk.csv"),
    row.names = FALSE
  )
  analysis_crosswalk <- dplyr::filter(crosswalk, .data$in_main_analysis)
  counts <- as.data.frame(table(
    assigned_dist_group = analysis_crosswalk$assigned_dist_group,
    realized_dist_group = analysis_crosswalk$realized_dist_group
  ))
  names(counts)[3L] <- "clusters"
  utils::write.csv(
    counts,
    file.path(output_dir, "distance-switch-counts.csv"),
    row.names = FALSE
  )
  manifest <- data.frame(
    distance_specification = specification,
    clusters = nrow(analysis_crosswalk),
    switched_clusters = sum(analysis_crosswalk$switched),
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  utils::write.csv(
    manifest,
    file.path(output_dir, "distance-manifest.csv"),
    row.names = FALSE
  )
  file.path(output_dir, c(
    "distance-crosswalk.csv", "distance-switch-counts.csv",
    "distance-manifest.csv"
  ))
}
