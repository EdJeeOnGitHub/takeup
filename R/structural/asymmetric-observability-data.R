# Build the linked peer-response sufficient statistics used by the asymmetric
# observability model. This file deliberately uses base R so it can run before
# the main sampling script attaches tidyverse packages.

main_core_empty_peer_response_data <- function() {
  list(
    core_num_peer_response_rows = 0L,
    core_peer_response_cluster_id = integer(),
    core_peer_true_takeup = integer(),
    core_peer_total = integer(),
    core_peer_recognized = integer(),
    core_peer_report_count = matrix(integer(), nrow = 0L, ncol = 3L),
    core_peer_link_audit = matrix(numeric(), nrow = 0L, ncol = 6L)
  )
}

main_core_prepare_peer_response_data <- function(
    sample_data,
    project_root = ".",
    write_audit = NULL) {
  required <- c(
    "analysis_data", "beliefs_obs_index", "obs_cluster_id",
    "num_beliefs_obs", "num_clusters", "cluster_assigned_treatment",
    "cluster_assigned_dist_group", "num_recognized", "num_knows_1ord"
  )
  missing_required <- setdiff(required, names(sample_data))
  if (length(missing_required)) {
    stop("Asymmetric observability data require: ",
         paste(missing_required, collapse = ", "), call. = FALSE)
  }
  analysis <- sample_data$analysis_data
  needed_analysis <- c("KEY.individ", "dewormed.any")
  if (!all(needed_analysis %in% names(analysis))) {
    stop("analysis_data lacks KEY.individ or dewormed.any.", call. = FALSE)
  }
  if (length(sample_data$beliefs_obs_index) != sample_data$num_beliefs_obs) {
    stop("Belief row indices are inconsistent.", call. = FALSE)
  }
  beliefs_cluster_index <- sample_data$obs_cluster_id[
    sample_data$beliefs_obs_index
  ]
  if (length(beliefs_cluster_index) != sample_data$num_beliefs_obs ||
      any(!beliefs_cluster_index %in% seq_len(sample_data$num_clusters))) {
    stop("Could not derive valid cluster indices for belief respondents.",
         call. = FALSE)
  }

  peer_path <- file.path(
    project_root, "data", "clean-data",
    "clean-endline-know-table-data-long.rds"
  )
  if (!file.exists(peer_path)) {
    stop("Missing asymmetric-observability input: ",
         peer_path,
         call. = FALSE)
  }
  peer <- readRDS(peer_path)
  peer_required <- c(
    "know.table.type", "recognized", "dewormed", "KEY.individ",
    "cluster.id", "sms.treatment", "KEY.individ.other.1",
    "actual.other.dewormed.any.1"
  )
  if (!all(peer_required %in% names(peer))) {
    stop("Malformed clean Table-A peer data.", call. = FALSE)
  }
  peer <- peer[as.character(peer$know.table.type) == "table.A", ]
  report_text <- as.character(peer$dewormed)
  table_a <- data.frame(
    respondent_key = peer$KEY.individ,
    original_cluster_id = peer$cluster.id,
    sms_treatment = as.character(peer$sms.treatment),
    peer_key = peer$KEY.individ.other.1,
    recognized = as.integer(as.character(peer$recognized) == "yes"),
    report_raw = ifelse(
      report_text == "yes", 1L,
      ifelse(report_text == "no", 0L,
             ifelse(report_text == "don't know", 98L, NA_integer_))
    ),
    true_takeup = as.integer(peer$actual.other.dewormed.any.1),
    stringsAsFactors = FALSE
  )

  # beliefs_obs_index is a cluster anchor, not a respondent index. Recover the
  # original cluster label attached to each active Stan cluster, then rebuild
  # the FOB/drop respondent sample from the endline microdata.
  if (!"cluster.id" %in% names(analysis)) {
    stop("analysis_data lacks the original cluster.id.", call. = FALSE)
  }
  active_original_cluster <- vapply(seq_len(sample_data$num_clusters), function(cluster) {
    value <- unique(analysis$cluster.id[sample_data$obs_cluster_id == cluster])
    value <- value[!is.na(value)]
    if (length(value) != 1L) {
      stop("Active Stan clusters do not map one-to-one to original clusters.",
           call. = FALSE)
    }
    as.character(value)
  }, character(1L))
  table_a$cluster_id <- match(
    as.character(table_a$original_cluster_id), active_original_cluster
  )
  table_a <- table_a[
    table_a$sms_treatment == "sms.control" & !is.na(table_a$cluster_id),
  ]

  respondent_recognized <- rowsum(
    table_a$recognized, table_a$respondent_key, reorder = FALSE
  )[, 1L]
  included_respondents <- names(respondent_recognized)[respondent_recognized > 0]
  table_a <- table_a[table_a$respondent_key %in% included_respondents, ]

  if (nrow(table_a) != sample_data$num_beliefs_obs *
      sample_data$know_table_A_sample_size) {
    stop("Peer rows do not reproduce the structural Table-A sample size: got ",
         nrow(table_a), ", expected ",
         sample_data$num_beliefs_obs * sample_data$know_table_A_sample_size,
         ".",
         call. = FALSE)
  }
  if (any(!table_a$recognized %in% 0:1)) {
    stop("Recognition must be coded 0/1.", call. = FALSE)
  }
  recognized_report <- table_a$report_raw[table_a$recognized == 1L]
  if (any(!recognized_report %in% c(0L, 1L, 98L))) {
    stop("Recognized-peer reports must be Yes, No, or Don't Know.",
         call. = FALSE)
  }
  definite <- as.integer(table_a$recognized == 1L &
                         table_a$report_raw %in% c(0L, 1L))
  respondent_key <- interaction(
    table_a$cluster_id, table_a$respondent_key, drop = TRUE
  )
  reproduced <- cbind(
    cluster = as.integer(rowsum(table_a$cluster_id, respondent_key)[, 1L] /
      rowsum(rep.int(1L, nrow(table_a)), respondent_key)[, 1L]),
    recognized = as.integer(rowsum(table_a$recognized, respondent_key)[, 1L]),
    definite = as.integer(rowsum(definite, respondent_key)[, 1L])
  )
  frozen <- cbind(
    cluster = beliefs_cluster_index,
    recognized = sample_data$num_recognized,
    definite = sample_data$num_knows_1ord
  )
  storage.mode(reproduced) <- "integer"
  storage.mode(frozen) <- "integer"
  row_order <- function(x) do.call(order, as.data.frame(x))
  if (!identical(unname(reproduced[row_order(reproduced), , drop = FALSE]),
                 unname(frozen[row_order(frozen), , drop = FALSE]))) {
    stop("Peer records do not reproduce the frozen belief sufficient statistics.",
         call. = FALSE)
  }

  truth_raw <- table_a$true_takeup
  if (is.logical(truth_raw)) {
    truth <- as.integer(truth_raw)
  } else if (is.factor(truth_raw)) {
    truth_text <- tolower(as.character(truth_raw))
    truth <- ifelse(truth_text %in% c("true", "yes", "1"), 1L,
                    ifelse(truth_text %in% c("false", "no", "0"), 0L, NA_integer_))
  } else {
    truth <- suppressWarnings(as.integer(truth_raw))
  }
  truth[!truth %in% 0:1] <- NA_integer_
  table_a$true_takeup <- truth
  table_a$linked_truth <- !is.na(truth)

  treatment <- sample_data$cluster_assigned_treatment[table_a$cluster_id]
  dist_group <- sample_data$cluster_assigned_dist_group[table_a$cluster_id]
  audit_key <- paste(treatment, dist_group, sep = "\r")
  audit_split <- split(seq_len(nrow(table_a)), audit_key)
  audit <- do.call(rbind, lapply(audit_split, function(index) {
    c(
      treatment = treatment[index[1L]],
      distance_group = dist_group[index[1L]],
      peer_rows = length(index),
      linked_rows = sum(table_a$linked_truth[index]),
      linked_share = mean(table_a$linked_truth[index]),
      recognized_share = mean(table_a$recognized[index])
    )
  }))
  audit <- unname(audit)
  colnames(audit) <- c(
    "treatment", "distance_group", "peer_rows", "linked_rows",
    "linked_share", "recognized_share"
  )

  linked <- table_a[table_a$linked_truth, ]
  aggregate_key <- paste(linked$cluster_id, linked$true_takeup, sep = "\r")
  aggregate_split <- split(seq_len(nrow(linked)), aggregate_key)
  aggregate_rows <- lapply(aggregate_split, function(index) {
    report <- linked$report_raw[index]
    recognized <- linked$recognized[index]
    c(
      cluster_id = linked$cluster_id[index[1L]],
      true_takeup = linked$true_takeup[index[1L]],
      total = length(index),
      recognized = sum(recognized),
      report_yes = sum(recognized == 1L & report == 1L),
      report_no = sum(recognized == 1L & report == 0L),
      report_dk = sum(recognized == 1L & report == 98L)
    )
  })
  aggregate <- do.call(rbind, aggregate_rows)
  storage.mode(aggregate) <- "integer"
  if (any(aggregate[, "recognized"] != rowSums(
      aggregate[, c("report_yes", "report_no", "report_dk"), drop = FALSE]
  ))) {
    stop("Conditional report counts do not partition recognized peers.",
         call. = FALSE)
  }

  if (!is.null(write_audit)) {
    dir.create(dirname(write_audit), recursive = TRUE, showWarnings = FALSE)
    write.csv(as.data.frame(audit), write_audit, row.names = FALSE)
  }
  list(
    core_num_peer_response_rows = nrow(aggregate),
    core_peer_response_cluster_id = aggregate[, "cluster_id"],
    core_peer_true_takeup = aggregate[, "true_takeup"],
    core_peer_total = aggregate[, "total"],
    core_peer_recognized = aggregate[, "recognized"],
    core_peer_report_count = aggregate[, c(
      "report_yes", "report_no", "report_dk"
    ), drop = FALSE],
    core_peer_link_audit = audit
  )
}
