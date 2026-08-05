# Shared draw-level contrasts for multiplier robustness summaries. Columns must
# follow the canonical Control, Ink, Calendar, Bracelet treatment order.

main_core_multiplier_contrasts <- function(multiplier_matrix) {
  stopifnot(is.matrix(multiplier_matrix), ncol(multiplier_matrix) == 4L)
  colnames(multiplier_matrix) <- c("Control", "Ink", "Calendar", "Bracelet")
  cbind(
    `No Signal - Any Signal` =
      (multiplier_matrix[, "Control"] + multiplier_matrix[, "Calendar"] -
       multiplier_matrix[, "Ink"] - multiplier_matrix[, "Bracelet"]) / 2,
    `Control - Bracelet` =
      multiplier_matrix[, "Control"] - multiplier_matrix[, "Bracelet"],
    `Calendar - Bracelet` =
      multiplier_matrix[, "Calendar"] - multiplier_matrix[, "Bracelet"],
    `Control - Ink` =
      multiplier_matrix[, "Control"] - multiplier_matrix[, "Ink"]
  )
}

main_core_summarize_contrast <- function(value) {
  finite <- value[is.finite(value)]
  if (length(finite) == 0L) {
    return(c(
      median = NA_real_, lower = NA_real_, upper = NA_real_,
      share_positive = NA_real_, finite_share = 0
    ))
  }
  c(
    median = median(finite),
    lower = unname(quantile(finite, 0.025)),
    upper = unname(quantile(finite, 0.975)),
    share_positive = mean(finite > 0),
    finite_share = length(finite) / length(value)
  )
}
