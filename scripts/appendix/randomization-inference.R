#!/usr/bin/env Rscript

# Conditional randomization inference for the four item arms.  This script uses
# base R deliberately: it can run on a login node without restoring the full
# historical package library.

args <- commandArgs(trailingOnly = TRUE)
source("R/distance/spec.R")
arg_value <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}

B <- as.integer(arg_value("--permutations", "99999"))
seed <- as.integer(arg_value("--seed", "20260812"))
cores <- as.integer(arg_value("--cores", "4"))
out_dir <- arg_value("--output-dir", "ref-reports/randomization-inference")
appendix_dir <- arg_value("--appendix-dir", "appendix/structural-robustness")
if (is.na(B) || B < 99L) stop("--permutations must be at least 99")
if (is.na(seed)) stop("--seed must be an integer")
if (is.na(cores) || cores < 1L) stop("--cores must be a positive integer")

dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(appendix_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(appendix_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
clean_id <- function(x) as.integer(trimws(as.character(x)))
arms <- c("control", "calendar", "ink", "bracelet")
arm_labels <- c(control = "Control", calendar = "Calendar", ink = "Ink", bracelet = "Bracelet")

message("Loading canonical take-up and observability samples")
takeup <- read.csv("temp-data/analysis-cluster-recentered-covariate-data.csv",
                   check.names = FALSE)
takeup <- data.frame(
  person_id = takeup$KEY.individ,
  cluster_id = clean_id(takeup$cluster.id.x),
  county = factor(takeup$county),
  arm = factor(takeup$assigned.treatment, levels = arms),
  distance_km = takeup$cluster.dist.to.pot / 1000,
  y = as.numeric(as.logical(takeup$dewormed)),
  female = as.numeric(takeup$female),
  age = takeup$age.census,
  mu_d = takeup$mu_d
)

distance_groups <- takeup_distance_crosswalk() |>
  dplyr::filter(.data$in_main_analysis) |>
  dplyr::transmute(
    cluster_id = .data$cluster.id,
    original_far = as.integer(.data$assigned_dist_group == "far"),
    final_far = as.integer(.data$realized_dist_group == "far")
  ) |>
  as.data.frame()
takeup <- merge(takeup, distance_groups, by = "cluster_id", all.x = TRUE, sort = FALSE)

endline <- readRDS("data/clean-data/clean-endline-data.rds")
knowledge <- readRDS("data/clean-data/clean-endline-know-table-data.rds")
knowledge <- knowledge[as.character(knowledge$know.table.type) == "table.A" &
                         knowledge$obs_know_person > 0, ]
obs <- merge(endline, knowledge, by = "KEY.individ", all = FALSE, sort = FALSE)
obs <- obs[as.character(obs$sms.treatment) == "sms.control", ]

full <- readRDS("data/clean-data/full-takeup-data.rds")
person_cov <- unique(data.frame(
  person_id = full$KEY.individ,
  female_cov = as.numeric(full$gender == "female"),
  age_cov = full$age.census
))
cluster_cov <- unique(takeup[c("cluster_id", "mu_d", "distance_km")])
obs <- merge(obs, person_cov, by.x = "KEY.individ", by.y = "person_id", all.x = TRUE, sort = FALSE)
obs$female <- obs$female_cov
obs$age <- obs$age_cov
obs$cluster_id <- clean_id(obs$cluster.id)
obs <- merge(obs, cluster_cov, by = "cluster_id", all.x = TRUE, sort = FALSE)
obs <- merge(obs, distance_groups, by = "cluster_id", all.x = TRUE, sort = FALSE)
obs$county <- factor(obs$county)
obs$arm <- factor(obs$assigned.treatment, levels = arms)
obs$definite <- obs$knows_other_dewormed / obs$obs_know_person
obs$correct <- obs$pct_correct_classification_yesnodk
obs$own_observed <- obs$thinks_other_knows / obs$obs_know_person

assert(nrow(takeup) == 9805L, "Take-up sample is not 9,805 adults")
assert(length(unique(takeup$cluster_id)) == 144L, "Take-up sample is not 144 communities")
switch_count <- sum((takeup$original_far != takeup$final_far) /
                    ave(rep(1, nrow(takeup)), takeup$cluster_id, FUN = sum))
assert(abs(switch_count - 26) < 1e-8,
       "Original/final distance classifications do not contain 26 community switches")
assert(!anyNA(takeup), "Take-up analysis fields contain missing values")

cluster_frame <- unique(takeup[c("cluster_id", "county", "arm", "final_far", "original_far", "distance_km")])
cluster_frame <- cluster_frame[order(cluster_frame$cluster_id), ]
cluster_frame$stratum <- interaction(cluster_frame$county, cluster_frame$original_far, drop = TRUE)
assert(nrow(cluster_frame) == 144L, "Cluster fields are not constant within community")

stratum_counts <- as.data.frame.matrix(table(cluster_frame$stratum, cluster_frame$arm))
stratum_counts$stratum <- rownames(stratum_counts)
rownames(stratum_counts) <- NULL
write.csv(stratum_counts, file.path(out_dir, "data", "assignment-stratum-counts.csv"), row.names = FALSE)

set.seed(seed)
observed_assignment <- as.character(cluster_frame$arm)
permutations <- matrix(NA_integer_, nrow = B, ncol = nrow(cluster_frame))
for (b in seq_len(B)) {
  perm_arm <- observed_assignment
  for (s in levels(cluster_frame$stratum)) {
    ii <- which(cluster_frame$stratum == s)
    perm_arm[ii] <- sample(observed_assignment[ii], length(ii), replace = FALSE)
  }
  permutations[b, ] <- match(perm_arm, arms)
}

# Build cluster sufficient statistics for X = [arm intercepts, arm slopes,
# female, age, expected distance, county indicators].  Treatment regressors
# are constant within community, so these statistics exactly recover the
# individual-level OLS fit and cluster sandwich without retaining X each draw.
prepare_problem <- function(data, outcome, moderator, outcome_label, moderator_label) {
  keep <- complete.cases(data[, c(outcome, moderator, "female", "age", "mu_d", "county",
                                  "cluster_id", "arm")])
  d <- data[keep, ]
  d$y_work <- as.numeric(d[[outcome]])
  d$m_work <- as.numeric(d[[moderator]])
  # Arm indicators supply the intercept, so omit the reference county to keep
  # the saturated design full rank (the fixed-effect fit is otherwise equal).
  county_mm <- model.matrix(~ county, d)[, -1, drop = FALSE]
  W <- cbind(female = d$female, age = d$age, mu_d = d$mu_d, county_mm)
  cluster_ids <- cluster_frame$cluster_id
  K <- length(cluster_ids)
  q <- ncol(W)
  n_g <- numeric(K); sy_g <- numeric(K); sy2_g <- numeric(K)
  sw_g <- matrix(0, K, q); swy_g <- matrix(0, K, q)
  sww_g <- array(0, c(q, q, K)); moderator_g <- numeric(K)
  for (g in seq_len(K)) {
    ii <- which(d$cluster_id == cluster_ids[g])
    if (!length(ii)) next
    wg <- W[ii, , drop = FALSE]; yg <- d$y_work[ii]
    n_g[g] <- length(ii); sy_g[g] <- sum(yg); sy2_g[g] <- sum(yg^2)
    sw_g[g, ] <- colSums(wg); swy_g[g, ] <- crossprod(wg, yg)
    sww_g[, , g] <- crossprod(wg)
    moderator_g[g] <- unique(d$m_work[ii])
    assert(length(unique(d$m_work[ii])) == 1L, "Moderator varies within community")
  }
  list(outcome = outcome_label, moderator = moderator_label, n = sum(n_g), G = sum(n_g > 0),
       p = 8L + q, cluster_ids = cluster_ids, n_g = n_g, sy_g = sy_g,
       sy2_g = sy2_g, sw_g = sw_g, swy_g = swy_g, sww_g = sww_g,
       moderator_g = moderator_g)
}

fit_assignment <- function(problem, arm_index) {
  p <- problem$p; qn <- p - 8L
  bread <- matrix(0, p, p); xy <- numeric(p)
  cluster_xx <- vector("list", length(problem$cluster_ids))
  cluster_xy <- vector("list", length(problem$cluster_ids))
  for (g in seq_along(problem$cluster_ids)) {
    if (problem$n_g[g] == 0) next
    q <- numeric(8L); a <- arm_index[g]
    q[a] <- 1; q[4L + a] <- problem$moderator_g[g]
    xx <- rbind(cbind(problem$n_g[g] * tcrossprod(q), tcrossprod(q, problem$sw_g[g, ])),
                cbind(tcrossprod(problem$sw_g[g, ], q), problem$sww_g[, , g]))
    xyg <- c(q * problem$sy_g[g], problem$swy_g[g, ])
    bread <- bread + xx; xy <- xy + xyg
    cluster_xx[[g]] <- xx; cluster_xy[[g]] <- xyg
  }
  beta <- tryCatch(solve(bread, xy), error = function(e) qr.solve(bread, xy, tol = 1e-10))
  invbread <- tryCatch(solve(bread), error = function(e) qr.solve(bread, diag(p), tol = 1e-10))
  meat <- matrix(0, p, p)
  for (g in seq_along(cluster_xx)) {
    if (is.null(cluster_xx[[g]])) next
    score <- cluster_xy[[g]] - cluster_xx[[g]] %*% beta
    meat <- meat + tcrossprod(score)
  }
  correction <- problem$G / (problem$G - 1) * (problem$n - 1) / (problem$n - p)
  vcov <- correction * invbread %*% meat %*% invbread
  list(beta = beta, vcov = vcov)
}

contrast_matrix <- rbind(
  any_signal = c(-0.5, -0.5, 0.5, 0.5),
  bracelet_calendar = c(0, -1, 0, 1)
)

evaluate_problem <- function(problem) {
  C <- matrix(0, 2, problem$p)
  C[, 5:8] <- contrast_matrix
  observed <- fit_assignment(problem, match(observed_assignment, arms))
  est <- as.numeric(C %*% observed$beta)
  V <- C %*% observed$vcov %*% t(C)
  se <- sqrt(pmax(diag(V), 0))
  t_obs <- est / se
  sharp_draws <- matrix(NA_real_, B, 2)
  t_draws <- matrix(NA_real_, B, 2)
  omnibus_draw <- numeric(B)

  H <- matrix(0, 3, problem$p)
  H[1, 5:8] <- c(-1, 1, 0, 0)
  H[2, 5:8] <- c(-1, 0, 1, 0)
  H[3, 5:8] <- c(-1, 0, 0, 1)
  hobs <- H %*% observed$beta
  HVH <- H %*% observed$vcov %*% t(H)
  omnibus_obs <- as.numeric(t(hobs) %*% qr.solve(HVH, hobs))

  for (b in seq_len(B)) {
    fit <- fit_assignment(problem, permutations[b, ])
    eb <- as.numeric(C %*% fit$beta)
    seb <- sqrt(pmax(diag(C %*% fit$vcov %*% t(C)), 0))
    sharp_draws[b, ] <- eb
    t_draws[b, ] <- eb / seb
    hb <- H %*% fit$beta
    hvh <- H %*% fit$vcov %*% t(H)
    omnibus_draw[b] <- as.numeric(t(hb) %*% qr.solve(hvh, hb))
  }
  result <- data.frame(
    outcome = problem$outcome, distance_specification = problem$moderator,
    contrast = c("Any Signal - No Signal", "Bracelet - Calendar"),
    estimate = est, std_error = se, conf_low = est - 1.96 * se,
    conf_high = est + 1.96 * se,
    fisher_p = (1 + colSums(abs(sharp_draws) >= rep(abs(est), each = B))) / (B + 1),
    studentized_p = (1 + colSums(abs(t_draws) >= rep(abs(t_obs), each = B))) / (B + 1),
    stringsAsFactors = FALSE
  )
  omnibus <- data.frame(
    outcome = problem$outcome, distance_specification = problem$moderator,
    statistic = omnibus_obs,
    randomization_p = (1 + sum(omnibus_draw >= omnibus_obs)) / (B + 1)
  )
  list(summary = result, t_draws = t_draws, t_obs = t_obs, omnibus = omnibus)
}

problems <- list(
  prepare_problem(takeup, "y", "original_far", "Administrative take-up", "Original assigned Close/Far"),
  prepare_problem(takeup, "y", "final_far", "Administrative take-up", "Final realized Close/Far"),
  prepare_problem(takeup, "y", "distance_km", "Administrative take-up", "Realized distance (per km)"),
  prepare_problem(obs, "definite", "original_far", "Definite peer classification", "Original assigned Close/Far"),
  prepare_problem(obs, "definite", "final_far", "Definite peer classification", "Final realized Close/Far"),
  prepare_problem(obs, "definite", "distance_km", "Definite peer classification", "Realized distance (per km)"),
  prepare_problem(obs, "correct", "original_far", "Correct classification", "Original assigned Close/Far"),
  prepare_problem(obs, "correct", "final_far", "Correct classification", "Final realized Close/Far"),
  prepare_problem(obs, "correct", "distance_km", "Correct classification", "Realized distance (per km)"),
  prepare_problem(obs, "own_observed", "original_far", "Perceived own observability", "Original assigned Close/Far"),
  prepare_problem(obs, "own_observed", "final_far", "Perceived own observability", "Final realized Close/Far"),
  prepare_problem(obs, "own_observed", "distance_km", "Perceived own observability", "Realized distance (per km)")
)

message("Running ", B, " conditional arm reallocations for ", length(problems),
        " models on ", cores, " core(s)")
run_problem <- function(i) {
  message("  ", i, "/", length(problems), ": ", problems[[i]]$outcome,
          " — ", problems[[i]]$moderator)
  evaluate_problem(problems[[i]])
}
if (.Platform$OS.type == "unix" && cores > 1L) {
  evaluated <- parallel::mclapply(seq_along(problems), run_problem,
                                  mc.cores = min(cores, length(problems)),
                                  mc.preschedule = TRUE)
} else {
  evaluated <- lapply(seq_along(problems), run_problem)
}

summary <- do.call(rbind, lapply(evaluated, `[[`, "summary"))
omnibus <- do.call(rbind, lapply(evaluated, `[[`, "omnibus"))

romano_wolf <- function(indices) {
  obs <- unlist(lapply(evaluated[indices], `[[`, "t_obs"))
  draws <- do.call(cbind, lapply(evaluated[indices], `[[`, "t_draws"))
  ord <- order(abs(obs), decreasing = TRUE)
  adjusted <- numeric(length(obs)); running <- 0
  for (k in seq_along(ord)) {
    remaining <- ord[k:length(ord)]
    p <- (1 + sum(apply(abs(draws[, remaining, drop = FALSE]), 1, max) >= abs(obs[ord[k]]))) / (B + 1)
    running <- max(running, p); adjusted[ord[k]] <- running
  }
  adjusted
}

# Problems 1 and 4/7/10 are primary original-assignment families.  Realized
# binary and continuous models form separate secondary families by domain.
takeup_primary <- 1L
obs_primary <- c(4L, 7L, 10L)
takeup_secondary <- c(2L, 3L)
obs_secondary <- c(5L, 6L, 8L, 9L, 11L, 12L)

summary$romano_wolf_p <- NA_real_
row_offset <- function(problem_indices) unlist(lapply(problem_indices, function(i) (2L * i - 1L):(2L * i)))
summary$romano_wolf_p[row_offset(takeup_primary)] <- romano_wolf(takeup_primary)
summary$romano_wolf_p[row_offset(obs_primary)] <- romano_wolf(obs_primary)
summary$romano_wolf_p[row_offset(takeup_secondary)] <- romano_wolf(takeup_secondary)
summary$romano_wolf_p[row_offset(obs_secondary)] <- romano_wolf(obs_secondary)
summary$mcse_studentized <- sqrt(summary$studentized_p * (1 - summary$studentized_p) / (B + 1))

write.csv(summary, file.path(out_dir, "data", "ri-results.csv"), row.names = FALSE)
write.csv(omnibus, file.path(out_dir, "data", "ri-omnibus-results.csv"), row.names = FALSE)

sample_audit <- data.frame(
  outcome = vapply(problems[c(1, 4, 7, 10)], `[[`, "", "outcome"),
  observations = vapply(problems[c(1, 4, 7, 10)], `[[`, 0, "n"),
  communities = vapply(problems[c(1, 4, 7, 10)], `[[`, 0, "G")
)
write.csv(sample_audit, file.path(out_dir, "data", "sample-audit.csv"), row.names = FALSE)

permutation_audit <- data.frame(
  permutations = B, seed = seed, cores = cores, communities = nrow(cluster_frame),
  strata = nlevels(cluster_frame$stratum), switches = 26,
  stratum_definition = "county_x_original_assigned_close_far",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write.csv(permutation_audit, file.path(out_dir, "data", "permutation-audit.csv"), row.names = FALSE)

fmt <- function(x) sprintf("%.3f", x)
tex_escape <- function(x) gsub("%", "\\\\%", x, fixed = TRUE)
make_table <- function(path) {
  con <- file(path, "w"); on.exit(close(con))
  writeLines(c("\\begin{tabular}{llrrrrrr}", "\\toprule",
    "Outcome & Contrast & Estimate & 95\\% CI & Fisher $p$ & Stud. RI $p$ & RW $p$ & $N$ \\\\",
    "\\midrule"), con)
  specs <- unique(summary$distance_specification)
  for (s in specs) {
    writeLines(paste0("\\multicolumn{8}{l}{\\textit{", tex_escape(s), "}} \\\\"), con)
    ss <- summary[summary$distance_specification == s, ]
    for (i in seq_len(nrow(ss))) {
      n_i <- problems[[which(vapply(problems, function(z) z$outcome == ss$outcome[i] && z$moderator == s, FALSE))]]$n
      line <- paste(tex_escape(ss$outcome[i]), tex_escape(ss$contrast[i]), fmt(ss$estimate[i]),
                    paste0("[", fmt(ss$conf_low[i]), ", ", fmt(ss$conf_high[i]), "]"),
                    fmt(ss$fisher_p[i]), fmt(ss$studentized_p[i]), fmt(ss$romano_wolf_p[i]),
                    format(n_i, big.mark = ","), sep = " & ")
      writeLines(paste0(line, " \\\\"), con)
    }
    if (s != tail(specs, 1)) writeLines("\\addlinespace", con)
  }
  writeLines(c("\\bottomrule", "\\end{tabular}"), con)
}

table_path <- file.path(out_dir, "tables", "randomization-inference.tex")
make_table(table_path)
file.copy(table_path, file.path(appendix_dir, "tables", "randomization-inference.tex"), overwrite = TRUE)

# Community-level diagnostic figures use base graphics for portability.
cluster_plot <- aggregate(takeup$y, list(cluster_id = takeup$cluster_id), mean)
names(cluster_plot)[2] <- "outcome"
cluster_plot <- merge(cluster_plot, cluster_frame, by = "cluster_id")
plot_one <- function(path, distance, xlab) {
  pdf(path, width = 8, height = 5)
  cols <- c("#333333", "#D89000", "#0072B2", "#CC79A7")
  plot(cluster_plot[[distance]], cluster_plot$outcome, type = "n", xlab = xlab,
       ylab = "Community mean take-up", ylim = c(0, 1))
  for (a in seq_along(arms)) {
    ii <- cluster_plot$arm == arms[a]
    points(jitter(cluster_plot[[distance]][ii], amount = if (distance == "distance_km") 0 else .035),
           cluster_plot$outcome[ii], pch = 16, col = adjustcolor(cols[a], .65))
    fit <- lm(outcome ~ cluster_plot[[distance]], data = cluster_plot, subset = ii,
              weights = ave(rep(1, nrow(takeup)), takeup$cluster_id, FUN = length)[match(cluster_plot$cluster_id, takeup$cluster_id)])
    abline(fit, col = cols[a], lwd = 2)
  }
  legend("topright", legend = arm_labels[arms], col = cols, pch = 16, lty = 1, bty = "n")
  dev.off()
}
plot_one(file.path(out_dir, "figures", "takeup-original-assignment.pdf"), "original_far", "Original assigned distance (Close/Far)")
plot_one(file.path(out_dir, "figures", "takeup-realized-distance.pdf"), "distance_km", "Realized centroid-to-PoT distance (km)")
file.copy(file.path(out_dir, "figures", "takeup-original-assignment.pdf"),
          file.path(appendix_dir, "figures", "ri-takeup-original-assignment.pdf"), overwrite = TRUE)
file.copy(file.path(out_dir, "figures", "takeup-realized-distance.pdf"),
          file.path(appendix_dir, "figures", "ri-takeup-realized-distance.pdf"), overwrite = TRUE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "data", "session-info.txt"))
message("Finished. Results written to ", out_dir)
