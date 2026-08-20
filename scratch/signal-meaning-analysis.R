#!/usr/bin/env Rscript

# Analyze the coded endline item-meaning responses: does the motivation /
# character content of a public signal's meaning differ between randomized
# Close and Far communities?  Uses the blind two-pass LLM coding produced from
# ref-reports/signal-meaning-coding/codebook.md.  Base R only.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}

B_weights <- as.integer(arg_value("--weight-draws", "1000"))
B_perm <- as.integer(arg_value("--permutations", "9999"))
seed <- as.integer(arg_value("--seed", "20260814"))
in_dir <- arg_value("--input-dir", "ref-reports/signal-meaning-coding")
appendix_dir <- arg_value("--appendix-dir", "appendix/structural-robustness")

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
labels <- c("M", "P", "C", "V", "O")

message("Loading coded responses and respondent data")
respondent <- readRDS(file.path(in_dir, "respondent-level.rds"))
final <- read.csv(file.path(in_dir, "coded-final.csv"))
passA <- read.csv(file.path(in_dir, "coded-passA.csv"))
passB <- read.csv(file.path(in_dir, "coded-passB.csv"))
dedup <- read.csv(file.path(in_dir, "responses-dedup.csv"))

for (d in list(final, passA, passB)) {
  assert(identical(sort(d$response_id), sort(dedup$response_id)),
         "A coding file does not cover exactly the dedup response ids")
  assert(all(rowSums(d[labels]) >= 1), "A coded row has no label")
  assert(all(d$O == 0 | rowSums(d[labels]) == 1), "O must be exclusive")
}

# ---------------------------------------------------------------------------
# Inter-coder agreement between the two independent passes.
passA <- passA[match(dedup$response_id, passA$response_id), ]
passB <- passB[match(dedup$response_id, passB$response_id), ]
kappa_stat <- function(a, b) {
  po <- mean(a == b)
  pe <- mean(a) * mean(b) + mean(1 - a) * mean(1 - b)
  if (abs(1 - pe) < 1e-12) return(NA_real_)
  (po - pe) / (1 - pe)
}
agreement <- data.frame(
  label = labels,
  share_passA = sapply(labels, function(l) mean(passA[[l]])),
  share_passB = sapply(labels, function(l) mean(passB[[l]])),
  agree = sapply(labels, function(l) mean(passA[[l]] == passB[[l]])),
  kappa = sapply(labels, function(l) kappa_stat(passA[[l]], passB[[l]]))
)
exact_match <- mean(rowSums(passA[labels] == passB[labels]) == length(labels))
write.csv(agreement, file.path(in_dir, "coder-agreement.csv"), row.names = FALSE)
message(sprintf("Exact five-label agreement: %.3f", exact_match))
print(agreement, row.names = FALSE)

# ---------------------------------------------------------------------------
# Merge final codes onto respondents.
final <- final[match(dedup$response_id, final$response_id), ]
code_row <- match(respondent$response_id, final$response_id)
for (l in labels) respondent[[paste0("code_", l)]] <- final[[l]][code_row]

respondent$any_M <- ifelse(respondent$has_response == 1L, respondent$code_M, NA)
respondent$P_only <- ifelse(
  respondent$has_response == 1L,
  as.integer(respondent$code_P == 1 & respondent$code_M == 0 &
               respondent$code_C == 0 & respondent$code_V == 0),
  NA)
respondent$any_V <- ifelse(respondent$has_response == 1L, respondent$code_V, NA)
respondent$county <- factor(respondent$county)

# ---------------------------------------------------------------------------
# Estimation helpers.
design <- function(d, far_var) {
  county_mm <- if (nlevels(droplevels(d$county)) > 1)
    model.matrix(~ county, droplevels(d))[, -1, drop = FALSE] else NULL
  item_mm <- if (length(unique(d$item)) > 1)
    model.matrix(~ item, d)[, -1, drop = FALSE] else NULL
  X <- cbind(intercept = 1, far = d[[far_var]], female = d$female,
             age = d$age, mu_d = d$mu_d)
  if (!is.null(county_mm)) X <- cbind(X, county_mm)
  if (!is.null(item_mm)) X <- cbind(X, item_mm)
  X
}

fit_wls <- function(X, y, w) {
  XtWX <- crossprod(X, X * w)
  beta <- tryCatch(solve(XtWX, crossprod(X, y * w)),
                   error = function(e) qr.solve(XtWX, crossprod(X, y * w), tol = 1e-10))
  drop(beta)
}

cr1_se <- function(X, y, beta, cluster) {
  resid <- y - drop(X %*% beta)
  bread <- tryCatch(solve(crossprod(X)),
                    error = function(e) qr.solve(crossprod(X), diag(ncol(X)), tol = 1e-10))
  ids <- unique(cluster)
  meat <- matrix(0, ncol(X), ncol(X))
  for (g in ids) {
    ii <- cluster == g
    score <- crossprod(X[ii, , drop = FALSE], resid[ii])
    meat <- meat + tcrossprod(score)
  }
  G <- length(ids); n <- nrow(X); p <- ncol(X)
  correction <- G / (G - 1) * (n - 1) / (n - p)
  vcov <- correction * bread %*% meat %*% bread
  sqrt(pmax(diag(vcov), 0))
}

analyze <- function(d, outcome, far_var, label, permute = TRUE) {
  keep <- complete.cases(d[, c(outcome, far_var, "female", "age", "mu_d",
                               "county", "cluster_id")])
  d <- d[keep, ]
  X <- design(d, far_var)
  y <- as.numeric(d[[outcome]])
  beta <- fit_wls(X, y, rep(1, nrow(d)))
  se <- cr1_se(X, y, beta, d$cluster_id)
  est <- beta[["far"]]
  t_obs <- est / se[which(colnames(X) == "far")]
  close_mean <- mean(y[d[[far_var]] == 0])

  # County-stratified exponential community weights (normalized within county
  # to sum to the county's number of communities), percentile interval.
  cframe <- unique(d[c("cluster_id", "county")])
  cframe <- cframe[order(cframe$cluster_id), ]
  cindex <- match(d$cluster_id, cframe$cluster_id)
  draws <- rep(NA_real_, B_weights)
  for (b in seq_len(B_weights)) {
    w_raw <- rexp(nrow(cframe))
    for (cty in unique(cframe$county)) {
      ii <- cframe$county == cty
      w_raw[ii] <- w_raw[ii] * sum(ii) / sum(w_raw[ii])
    }
    draws[b] <- tryCatch(fit_wls(X, y, w_raw[cindex])[["far"]],
                         error = function(e) NA_real_)
  }
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

  # Conditional randomization inference: permute the original Close/Far label
  # across this sample's communities within county, preserving observed counts.
  fisher_p <- NA_real_; stud_p <- NA_real_
  if (permute && far_var == "original_far") {
    far_g <- tapply(d[[far_var]], d$cluster_id, function(v) v[1])
    cty_g <- tapply(as.character(d$county), d$cluster_id, function(v) v[1])
    g_ids <- as.integer(names(far_g))
    est_perm <- numeric(B_perm); t_perm <- numeric(B_perm)
    far_col <- which(colnames(X) == "far")
    for (b in seq_len(B_perm)) {
      perm_far_g <- far_g
      for (cty in unique(cty_g)) {
        ii <- which(cty_g == cty)
        perm_far_g[ii] <- sample(far_g[ii], length(ii), replace = FALSE)
      }
      Xp <- X
      Xp[, far_col] <- perm_far_g[match(d$cluster_id, g_ids)]
      bp <- fit_wls(Xp, y, rep(1, nrow(d)))
      sp <- cr1_se(Xp, y, bp, d$cluster_id)
      est_perm[b] <- bp[["far"]]
      t_perm[b] <- bp[["far"]] / sp[far_col]
    }
    fisher_p <- (1 + sum(abs(est_perm) >= abs(est) - 1e-12)) / (B_perm + 1)
    stud_p <- (1 + sum(abs(t_perm) >= abs(t_obs) - 1e-12)) / (B_perm + 1)
  }

  data.frame(sample = label, outcome = outcome, far_var = far_var,
             close_mean = close_mean, estimate = est,
             ci_lower = ci[1], ci_upper = ci[2],
             fisher_p = fisher_p, studentized_p = stud_p,
             n = nrow(d), clusters = length(unique(d$cluster_id)))
}

# ---------------------------------------------------------------------------
set.seed(seed)
resp_seen <- respondent[respondent$has_response == 1L, ]
resp_seen_sms0 <- resp_seen[resp_seen$sms_treatment == "sms.control", ]

samples <- list(
  list(d = resp_seen_sms0[resp_seen_sms0$item == "bracelet", ],
       outcome = "any_M", far = "original_far", label = "Bracelet (primary)"),
  list(d = resp_seen[resp_seen$item == "bracelet", ],
       outcome = "any_M", far = "original_far", label = "Bracelet, incl. SMS"),
  list(d = resp_seen_sms0[resp_seen_sms0$item == "bracelet", ],
       outcome = "any_M", far = "final_far", label = "Bracelet, realized Far"),
  list(d = resp_seen_sms0[resp_seen_sms0$item == "ink", ],
       outcome = "any_M", far = "original_far", label = "Ink"),
  list(d = resp_seen_sms0[resp_seen_sms0$item %in% c("bracelet", "ink"), ],
       outcome = "any_M", far = "original_far", label = "Pooled signals"),
  list(d = resp_seen_sms0[resp_seen_sms0$item == "calendar", ],
       outcome = "any_M", far = "original_far", label = "Calendar (placebo)"),
  list(d = resp_seen_sms0[resp_seen_sms0$item == "calendar", ],
       outcome = "any_V", far = "original_far", label = "Calendar, private value"),
  list(d = resp_seen_sms0[resp_seen_sms0$item == "bracelet", ],
       outcome = "P_only", far = "original_far", label = "Bracelet, marker only")
)

# Unconditional sensitivity: every arm respondent, non-response coded 0.
uncond <- respondent[respondent$item == "bracelet" &
                       respondent$sms_treatment == "sms.control", ]
uncond$any_M_uncond <- ifelse(uncond$has_response == 1L, uncond$code_M, 0)
samples <- c(samples, list(
  list(d = transform(uncond, any_M = any_M_uncond),
       outcome = "any_M", far = "original_far", label = "Bracelet, unconditional")))

results <- do.call(rbind, lapply(samples, function(s)
  analyze(s$d, s$outcome, s$far, s$label)))

# Selection check: does reporting having seen the item respond to distance?
selection <- do.call(rbind, lapply(c("bracelet", "ink", "calendar"), function(itm) {
  d <- respondent[respondent$item == itm &
                    respondent$sms_treatment == "sms.control", ]
  d$seen01 <- as.numeric(d$seen_item)
  analyze(d, "seen01", "original_far", paste0("Seen ", itm), permute = FALSE)
}))

write.csv(rbind(results, selection),
          file.path(in_dir, "meaning-effects.csv"), row.names = FALSE)
print(rbind(results, selection), row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# Category shares by item and original assignment.
shares <- do.call(rbind, lapply(c("bracelet", "ink", "calendar"), function(itm) {
  do.call(rbind, lapply(0:1, function(far) {
    d <- resp_seen[resp_seen$item == itm & resp_seen$original_far == far, ]
    out <- data.frame(item = itm, original_far = far, n = nrow(d))
    for (l in labels) out[[l]] <- mean(d[[paste0("code_", l)]])
    out
  }))
}))
write.csv(shares, file.path(in_dir, "meaning-shares.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# LaTeX tables.
fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")
tex_escape <- function(x) gsub("([&%$#_])", "\\\\\\1", x)

shares_tex <- c(
  "\\begin{tabular}{llcccccc}", "\\toprule",
  "Item & Assignment & $N$ & M & P & C & V & O \\\\", "\\midrule")
for (i in seq_len(nrow(shares))) {
  s <- shares[i, ]
  shares_tex <- c(shares_tex, sprintf(
    "%s & %s & %d & %s & %s & %s & %s & %s \\\\",
    if (s$original_far == 0) tools::toTitleCase(s$item) else "",
    ifelse(s$original_far == 1, "Far", "Close"), s$n,
    fmt(s$M), fmt(s$P), fmt(s$C), fmt(s$V), fmt(s$O)))
  if (s$original_far == 1 && i < nrow(shares))
    shares_tex <- c(shares_tex, "\\addlinespace")
}
shares_tex <- c(shares_tex, "\\bottomrule", "\\end{tabular}")
writeLines(shares_tex, file.path(in_dir, "signal-meaning-shares.tex"))

effects_tex <- c(
  "\\begin{tabular}{lccccccc}", "\\toprule",
  "Sample & Outcome & Close mean & Far $-$ Close & 95\\% interval & Fisher $p$ & Stud.\\ RI $p$ & $N$ \\\\",
  "\\midrule")
pretty_outcome <- c(any_M = "Motivation (M)", any_V = "Private value (V)",
                    P_only = "Marker only", seen01 = "Seen item")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  effects_tex <- c(effects_tex, sprintf(
    "%s & %s & %s & %s & [%s, %s] & %s & %s & %s \\\\",
    tex_escape(r$sample), pretty_outcome[[r$outcome]],
    fmt(r$close_mean), fmt(r$estimate), fmt(r$ci_lower), fmt(r$ci_upper),
    ifelse(is.na(r$fisher_p), "--", fmt(r$fisher_p)),
    ifelse(is.na(r$studentized_p), "--", fmt(r$studentized_p)),
    format(r$n, big.mark = ",")))
}
effects_tex <- c(effects_tex, "\\addlinespace")
for (i in seq_len(nrow(selection))) {
  r <- selection[i, ]
  effects_tex <- c(effects_tex, sprintf(
    "%s & %s & %s & %s & [%s, %s] & -- & -- & %s \\\\",
    tex_escape(r$sample), "Seen item",
    fmt(r$close_mean), fmt(r$estimate), fmt(r$ci_lower), fmt(r$ci_upper),
    format(r$n, big.mark = ",")))
}
effects_tex <- c(effects_tex, "\\bottomrule", "\\end{tabular}")
writeLines(effects_tex, file.path(in_dir, "signal-meaning-effects.tex"))

# Example quotes: most common coded examples per category (bracelet arm).
final_with_text <- merge(final, dedup, by = "response_id")
quotes_tex <- c("\\begin{tabular}{lp{0.62\\linewidth}r}", "\\toprule",
                "Category & Example responses (verbatim) & $N$ \\\\", "\\midrule")
for (l in labels) {
  ex <- final_with_text[final_with_text[[l]] == 1 &
                          final_with_text$item == "bracelet", ]
  ex <- ex[order(-ex$n_respondents), ]
  ex <- head(ex, 3)
  if (!nrow(ex)) next
  for (j in seq_len(nrow(ex))) {
    quotes_tex <- c(quotes_tex, sprintf(
      "%s & ``%s'' & %d \\\\",
      if (j == 1) l else "", tex_escape(ex$text[j]), ex$n_respondents[j]))
  }
  if (l != labels[length(labels)]) quotes_tex <- c(quotes_tex, "\\addlinespace")
}
quotes_tex <- c(quotes_tex, "\\bottomrule", "\\end{tabular}")
writeLines(quotes_tex, file.path(in_dir, "signal-meaning-quotes.tex"))

# This exercise is archived rather than paper-facing: the meaning question
# elicits the signal's denotation, so its near-zero motivation base rate does
# not discipline the inference channel.  Outputs stay in ref-reports only.
message("Tables written to ", in_dir, " (archived exercise; not copied to the bundle)")
