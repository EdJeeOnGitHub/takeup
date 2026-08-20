#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
})

fit_version <- "104"
model_name <- "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
arms <- c("Control", "Ink", "Calendar", "Bracelet")
num_clusters <- 144L
read_chain <- function(chain) {
  read_selected <- function(path, variables) {
    fread(
      cmd = sprintf("awk '!/^#/' %s", shQuote(path)),
      select = variables,
      showProgress = FALSE
    )
  }
  variables <- c(
    "base_mu_rep",
    unlist(lapply(seq_along(arms), function(arm) {
      sprintf("cluster_mu_rep.16.%d.%d", seq_len(num_clusters), arm)
    }))
  )
  fit <- read_selected(
    sprintf(
      "data/stan_analysis_data/dist_fit%s_%s-%d.csv",
      fit_version, model_name, chain
    ),
    variables
  )

  visibility <- sapply(seq_along(arms), function(arm) {
    cols <- sprintf(
      "cluster_mu_rep.16.%d.%d",
      seq_len(num_clusters), arm
    )
    rowMeans(fit[, ..cols]) / fit[["base_mu_rep"]]
  })
  colnames(visibility) <- arms

  likelihood_names <- c(
    "sensitivity_log_lik_takeup",
    "sensitivity_log_lik_beliefs",
    "sensitivity_log_lik_wtp"
  )
  likelihood <- read_selected(
    sprintf(
      "temp-data/bayesian-sensitivity/dist_fit%s_%s_sensitivity-%d.csv",
      fit_version, model_name, chain
    ),
    likelihood_names
  ) |> as.matrix()

  list(objects = visibility, likelihood = likelihood)
}

chain_results <- parallel::mclapply(1:4, read_chain, mc.cores = 4L)
failed <- vapply(chain_results, inherits, logical(1), "try-error")
if (any(failed)) {
  stop(paste(unlist(chain_results[failed]), collapse = "\n"))
}
objects <- do.call(rbind, lapply(chain_results, `[[`, "objects"))
likelihood <- do.call(rbind, lapply(chain_results, `[[`, "likelihood"))
colnames(likelihood) <- c("Takeup", "Beliefs", "WTP")
sensitivity <- cov(objects, likelihood)
sensitivity <- sweep(sensitivity, 1, apply(objects, 2, sd), "/")

posterior_mean <- colMeans(objects)
posterior_variance <- apply(objects, 2, var)
centered_squared <- sweep(objects, 2, posterior_mean, "-")^2
log_sd_sensitivity <- cov(centered_squared, likelihood)
log_sd_sensitivity <- sweep(
  log_sd_sensitivity,
  1,
  2 * posterior_variance,
  "/"
)

output <- data.frame(
  object = "Visibility",
  row = rownames(sensitivity),
  setNames(as.data.frame(sensitivity), paste0("mean_", colnames(sensitivity))),
  setNames(
    as.data.frame(log_sd_sensitivity),
    paste0("log_sd_", colnames(log_sd_sensitivity))
  ),
  check.names = FALSE
)
write_csv(
  output,
  "temp-data/bayesian-sensitivity/visibility-sensitivity-1500m.csv"
)
