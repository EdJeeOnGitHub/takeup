library(tidyverse)
library(nleqslv)
library(microbenchmark)
library(owent)

source("R/structural/legacy-utils.R")
source("optim/optim-functions.R")

# Realistic synthetic inputs: 1252 village-PoT pairs, unbounded (bounds = c(-Inf, Inf))
set.seed(42)
N       <- 1252
b       <- rnorm(N, mean = -0.5, sd = 0.3)
mu_rep  <- runif(N, 0.1, 1.5)
u_sd    <- 0.8
total_error_sd <- sqrt(u_sd^2 + 1)
bounds  <- c(-Inf, Inf)

# ── Build the per-element fixed-point functions (shared setup cost) ──────────
v_fs <- map2(b, mu_rep,
    ~generate_v_cutoff_fixedpoint(
        b = .x, mu = .y,
        total_error_sd = total_error_sd,
        u_sd = u_sd,
        bounds = NULL   # unbounded
    ))

# ── Approach 1: current — nleqslv per element ────────────────────────────────
solver_nleqslv <- function() {
    fits <- map2(v_fs, b, ~nleqslv(x = -.y, fn = .x))
    map_dbl(fits, "x")
}

# ── Approach 2: uniroot per element ─────────────────────────────────────────
# Fixed-point eq: v + b + mu*delta(v) = 0, which is monotone in v.
# Bracket: fn(-10) and fn(10) should always straddle the root for our params.
solver_uniroot <- function() {
    map2_dbl(v_fs, b, ~{
        tryCatch(
            uniroot(.x, lower = -10, upper = 10, tol = 1e-8)$root,
            error = function(e) -.y   # fallback to initial guess on failure
        )
    })
}

# ── Correctness check first ──────────────────────────────────────────────────
v1 <- solver_nleqslv()
v2 <- solver_uniroot()
cat("Max abs diff (nleqslv vs uniroot):", max(abs(v1 - v2), na.rm = TRUE), "\n")
cat("Failures in uniroot:", sum(is.na(v2)), "\n\n")

# ── Benchmark ────────────────────────────────────────────────────────────────
mb <- microbenchmark(
    nleqslv = solver_nleqslv(),
    uniroot  = solver_uniroot(),
    times    = 20
)
print(mb)
cat("\nPer draw cost (ms):\n")
print(summary(mb)[, c("expr", "median", "mean")])
cat("\nProjected time for 200 draws (seconds):\n")
summ <- summary(mb)
cat("  nleqslv:", round(summ$median[summ$expr == "nleqslv"] * 200 / 1e9, 1), "s\n")
cat("  uniroot: ", round(summ$median[summ$expr == "uniroot"]  * 200 / 1e9, 1), "s\n")
