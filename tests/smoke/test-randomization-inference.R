#!/usr/bin/env Rscript

# End-to-end smoke and regression test. Outputs are confined to tempdir().
root <- normalizePath(".")
out <- file.path(tempdir(), "takeup-ri-test")
app <- file.path(tempdir(), "takeup-ri-test-appendix")
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("scripts/appendix/randomization-inference.R", "--permutations=99", "--seed=20260812",
    paste0("--output-dir=", out), paste0("--appendix-dir=", app))
)
stopifnot(status == 0L)

r <- read.csv(file.path(out, "data", "ri-results.csv"), check.names = FALSE)
a <- read.csv(file.path(out, "data", "sample-audit.csv"), check.names = FALSE)
s <- read.csv(file.path(out, "data", "assignment-stratum-counts.csv"), check.names = FALSE)
p <- read.csv(file.path(out, "data", "permutation-audit.csv"), check.names = FALSE)

stopifnot(nrow(r) == 24L, nrow(a) == 4L, nrow(s) == 6L)
stopifnot(a$observations[1] == 9805L, a$communities[1] == 144L)
stopifnot(all(a$communities == 144L), p$switches == 26L, p$permutations == 99L)
stopifnot(all(rowSums(s[c("control", "calendar", "ink", "bracelet")]) > 0))
stopifnot(all(r$fisher_p > 0 & r$fisher_p <= 1))
stopifnot(all(r$studentized_p > 0 & r$studentized_p <= 1))
stopifnot(all(r$romano_wolf_p >= r$studentized_p - 1e-12))

primary <- r[r$outcome == "Administrative take-up" &
             r$distance_specification == "Original assigned Close/Far", ]
stopifnot(abs(primary$estimate[primary$contrast == "Any Signal - No Signal"] -
              0.123569502770411) < 1e-10)
stopifnot(abs(primary$estimate[primary$contrast == "Bracelet - Calendar"] -
              0.0746489572930473) < 1e-10)

required <- c(
  file.path(out, "tables", "randomization-inference.tex"),
  file.path(out, "figures", "takeup-original-assignment.pdf"),
  file.path(out, "figures", "takeup-realized-distance.pdf"),
  file.path(app, "tables", "randomization-inference.tex")
)
stopifnot(all(file.exists(required)), all(file.info(required)$size > 0))
cat("randomization-inference smoke test passed\n")
