takeup_project_root <- function() {
  normalizePath(Sys.getenv("TAKEUP_PROJECT_ROOT", "."), mustWork = TRUE)
}

takeup_build_root <- function() {
  file.path(takeup_project_root(), Sys.getenv("TAKEUP_BUILD_ROOT", "build"))
}

takeup_run <- function(command, args = character(), env = character(),
                       wd = takeup_project_root()) {
  status <- withr::with_dir(
    wd,
    system2(command, args, env = env, stdout = "", stderr = "")
  )
  if (!identical(status, 0L)) {
    stop("Command failed (", status, "): ", command, " ",
         paste(args, collapse = " "), call. = FALSE)
  }
  invisible(status)
}

takeup_write_build_manifest <- function(specification, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  git_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(...) NA_character_
  )
  manifest <- data.frame(
    distance_specification = specification,
    git_commit = git_commit[1L],
    r_version = R.version.string,
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    threads = Sys.getenv("TAKEUP_THREADS", "8")
  )
  path <- file.path(output_dir, "build-manifest.csv")
  utils::write.csv(manifest, path, row.names = FALSE)
  path
}

takeup_prepare_run_dir <- function(specification) {
  root <- takeup_project_root()
  run_dir <- file.path(takeup_build_root(), "work", specification)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  legacy_links <- c(
    "analysis_util.R", "balance-functions.R", "balance.R",
    "clean-analysis-util.R", "dist_structural_util.R", "scratch"
  )
  for (entry in legacy_links) {
    destination <- file.path(run_dir, entry)
    if (nzchar(Sys.readlink(destination))) unlink(destination)
  }
  source_entries <- c(
    "R", "scripts", "rct-design-fieldwork", "multilvlr",
    "simulate-treatment-assignment", "stan_models", "data",
    "takeup.Rproj", "renv.lock"
  )
  for (entry in source_entries) {
    source <- file.path(root, entry)
    destination <- file.path(run_dir, entry)
    if (!file.exists(source)) next
    link_target <- Sys.readlink(destination)
    is_link <- !is.na(link_target) && nzchar(link_target)
    if (is_link &&
        !identical(normalizePath(destination), normalizePath(source))) {
      unlink(destination)
    }
    if (file.exists(destination)) next
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.symlink(source, destination)) {
      stop("Could not link build input: ", source, call. = FALSE)
    }
  }
  dir.create(file.path(run_dir, "temp-data", "tidy-rf-tes"),
             recursive = TRUE, showWarnings = FALSE)
  carry_inputs <- c(
    "analysis-cluster-covariate-data.csv",
    "balance-cts-dist-ri-fe.rds",
    "balance-ri-fe.rds"
  )
  for (filename in carry_inputs) {
    source <- file.path(root, "temp-data", filename)
    destination <- file.path(run_dir, "temp-data", filename)
    if (file.exists(source) && !file.exists(destination)) {
      file.copy(source, destination)
    }
  }
  dir.create(file.path(run_dir, "presentations", "rf-tables", "main-specs"),
             recursive = TRUE, showWarnings = FALSE)
  presentation_inputs <- c("balance-tables.Rmd", "plot-fn.R")
  for (filename in presentation_inputs) {
    source <- file.path(root, "presentations", filename)
    destination <- file.path(run_dir, "presentations", filename)
    if (file.exists(source)) {
      file.copy(source, destination, overwrite = TRUE)
    }
  }
  run_dir
}

takeup_run_reduced_form <- function(specification, output_dir, context_files,
                                    bootstrap_draws,
                                    dependencies = NULL) {
  table_dir <- file.path(output_dir, "presentations", "rf-tables", "main-specs")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(output_dir, "reduced-form.log")
  env <- c(
    paste0("TAKEUP_DISTANCE_SPEC=", specification),
    paste0("TAKEUP_BUILD_OUTPUT=", normalizePath(output_dir, mustWork = FALSE))
  )
  run_dir <- takeup_prepare_run_dir(specification)
  context_path <- normalizePath(
    context_files[grepl("analysis-context[.]rds$", context_files)],
    mustWork = TRUE
  )
  status <- withr::with_dir(run_dir, system2(
      "Rscript",
      c("--no-save", "--no-restore", "scripts/reduced-form/bootstrap.R",
        paste0("--distance-definition=", specification),
        paste0("--bootstrap-draws=", bootstrap_draws),
        paste0("--context-path=", context_path),
        paste0("--table-output-path=", table_dir)),
      env = env,
      stdout = log_path,
      stderr = log_path
    ))
  if (status != 0L) stop("Reduced-form build failed; see ", log_path, call. = FALSE)
  generated <- list.files(table_dir, pattern = "[.]tex$", full.names = TRUE)
  if (!length(generated) || any(file.info(generated)$size == 0L)) {
    stop("Reduced-form build produced no valid tables.", call. = FALSE)
  }
  marker <- file.path(output_dir, "reduced-form.complete")
  utils::write.csv(data.frame(
    distance_definition = specification, bootstrap_draws = bootstrap_draws,
    generated_tables = length(generated),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ), marker, row.names = FALSE)
  marker
}

takeup_run_balance <- function(specification, output_dir,
                               dependencies = NULL) {
  balance_dir <- file.path(output_dir, "balance")
  dir.create(balance_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(output_dir, "balance.log")
  run_dir <- takeup_prepare_run_dir(specification)
  status <- withr::with_dir(run_dir, system2(
      "Rscript",
      c("--no-save", "--no-restore", "scripts/balance/run.R", "--main", "--attrition",
        "--monitored-attrition", "--sms", paste0("--output-path=", balance_dir),
        paste0("--distance-definition=", specification)),
      env = paste0("TAKEUP_DISTANCE_SPEC=", specification),
      stdout = log_path,
      stderr = log_path
    ))
  if (status != 0L) stop("Balance build failed; see ", log_path, call. = FALSE)
  generated_balance_files <- list.files(balance_dir, full.names = TRUE)
  if (length(generated_balance_files)) {
    file.copy(generated_balance_files, file.path(run_dir, "temp-data"),
              overwrite = TRUE)
  }
  marker <- file.path(output_dir, "balance.complete")
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), marker)
  marker
}

takeup_run_balance_section <- function(specification, output_dir, section,
                                       context_files,
                                       ri_draws,
                                       dependencies = NULL) {
  valid_sections <- c(
    "main", "orig", "fit-ri", "attrition", "monitored-attrition", "sms"
  )
  if (!section %in% valid_sections) {
    stop("Unknown balance section: ", section, call. = FALSE)
  }
  balance_dir <- file.path(output_dir, "balance")
  dir.create(balance_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(output_dir, paste0("balance-", section, ".log"))
  run_dir <- takeup_prepare_run_dir(specification)
  context_path <- normalizePath(
    context_files[grepl("analysis-context[.]rds$", context_files)],
    mustWork = TRUE
  )
  status <- withr::with_dir(run_dir, system2(
    "Rscript",
    c("--no-save", "--no-restore", "scripts/balance/run.R", paste0("--", section),
      paste0("--output-path=", balance_dir),
      paste0("--ri-draws=", ri_draws),
      paste0("--context-path=", context_path),
      paste0("--distance-definition=", specification)),
    env = paste0("TAKEUP_DISTANCE_SPEC=", specification),
    stdout = log_path, stderr = log_path
  ))
  if (status != 0L) {
    stop("Balance section '", section, "' failed; see ", log_path,
         call. = FALSE)
  }
  expected <- c(
    main = "main-balance-samples.rds",
    orig = "cluster_balance_fits.rds",
    `fit-ri` = "balance-ri-fe.rds",
    attrition = "attrition_treat_comp_df.csv",
    `monitored-attrition` = "monitoring_attrition_comp_df.csv",
    sms = "sms_enrollment_balance.rds"
  )[[section]]
  expected_path <- file.path(balance_dir, expected)
  if (!file.exists(expected_path) || file.info(expected_path)$size == 0L) {
    stop("Balance section '", section, "' did not produce ", expected,
         ".", call. = FALSE)
  }
  generated_balance_files <- list.files(balance_dir, full.names = TRUE)
  if (length(generated_balance_files)) {
    file.copy(generated_balance_files, file.path(run_dir, "temp-data"),
              overwrite = TRUE)
  }
  marker <- file.path(output_dir, paste0("balance-", section, ".complete"))
  utils::write.csv(data.frame(
    distance_definition = specification, section = section,
    ri_draws = ri_draws, required_output = expected,
    required_output_sha256 = unname(tools::md5sum(expected_path)),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ), marker, row.names = FALSE)
  marker
}

takeup_render_balance_tables <- function(specification, output_dir,
                                         balance_dependency,
                                         render_dependencies = NULL) {
  balance_dependency
  render_dependencies
  run_dir <- takeup_prepare_run_dir(specification)
  table_dir <- file.path(output_dir, "presentations", "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(output_dir, "balance-tables.log")
  expression <- paste0(
    "rmarkdown::render('presentations/balance-tables.Rmd', ",
    "output_format='md_document', output_file='balance-tables.md', ",
    "output_dir='", normalizePath(file.path(run_dir, "rendered"),
                                  mustWork = FALSE), "', knit_root_dir='",
    normalizePath(run_dir, mustWork = FALSE), "', params=list(",
    "output_path='temp-data', table_output_path='",
    normalizePath(table_dir, mustWork = FALSE), "', cache=FALSE), ",
    "quiet=TRUE)"
  )
  dir.create(file.path(run_dir, "rendered"), recursive = TRUE,
             showWarnings = FALSE)
  status <- withr::with_dir(run_dir, system2(
    "Rscript", c("--no-save", "--no-restore", "-e", shQuote(expression)),
    env = paste0("TAKEUP_DISTANCE_SPEC=", specification),
    stdout = log_path, stderr = log_path
  ))
  if (status != 0L) stop("Balance table render failed; see ", log_path,
                         call. = FALSE)
  marker <- file.path(output_dir, "balance-tables.complete")
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), marker)
  marker
}

takeup_saved_draws <- function(specification) {
  deposited_root <- Sys.getenv(
    paste0("TAKEUP_", toupper(specification), "_DRAWS"),
    file.path("replication", "structural-draws", specification)
  )
  generated_root <- file.path("build", "structural-fit", specification)
  root <- if (dir.exists(deposited_root)) deposited_root else generated_root
  files <- list.files(
    root, full.names = TRUE,
    pattern = "^dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-[1-4][.]csv$"
  )
  if (length(files) != 4L) {
    stop(
      "Expected four saved slim-chain CSVs in ", root,
      ". Run `make structural-fit DISTANCE_SPEC=", specification,
      "` or deposit the replication draws.", call. = FALSE
    )
  }
  files
}

takeup_run_compact_gq <- function(specification, output_dir, draws) {
  gq_dir <- file.path(output_dir, "structural", "compact-gq")
  dir.create(gq_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(output_dir, "structural-postprocess.log")
  args <- c(
    "--no-save", "--no-restore", "scripts/structural/generate-compact-gq.R",
    paste0("--workspace=", Sys.getenv(
      "TAKEUP_STRUCTURAL_WORKSPACE",
      "data/stan_analysis_data/dist_fit104.RData"
    )),
    paste0("--fit-csvs=", paste(draws, collapse = ",")),
    paste0("--output-path=", gq_dir),
    paste0("--distance-definition=", specification)
  )
  status <- system2(
    "Rscript", args,
    env = paste0("TAKEUP_DISTANCE_SPEC=", specification),
    stdout = log_path, stderr = log_path
  )
  if (status != 0L) stop("Compact GQ failed; see ", log_path, call. = FALSE)
  files <- Sys.glob(file.path(gq_dir, "*.csv"))
  if (!length(files)) stop("Compact GQ produced no CSV files.", call. = FALSE)
  files
}

takeup_compare_distance_outputs <- function(assigned_dir, realized_dir, output_dir,
                                            dependencies = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  relative_files <- function(root) {
    paths <- list.files(root, recursive = TRUE, full.names = TRUE,
                        pattern = "\\.(csv|tex)$")
    sub(paste0("^", normalizePath(root), "/?"), "", normalizePath(paths))
  }
  a <- relative_files(assigned_dir)
  r <- relative_files(realized_dir)
  registry <- merge(data.frame(path = a, assigned = TRUE),
                    data.frame(path = r, realized = TRUE), all = TRUE)
  registry$assigned[is.na(registry$assigned)] <- FALSE
  registry$realized[is.na(registry$realized)] <- FALSE
  registry$identical <- vapply(registry$path, function(path) {
    if (!file.exists(file.path(assigned_dir, path)) ||
        !file.exists(file.path(realized_dir, path))) return(FALSE)
    identical(readBin(file.path(assigned_dir, path), "raw",
                      file.info(file.path(assigned_dir, path))$size),
              readBin(file.path(realized_dir, path), "raw",
                      file.info(file.path(realized_dir, path))$size))
  }, logical(1))
  path <- file.path(output_dir, "artifact-comparison.csv")
  utils::write.csv(registry, path, row.names = FALSE)
  path
}
