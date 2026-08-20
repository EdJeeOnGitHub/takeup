library(targets)

source("R/distance-spec.R")
source("R/pipeline.R")

tar_option_set(
  packages = c("dplyr", "withr"),
  error = "abridge",
  memory = "transient",
  garbage_collection = TRUE
)

requested_specs <- Sys.getenv("TAKEUP_BUILD_SPECS", "assigned")
specifications <- if (identical(requested_specs, "both")) {
  c("assigned", "realized")
} else {
  takeup_distance_spec(requested_specs)
}
requested_balance_sections <- Sys.getenv("TAKEUP_BALANCE_SECTIONS", "all")
balance_sections <- if (identical(requested_balance_sections, "all")) {
  c("main", "orig", "fit-ri", "attrition", "monitored-attrition", "sms")
} else {
  strsplit(requested_balance_sections, ",", fixed = TRUE)[[1L]]
}

list(
  tar_target(
    distance_input_files,
    c(
      "data/rct_targetable_schools_2.0.rds",
      "data/takeup_processed_cluster_strat.rds",
      "data/clean-data/monitored-nosms-takeup-data.rds"
    ),
    format = "file"
  ),
  tar_target(
    reduced_form_source_files,
    c(
      "R/distance-spec.R", "scratch/reduced-form-bootstrap.R",
      "scratch/reduced-form-setup.R", "scratch/reduced-form-functions.R",
      "analysis_util.R", "clean-analysis-util.R", "dist_structural_util.R",
      "rct-design-fieldwork/takeup_rct_assign_clusters.R",
      "multilvlr/multilvlr_util.R"
    ),
    format = "file"
  ),
  tar_target(
    balance_source_files,
    c(
      "R/distance-spec.R", "scratch/reduced-form-setup.R",
      "scratch/reduced-form-functions.R", "analysis_util.R",
      "clean-analysis-util.R",
      "rct-design-fieldwork/takeup_rct_assign_clusters.R",
      "balance.R", "balance-functions.R"
    ),
    format = "file"
  ),
  tar_target(
    balance_table_source_files,
    c("presentations/balance-tables.Rmd", "presentations/plot-fn.R"),
    format = "file"
  ),
  tar_target(distance_specification, specifications, iteration = "vector"),
  tar_target(balance_section, balance_sections, iteration = "vector"),
  tar_target(
    distance_crosswalk,
    {
      distance_input_files
      takeup_distance_crosswalk()
    }
  ),
  tar_target(
    distance_audit,
    takeup_write_distance_audit(
      distance_crosswalk,
      file.path(takeup_build_root(), distance_specification, "audit"),
      distance_specification
    ),
    pattern = map(distance_specification),
    format = "file"
  ),
  tar_target(
    build_manifest,
    takeup_write_build_manifest(
      distance_specification,
      file.path(takeup_build_root(), distance_specification)
    ),
    pattern = map(distance_specification),
    format = "file"
  ),
  tar_target(
    reduced_form,
    takeup_run_reduced_form(
      distance_specification,
      file.path(takeup_build_root(), distance_specification),
      dependencies = reduced_form_source_files
    ),
    pattern = map(distance_specification),
    format = "file"
  ),
  tar_target(
    balance,
    takeup_run_balance_section(
      distance_specification,
      file.path(takeup_build_root(), distance_specification),
      balance_section,
      dependencies = balance_source_files
    ),
    pattern = cross(distance_specification, balance_section),
    format = "file"
  ),
  tar_target(
    balance_tables,
    takeup_render_balance_tables(
      distance_specification,
      file.path(takeup_build_root(), distance_specification),
      balance_dependency = balance,
      render_dependencies = balance_table_source_files
    ),
    pattern = map(distance_specification),
    format = "file"
  ),
  tar_target(
    saved_structural_draws,
    takeup_saved_draws(distance_specification),
    pattern = map(distance_specification),
    format = "file"
  ),
  tar_target(
    compact_gq,
    takeup_run_compact_gq(
      distance_specification,
      file.path(takeup_build_root(), distance_specification),
      saved_structural_draws
    ),
    pattern = map(distance_specification),
    format = "file"
  ),
  tar_target(
    distance_comparison,
    takeup_compare_distance_outputs(
      file.path(takeup_build_root(), "assigned"),
      file.path(takeup_build_root(), "realized"),
      file.path(takeup_build_root(), "comparison"),
      dependencies = c(reduced_form, balance)
    ),
    format = "file"
  )
)
