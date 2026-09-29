.libPaths(c(.libPaths(), '/home/ed/R/x86_64-pc-linux-gnu-library/4.5'))
Sys.setenv(TAKEUP_THREADS=1,OPENBLAS_NUM_THREADS=1,TZ='UTC')
options(dplyr.summarise.inform=FALSE)
fixest::setFixest_nthreads(1)
# Replay exactly the production script, replacing only the context load with
# an explicitly historical reconstruction. Never relax production validation.
expressions <- parse(Sys.getenv('TAKEUP_AUDIT_SCRIPT', 'scripts/reduced-form/bootstrap.R'))
replaced <- FALSE
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
      identical(expr[[2]], as.name('analysis_context'))) {
    analysis_context <- takeup_get_analysis_context(Sys.getenv("TAKEUP_ANALYSIS_CONTEXT", "/home/ed/projects/takeup/build/introduction-pvalue-audit-20260929/corrected-context.rds"))
    if (Sys.getenv('TAKEUP_AUDIT_VERSION') != 'corrected') {
    expected100 <- readr::read_csv(Sys.getenv('TAKEUP_AUDIT_EXPECTED', '/home/ed/projects/takeup/build/introduction-pvalue-audit-20260929/expected-distance-100-seeds.csv'),show_col_types=FALSE) |> dplyr::filter(!is.na(cluster.id))
    ordered_ids <- sort(expected100$cluster.id)
    cov <- analysis_context$data$cov_analysis_data
    historical_dist <- expected100$dist[match(ordered_ids[cov$cluster_id],expected100$cluster.id)]
    stopifnot(!anyNA(historical_dist))
    cov$clust_expected_dist <- historical_dist
    cov$standard_clust_expected_dist <- historical_dist / analysis_context$config$distance_sd
    cov$mu_d <- cov$standard_clust_expected_dist
    end <- analysis_context$data$endline_data
    end$mu_d <- cov$mu_d[match(as.integer(as.character(end$cluster.id)),cov$cluster.id.x)]
    end$standard_clust_expected_dist <- end$mu_d
    stopifnot(!anyNA(end$mu_d))
    analysis_context$data$cov_analysis_data <- cov
    analysis_context$data$endline_data <- end
    # Other frames join expected distances by original ID and were unaffected
    # by the rank-join bug, but use the historical simulation precision here.
    ad <- analysis_context$data$analysis_data
    ad$clust_expected_dist <- expected100$dist[match(as.integer(as.character(ad$cluster.id)),expected100$cluster.id)]
    ad$standard_clust_expected_dist <- ad$clust_expected_dist / analysis_context$config$distance_sd
    analysis_context$data$analysis_data <- ad
    analysis_context$data$cluster_expected_dist_df <- expected100 |>
      dplyr::transmute(cluster.id=factor(cluster.id),clust_expected_dist=dist)
    message('AUDIT ONLY: reconstructing the historical rank join; expected-distance input: ', Sys.getenv('TAKEUP_AUDIT_EXPECTED', 'reconstructed 100-draw control'))
    }
    replaced <- TRUE
  } else if (is.call(expr) && identical(expr[[1]], as.name("takeup_context_into_environment"))) {
    stopifnot(replaced)
    list2env(analysis_context$data, envir=.GlobalEnv)
  } else {
    eval(expr, envir=.GlobalEnv)
    if (is.call(expr) && identical(expr[[1]], as.name('source')) &&
        grepl('functions.R', paste(deparse(expr), collapse=''), fixed=TRUE)) {
      audit_original_wrapper <- wrapper_function
      wrapper_function <- function(...) {
        args <- list(...)
        path <- args$tidy_summ_path
        tex <- file.path(params$table_output_path,paste0(args$table_name,'.tex'))
        if (file.exists(path) && file.exists(tex)) {
          message('AUDIT CHECKPOINT: using completed ', path)
          return(list(tidy_summary=readr::read_csv(path,show_col_types=FALSE),
            default_tbl=paste(readLines(tex),collapse='\n')))
        }
        do.call(audit_original_wrapper,args)
      }
    }
  }
}
stopifnot(replaced)
