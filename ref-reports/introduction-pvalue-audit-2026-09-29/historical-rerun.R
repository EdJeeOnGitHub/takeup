.libPaths(c(.libPaths(), '/home/ed/R/x86_64-pc-linux-gnu-library/4.5'))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr);library(stringr);library(fixest)})
source('R/reduced-form/functions.R');setFixest_nthreads(1);options(dplyr.summarise.inform=FALSE)
Sys.setenv(TAKEUP_THREADS=1)
ctx<-readRDS('build/introduction-pvalue-audit-20260929/fresh-context.rds')
models<-readRDS('build/introduction-pvalue-audit-20260929/model-frames.rds')
expected<-read_csv('build/introduction-pvalue-audit-20260929/expected-distance-100-seeds.csv',show_col_types=FALSE)
# Preserve the historical erroneous rank join to isolate simulation precision.
original_ids<-sort(unique(expected$cluster.id))
lookup<-ctx$data$cov_analysis_data %>% distinct(cluster.id.x,cluster_id) %>% mutate(mu_d=expected$dist[match(original_ids[cluster_id],expected$cluster.id)]/ctx$config$distance_sd)
for (name in names(models)) {
 m<-models[[name]];d<-m$old
 ids<-if(name=='takeup') d$cluster.id.x else as.numeric(as.character(d$cluster.id))
 d$mu_d<-lookup$mu_d[match(ids,lookup$cluster.id.x)]
 fit<-actual_bayesian_bs_fit('realised fit',m$f,d)
 results<-bind_rows(clean_te_draws(fit$bs_fit),clean_signal_draws(fit$bs_fit)) %>% filter(assigned_treatment!='no signal')
 write_csv(results,paste0('build/introduction-pvalue-audit-20260929/',name,'-100-seeds-points.csv'))
 f<-enable_fast_discrete_wls(m$f,m$outcome,c('female','age.census','mu_d'))
 bs<-bootstrap_regression_draws(f,d,500)
 results<-bind_rows(lapply(list(clean_te_draws,clean_signal_draws),function(clean){add_summ_stats(clean(bs),clean(fit$bs_fit) %>% transmute(assigned_treatment,assigned_dist_group,realised_pred=estimate))})) %>% filter(assigned_treatment!='no signal')
 write_csv(results,paste0('build/introduction-pvalue-audit-20260929/',name,'-100-seeds-results.csv'))
}
