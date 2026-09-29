.libPaths(c(.libPaths(), '/home/ed/R/x86_64-pc-linux-gnu-library/4.5'))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr);library(stringr);library(forcats);library(fixest)})
source('R/reduced-form/functions.R')
options(dplyr.summarise.inform=FALSE,fixest_notes=FALSE)
setFixest_nthreads(1)
Sys.setenv(TAKEUP_THREADS=1,OPENBLAS_NUM_THREADS=1)
ctx <- readRDS('build/introduction-pvalue-audit-20260929/fresh-context.rds')
expected <- read_csv('data/cluster_expected_dist.csv',show_col_types=FALSE)
correct_cov <- ctx$data$cov_analysis_data
correct_cov$mu_d <- expected$dist[match(correct_cov$cluster.id.x,expected$cluster.id)] / ctx$config$distance_sd
stopifnot(!anyNA(correct_cov$mu_d))
correct_end <- ctx$data$endline_data
correct_end$mu_d <- expected$dist[match(as.numeric(as.character(correct_end$cluster.id)),expected$cluster.id)] / ctx$config$distance_sd
stopifnot(!anyNA(correct_end$mu_d))
f_take <- function(data,weights) feols(dewormed ~0+assigned_treatment*assigned_dist_group+female+age.census+mu_d|county,data=data,weights=weights,notes=FALSE)
f_obs <- function(data,weights) feols(prop_knows ~assigned_treatment+assigned_dist_group+i(assigned_treatment,assigned_dist_group,'control')+female+age.census+mu_d|county,data=data,weights=weights,notes=FALSE)
f_pred <- function(data,weights) feols(dworm_frac ~0+assigned_treatment+assigned_dist_group+i(assigned_treatment,assigned_dist_group,'control')+female+age.census+mu_d|county,data=data,weights=weights,notes=FALSE)
obs_frame <- function(x) x %>% left_join(filter(ctx$data$summ_endline_know_table,know.table.type=='table.A'),by='KEY.individ') %>% filter(sms.treatment=='sms.control',obs_know_person>0) %>% mutate(assigned_treatment=assigned.treatment,assigned_dist_group=dist.pot.group,prop_knows=knows_other_dewormed/obs_know_person)
pred_frame <- function(x) x %>% mutate(assigned_treatment=as_factor(assigned.treatment),assigned_dist_group=as_factor(dist.pot.group),cluster_id=as_factor(cluster.id),dworm_frac=dworm_rate/10)
models <- list(takeup=list(f=f_take,outcome='dewormed',old=ctx$data$cov_analysis_data,new=correct_cov),observability=list(f=f_obs,outcome='prop_knows',old=obs_frame(ctx$data$endline_data),new=obs_frame(correct_end)),predicted=list(f=f_pred,outcome='dworm_frac',old=pred_frame(ctx$data$endline_data),new=pred_frame(correct_end)))
saveRDS(models,'build/introduction-pvalue-audit-20260929/model-frames.rds')
for (name in names(models)) for (version in c('old','new')) {
  m<-models[[name]]; d<-m[[version]]; f<-enable_fast_discrete_wls(m$f,m$outcome,c('female','age.census','mu_d'))
  cat('\nRUN',name,version,'rows',nrow(d),'clusters',n_distinct(d$cluster.id),'\n')
  actual<-actual_bayesian_bs_fit('realised fit',f,d)
  bs<-bootstrap_regression_draws(f,d,500)
  results<-bind_rows(lapply(list(clean_te_draws,clean_signal_draws),function(clean){add_summ_stats(clean(bs),clean(actual$bs_fit) %>% transmute(assigned_treatment,assigned_dist_group,realised_pred=estimate))})) %>% filter(assigned_treatment!='no signal')
  write_csv(results,paste0('build/introduction-pvalue-audit-20260929/',name,'-',version,'-results.csv'))
  saveRDS(bs,paste0('build/introduction-pvalue-audit-20260929/',name,'-',version,'-draws.rds'))
  print(results %>% filter(assigned_treatment %in% c('calendar','bracelet','control','bracelet - calendar'),assigned_dist_group %in% c('combined','far - close')))
}
