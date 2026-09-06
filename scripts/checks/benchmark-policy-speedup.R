#!/usr/bin/env Rscript
# Bounded throughput pilot on allocated compute, at unchanged economic inputs.
args <- commandArgs(TRUE); source('R/policy/bootstrap.R')
reference <- normalizePath(policy_option_value(args,'--reference-code'))
output <- policy_option_value(args,'--output-path')
source_root <- policy_option_value(args,'--source-root','/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-model-robustness-1451')
model <- policy_option_value(args,'--model','benchmark')
n <- as.integer(policy_option_value(args,'--draws','32'))
stopifnot(n>=8L,n<=64L)
candidate <- getwd();dir.create(output,recursive=TRUE,showWarnings=FALSE); output <- normalizePath(output)
base <- '/project/akaring/takeup-data'
target <- paste0(base,'/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv')
object <- readRDS(file.path(source_root,model,'policy-model-parameters.rds'))
object$parameters <- object$parameters[seq_len(n),,drop=FALSE]
if (!is.null(object$cluster_shock)) object$cluster_shock <- object$cluster_shock[seq_len(n),,drop=FALSE]
run <- function(code,script,opts,log) {
 old <- setwd(code);on.exit(setwd(old)); start <- proc.time()[[3L]]
 status <- system2('/usr/bin/time',c('-f','%M','-o',shQuote(paste0(log,'.rss')),shQuote(file.path(R.home('bin'),'Rscript')),'--vanilla',shQuote(script),shQuote(opts)),stdout=log,stderr=log)
 if(status)stop('Failed: ',log)
 proc.time()[[3L]]-start
}
rows <- list()
for (variant in c('reference','serial-auto','workers-1','workers-4','workers-8')) {
 dest <- file.path(output,variant);dir.create(dest)
 saveRDS(object,file.path(dest,'policy-model-parameters.rds'),compress=FALSE)
 code <- if(variant=='reference')reference else candidate
 run(code,'scripts/policy/predict-model-robustness.R',c(paste0('--parameter-rds=',dest,'/policy-model-parameters.rds'),
 paste0('--output-path=',dest),paste0('--distance-data=',base,'/optim/data/full-many-pots-experiment.rds'),
 '--num-cores=8','--distance-cap=3500'),file.path(dest,'predict.log'))
 for (repeat_id in 1:2) {
  # Fresh assignments, same loaded-on-demand cache on each process invocation.
  if(repeat_id==2)file.rename(file.path(dest,'allocations'),file.path(dest,'allocations-repeat1'))
  opts <- c(paste0('--input-path=',dest),paste0('--target-csv=',target),'--scenario-id=1',
            paste0('--num-replicates=',n),'--solver=gurobi','--time-limit=120')
  if(variant!='reference')opts <- c(opts,paste0('--num-cores=',if(variant=='serial-auto')1 else sub('workers-','',variant)),
    '--draw-batch-size=4',paste0('--solver-threads=',if(variant=='serial-auto')'auto' else '1'),
    paste0('--solver-seed=',if(variant=='serial-auto')'auto' else '0'))
  seconds <- run(code,'scripts/policy/optimize-cluster-bootstrap.R',opts,file.path(dest,paste0('optimize-',repeat_id,'.log')))
  rows[[length(rows)+1L]]<-data.frame(variant=variant,repeat_id=repeat_id,draws=n,seconds=seconds,draws_per_minute=60*n/seconds,
    peak_rss_kb=as.numeric(readLines(file.path(dest,paste0("optimize-",repeat_id,".log.rss")))) )
  write.csv(do.call(rbind,rows),file.path(output,'performance.csv'),row.names=FALSE)
 }
 if(variant %in% c('serial-auto','workers-4','workers-8')) {
  ref <- file.path(output,if(variant=='serial-auto')'reference' else 'workers-1')
  run(candidate,'scripts/checks/check-policy-speedup-equivalence.R',c(paste0('--reference=',ref),
   paste0('--candidate=',dest),paste0('--output=',dest,'/comparison.csv')),file.path(dest,'comparison.log'))
 }
}
