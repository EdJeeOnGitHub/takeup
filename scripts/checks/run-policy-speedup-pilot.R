#!/usr/bin/env Rscript
# Run in an allocated Midway session. Never writes to production caches.
args <- commandArgs(TRUE)
source('R/policy/bootstrap.R')
reference_code <- normalizePath(policy_option_value(args,'--reference-code'),mustWork=TRUE)
output <- policy_option_value(args,'--output-path')
source_root <- policy_option_value(args,'--source-root','/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-model-robustness-1451')
data_root <- '/project/akaring/takeup-data'
distance_data <- file.path(data_root,'optim/data/full-many-pots-experiment.rds')
target <- file.path(data_root,'optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv')
candidate_code <- getwd()
dir.create(output,recursive=TRUE,showWarnings=FALSE); output <- normalizePath(output)
run <- function(code,script,options,log) {
  old <- setwd(code); on.exit(setwd(old))
  started <- proc.time()[[3L]]
  status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(script),shQuote(options)),stdout=log,stderr=log)
  elapsed <- proc.time()[[3L]]-started
  if (status!=0) stop('Command failed: ',log)
  elapsed
}
models <- c('benchmark','cluster-shock','full-information','tight-multinomial')
results <- list(); comparisons <- list()
for (model in models) {
  object <- readRDS(file.path(source_root,model,'policy-model-parameters.rds'))
  indices <- 1:2
  if (model=='tight-multinomial') {
    st <- read.csv(file.path(source_root,model,'allocations/suppress-reputation/status.csv'))
    ids <- head(st$draw[st$status=='target_infeasible'],2)
    if (length(ids)) indices <- unique(c(match(ids,object$parameters$draw),1L))[1:2]
  }
  object$parameters <- object$parameters[indices,,drop=FALSE]
  if (!is.null(object$cluster_shock)) object$cluster_shock <- object$cluster_shock[indices,,drop=FALSE]
  for (variant in c('reference','serial-auto','serial','parallel')) {
    dest <- file.path(output,model,variant);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    saveRDS(object,file.path(dest,'policy-model-parameters.rds'),compress=FALSE)
    write.csv(object$parameters,file.path(dest,'selected-parameters.csv'),row.names=FALSE)
    code <- if (variant=='reference') reference_code else candidate_code
    timing <- run(code,'scripts/policy/predict-model-robustness.R',c(
      paste0('--parameter-rds=',dest,'/policy-model-parameters.rds'),paste0('--distance-data=',distance_data),
      paste0('--output-path=',dest),'--distance-cap=3500','--num-cores=2',
      paste0('--household-workspace=',data_root,'/data/stan_analysis_data/dist_fit104.RData')),
      file.path(dest,'predict.log'))
    for (scenario in 1:5) {
      opts <- c(paste0('--input-path=',dest),paste0('--target-csv=',target),paste0('--scenario-id=',scenario),
                '--num-replicates=2','--solver=gurobi','--time-limit=120')
      if (variant!='reference') opts <- c(opts,paste0('--num-cores=',if(variant=='parallel')2 else 1),
        '--draw-batch-size=1',paste0('--solver-threads=',if(variant=='serial-auto')'auto' else '1'),
        paste0('--solver-seed=',if(variant=='serial-auto')'auto' else '0'))
      elapsed <- run(code,'scripts/policy/optimize-cluster-bootstrap.R',opts,file.path(dest,paste0('optimize-',scenario,'.log')))
      results[[length(results)+1L]] <- data.frame(model=model,variant=variant,scenario=scenario,seconds=elapsed)
    }
    if (variant=='reference') {
      logs <- list.files(file.path(dest,'allocations'),pattern='^reference-.*[.]log$',recursive=TRUE,full.names=TRUE)
      saved <- list.files(file.path(dest,'allocations'),pattern='^replicate-.*[.]rds$',recursive=TRUE,full.names=TRUE)
      complete <- sum(vapply(saved,function(p)readRDS(p)$status$status=='complete',logical(1)))
      stopifnot(length(logs)==complete,all(vapply(logs,function(p)any(grepl('^Optimal solution found',readLines(p))),logical(1))))
    }
    if (variant!='reference') {
      previous <- switch(variant,'serial-auto'='reference','serial'='serial-auto','parallel'='serial')
      comparison <- file.path(output,model,paste0(previous,'-vs-',variant,'.csv'))
      # Keep collecting all variants if different solver thread settings choose ties.
      tryCatch(run(candidate_code,'scripts/checks/check-policy-speedup-equivalence.R',c(
        paste0('--reference=',output,'/',model,'/',previous),paste0('--candidate=',dest),paste0('--output=',comparison)),
        file.path(dest,'comparison.log')),error=function(e)message(conditionMessage(e)))
      comparisons[[length(comparisons)+1L]] <- data.frame(model=model,comparison=paste0(previous,'-vs-',variant),
        pass=file.exists(comparison) && all(read.csv(comparison)$pass))
    }
  }
}
write.csv(do.call(rbind,results),file.path(output,'timings.csv'),row.names=FALSE)
write.csv(do.call(rbind,comparisons),file.path(output,'comparisons.csv'),row.names=FALSE)
if (!all(do.call(rbind,comparisons)$pass)) stop('Some matched comparisons failed; inspect comparisons.csv')
