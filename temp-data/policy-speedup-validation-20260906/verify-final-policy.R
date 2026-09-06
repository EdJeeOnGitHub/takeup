source('R/policy/bootstrap.R')
base <- '/project/akaring/takeup-data/scratch/policy-speedup-20260906'
source_root <- file.path(base,'equivalence')
output <- file.path(base,'final-code-verification')
target <- '/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv'
run <- function(script,args,log) {
 status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',script,shQuote(args)),stdout=log,stderr=log)
 if(status)stop(log)
}
for (model in c('benchmark','cluster-shock','full-information','tight-multinomial')) {
 ref <- file.path(source_root,model,'reference')
 for (variant in c('default','parallel')) {
  dest <- file.path(output,model,variant);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  if(variant=='default') {
   file.copy(file.path(ref,'policy-model-parameters.rds'),dest)
   run('scripts/policy/predict-model-robustness.R',c(paste0('--parameter-rds=',dest,'/policy-model-parameters.rds'),
    paste0('--output-path=',dest),'--distance-data=/project/akaring/takeup-data/optim/data/full-many-pots-experiment.rds',
    '--household-workspace=/project/akaring/takeup-data/data/stan_analysis_data/dist_fit104.RData','--num-cores=2'),file.path(dest,'predict.log'))
  } else file.copy(list.files(file.path(output,model,'default'),pattern='^policy-',full.names=TRUE),dest)
  for (scenario in 1:5) {
   opts <- c(paste0('--input-path=',dest),paste0('--target-csv=',target),paste0('--scenario-id=',scenario),'--solver=gurobi', '--time-limit=120','--num-replicates=2')
   if(variant=='parallel')opts<-c(opts,'--num-cores=8','--solver-threads=1','--solver-seed=0','--draw-batch-size=1')
   run('scripts/policy/optimize-cluster-bootstrap.R',opts,file.path(dest,paste0('optimize-',scenario,'.log')))
  }
  other <- if(variant=='default')ref else file.path(source_root,model,'serial')
  run('scripts/checks/check-policy-speedup-equivalence.R',c(paste0('--reference=',other),paste0('--candidate=',dest),
    paste0('--output=',dest,'/comparison.csv')),file.path(dest,'comparison.log'))
 }
}
