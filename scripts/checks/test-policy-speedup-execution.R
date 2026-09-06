#!/usr/bin/env Rscript
# Small synthetic execution-contract tests; these are not empirical results.
args <- commandArgs(TRUE)
source('R/policy/bootstrap.R'); source('R/policy/cost-sensitivity.R')
solver <- policy_option_value(args,'--solver',if(nzchar(Sys.which('glpsol')))'glpk' else 'gurobi')
root <- policy_option_value(args,'--output-path',tempfile('policy-execution-'))
dir.create(root,recursive=TRUE,showWarnings=FALSE); root <- normalizePath(root)
edges <- data.frame(village_i=c(1L,1L,2L,2L),pot_j=c(1L,2L,1L,3L),
                    distance=c(100,200,100,300),distance_km=c(.1,.2,.1,.3))
# Exact LP text preservation, including population-weighted coefficients.
a <- tempfile(); b <- tempfile()
for (target in c(.8,1.1)) {
  policy_cost_write_lp(a,edges,rep(0,4),rep(1,3),c(.7,.5,.6,.4),target)
  policy_cost_write_lp(b,edges,rep(0,4),rep(1,3),c(.7,.5,.6,.4),target,prepared=policy_cost_prepare(edges))
  stopifnot(identical(readLines(a),readLines(b)))
}
demand <- rbind(rep(c(.7,.5,.6,.4),5),rep(c(.2,.1,.2,.1),5),rep(NA_real_,20))
target <- file.path(root,'target.csv'); write.csv(data.frame(social_welfare=1),target,row.names=FALSE)
create <- function(name) {
  p <- file.path(root,name);dir.create(p)
  saveRDS(edges,file.path(p,'policy-feasible-edges.rds'))
  saveRDS(demand,file.path(p,'policy-edge-demand-matrix.rds'))
  write.csv(data.frame(draw=1:3,replicate=1:3),file.path(p,'policy-edge-demand-draw-map.csv'),row.names=FALSE)
  p
}
run <- function(path,cores=1L,extra=character(),success=TRUE) {
  code <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/policy/optimize-cluster-bootstrap.R',
    shQuote(c(paste0('--input-path=',path),paste0('--target-csv=',target),'--scenario-id=1',
    paste0('--solver=',solver),paste0('--num-cores=',cores),'--draw-batch-size=1','--solver-threads=1','--solver-seed=0',extra))),
    stdout=file.path(path,'test.log'),stderr=file.path(path,'test.log'))
  stopifnot(if(success)code==0L else code!=0L)
}
serial <- create('serial'); parallel <- create('parallel')
run(serial); run(parallel,2)
check <- function(r,c) {
  status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/checks/check-policy-speedup-equivalence.R',
    shQuote(c(paste0('--reference=',r),paste0('--candidate=',c),paste0('--output=',root,'/comparison.csv')))))
  stopifnot(status==0L)
}
check(serial,parallel)
# Simulate interrupted completion by deleting one atomic draw result, then resume.
unlink(file.path(parallel,'allocations/control/replicate-0001.rds'))
run(parallel,2); check(serial,parallel)
# Changed settings and corrupted resumable results must fail visibly.
run(parallel,2,'--time-limit=12',FALSE)
p <- file.path(parallel,'allocations/control/replicate-0001.rds')
saved <- readRDS(p); saved$contract_hash <- 'corrupt'; saveRDS(saved,p)
run(parallel,2,success=FALSE)
st <- read.csv(file.path(parallel,'allocations/control/status.csv'))
stopifnot(st$status[st$draw==1]=='failed')
# Explicit format selection prevents silently preferring a stale cache.
saveRDS(data.frame(),file.path(serial,'policy-demand-curves.rds'))
run(serial,success=FALSE)
cat('PASS: LP identity, serial/parallel, infeasible/undefined, restart, stale inputs, worker failure; ',root,'\n')
