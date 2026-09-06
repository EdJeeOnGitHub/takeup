#!/usr/bin/env Rscript
# Terminate a real draw run and verify recovery against an uninterrupted run.
source('R/policy/bootstrap.R')
args <- commandArgs(TRUE)
root <- policy_option_value(args, '--output-path', tempfile('policy-interruption-'))
solver <- policy_option_value(args, '--solver', if (nzchar(Sys.which('glpsol'))) 'glpk' else 'gurobi')
dir.create(root, recursive=TRUE, showWarnings=FALSE)
root <- normalizePath(root)
for(name in c('reference','resumed')) {
 p<-file.path(root,name);dir.create(p,showWarnings=FALSE)
 e<-data.frame(village_i=c(1L,1L,2L,2L),pot_j=c(1L,2L,1L,3L),distance=c(100,200,100,300),distance_km=c(.1,.2,.1,.3))
 saveRDS(e,file.path(p,'policy-feasible-edges.rds'))
 saveRDS(matrix(rep(c(.7,.5,.6,.4),5*200),nrow=200,byrow=TRUE),file.path(p,'policy-edge-demand-matrix.rds'))
 write.csv(data.frame(draw=1:200,replicate=1:200),file.path(p,'policy-edge-demand-draw-map.csv'),row.names=FALSE)
}
target<-file.path(root,'target.csv');write.csv(data.frame(social_welfare=1),target,row.names=FALSE)
opts<-function(p)c('--vanilla','scripts/policy/optimize-cluster-bootstrap.R',paste0('--input-path=',root,'/',p),paste0('--target-csv=',target),
 '--scenario-id=1',paste0('--solver=',solver),'--solver-threads=1','--solver-seed=0','--num-cores=2','--draw-batch-size=4')
status<-system2('timeout',c('1s',file.path(R.home('bin'),'Rscript'),opts('resumed')),stdout=file.path(root,'interrupted.log'),stderr=file.path(root,'interrupted.log'))
stopifnot(status==124)
partial<-length(list.files(file.path(root,'resumed/allocations/control'),pattern='^replicate-.*rds$'))
stopifnot(partial>0,partial<200)
for(p in c('reference','resumed'))stopifnot(system2(file.path(R.home('bin'),'Rscript'),opts(p),stdout=file.path(root,paste0(p,'.log')),stderr=file.path(root,paste0(p,'.log')))==0)
stopifnot(system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/checks/check-policy-speedup-equivalence.R',paste0('--reference=',root,'/reference'),paste0('--candidate=',root,'/resumed'),paste0('--output=',root,'/comparison.csv')))==0)
cat('Actual SIGTERM interruption after',partial,'draws: resumed 200/200 identical allocations\n')
