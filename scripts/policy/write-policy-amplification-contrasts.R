#!/usr/bin/env Rscript
source("R/policy/bootstrap.R")
args <- commandArgs(TRUE)
path <- policy_option_value(args,"--report-path")
x <- read.csv(file.path(path,"all-model-cap-replicates.csv"))
rows <- list(); summaries <- list()
for (panel in split(x,paste(x$model_id,x$cap_m))) {
  control <- panel[panel$scenario=="control", ]
  get <- function(scenario) {
    z<-panel[panel$scenario==scenario, ];z<-z[match(control$draw,z$draw), ]
    stopifnot(identical(z$draw,control$draw));z
  }
  bracelet<-get("bracelet");sc<-get("static-control");sb<-get("static-bracelet")
  valid<-control$status=="complete" & bracelet$status=="complete" & sc$status=="complete" & sb$status=="complete"
  endogenous<-control$n_pot-bracelet$n_pot;static<-sc$n_pot-sb$n_pot
  delta<-endogenous-static
  z<-data.frame(model_id=control$model_id,cap_m=control$cap_m,draw=control$draw,replicate=control$replicate,
    all_four_targets_met=valid,endogenous_sites_saved=endogenous,static_sites_saved=static,
    amplification_sites_saved_including_fallbacks=delta,amplification_sites_saved_at_preserved_targets=ifelse(valid,delta,NA_real_),
    any_undefined=control$status=="equilibrium_undefined" | bracelet$status=="equilibrium_undefined" | sc$status=="equilibrium_undefined" | sb$status=="equilibrium_undefined",
    any_infeasible=control$status=="target_infeasible" | bracelet$status=="target_infeasible" | sc$status=="target_infeasible" | sb$status=="target_infeasible")
  rows[[length(rows)+1L]]<-z
  q<-if(any(valid))as.numeric(quantile(delta[valid],c(.025,.5,.975)))else rep(NA_real_,3)
  summaries[[length(summaries)+1L]]<-data.frame(model_id=z$model_id[1],cap_m=z$cap_m[1],draws=nrow(z),all_four_targets_met=sum(valid),undefined=sum(z$any_undefined),infeasible=sum(z$any_infeasible),low=q[1],median=q[2],high=q[3],probability_positive=if(any(valid))mean(delta[valid]>0)else NA_real_)
}
write.csv(do.call(rbind,rows),file.path(path,"policy-amplification-paired-contrasts.csv"),row.names=FALSE)
write.csv(do.call(rbind,summaries),file.path(path,"policy-amplification-summary.csv"),row.names=FALSE)
