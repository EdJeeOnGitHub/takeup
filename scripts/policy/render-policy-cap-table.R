#!/usr/bin/env Rscript
source("R/policy/bootstrap.R")
args <- commandArgs(TRUE)
path <- policy_option_value(args,"--report-path")
x <- read.csv(file.path(path,"all-model-cap-summary.csv"))
x <- x[x$model_id=="benchmark", ]
stopifnot(nrow(x)==20L,!anyDuplicated(x[c("cap_m","scenario")]))
cell <- function(row) {
  value<-if(is.finite(row$sites_median))sprintf("%.0f\\\\{}[%.0f, %.0f]",row$sites_median,row$sites_low,row$sites_high)else "--"
  if(row$target_infeasible>0)value<-paste0(value,sprintf("\\\\{\\scriptsize %.1f\\%% infeasible}",100*row$target_infeasible/row$draws))
  if(row$equilibrium_undefined>0)value<-paste0(value,sprintf("\\\\{\\scriptsize %.1f\\%% undefined}",100*row$equilibrium_undefined/row$draws))
  paste0("\\makecell[c]{",value,"}")
}
lines<-c("\\begin{tabular}{lccccc}","\\toprule",
 "Cap & Control & Bracelet & \\makecell{Control returns\\\\at 0.5 km} & \\makecell{Bracelet returns\\\\at 0.5 km} & \\makecell{No social\\\\image} \\\\","\\midrule")
for(cap in c(3500,4500,5500,10000)) {
  values<-vapply(policy_scenarios$scenario,function(s)cell(x[x$cap_m==cap & x$scenario==s, ]),character(1))
  lines<-c(lines,paste0(sprintf("%.1f km",cap/1000)," & ",paste(values,collapse=" & ")," \\\\"))
}
lines<-c(lines,"\\bottomrule","\\end{tabular}")
writeLines(lines,file.path(path,"tables/benchmark-distance-caps.tex"))
