#!/usr/bin/env Rscript
# Compare assignments and economics, ignoring timing and serialization metadata.
args <- commandArgs(TRUE)
source('R/policy/bootstrap.R')
reference <- policy_option_value(args, '--reference')
candidate <- policy_option_value(args, '--candidate')
output <- policy_option_value(args, '--output', 'policy-speedup-comparison.csv')
if (is.null(reference) || is.null(candidate)) stop('--reference and --candidate are required')
near <- function(a,b,tol=1e-8) isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE))
files <- function(root) sort(list.files(file.path(root,'allocations'), pattern='^replicate-[0-9]+[.]rds$',recursive=TRUE))
rf <- files(reference); cf <- files(candidate)
if (!length(rf) || !identical(rf,cf)) stop('Draw/scenario file inventories differ or are empty')
edges <- readRDS(file.path(reference,'policy-feasible-edges.rds'))
stopifnot(identical(edges,readRDS(file.path(candidate,'policy-feasible-edges.rds'))))
read_predictions <- function(root) {
  path <- file.path(root,'policy-edge-demand-matrix.rds')
  if (file.exists(path)) return(list(values=readRDS(path),map=read.csv(file.path(root,'policy-edge-demand-draw-map.csv'))))
  x <- readRDS(file.path(root,'policy-demand-curves.rds'))
  draws <- sort(unique(x$draw))
  values <- do.call(rbind,lapply(draws,function(d) unlist(lapply(1:5,function(s) {
    z <- x[x$draw==d & x$scenario_id==s,]
    if ('village_i' %in% names(z) && all(!is.na(z$village_i))) {
      stopifnot(identical(as.integer(z$village_i),as.integer(edges$village_i)),near(z$distance,edges$distance,0))
      z$demand
    } else z$demand[match(edges$distance,z$distance)]
  }),use.names=FALSE)))
  list(values=values,map=data.frame(draw=draws,replicate=vapply(draws,function(d)x$replicate[match(d,x$draw)],numeric(1))))
}
rp <- read_predictions(reference); cp <- read_predictions(candidate)
stopifnot(identical(rp$map$draw,cp$map$draw),near(rp$map$replicate,cp$map$replicate,0),
          identical(is.na(rp$values),is.na(cp$values)),near(rp$values,cp$values,1e-12))
validate <- function(x, predictions) {
  a <- x$allocation; s <- x$status
  demand <- predictions$values[match(s$draw,predictions$map$draw),(s$scenario_id-1L)*nrow(edges)+seq_len(nrow(edges))]
  if (s$status=='equilibrium_undefined') return(nrow(a)==0 && any(!is.finite(demand)))
  k <- match(paste(a$village_i,a$pot_j),paste(edges$village_i,edges$pot_j))
  if (anyNA(k) || anyDuplicated(a$village_i) || !setequal(a$village_i,edges$village_i)) return(FALSE)
  if (!near(a$distance,edges$distance[k],0) || !near(a$demand,demand[k],1e-12)) return(FALSE)
  achieved <- sum(a$demand)
  if (!near(achieved,s$achieved_welfare) || length(unique(a$pot_j))!=s$n_pot ||
      !near(mean(a$demand),s$mean_demand) || !near(mean(a$distance),s$mean_distance)) return(FALSE)
  if (s$status=='target_infeasible') {
    maximum <- sum(vapply(split(demand,edges$village_i),max,numeric(1)))
    return(maximum+1e-5 < s$target_welfare && near(achieved,maximum))
  }
  s$status=='complete' && achieved+1e-5 >= s$target_welfare
}
rows <- lapply(rf,function(f) {
  r <- readRDS(file.path(reference,'allocations',f)); c <- readRDS(file.path(candidate,'allocations',f))
  status_fields <- setdiff(names(r$status),c('elapsed_seconds','error'))
  summaries <- near(r$status[,status_fields,drop=FALSE],c$status[,status_fields,drop=FALSE])
  assignments <- identical(as.integer(r$allocation$pot_j),as.integer(c$allocation$pot_j)) &&
                 identical(as.integer(r$allocation$village_i),as.integer(c$allocation$village_i))
  economic_fields <- intersect(c('population','expected_takers','signal_cost','travel_cost','demand','distance'),names(r$allocation))
  economics <- all(economic_fields %in% names(c$allocation)) && near(r$allocation[,economic_fields,drop=FALSE],c$allocation[,economic_fields,drop=FALSE])
  # Candidate retains independently parsed optimality evidence. Old reference
  # solver logs are verified by the pilot driver when running frozen code.
  certified <- c$status$status!='complete' || (isTRUE(c$solver$optimal) && near(c$solver$objective,c$status$n_pot))
  valid <- validate(r,rp) && validate(c,cp)
  data.frame(file=f,assignments_equal=assignments,summaries_equal=summaries,
             economics_equal=economics,valid=valid,candidate_optimal=certified,
             pass=assignments && summaries && economics && valid && certified)
})
report <- do.call(rbind,rows)
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
write.csv(report,output,row.names=FALSE)
if (!all(report$pass)) stop('Equivalence failed: ',paste(report$file[!report$pass],collapse=', '))
cat('PASS:',nrow(report),'matched allocations; predictions agree\n')
