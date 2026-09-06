#!/usr/bin/env Rscript
source("R/policy/bootstrap.R")
args <- commandArgs(TRUE)
path <- policy_option_value(args, "--input-path")
parameters <- read.csv(file.path(path, "policy-bootstrap-parameters.csv"))
parameters$model_id <- "cluster-weighted"
parameters$model_label <- "Exponential cluster-weighted modes"
parameters$model_family <- "gaussian"
parameters$chain <- 1L
parameters$iteration <- parameters$replicate
write.csv(parameters, file.path(path, "policy-model-parameters.csv"), row.names = FALSE)
saveRDS(list(parameters = parameters, cluster_shock = NULL, cluster_external_id = NULL),
        file.path(path, "policy-model-parameters.rds"), compress = FALSE)
write.csv(data.frame(model_id = "cluster-weighted", model_label = "Exponential cluster-weighted modes",
                    model_family = "gaussian", draws = nrow(parameters), chains = 1L),
          file.path(path, "policy-model-parameter-status.csv"), row.names = FALSE)
