#!/usr/bin/env Rscript
# Exercise the actual Stan derivative, including legacy nesting, in compiled code.
library(cmdstanr)
root <- normalizePath("stan_models")
stan <- paste0('functions {\n#include struct_section_functions.stan\n#include core_asymmetric_observability_functions.stan\n}\ngenerated quantities {\n real max_error = 0;\n for (truth in 1:2) {\n  for (j in 1:20) {\n   real d = -1 + j * 0.15;\n   real h = 1e-5;\n   real ds = -0.7;\n   real acs = 0.8;\n   real p = inv_logit(0.3 + ds*d);\n   real a = inv_logit(-0.2 + acs*d);\n   vector[3] analytic = core_two_stage_report_row_derivative(p, ds, a, acs, truth);\n   vector[3] numerical = (core_two_stage_report_row(inv_logit(0.3+ds*(d+h)), inv_logit(-0.2+acs*(d+h)),truth) - core_two_stage_report_row(inv_logit(0.3+ds*(d-h)), inv_logit(-0.2+acs*(d-h)),truth))/(2*h);\n   max_error = fmax(max_error,max(abs(analytic-numerical)));\n   max_error = fmax(max_error,max(abs(core_two_stage_report_row_derivative(p,ds,a,0,truth)-core_two_stage_report_row_derivative(p,ds,a,truth))));\n   max_error = fmax(max_error,abs(sum(analytic)));\n  }\n }\n if (max_error > 1e-8) reject("Accuracy derivative test failed: ",max_error);\n}\n')
f <- write_stan_file(stan,dir=tempdir())
m <- cmdstan_model(f,include_paths=root,quiet=TRUE)
r <- m$sample(data=list(),fixed_param=TRUE,chains=1,iter_sampling=1,iter_warmup=0,refresh=0)
stopifnot(r$summary("max_error")$mean < 1e-8)
cat("Compiled accuracy-gradient derivative and legacy nesting: OK\n")
