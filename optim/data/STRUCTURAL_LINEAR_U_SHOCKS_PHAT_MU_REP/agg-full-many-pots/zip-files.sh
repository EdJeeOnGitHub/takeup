#!/bin/bash

# Create zips for posterior files
zip -m posterior_files.zip posterior-clean-summ-optim.csv posterior-suppress-rep-summ-optim.csv

# Create zips for pred-demand-dist-fit86- files grouped by their parameters
zip -m pred_demand_dist_fit86_bracelet.zip pred-demand-dist-fit86-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-delta0-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-delta1250-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-delta2500-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-mu0-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-mu1250-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-mu2500-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-no-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-static-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-suppress-rep-cutoff-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv

zip -m pred_demand_dist_fit86_calendar.zip pred-demand-dist-fit86-cutoff-b-calendar-mu-calendar-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-static-cutoff-b-calendar-mu-calendar-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-suppress-rep-cutoff-b-calendar-mu-calendar-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv

zip -m pred_demand_dist_fit86_control.zip pred-demand-dist-fit86-cutoff-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-cutoff-b-control-mu-calendar-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-cutoff-b-control-mu-ink-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-no-cutoff-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-no-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-static-cutoff-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-static-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-suppress-rep-cutoff-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-suppress-rep-cutoff-b-control-mu-calendar-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-suppress-rep-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-suppress-rep-cutoff-b-control-mu-ink-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv

zip -m pred_demand_dist_fit86_ink.zip pred-demand-dist-fit86-cutoff-b-ink-mu-ink-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-static-cutoff-b-ink-mu-ink-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    pred-demand-dist-fit86-suppress-rep-cutoff-b-ink-mu-ink-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv

# Create zips for remaining .zip files
zip -m pred_demand_dist_fit86.zip pred-demand-dist-fit86.zip
zip -m pred_demand_dist_fit86_2.zip pred-demand-dist-fit86-2.zip

# Create zips for summ- files
zip -m summ_agg_files.zip summ-agg-identity-experiment-target-constraint.csv \
    summ-agg-log-experiment-target-constraint.csv \
    summ-agg-log-target-closest-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    summ-agg-log-target-rep-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv \
    summ-agg-log-target-suppress-rep-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv

# Create zips for rep-target and target-rep files
zip -m rep_target_files.zip rep-target-rep-agg-log-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-experimental-control-allocation-data.rds \
    target-rep-agg-identity-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-experimental-control-allocation-data.rds

