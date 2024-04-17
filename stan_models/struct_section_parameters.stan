
#include dist_parameters_sec.stan
#include wtp_parameters.stan
#include beliefs_parameters_sec.stan
  
  // Levels: control ink calendar bracelet
  real beta_intercept;
  real beta_ink_effect;
  real<lower = (use_private_incentive_restrictions ? 0 : negative_infinity())> beta_calendar_effect;
  real<lower = (use_private_incentive_restrictions ? 0 : negative_infinity())> beta_bracelet_effect;
  
  vector[use_cluster_fixed_effect ? num_clusters : 0] structural_beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_fixed_effect ? 1 : 0] structural_beta_cluster_sd;
  
  vector[use_county_effects ? num_counties : 0] structural_beta_county_raw;
  row_vector<lower = 0>[use_county_effects ? 1 : 0] structural_beta_county_sd;
  
  // Reputational Returns
  
  real<lower = 0> base_mu_rep;
  // vector<lower = 0>[use_cluster_effects ? num_clusters : 1] raw_u_sd;
  vector<lower = 0>[1] raw_u_sd;
  vector<lower=0>[use_cluster_effects ? num_clusters : 0] raw_cluster_sd_tilde;
  
  // Linear Parametric Cost
  
  vector[num_dist_param] dist_beta_v; // Linear distance*treatment effects
  
  matrix[use_param_dist_cluster_effects ? num_clusters : 0, num_dist_param] dist_beta_cluster_raw;
  row_vector<lower = 0>[use_param_dist_cluster_effects ? num_dist_param : 0] dist_beta_cluster_sd;
  
  matrix[use_param_dist_county_effects ? num_counties : 0, num_dist_param] dist_beta_county_raw;
  row_vector<lower = 0>[use_param_dist_county_effects ? num_dist_param : 0] dist_beta_county_sd;
  
  // Quadratic Cost Model
  
  vector<lower = 0>[num_dist_param_quadratic] dist_quadratic_beta_v; // Quadratic distance*treatment effects
  
  // WTP valuation parameters
  real<lower = 0> wtp_value_utility;