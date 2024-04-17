
  #include base_data_sec.stan
  #include takeup_data_sec.stan
  #include wtp_data.stan
  #include beliefs_data_sec.stan
  #include dist_data_sec.stan

  int MIN_COST_MODEL_TYPE_VALUE;
  int MAX_COST_MODEL_TYPE_VALUE;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_LINEAR;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_QUADRATIC;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_DISCRETE;

  int<lower= 1, upper = 2> BELIEFS_ORDER;  
  int<lower = 0, upper = 1> use_wtp_model;
  int<lower = 0, upper = 1> use_homoskedastic_shocks;
  
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> use_cost_model;
  int<lower = 0, upper = 1> suppress_reputation;
  int<lower = 0, upper = 1> use_u_in_delta;
  
  int<lower = 0, upper = 1> use_param_dist_cluster_effects; // These are used for parameteric (linear, quadratic) distance cost models only
  int<lower = 0, upper = 1> use_param_dist_county_effects;

  int<lower = 0, upper = 1> lnorm_wtp_value_utility_prior;

  // 0 => exponential, 1 => log, 2 => linear, 3 => reserved for optim R
  // 4 => \hat{p}
  int<lower = 0, upper = 4> mu_rep_type;   
  // Rate of Change
  
  int<lower = 1, upper = num_treatments> roc_compare_treatment_id_left;
  int<lower = 1, upper = num_treatments> roc_compare_treatment_id_right;
  int<lower = 0> num_roc_distances;
  vector<lower = 0>[num_roc_distances] roc_distances; // These should be standardized already
 
  // Optim Prediction
  int<lower=0> num_B_treatments;
  int<lower=0> num_mu_treatments;
  int<lower=0> num_optim_distances;
  vector[num_optim_distances] optim_distances; // These should be standardized already too
  int<lower = 0, upper = 1> USE_MAP_IN_OPTIM;
  int<lower = 0, upper = 1> GEN_OPTIM;

  real<lower=0> raw_cluster_sd_sd; // If cluster level shock, SD of that shock SD
  int<lower=0, upper=1> use_cluster_fixed_effect;

  // Sim Delta
  
  int<lower = 0> num_sim_delta_w;
  vector[num_sim_delta_w] sim_delta_w;
  
  // Hyperparam
  real<lower = 0> mu_rep_sd;