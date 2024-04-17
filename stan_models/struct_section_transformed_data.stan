
#include base_transformed_data.stan 
#include wtp_transformed_data.stan
#include takeup_transformed_data.stan
#include beliefs_transformed_data.stan

  array[1] real dummy_xr = { 1.0 }; 
  array[1] int dummy_xi = { 1 }; 
  
  int<lower = 0, upper = num_treatments> num_treatment_shocks = num_treatments;
  
  int num_dist_param = 1;
  int num_dist_param_quadratic = use_cost_model == COST_MODEL_TYPE_PARAM_QUADRATIC ? 1 : 0; 
  
  if (num_age_groups > 1) {
    reject("Age groups not suported in structural model.");
  }