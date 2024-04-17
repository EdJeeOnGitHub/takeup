
#include dist_transformed_parameters.stan
#include wtp_transformed_parameters.stan
#include beliefs_transformed_parameters.stan

  matrix[num_clusters, num_treatments] centered_cluster_beta_beliefs;
  matrix[num_clusters, num_treatments] centered_cluster_dist_beta_beliefs;

  vector[num_dist_group_treatments] beta;
  
  vector[num_dist_group_treatments] structural_treatment_effect;
  vector[num_clusters] structural_cluster_benefit_cost;
  vector[num_clusters] structural_beta_cluster = rep_vector(0, num_clusters);
  vector[num_counties] structural_beta_county = rep_vector(0, num_counties);
 
  vector[suppress_reputation ? 0 : num_clusters] obs_cluster_mu_rep;
  
  vector[num_clusters] structural_cluster_obs_v = rep_vector(0, num_clusters);
  row_vector<lower = 0, upper = 1>[fit_model_to_data ? num_clusters : 0] structural_cluster_takeup_prob;
  
  vector[num_dist_group_treatments] linear_dist_cost = rep_vector(0, num_dist_group_treatments);
  vector[num_dist_group_treatments] quadratic_dist_cost = rep_vector(0, num_dist_group_treatments);
  matrix[num_clusters, num_dist_group_treatments] cluster_linear_dist_cost = rep_matrix(0, num_clusters, num_dist_group_treatments);
  matrix[num_clusters, num_dist_group_treatments] cluster_quadratic_dist_cost = rep_matrix(0, num_clusters, num_dist_group_treatments);
  vector[num_clusters] cluster_dist_cost = rep_vector(0, num_clusters);
  
  vector<lower = 0>[use_cluster_effects ?  num_clusters : num_treatments] u_sd;
  vector<lower = 0>[use_cluster_effects ? num_clusters : num_treatments] total_error_sd;
  vector<lower=0>[use_cluster_effects ? num_clusters : 0] raw_cluster_sd;


  if (BELIEFS_ORDER == 1) {
    centered_cluster_beta_beliefs = centered_cluster_beta_1ord;
    centered_cluster_dist_beta_beliefs = centered_cluster_dist_beta_1ord;
  }
  if (BELIEFS_ORDER == 2) {
    centered_cluster_beta_beliefs = centered_cluster_beta_2ord;
    centered_cluster_dist_beta_beliefs = centered_cluster_dist_beta_2ord;
  }

  raw_cluster_sd = raw_cluster_sd_tilde .* raw_cluster_sd_sd;
  if (use_cluster_effects) {
    u_sd = raw_u_sd[1] + raw_cluster_sd;
  }
  if (use_homoskedastic_shocks && use_cluster_effects == 0) {
    u_sd = rep_vector(raw_u_sd[1], num_treatments);
  } 
  if (use_homoskedastic_shocks == 0 && use_cluster_effects == 0) {
    u_sd = raw_u_sd;
  }
    
  total_error_sd = sqrt(1 + square(u_sd));
  
  for (dist_index in 1:num_discrete_dist) {
    if (dist_index > 1) {
      beta[(num_treatments + 1):] = rep_vector(0, num_treatments); 
    } else if (use_wtp_model) { 
      beta[1:2] = [ beta_intercept, beta_ink_effect ]';
      beta[CALENDAR_TREATMENT_INDEX] = beta_bracelet_effect + wtp_value_utility * hyper_wtp_mu;
      beta[BRACELET_TREATMENT_INDEX] = beta_bracelet_effect;
    } else {
      beta[1:num_treatments] = [ beta_intercept, beta_ink_effect, beta_calendar_effect, beta_bracelet_effect ]';
    }
  }
 
  structural_treatment_effect = treatment_map_design_matrix * beta;
  
  // Levels: control ink calendar bracelet
 
  if (!suppress_reputation) { 
    obs_cluster_mu_rep = calculate_mu_rep(
      cluster_incentive_treatment_id, cluster_standard_dist, 
      base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_beliefs, centered_cluster_dist_beta_beliefs,
      mu_rep_type);
  }

  linear_dist_cost = rep_vector(dist_beta_v[1], num_dist_group_treatments);
  
  if (use_param_dist_cluster_effects) {
    cluster_linear_dist_cost += rep_matrix(dist_beta_cluster_raw[, 1] * dist_beta_cluster_sd[1], num_dist_group_treatments);
  }
  
  if (use_param_dist_county_effects) {
    cluster_linear_dist_cost += rep_matrix((dist_beta_county_raw[, 1] * dist_beta_county_sd[1])[cluster_county_id], num_dist_group_treatments);
  }
  
  if (num_dist_param_quadratic > 0) {
    quadratic_dist_cost = rep_vector(dist_quadratic_beta_v[1], num_dist_group_treatments);
    
    if (use_param_dist_cluster_effects) {
      // TODO
    }
    
    if (use_param_dist_county_effects) {
      // TODO
    }
  }
  
  cluster_linear_dist_cost += rep_matrix(linear_dist_cost', num_clusters);
  cluster_quadratic_dist_cost += rep_matrix(quadratic_dist_cost', num_clusters);
  
    
  cluster_dist_cost = param_dist_cost(
                                      cluster_standard_dist, 
                                      to_vector(cluster_linear_dist_cost)[long_cluster_by_treatment_index],
                                      to_vector(cluster_quadratic_dist_cost)[long_cluster_by_treatment_index]);
  
  structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_dist_group_treatment] - cluster_dist_cost;

  // If we're using cluster effects, we want u_{ic} = e_{i} + \gamma_c
  // i.e. latent shock has common cluster component.
  // This code is adding another cluster on top of that directly however
  // We're double counting the cluster effect: cluster level latent shock +
  // cluster effect directly on beta 
  if (use_cluster_fixed_effect) {
    //use_cluster_effects
    // for now turn off beta level shock
    structural_beta_cluster = structural_beta_cluster_raw * structural_beta_cluster_sd[1];
    structural_cluster_benefit_cost += structural_beta_cluster;
  }
  
  if (use_county_effects) {
    structural_beta_county = structural_beta_county_raw * structural_beta_county_sd[1];
    
    structural_cluster_benefit_cost += structural_beta_county[cluster_county_id];
  }

  if (fit_model_to_data) {
    if (suppress_reputation) {
      structural_cluster_obs_v = - structural_cluster_benefit_cost;
    } else {
      if (multithreaded) {
        if (use_cluster_effects) {
          structural_cluster_obs_v = map_find_fixedpoint_solution(
            structural_cluster_benefit_cost, 
            obs_cluster_mu_rep,
            total_error_sd, // indexed by cluster_id
            u_sd, // indexed by cluster_id
            
            use_u_in_delta,
            alg_sol_rel_tol, // 1e-10,
            alg_sol_f_tol, // 1e-5,
            alg_sol_max_steps
          ); 
        } else {
          structural_cluster_obs_v = map_find_fixedpoint_solution(
            structural_cluster_benefit_cost, 
            obs_cluster_mu_rep,
            total_error_sd[cluster_incentive_treatment_id],
            u_sd[cluster_incentive_treatment_id],
            
            use_u_in_delta,
            alg_sol_rel_tol, // 1e-10,
            alg_sol_f_tol, // 1e-5,
            alg_sol_max_steps
          ); 
        }
      } else {
        for (cluster_index in 1:num_clusters) {
          if (use_cluster_effects) {
            structural_cluster_obs_v[cluster_index] = find_fixedpoint_solution(
              structural_cluster_benefit_cost[cluster_index],
              obs_cluster_mu_rep[cluster_index],
              total_error_sd[cluster_index],
              u_sd[cluster_index],

              use_u_in_delta,
              alg_sol_rel_tol, // 1e-10,
              alg_sol_f_tol, // 1e-5,
              alg_sol_max_steps
            );
          } else {
            structural_cluster_obs_v[cluster_index] = find_fixedpoint_solution(
              structural_cluster_benefit_cost[cluster_index],
              obs_cluster_mu_rep[cluster_index],
              total_error_sd[cluster_incentive_treatment_id[cluster_index]],
              u_sd[cluster_incentive_treatment_id[cluster_index]],

              use_u_in_delta,
              alg_sol_rel_tol, // 1e-10,
              alg_sol_f_tol, // 1e-5,
              alg_sol_max_steps
            );
          }
        }
      }
    }
    if (use_cluster_effects) {
      structural_cluster_takeup_prob = Phi_approx(- structural_cluster_obs_v ./ total_error_sd)';
    } else {
      // Index clusters by their treatment (in case error sd varying w/ treatment)
      structural_cluster_takeup_prob = Phi_approx(- structural_cluster_obs_v ./ total_error_sd[cluster_incentive_treatment_id])';
    }
  }