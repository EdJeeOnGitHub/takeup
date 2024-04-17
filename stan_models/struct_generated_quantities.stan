#include dist_generated_quantities.stan
#include beliefs_generated_quantities.stan

  matrix[num_clusters, num_dist_group_treatments] structural_cluster_benefit = 
        rep_matrix(structural_treatment_effect', num_clusters) + 
        rep_matrix((structural_beta_cluster + structural_beta_county[cluster_county_id]), num_dist_group_treatments);
  
  array[num_dist_group_treatments] vector[num_clusters] cluster_cf_benefit_cost; 
  
  array[num_dist_group_treatments, num_treatments] vector[num_clusters] cluster_cf_cutoff; 
  
  vector[generate_rep ? num_clusters : 0] cluster_rep_benefit_cost; 
  vector[generate_rep ? num_clusters : 0] cluster_rep_cutoff; 
  
  array[num_clusters] vector[num_sim_delta_w] sim_delta;
  
  real wtp_travel_dist = dist_beta_v[1] / wtp_value_utility;
  real calendar_preference_in_dist = hyper_wtp_mu / wtp_travel_dist;
  
#include wtp_generated_quantities.stan
#include takeup_struct_cv.stan
#include takeup_struct_quantities.stan

  {
    int treatment_cluster_pos = 1;
    int cluster_treatment_cf_pos = 1;
    
    for (treatment_index in 1:num_dist_group_treatments) {
      vector[num_clusters] treatment_dist_cost;
     
      int curr_assigned_treatment = cluster_treatment_map[treatment_index, 1];
      int curr_assigned_dist_group = cluster_treatment_map[treatment_index, 2];
      
      treatment_dist_cost = param_dist_cost(
                                            all_cluster_standard_dist[, curr_assigned_dist_group],
                                            cluster_linear_dist_cost[, treatment_index],
                                            cluster_quadratic_dist_cost[, treatment_index]);
                                                                 
      cluster_cf_benefit_cost[treatment_index] = structural_cluster_benefit[, treatment_index] - treatment_dist_cost;
      
      for (mu_treatment_index in 1:num_treatments) {
        vector[num_clusters] curr_cluster_mu_rep = calculate_mu_rep(
          { mu_treatment_index }, all_cluster_standard_dist[, curr_assigned_dist_group], 
          base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_beliefs, centered_cluster_dist_beta_beliefs,
          mu_rep_type);
        
        if (multithreaded) {
          if (use_cluster_effects) {
            cluster_cf_cutoff[treatment_index, mu_treatment_index] = map_find_fixedpoint_solution(
              cluster_cf_benefit_cost[treatment_index],
              curr_cluster_mu_rep,
              total_error_sd,
              u_sd,
              
              use_u_in_delta, 
              alg_sol_rel_tol, 
              alg_sol_f_tol, 
              alg_sol_max_steps
            );
          } else {
            cluster_cf_cutoff[treatment_index, mu_treatment_index] = map_find_fixedpoint_solution(
              cluster_cf_benefit_cost[treatment_index],
              curr_cluster_mu_rep,
              rep_vector(total_error_sd[curr_assigned_treatment], num_clusters),
              rep_vector(u_sd[curr_assigned_treatment], num_clusters),
              
              use_u_in_delta, 
              alg_sol_rel_tol, 
              alg_sol_f_tol, 
              alg_sol_max_steps
            );
          }
        } else {
          for (cluster_index in 1:num_clusters) {
            if (use_cluster_effects) {
              cluster_cf_cutoff[treatment_index, mu_treatment_index, cluster_index] = find_fixedpoint_solution(
                cluster_cf_benefit_cost[treatment_index, cluster_index],
                curr_cluster_mu_rep[cluster_index],
                total_error_sd[cluster_index],
                u_sd[cluster_index],
                
                use_u_in_delta, 
                alg_sol_rel_tol, 
                alg_sol_f_tol, 
                alg_sol_max_steps
              ); 

            } else {
              cluster_cf_cutoff[treatment_index, mu_treatment_index, cluster_index] = find_fixedpoint_solution(
                cluster_cf_benefit_cost[treatment_index, cluster_index],
                curr_cluster_mu_rep[cluster_index],
                total_error_sd[curr_assigned_treatment],
                u_sd[curr_assigned_treatment],
                
                use_u_in_delta, 
                alg_sol_rel_tol, 
                alg_sol_f_tol, 
                alg_sol_max_steps
              ); 
            }
          }
        }
      }
    }
  }
  
  // Replicated Data  
  if (generate_rep) {
    for (cluster_index in 1:num_clusters) {
      vector[num_dist_group_treatments] rep_beta_cluster = rep_vector(0, num_dist_group_treatments);
        
      vector[num_dist_group_treatments] rep_beta_county = rep_vector(0, num_dist_group_treatments);
      
      int curr_assigned_treatment = cluster_incentive_treatment_id[cluster_index];
      
      real rep_dist_cost;
      
      real rep_linear_dist_cost = linear_dist_cost[cluster_assigned_dist_group_treatment[cluster_index]];
      real rep_quadratic_dist_cost = quadratic_dist_cost[cluster_assigned_dist_group_treatment[cluster_index]];

    rep_beta_cluster[1:num_treatments] =  rep_vector(0, num_treatments);
     if (use_cluster_fixed_effect) {
      rep_beta_cluster[1:num_treatments] =  to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_cluster_sd'));
     } 
      rep_beta_cluster[(num_treatments + 2):num_dist_group_treatments] = rep_beta_cluster[2:num_treatments];
      
      rep_beta_county[1:num_treatments] = use_county_effects ? to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_county_sd')) : rep_vector(0, num_treatments);
      rep_beta_county[(num_treatments + 2):num_dist_group_treatments] = rep_beta_county[2:num_treatments];
        
      if (use_param_dist_cluster_effects) {
        rep_linear_dist_cost += normal_rng(0, dist_beta_cluster_sd[1]); 
      }
      
      if (use_param_dist_county_effects) {
        rep_linear_dist_cost += normal_rng(0, dist_beta_county_sd[1]); 
      }
        
      rep_dist_cost = param_dist_cost(
                                      [ cluster_standard_dist[cluster_index] ]', 
                                      [ rep_linear_dist_cost ]',
                                      [ rep_quadratic_dist_cost ]')[1]; 
      
      cluster_rep_benefit_cost[cluster_index] = structural_treatment_effect[cluster_assigned_dist_group_treatment[cluster_index]] 
        + treatment_map_design_matrix[cluster_assigned_dist_group_treatment[cluster_index]] * (rep_beta_cluster + rep_beta_county) - rep_dist_cost;
        
      cluster_rep_cutoff[cluster_index] = find_fixedpoint_solution(
          cluster_rep_benefit_cost[cluster_index],
          obs_cluster_mu_rep[cluster_index],
          total_error_sd[curr_assigned_treatment],
          u_sd[curr_assigned_treatment],
          
          use_u_in_delta,
          alg_sol_rel_tol,
          alg_sol_f_tol,
          alg_sol_max_steps
        );
    }
  }

  if (use_cluster_effects) {
    for (cluster_idx in 1:num_clusters) {
      for (sim_delta_index in 1:num_sim_delta_w) {
        sim_delta[cluster_idx, sim_delta_index] = expected_delta(
          sim_delta_w[sim_delta_index], 
          total_error_sd[cluster_idx], 
          u_sd[cluster_idx], 
          dummy_xr, 
          dummy_xi);
      }
    }
  } else {
    for (sim_delta_index in 1:num_sim_delta_w) {
      sim_delta[1:num_clusters, sim_delta_index] = rep_array(
        expected_delta(
          sim_delta_w[sim_delta_index], 
          total_error_sd[1], 
          u_sd[1], 
          dummy_xr, 
          dummy_xi), 
          num_clusters
      );
    }
  } 

  if (GEN_OPTIM) {
    #include takeup_optim_quantities.stan
  }
