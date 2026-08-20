// Sampling-only entry point for the individual fixed-point robustness model.
//
// The shared transformed-parameter block creates many large deterministic
// individual-by-treatment arrays. They are required while evaluating the
// likelihood, but they need not be written for every retained MCMC draw.
// Keeping the block local to `model` preserves the original log density while
// limiting the CSV output to sampler diagnostics and sampled parameters.

functions {
  #include struct_section_functions.stan
}

data {
  #include struct_section_data.stan
}

transformed data {
  #include struct_section_transformed_data.stan
}

parameters {
  #include struct_section_parameters.stan
}

model {
  {
    #include struct_section_model_locals.stan
    #include struct_section_model.stan
  }
}

generated quantities {
}
