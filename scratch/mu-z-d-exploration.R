library(tidyverse)

fit_version = 66

struct_param_draws = read_csv(
    file.path(
        script_options$input_path, 
        str_interp(
            "param_posterior_draws_dist_fit${fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.csv"
        )
    )
)