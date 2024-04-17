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
transformed parameters {
  #include struct_section_transformed_parameters.stan
}
model {
  #include struct_section_model.stan
}
generated quantities {

}
