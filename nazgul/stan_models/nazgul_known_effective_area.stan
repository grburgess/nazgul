functions {
#include functions.stan
}

data {
#include data.stan
}

transformed data {
#include transform_data.stan
#include standard_trans_data_ops.stan
}

parameters {
#include standard_parameters.stan
}

transformed parameters {
#include standard_transformed_parameters.stan
#include standard_transformed_parameter_ops.stan

}

model {
#include standard_model_dec.stan
#include known_effective_area_model.stan
#include standard_priors.stan  
#include standard_likelihood.stan
   
}

generated quantities {
  real grb_phi = atan2(grb_xyz[2], grb_xyz[1]);
  real grb_theta = -( acos(grb_xyz[3]) - 0.5*pi());
}
