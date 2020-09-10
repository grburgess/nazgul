functions {
#include functions.stan
}

data {


  int N_detectors; // number of detectors used
  int N_time_bins[N_detectors]; // the number of time bins in each detector
  int max_N_time_bins; // the max number of time bins

  int counts[N_detectors, max_N_time_bins]; // matrix of counts per detector

  vector[max_N_time_bins] time[N_detectors]; // matrix of time bin midpoints per detector
  vector[max_N_time_bins] exposure[N_detectors]; // matrix of exposures per detector


  vector[3] sc_pos[N_detectors]; // the 3D vector for each space craft

  int<lower=1> k; // number of FFs

  int grainsize;



}
transformed data {

  row_vector[k] omega1; // this weird MC integration thing.
  row_vector[k] omega2; // this weird MC integration thing.



  // for the non-delayed LC, let's go ahead and compute the fucking matrices

  for (i in 1:k) {

    omega1[i] = normal_rng(0, 1);
    omega2[i] = normal_rng(0, 1);
  }

}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  real log_bkg;
  real log_amplitude; // independent amplitude1 of LC 1; probably do not need right now...

  real log_scale;
  
  real<lower=0> bw;



}

transformed parameters {

  real bkg = exp(log_bkg);
  real amplitude = exp(log_amplitude);

  real scale = exp(log_scale) * inv_sqrt(k);
 

}

model {

  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();

  log_scale ~ std_normal();

  log_bkg ~ normal(log(500), 1);
  bw ~ normal(1, 1);


  log_amplitude ~ std_normal();



  //  log_duration ~ normal(1,.2);
  //tstart ~ normal(1,5);

  target += reduce_sum(partial_log_like, counts[1], grainsize,
                       time[1], exposure[1],
                       omega1, omega2, beta1, beta2,bw,
                       0., bkg, scale, amplitude);





}

generated quantities {

  vector[k] omega[2];

  omega[1, :] = omega1';
  omega[2, :] = omega2';


}
