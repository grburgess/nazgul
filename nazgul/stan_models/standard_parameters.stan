
vector[k] beta1; // the amplitude along the cos basis
vector[k] beta2; // the amplitude along the sin basis

row_vector[k] omega_var[2]; // this weird MC integration thing.
  
vector[N_detectors]  log_bkg;
vector[N_detectors-1] log_amplitude; // independent amplitude1 of LC 1; probably do not need right now...

positive_ordered [2] raw_scale;

real<lower=0, upper=1> range1_raw;
real<lower=0, upper=range1_raw> range2_raw;
  

unit_vector[3] grb_xyz; // GRB cartesian location
