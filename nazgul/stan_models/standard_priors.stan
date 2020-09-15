// priors

beta1 ~ std_normal();
beta2 ~ std_normal();

raw_scale ~ normal(1,1);

range1_raw ~ lognormal(0, .2);
range2_raw ~ lognormal(log(1e-2), .2);
  
omega_var[1] ~ std_normal();
omega_var[2] ~ std_normal();


log_bkg ~ normal(log(500), log(100));  


