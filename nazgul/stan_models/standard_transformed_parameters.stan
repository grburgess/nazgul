
vector[N_detectors] bkg = exp(log_bkg);


vector[2] scale = raw_scale * inv_sqrt(k);

vector[2] range;
vector[2] bw;
  
row_vector[k] omega[2]; // this weird MC integration thing.
  
vector[N_detectors-1] dt;
