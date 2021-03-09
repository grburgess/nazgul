
real partial_log_like_bw_multi_scale(int[] counts_slice, int start, int end, vector time, vector exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real dt, real bkg, real scale1, real scale2, real amplitude, int k) {

  int N = end - start +1;

  vector[N] time_slice = time[start:end] - dt;
  vector[N] expected_counts_log;
  vector[N] expected_counts;

  matrix [N,k] tw1 = time_slice * omega1;
  matrix [N,k] tw2 = time_slice * omega2;

  
  expected_counts_log = ((scale1 * cos(tw1) + scale2 * cos(tw2)  ) * beta1) + ((scale1 * sin(tw1) + scale2 * sin(tw2)  ) * beta2);

  expected_counts = exposure[start:end] .* (exp(expected_counts_log) * amplitude + bkg);


  return poisson_lpmf(counts_slice | expected_counts);

}

real partial_log_like_bw_multi_scale_fast(int[] counts_slice, int start, int end, vector time, vector exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real dt, real bkg, real scale1, real scale2, real amplitude, int k) {

  int N = end - start +1;

  vector[N] time_slice = time[start:end] - dt;
  
  matrix [N,k] tw1 = time_slice * omega1;
  matrix [N,k] tw2 = time_slice * omega2;


  return poisson_lpmf(counts_slice | exposure[start:end] .* (exp(((scale1 * cos(tw1) + scale2 * cos(tw2)  ) * beta1) + ((scale1 * sin(tw1) + scale2 * sin(tw2)  ) * beta2)) * amplitude + bkg));

}

real partial_total_like(int[] detector_slice, int start, int end, int[,] counts,
			vector[] time, vector[] exposure,
			row_vector omega1, row_vector omega2, vector beta1, vector beta2,
			vector dt, vector bkg, real scale1, real scale2, vector amplitude,
			int k, int[] gs_array, int[] N_time_bins) {

  real lp = 0.;
  int num_slice_terms = end - start +1;
  
  for (m in 1:num_slice_terms) {

    int n = detector_slice[m];
    
    lp += reduce_sum(partial_log_like_bw_multi_scale_fast, counts[n,:N_time_bins[n]], gs_array[n],
                         time[n,:N_time_bins[n]], exposure[n,:N_time_bins[n]],
                         omega1, omega2, beta1, beta2,
                         dt[n], bkg[n], scale1, scale2, amplitude[n], k);



  }

  return lp;

}



real c() {
  // the speed of light 
  return 299792.46; // km/s

}


  



real angular_separation(vector grb_xyz, vector sc_pointing_norm) {

  real tmp;

  tmp = dot_product(grb_xyz, sc_pointing_norm) ;

  if (tmp >0) {
  
    return tmp;  
    }

  else {

    return 0;
  }

}



real calculate_horizon_angle(vector sc_position) {

  real earth_radius = 6371.0;
  real altitude = sqrt(sum((sc_position .* sc_position)));
  real horizon_angle = 0.5 * pi() - acos(earth_radius/altitude);

  return horizon_angle;


  

}


vector norm_vector(vector sc_position) {

  real norm = sqrt(sum((sc_position .* sc_position)));

  return sc_position/norm;

}


real earth_occulation(real horizon_angle, vector sc_position, vector grb_xyz) {

  real tmp;

  real angle = acos(dot_product(grb_xyz, -sc_position));

  if (angle < horizon_angle) {

    return 0.;
    }

  else {

    return 1.;

  }

  

}




/* real time_delay( vector grb_xyz, vector sc_pos1, vector sc_pos2) { */
/*   // compute the time delay between the signals recieved from the GRB */

  
/*   real t1 = dot_product(grb_xyz, sc_pos1); */
/*   real t2 = dot_product(grb_xyz, sc_pos2); */


/*   return (t1 - t2) * inv(c()); */


/* } */


/* real time_delay( vector grb_xyz, vector sc_pos1, vector sc_pos2) { */
/*   // compute the time delay between the signals recieved from the GRB */

/*   return dot_product(grb_xyz, sc_pos1 - sc_pos2) * inv(c()); */

/* } */

real time_delay( vector grb_xyz, vector sc_pos_diff) {
  // compute the time delay between the signals recieved from the GRB

  return dot_product(grb_xyz, sc_pos_diff) * inv(c());

}
