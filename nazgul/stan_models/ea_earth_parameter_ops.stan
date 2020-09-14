for (n in 1:N_detectors) {

  ang_sep[n]= angular_separation(grb_xyz, sc_pointing_norm[n]);
  earth_occulted[n]= earth_occulation(horizon_angle[n], sc_pos_norm[n], grb_xyz);
  
 }
