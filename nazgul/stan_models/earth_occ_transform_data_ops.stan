  for (n in 1:N_detectors) {

    horizon_angle[n] = calculate_horizon_angle(sc_pos[n]);
    sc_pos_norm[n] = norm_vector(sc_pos[n]);
    
  }
