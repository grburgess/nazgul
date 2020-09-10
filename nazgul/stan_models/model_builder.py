from .stan_builder import DataBlock, ParametersBlock, TransDataBlock, ModelBlock, TransParametersBlock, StanGenerator, FunctionsBlock


class Nazgul(object):

    def __init__(self, assume_effective_area=True, earth_occultation=False, angular_dependence=False):

        self._assume_effective_area = assume_effective_area
        self._earth_occultation = earth_occultation
        self._angular_dependence = angular_dependence


    def _build_code(self):


        self._functions_block = FunctionsBlock()
        self._functions_block.add_include("functions.stan")


        self._data_block = DataBlock(data_size"N_detectors")
        self._data_block.add_vector_data(name="N_time_bins",stan_type="int")
        self._data_block.add_data(name="max_N_time_bins" stan_type="int")
        self._data_block.add_vector_data(name="counts",stan_type="int", size="max_N_time_bins", n_vectors="N_detectors")
        
        


_functions_block = "functions {\n"
_functions_block += "#include functions.stan\n"
_functions_block += "}\n"

_data_block ="\n"
_data_block += "data {\n"
_data_block += "\tint N_detectors; // number of detectors used\n"
_data_block += "\tint N_time_bins[N_detectors]; // the number of time bins in each detector\n"
_data_block += "\tint max_N_time_bins; // the max number of time bins\n"
_data_block += "\tint counts[N_detectors, max_N_time_bins]; // matrix of counts per detector\n"
_data_block += "\tvector[max_N_time_bins] time[N_detectors]; // matrix of time bin midpoints per detector\n"
_data_block += "\tvector[max_N_time_bins] exposure[N_detectors]; // matrix of exposures per detector\n"
_data_block += "\tvector[3] sc_pos[N_detectors]; // the 3D vector for each space craft\n"
  data_block += "vector[3] sc_pointing_norm[N_detectors]; // the 3D vector for each space craft\n"
_data_block += "\tint<lower=1> k; // number of FFs\n"
_data_block += "\tint grainsize[N_detectors];\n"
_data_block += "}\n"



_transformed_data_block = "transformed data {\n"

_transformed_data_block += "\tvector[max_N_time_bins] log_exposure[N_detectors] = log(exposure);\n"
  
_transformed_data_block += "\treal inv_sqrt_k = inv_sqrt(k);\n"
_transformed_data_block += "\tvector[N_detectors] horizon_angle;\n"
_transformed_data_block += "\tvector[3] sc_pos_norm[N_detectors];\n"
_transformed_data_block += "\treal max_range;\n"
_transformed_data_block += "\tvector[N_detectors] maxs;\n"
_transformed_data_block += "\tvector[N_detectors] mins;\n\n"

_transformed_data_block += "\tfor (n in 1:N_detectors){\n"
_transformed_data_block += "\t\tmaxs[n] = max(time[n]);\n"
_transformed_data_block += "\t\tmins[n] = min(time[n]);\n"
_transformed_data_block += "\t}\n\n"
_transformed_data_block += "\tmax_range = max(maxs) - min(mins);\n"
_transformed_data_block += "\tfor (n in 1:N_detectors) {\n"
_transformed_data_block += "\t\thorizon_angle[n] = calculate_horizon_angle(sc_pos[n]);\n"
_transformed_data_block += "\t\tsc_pos_norm[n] = norm_vector(sc_pos[n]);\n"
_transformed_data_block += "\t}\n"
_transformed_data_block += "}\n"


_parameters_block = "parameters {\n"

_parameters_block += "\tvector[k] beta1; // the amplitude along the cos basis\n"
_parameters_block += "\tvector[k] beta2; // the amplitude along the sin basis\n"

_parameters_block += "\trow_vector[k] omega_var[2]; // this weird MC integration thing.\n"
  



_parameters_block += "\tvector[N_detectors]  log_bkg;\n"

_parameters_block += "\tpositive_ordered [2] raw_scale;\n"
_parameters_block += "\t//real<lower=0> log_scale;\n"
_parameters_block += "\t//positive_ordered[2] bw;\n"
_parameters_block += "\t//ordered[2] log_bw;\n"

_parameters_block += "\treal<lower=0, upper=1> range1_raw;\n"
_parameters_block += "\treal<lower=0, upper=range1_raw> range2_raw;\n"
  

_parameters_block += "\tunit_vector[3] grb_xyz; // GRB cartesian location\n"
_parameters_block += "}\n"


_parameters_block += "\tvector[N_detectors-1] log_amplitude; // independent amplitude1 of LC 1; probably do not need right now...\n"




