import os
import pkg_resources
from .stan_models import get_stan_model, list_stan_models


class Nazgul(object):

    def __init__(self, assume_effective_area=False, earth_occultation=False, angular_dependence=False, force_compile=False):

        self._assume_effective_area = assume_effective_area
        self._earth_occultation = earth_occultation
        self._angular_dependence = angular_dependence

        if assume_effective_area:

            if earth_occultation:

                if angular_dependence:

                    pass

            elif angular_dependence:

                pass

            else:

                self._model_name = "nazgul_known_effective_area"

        elif earth_occultation:

            if angular_dependence:

                self._model_name = "nazgul_ea_earth"

            else:

                self._model_name = "nazgul_earth_occ"

        elif angular_dependence:

            self._model_name = "nazgul_effective_area"

        else:

            self._model_name = "nazgul_standard"


        self._model = get_stan_model(f"{self._model_name}.stan")
            
        if force_compile:

            self._model.compile(force=True)


    def clean_model(self):

        model_binary = pkg_resources.resource_filename(
        "nazgul", os.path.join("stan_models", self._model_name)
    )


        os.remove(model_binary)
        
        
        
            
    @property
    def model(self):

        return self._model
