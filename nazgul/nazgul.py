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

                self._model = get_stan_model("nazgul_known_effective_area.stan")

        elif earth_occultation:

            if angular_dependence:

                self._model = get_stan_model("nazgul_ea_earth.stan")

            else:

                self._model = get_stan_model("nazgul_earth_occ.stan")

        elif angular_dependence:

            self._model = get_stan_model("nazgul_effective_area.stan")

        else:

            self._model = get_stan_model("nazgul_standard.stan")

        if force_compile:

            self._model.compile(force=True)

    @property
    def model(self):

        return self._model
