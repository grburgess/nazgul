import cmdstanpy
import arviz
import ipyvolume as ipv
import numpy as np
import pkg_resources
import os

_available_models = [
    "rff_ipn.stan",
    "nazgul_standard.stan",
    "nazgul_ea_earth.stan",
    "nazgul_effective_area.stan",
    "nazgul_earth_occ.stan",
]


def get_stan_model(stan_model, mpi=False, threads=True):

    assert (
        stan_model in _available_models
    ), f"{stan_model} is not in {','.join(_available_models)}"

    stan_file = pkg_resources.resource_filename(
        "nazgul", os.path.join("stan_models", stan_model)
    )

    cpp_options = {}

    if mpi:
        cpp_options["STAN_MPI"] = True

    if threads:

        cpp_options["STAN_THREADS"] = True

    model = cmdstanpy.CmdStanModel(
        stan_file=stan_file, cpp_options=cpp_options
    )

    return model


def list_stan_models():
    for m in _available_models:

        print(m)
