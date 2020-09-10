import cmdstanpy
import arviz
import ipyvolume as ipv
import numpy as np
import pkg_resources
import os

_available_models = [
    "rff_ipn.stan",
    "rff_omega_ipn.stan",
    "rff_omega_ipn_cos.stan",
    "rff_omega_ipn_earth.stan",
    "rff_bw_ipn.stan",
    "rff.stan",
    "rff_bw.stan",
    "rff_omega_ipn_dt.stan",
    "rff_omega.stan",
]


def get_stan_model(stan_model, mpi=False, threads=True):

    assert (
        stan_model in _available_models
    ), f"{stan_model} is not in {','.join(_available_models)}"

    stan_file = pkg_resources.resource_filename(
        "pyipn", os.path.join("stan_models", stan_model)
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


# def plot_stan_fit(fit, universe, cmap="Set1", color="blue"):

#     ar = arviz.from_cmdstanpy(fit)

#     raw_xyz = np.array(ar.posterior.grb_xyz).reshape(-1,
#                                                      ar.posterior.grb_xyz.shape[-1])

#     rad = universe.grb_radius

#     # scatter = rad + np.random.normal(0, rad * 0.05, size=len(raw_xyz))

#     universe.plot_all_annuli(cmap=cmap, lw=3, threeD=True)

#     xyz = rad * raw_xyz

#     ipv.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2],
#                 marker="sphere", color=color, size=0.7)


def list_stan_models():
    for m in _available_models:

        print(m)
