import arviz as av
import numpy as np
import numba as nb
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import astropy.units as u

import ligo.skymap.plot
from ligo.skymap.plot import reticle
from ligo.skymap.kde import Clustered2DSkyKDE
import astropy_healpix as ah
from astropy.wcs.utils import skycoord_to_pixel
from mocpy import MOC

from pyipn.rff import RFF, RFF_multiscale
from .utils.timing import compute_annulus_from_time_delay

from .io.plotting.projection import create_skw_dict
from .universe import Universe


class Fit(object):
    def __init__(self, *inference_data, universe_save=None, npix=2 ** 5, fast_open=False):
        """FIXME! briefly describe function

        :param inference_data:
        :param universe_save:
        :param npix:
        :returns: we love
        :rtype:

        """

        if len(inference_data) == 1:

            self._posterior = inference_data[0]

        else:

            self._posterior = inference_data[0]

            for idata in inference_data[1:]:

                self._posterior = av.concat(
                    self._posterior, idata, dim="chain")

        self._npix = npix

        try:
            ang_sep = self._posterior.posterior.ang_sep.stack(
                sample=("chain", "draw")
            ).values

            self._do_contour = True

        except:

            self._do_contour = False


    

        
        self._beta1 = self._posterior.posterior.beta1.stack(
            sample=("chain", "draw")
        ).values

        # set the number of samples.. flattened over chain

        self._n_samples = self._beta1.shape[-1]

        self._beta2 = self._posterior.posterior.beta2.stack(
            sample=("chain", "draw")
        ).values

        self._omega1 = self._posterior.posterior.omega.stack(
            sample=("chain", "draw")
        ).values[0]

        self._omega2 = self._posterior.posterior.omega.stack(
            sample=("chain", "draw")
        ).values[1]

        try:

            self._amplitude = self._posterior.posterior.amplitude.stack(
                sample=("chain", "draw")
            ).values

        except:

            self._amplitude = np.ones(self._n_samples)

        self._background = self._posterior.posterior.bkg.stack(
            sample=("chain", "draw")
        ).values

        self._scale = self._posterior.posterior.scale.stack(
            sample=("chain", "draw")
        ).values

        if self._scale.shape[0] == 2:

            self._multi_scale = True

        else:

            self._multi_scale = False

        try:

            self._dt = self._posterior.posterior.dt.stack(
                sample=("chain", "draw")
            ).values

            self._grb_theta = self._posterior.posterior.grb_theta.stack(
                sample=("chain", "draw")
            ).values
            self._grb_phi = self._posterior.posterior.grb_phi.stack(
                sample=("chain", "draw")
            ).values

            self._is_dt_fit = True

            self._n_dets = self._background.shape[0]

        except:

            self._is_dt_fit = False

            self._dt = None
            self._grb_theta = None
            self._grb_phi = None

            self._n_dets = 1

        self._use_bw = False

        try:

            self._bw1 = self._posterior.posterior.bw1.stack(
                sample=("chain", "draw")
            ).values

            self._bw2 = self._posterior.posterior.bw2.stack(
                sample=("chain", "draw")
            ).values

            self._multi_bw = True

        except:

            try:

                self._bw = self._posterior.posterior.bw.stack(
                    sample=("chain", "draw")
                ).values

                if self._bw.shape[0] == 2:

                    self._multi_bw = True

                else:

                    self._multi_bw = False

            except:

                self._bw = self._posterior.posterior.bw_out.stack(
                    sample=("chain", "draw")
                ).values

                self._use_bw = False
                self._multi_bw = True

        self.grb_color = "k"
        self._grb_style = "lrtb"

        self._has_universe = False

        if universe_save is not None:

            self._universe = Universe.from_save_file(universe_save)

            self._has_universe = True

        if self._is_dt_fit and self._n_dets > 2 and (not fast_open):

            self._build_moc_map()

        elif self._do_contour:

            self._build_moc_map()
            
    def _build_moc_map(self):

        pts = np.column_stack((self._grb_phi, self._grb_theta))

        self._kde_map = Clustered2DSkyKDE(pts, jobs=12)

        data = self._kde_map.as_healpix(top_nside=self._npix)

        self._uniq = data["UNIQ"]
        self._prob_density = data["PROBDENSITY"]

        level, ipix = ah.uniq_to_level_ipix(self._uniq)
        area = ah.nside_to_pixel_area(
            ah.level_to_nside(level)).to_value(u.steradian)
        self._prob = self._prob_density * area

    def _detector_check(self, det_number):

        assert det_number in range(self._n_dets)

    @classmethod
    def from_cmdstanpy(cls, fit):

        inference_data = av.from_cmdstanpy(fit)

        return cls(inference_data)

    @classmethod
    def from_netcdf(cls, *file_name, fast_open=False):

        if len(file_name) == 1:

            inference_data = [av.from_netcdf(file_name[0])]

        else:

            inference_data = [av.from_netcdf(f) for f in file_name]

        return cls(*inference_data, fast_open=fast_open)

    def set_universe(self, universe):

        self._universe = universe
        self._has_universe = True

    def expected_rate(self, time, detector):
        """FIXME! briefly describe function

        :param time:
        :param detector:
        :returns:
        :rtype:

        """

        self._detector_check(detector)

        if self._is_dt_fit and (detector > 0):

            dt = self._dt[detector - 1]

        else:

            dt = np.zeros(self._n_samples)

        if not self._is_dt_fit:

            amp = self._amplitude

        elif detector == 0:

            amp = np.ones(self._n_samples)

        else:

            amp = self._amplitude[detector - 1]

        if self._use_bw:

            bw = self._bw

        else:

            bw = np.ones(self._n_samples)

        if not self._multi_scale:
            out = _expeced_rate(
                time,
                self._omega1,
                self._omega2,
                self._beta1,
                self._beta2,
                bw,
                self._scale,
                amp,
                dt,
                self._n_samples,
            )
        else:

            out = _expeced_rate_multiscale(
                time,
                self._omega1,
                self._omega2,
                self._beta1,
                self._beta2,
                bw,
                self._scale,
                amp,
                dt,
                self._n_samples,
            )

        return out

    def _contour_two_detectors(
        self, levels=[0.68], colors=["green"], ax=None, **kwargs
    ):

        dt = self._dt[0]

        assert len(levels) == len(colors)

        dkey = list(self._universe.detectors.keys())

        d1 = self._universe.detectors[dkey[0]]
        d2 = self._universe.detectors[dkey[1]]

        for i, level in enumerate(levels):

            dt1, dt2 = av.hdi(dt, hdi_prob=level)

            compute_annulus_from_time_delay(
                dt1 * u.s, dt2 * u.s, d1, d2, color=colors[i], ax=ax, **kwargs
            )

    def _contour_more_detectors(
        self, levels=[0.68], colors=["green"], ax=None, **kwargs
    ):

        assert len(levels) == len(colors)

        mocs = [
            MOC.from_valued_healpix_cells(self._uniq, self._prob, cumul_to=c)
            for c in levels
        ]

        for moc, col in zip(mocs, colors[::-1]):
            skycoords = moc.get_boundaries()

            for sc in skycoords:

                x, y = skycoord_to_pixel(sc, ax.wcs)
                p = Path(np.vstack((x, y)).T)
                patch = PathPatch(p, color=col, **kwargs)
                ax.add_patch(patch)

    def _show_grb(self, ax):

        ax.plot(
            self._universe.grb.location.coord.ra.deg,
            self._universe.grb.location.coord.dec.deg,
            transform=ax.get_transform("icrs"),
            marker=reticle(inner=0.3, which=self._grb_style),
            markersize=20,
            markeredgewidth=1,
            color=self.grb_color,
        )

    def location_contour(
        self,
        levels=[0.68],
        colors=["green"],
        ax=None,
        projection="astro degrees mollweide",
        center=None,
        radius=None,
        show_grb=True,
        **kwargs
    ):
        """FIXME! briefly describe function

        :param levels:
        :param colors:
        :param ax:
        :param projection:
        :param center:
        :param radius:
        :param show_grb:
        :returns:
        :rtype:

        """

        if ax is None:

            skw_dict = create_skw_dict(projection, center, radius)

            fig, ax = plt.subplots(subplot_kw=skw_dict)

        else:

            fig = ax.get_figure()

        if self._n_dets == 2 and (not self._do_contour):

            self._contour_two_detectors(levels, colors, ax=ax, **kwargs)

        else:

            self._contour_more_detectors(levels, colors, ax=ax, **kwargs)

        if show_grb and self._has_universe:

            self._show_grb(ax)

        return fig

    def location_scatter(
        self,
        color="green",
        projection="astro degrees mollweide",
        center=None,
        radius=None,
        show_grb=True,
        ax=None,
        **kwargs
    ):

        if ax is None:

            skw_dict = create_skw_dict(projection, center, radius)

            fig, ax = plt.subplots(subplot_kw=skw_dict)

        else:

            fig = ax.get_figure()

        theta = np.rad2deg(self._grb_theta)
        phi = np.rad2deg(self._grb_phi)

        idx = phi <= 0

        phi[idx] += 360

        ax.scatter(
            phi, theta, transform=ax.get_transform("icrs"), color=color, **kwargs
        )

        if self._has_universe and show_grb:
            self._show_grb(ax)

        return fig

    def plot_light_curve_fit(self, detector, tstart, tstop, dt=0.2, thin=1, color='r', **kwargs):

        self._detector_check(detector)
        assert self._has_universe

        lc = self._universe.light_curves[
            list(self._universe.light_curves.keys())[detector]
        ]

        rate, edges, counts = lc.get_binned_light_curve(tstart, tstop, dt)

        exposure = np.diff(edges)

        mid_points = 0.5 * (edges[:-1] + edges[1:])

        fig, ax = plt.subplots()

        pred_rate = self.expected_rate(mid_points, detector)

        if self._is_dt_fit:

            bkg = self._background[detector]

        else:

            bkg = self._background

        for i in range(self._n_samples)[::thin]:

            ax.plot(mid_points, pred_rate[i] + bkg[i], color=color, alpha=0.05)

        ax.scatter(mid_points, rate, **kwargs)

        ax.set(xlabel="time (s)", ylabel="rate (cnts/s)")

        return fig

    def plot_light_curve_ppcs(self, detector, tstart, tstop, dt=0.2, levels=[99, 95, 68], colors=["r", "g", "b"], **kwargs):

        self._detector_check(detector)
        assert self._has_universe

        lc = self._universe.light_curves[
            list(self._universe.light_curves.keys())[detector]
        ]

        rate, edges, counts = lc.get_binned_light_curve(tstart, tstop, dt)

        exposure = np.diff(edges)

        mid_points = 0.5 * (edges[:-1] + edges[1:])

        fig, ax = plt.subplots()

        if self._is_dt_fit:

            bkg = self._background[detector]

        else:

            bkg = self._background

        # compute the PPC bounds

        ppcs = self._compute_ppcs(detector, tstart, tstop, dt)

        ppc_low = []
        ppc_high = []

        for level in levels:

            tmp_low = np.percentile(ppcs / exposure, 50. - level / 2., axis=0)
            tmp_high = np.percentile(ppcs / exposure, 50. + level / 2., axis=0)

            ppc_low.append(tmp_low)
            ppc_high.append(tmp_high)

        #colors = [light,mid,dark]

        for j, (lo, hi) in enumerate(zip(ppc_low, ppc_high)):

            for i in range(len(edges) - 1):
                ax.fill_between([edges[i], edges[i + 1]],
                                lo[i], hi[i], fc=colors[j], ec="none", lw=0)

        ax.scatter(mid_points, rate, **kwargs)

        ax.set(xlabel="time (s)", ylabel=r"rate (cnt s$^{-1}$)")

        return fig

    def _compute_ppcs(self, detector, tstart, tstop, dt):

        lc = self._universe.light_curves[
            list(self._universe.light_curves.keys())[detector]
        ]

        rate, edges, counts = lc.get_binned_light_curve(tstart, tstop, dt)

        mid_points = 0.5 * (edges[:-1] + edges[1:])

        exposure = np.diff(edges)

        pred_rate = self.expected_rate(mid_points, detector)

        if self._is_dt_fit:

            bkg = self._background[detector]

        else:

            bkg = self._background

        ppcs = ppc_generator(mid_points, exposure,
                             pred_rate, bkg, self._n_samples)

        return ppcs

    @property
    def beta1(self):

        return self._beta1

    @property
    def beta2(self):

        return self._beta2

    @property
    def omega1(self):

        return self._omega1

    @property
    def omega2(self):

        return self._omega2

    @property
    def bw(self):

        return self._bw

    @property
    def bw1(self):

        return self._bw1

    @property
    def bw2(self):

        return self._bw2

    @property
    def scale(self):
        return self._scale

    @property
    def dt(self):

        return self._dt

    @property
    def grb_theta(self):

        return self._grb_theta

    @property
    def grb_phi(self):

        return self._grb_phi

    @property
    def amplitude(self):

        return self._amplitude

    @property
    def background(self):

        return self._background

    @property
    def posterior(self):

        return self._posterior


@nb.njit(fastmath=True, cache=True)
def _expeced_rate(time, omega1, omega2, beta1, beta2, bw, scale, amplitude, dt, N):

    out = np.empty((N, len(time)))

    for n in range(N):

        out[n] = amplitude[n] * RFF(
            time - dt[n],
            omega1[:, n],
            omega2[:, n],
            beta1[:, n],
            beta2[:, n],
            bw[n],
            scale[n],
        )

    return out


@nb.njit(fastmath=True, cache=True)
def _expeced_rate_multiscale(
    time, omega1, omega2, beta1, beta2, bw, scale, amplitude, dt, N
):

    out = np.empty((N, len(time)))

    for n in range(N):

        out[n] = amplitude[n] * RFF_multiscale(
            time - dt[n],
            omega1[:, n],
            omega2[:, n],
            beta1[:, n],
            beta2[:, n],
            bw[n],
            scale[0, n],
            scale[1, n],
        )

    return out


@nb.njit(fastmath=True, cache=True)
def ppc_generator(time, exposure, rates, bkg, N):

    out = np.empty((N, len(time)))

    for n in range(N):

        for i in range(len(time)):

            out[n, i] = np.random.poisson((rates[n, i] + bkg[n]) * exposure[i])

    return out
