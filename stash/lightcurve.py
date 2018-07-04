
import numpy as np
import batman
import matplotlib.pyplot as plt
from copy import deepcopy

__all__ = ['LightCurve', 'transit_duration']


class LightCurve(object):
    """
    Container for light curves.
    """
    def __init__(self, times, fluxes):
        """
        Parameters
        ----------
        times : `~np.ndarray`
            Times
        fluxes : `~np.ndarray`
            Fluxes
        """
        self.times = times
        self.fluxes = fluxes

    def plot(self, *args, **kwargs):
        if kwargs.get('ax') is not None:
            ax = kwargs.pop('ax')
            ax.plot(self.times, self.fluxes, *args, **kwargs)
        else:
            plt.plot(self.times, self.fluxes, *args, **kwargs)

    def get_transit_model(self, init_params, yerr):
        """
        Fit a Mandel & Agol transit model to the light curve.

        Free parameters: inclination, midtransit time, limb-darkening parameters

        Parameters
        ----------
        init_params : list of length 4
            Initial fitting parameters for inclination, midtransit time, and
            two quadratic limb-darkening parmaeters

        yerr : float
            Uncertainty on flux measurements

        Returns
        -------
        lc : `~stash.LightCurve`
            Best-fit transit model
        """
        def transit_model(p):
            rp, inc, t0, u1, u2 = p

            trial_params = deepcopy(init_params)
            trial_params.rp = rp
            trial_params.inc = inc
            trial_params.t0 = t0
            trial_params.u = [u1, u2]

            m = batman.TransitModel(trial_params, self.times,
                                    exp_time=self.times[1]-self.times[0],
                                    supersample_factor=3)
            model = m.light_curve(trial_params)

            return model

        y = self.fluxes/self.fluxes.max()

        def chi2(p):
            return np.sum((transit_model(p) - y)**2 / yerr**2)

        from scipy.optimize import fmin_l_bfgs_b

        result = fmin_l_bfgs_b(chi2, [init_params.rp, init_params.inc,
                                      init_params.t0, init_params.u[0],
                                      init_params.u[1]],
                               approx_grad=True,
                               bounds=[[0, 1], [0, 90], [-0.5, 0.5], [-1, 1],
                                       [-1, 1]])[0]
        return LightCurve(self.times, transit_model(result))


def transit_duration(R_star, R_planet, orbital_period, semimajor_axis,
                     impact_parameter):
    """
    Compute transit duration from first through fourth contact given orbital
    and physical parameters.

    Parameters
    ----------
    R_star : `~astropy.units.Quantity`
        Stellar radius
    R_planet : `~astropy.units.Quantity`
        Planet radius
    orbital_period : `~astropy.units.Quantity`
        Orbital period
    semimajor_axis : `~astropy.units.Quantity`
        Orbital semimajor axis
    impact_parameter : float
        Impact parameter (-1, 1)

    Returns
    -------
    duration : `~astropy.units.Quantity`
        Transit duration from first through fourth contact.
    """
    inclination = np.arccos(impact_parameter / float(semimajor_axis/R_star))
    # Winn 2011 (Eqn 14)
    return orbital_period / np.pi * np.arcsin(float(R_star/semimajor_axis) *
                                              np.sqrt((1+float(R_planet/R_star))**2 -
                                                      impact_parameter**2) /
                                              np.sin(inclination))