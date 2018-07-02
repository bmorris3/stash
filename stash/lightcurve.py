
import numpy as np
import batman
import matplotlib.pyplot as plt
from copy import deepcopy

__all__ = ['LightCurve']


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

    def plot(self):
        plt.plot(self.times, self.fluxes)

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
        lc : `~hcts.LightCurve`
            Best-fit transit model
        """
        def transit_model(p):
            inc, t0, u1, u2 = p

            trial_params = deepcopy(init_params)
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

        result = fmin_l_bfgs_b(chi2, [init_params.inc, init_params.t0,
                                      init_params.u[0], init_params.u[1]],
                               approx_grad=True,
                               bounds=[[0, 90], [0, 0.5], [-1, 1], [-1, 1]])[0]
        return LightCurve(self.times, transit_model(result))