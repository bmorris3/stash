import numpy as np
import astropy.units as u
from astropy.constants import R_sun, R_earth
from . import generate_lightcurve as lc

__all__ = ['simulate_lightcurve']

from .lightcurve import LightCurve


def simulate_lightcurve(image, period, a, R_planet_physical=R_earth,
                        background=269, R_star_physical=R_sun):
    """
    Simulate a light curve of a planet with radius ``R_planet_physical`` with
    orbital period ``period``,  semimajor axis ``a``, and assuming the Sun has
    radius ``R_star_physical``.

    Parameters
    ----------
    image : `~numpy.ndarray`
        SDO HMI continuum image of the Sun
    period : `~astropy.units.Quantity`
        Orbital period of the planet
    a : `~astropy.units.Quantity`
        Semimajor axis of the planet in absolute units
    R_planet_physical : `~astropy.units.Quantity`
        Radius of the planet
    R_star_physical : `~astropy.units.Quantity`
        Radius of the star

    Returns
    -------
    lc : `~hcts.LightCurve`
        Simulated light curve of a planet transiting the star in ``image``.
    """

    diff = np.diff(np.sum(image, axis=1))
    left, right = np.argmax(diff), np.argmin(diff)

    star_width_pixels = right - left
    star_width = 2 * R_star_physical
    distance_per_pixel = star_width / star_width_pixels

    orbital_velocity_earth = (2*np.pi*a)/period
    time_per_step = (distance_per_pixel/orbital_velocity_earth).to(u.day).value # days

    radius_planet_pixels = (R_planet_physical / distance_per_pixel).value

    fluxes = lc.generate_lightcurve(image, R_planet_pixels=radius_planet_pixels,
                                  background=background)

    times = np.arange(0, len(fluxes)) * time_per_step
    times -= times.mean()

    return LightCurve(times, fluxes)
