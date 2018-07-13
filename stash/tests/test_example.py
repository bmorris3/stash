import numpy as np
import astropy.units as u
from batman import TransitParams
from astropy.constants import R_sun, R_earth

def test_simple():

    from ..lightcurve import LightCurve

    lc = LightCurve(np.arange(10), np.ones(10))

    assert len(lc.times) == 10


def test_lc():
    from ..wrapper import simulate_lightcurve

    n = 4096
    x, y = np.meshgrid(np.arange(n), np.arange(n))
    r = 950 * 2
    u1, u2 = 0.6, 0.2

    on_star = (x - n/2)**2 + (y - n/2)**2 <= r**2
    d_sq = (x - n/2)**2/r**2 + (y - n/2)**2/r**2
    mu = np.sqrt(1 - d_sq)

    limb_dark = lambda mu: (1 - u1 * (1 - mu) - u2 * (1 - mu)**2) / (1 - u1/3 - u2/6) / np.pi
    image = limb_dark(mu) / limb_dark(0)
    image[~on_star] = 0

    # Earth's orbit:
    period = 365.25 * u.day
    a = 1.0 * u.AU
    b = 0
    lc = simulate_lightcurve(image, period, a, b, sdo_hmi=False, background=0)


    params = TransitParams()
    params.rp = float(R_earth/R_sun)
    params.inc = 90
    params.t0 = 0
    params.u = [0, 0]
    params.limb_dark = 'quadratic'
    params.per = period.value
    params.ecc = 0
    params.w = 90
    params.a = float(a/R_sun)
    best_fit = lc.get_transit_model(params, 1e-6)

    np.testing.assert_allclose(lc.fluxes, best_fit.fluxes, rtol=1e-3)