import astropy.units as u
from astropy.constants import R_earth, R_sun
from sunpy.map import Map
from stash import simulate_lightcurve
import numpy as np
import matplotlib.pyplot as plt

# Load image from the `fetch_and_plot_example.py` script
image = Map('data/hmi_ic_45s_2013_05_13_15_56_15_tai_continuum.fits').data

# Set up planet parameters
orbital_period = 365 * u.day
semimajor_axis = 1 * u.AU
R_planet = R_earth
R_star = R_sun

# Simulate a light curve for that system, return a `LightCurve` object
lc = simulate_lightcurve(image, orbital_period, semimajor_axis,
                         R_planet, R_star)


# Plot the resulting light curve
lc.plot()

# Show me the plot!
plt.show()


