import astropy.units as u
from astropy.constants import R_earth, R_sun
from sunpy.map import Map
from stash import simulate_lightcurve
import matplotlib.pyplot as plt
import drms
from astropy.io import fits

# make a connection to the database
c = drms.Client()

# Use drms to search for an SDO HMI continuum intensity image closest to the time below.
# We query the FITS header keywords and the data arrays, or segments, separately:
keys, segments = c.query(
    'hmi.Ic_45s[2013.05.13]', key=drms.const.all, seg='continuum')

# Download the file(s)
url = 'http://jsoc.stanford.edu' + segments.continuum[0]
image = fits.getdata(url)

# Set up planet parameters
orbital_period = 365 * u.day
semimajor_axis = 1 * u.AU
impact_parameter = 0.2
R_planet = R_earth
R_star = R_sun

# Simulate a light curve for that system, return a `LightCurve` object
lc = simulate_lightcurve(image, orbital_period, semimajor_axis,
                         impact_parameter, R_planet, R_star)

# Plot the resulting light curve
lc.plot()

# Show me the plot!
plt.show()
