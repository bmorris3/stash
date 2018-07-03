import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
import astropy.units as u
from astropy.constants import R_earth, R_sun
from sunpy.map import Map
from stash import simulate_lightcurve
from glob import glob
from astropy.time import Time
import numpy as np

# Set the HMI observing cadence (should be smaller than convective timescale)
cadence = 5 * u.min

paths = glob('data/*continuum.fits')

# If you haven't already done the download...
if len(paths) < 10:

    # Download four hours worth of images from JSOC, at ``cadence``
    from sunpy.net import jsoc
    client = jsoc.JSOCClient()
    response = client.search(a.jsoc.Time('2013/5/13 13:00', '2013/5/13 17:00'),
                             a.jsoc.Series('hmi.Ic_45s'),
                             a.jsoc.Notify("brettmorris21@gmail.com"),
                             a.vso.Sample(cadence))

    requests = client.request_data(response)

    res = client.get_request(requests, path='data/.')
    res.wait(progress=True)

    paths = glob('data/*continuum.fits')

lcs = []
times = []

# Set up planet parameters
orbital_period = 2 * u.day
semimajor_axis = 0.03 * u.AU
impact_parameter = 0.2
R_planet = R_earth
R_star = R_sun

# Make a light curve for each image
for path in paths:
    map = Map(path)
    image = map.data
    time = Time(map.date).jd

    # Simulate a light curve for that system, return a `LightCurve` object
    lc = simulate_lightcurve(image, orbital_period, semimajor_axis,
                             impact_parameter, R_planet, R_star)

    lcs.append(lc)
    times.append(time)

mid_transit = np.median(times)

# Grab chunks of light curves from each image in ``cadence`` intervals
time_chunks = []
flux_chunks = []
plt.figure()
for lc, t in zip(lcs, times):
    condition = np.abs(lc.times - (mid_transit - t)) <= (cadence/2).to(u.day).value
    time_chunks.append(lc.times[condition])
    flux_chunks.append(lc.fluxes[condition])

# Concatenate the chunks, and sort them in time order
times = np.concatenate(time_chunks)
fluxes = np.concatenate(flux_chunks)
fluxes = fluxes[np.argsort(times)]
times = times[np.argsort(times)]

# Save resulting light curve to a text file
np.savetxt('data/lc.txt', np.vstack([times, fluxes]).T)

# Plot the results
plt.plot(times, fluxes)

# Show me the plot!
plt.show()

