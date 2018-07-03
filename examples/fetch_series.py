import matplotlib.pyplot as plt
from sunpy.net import jsoc, attrs as a
import astropy.units as u
from astropy.constants import R_earth, R_sun
from sunpy.map import Map
from stash import simulate_lightcurve, transit_duration
from glob import glob
from astropy.time import Time
import numpy as np

# Set the HMI observing cadence (should be smaller than convective timescale)
cadence = 5 * u.min
# Set up planet parameters
orbital_period = 365 * u.day
semimajor_axis = 1 * u.AU
impact_parameter = 0.2
R_planet = R_earth
R_star = R_sun

paths = glob('data/*continuum.fits')

# If you haven't already done the download...
if len(paths) < 10:
    # Download four hours worth of images from JSOC, at ``cadence``
    client = jsoc.JSOCClient()
    response = client.search(a.jsoc.Time('2013/5/13 12:00', '2013/5/13 23:59'),
                             a.jsoc.Series('hmi.Ic_45s'),
                             a.jsoc.Notify("brettmorris21@gmail.com"),
                             a.vso.Sample(cadence))
    print(response)
    requests = client.request_data(response)
    print(requests)
    res = client.get_request(requests, path='data/.', max_conn=10)
    res.wait(progress=True)

    paths = glob('data/*continuum.fits')

lcs = []
times = []

duration = transit_duration(R_star, R_planet, orbital_period, semimajor_axis,
                            impact_parameter)

# number of frames needed for transit simulation:
n_frames = int(float(duration/cadence))

# Make a light curve for each image
for path in paths[:n_frames]:
    map = Map(path)
    image = map.data
    time = Time(map.date).jd

    # Simulate a light curve for that system, return a `LightCurve` object
    lc = simulate_lightcurve(image, orbital_period, semimajor_axis,
                             impact_parameter, R_planet, R_star,
                             supersample_factor=1)

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
plt.xlabel('Time [d]')
plt.ylabel('Flux')
# Show me the plot!
plt.show()

