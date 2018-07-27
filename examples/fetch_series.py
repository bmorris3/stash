import matplotlib.pyplot as plt
from sunpy.net import jsoc, attrs as a
import astropy.units as u
from astropy.constants import R_earth, R_sun
from astropy.io import fits
from stash import simulate_lightcurve, transit_duration
from glob import glob
from astropy.time import Time
import numpy as np
from astropy.utils.console import ProgressBar
from datetime import datetime as dt_obj
import drms

# Set up planet parameters
orbital_period = 365 * u.day
semimajor_axis = 1 * u.AU
impact_parameter = 0.2
R_planet = 1 * R_earth
R_star = R_sun

# determine the transit duration
duration = transit_duration(R_star, R_planet, orbital_period, semimajor_axis, impact_parameter)

# make a connection to the database
c = drms.Client()

# query the FITS header keywords and the data arrays, or segments, separately:
# setting the HMI observing cadence to 315 seconds (should be smaller than convective turnover time, which is ~10 minutes)
keys, segments = c.query('hmi.Ic_45s[2013.05.13/1d@315s]', key=drms.const.all, seg='continuum')
    
# convert the time into a datetime object
def parse_tai_string(tstr,datetime=True):
    year   = int(tstr[:4])
    month  = int(tstr[5:7])
    day    = int(tstr[8:10])
    hour   = int(tstr[11:13])
    minute = int(tstr[14:16])
    if datetime: return dt_obj(year,month,day,hour,minute)
    else: return year,month,day,hour,minute
    
t_recs = np.array([parse_tai_string(keys.T_REC[i],datetime=True) for i in range(keys.T_REC.size)])    
    
lcs = []
times = []

# number of frames needed for transit simulation:
cadence = 315 * u.second
n_frames = int(float(duration/cadence))

# Make a command line progress bar that shows estimated runtime:
with ProgressBar(n_frames) as bar:

    # Make a light curve for each image
    for i in range(len(segments)):
        url = 'http://jsoc.stanford.edu' + segments.continuum[i]
        f = fits.open(url, memmap=False, lazy_load_hdus=True)
        f[1].verify('silentfix')
        image = f[1].data
        time = Time(t_recs[i], format='datetime', scale='tai').jd

        # Simulate a light curve for that system, return a `LightCurve` object
        lc = simulate_lightcurve(image, orbital_period, semimajor_axis,
                                 impact_parameter, R_planet, R_star,
                                 supersample_factor=4)

        lcs.append(lc)
        times.append(time)
        bar.update()

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
