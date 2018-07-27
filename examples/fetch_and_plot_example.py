import matplotlib.pyplot as plt
from sunpy.net import attrs as a
from sunpy.map import Map
import drms
from matplotlib import mlab, cm
import matplotlib.colors as colors
from astropy.io import fits

# make a connection to the database
c = drms.Client()

# Use drms to search for an SDO HMI continuum intensity image closest to the time below.
# We query the FITS header keywords and the data arrays, or segments, separately:
keys, segments = c.query(
    'hmi.Ic_45s[2013.05.13]', key=drms.const.all, seg='continuum')

# Download the file(s) found to the `data` directory
url = 'http://jsoc.stanford.edu' + segments.continuum[0]
photosphere_full_image = fits.getdata(url)

# Create a header (with WCS keywords)
header_hmi = dict(keys.iloc[0])
# Add a DATE-OBS keyword which seems to be required by sunpy.map.Map.
header_hmi['DATE-OBS'] = keys.DATE__OBS[0]
# Add HGLN_OBS keyword to avoid a warning in sunpy.map.Map.
header_hmi['HGLN_OBS'] = 0.0

# Set the color table
hmimag = plt.get_cmap('hinodesotintensity')

# Make the figure
m = Map(photosphere_full_image, header_hmi)

# and use sunpy's built-in features to plot it
m.peek()
