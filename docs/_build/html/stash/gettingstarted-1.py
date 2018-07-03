import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from sunpy.map import Map

# Use fido to search for an SDO HMI continuum intensity image
# between the two times below
result = Fido.search(a.Time('2013/5/13 15:55', '2013/5/13 15:55:30'),
                     a.Instrument('HMI'), a.vso.Physobs('intensity'))

# Download the file(s) found to the `data` directory
downloaded_files = Fido.fetch(result, path='data/.')

# Use Sunpy to open up the image
m = Map(downloaded_files[0])

# and use sunpy's built-in features to plot it
m.peek()

plt.show()