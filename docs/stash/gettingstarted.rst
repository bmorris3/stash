***************
Getting Started
***************


Fetching an image with Fido
---------------------------

For this demo, we'll be using sunpy and matplotlib. Let's do some imports::

    import matplotlib.pyplot as plt
    from sunpy.net import Fido, attrs as a
    from sunpy.map import Map

Great! Now we'll use Fido to search between two times for an image taken by
SDO/HMI in continuum intensity mode::

    # Use fido to search for an SDO HMI continuum intensity image
    # between the two times below
    result = Fido.search(a.Time('2013/5/13 15:55', '2013/5/13 15:55:30'),
                         a.Instrument('HMI'), a.vso.Physobs('intensity'))

The variable ``result`` stores the information about the images that match the
search criteria. if we feed it to ``Fido.fetch``, ``sunpy`` will download the
image::

    # Download the file(s) found to the `data` directory
    downloaded_files = Fido.fetch(result, path='data/.')

We can open the downloaded image using `~sunpy.map.Map`, like so::

    # Use Sunpy to open up the image
    m = Map(downloaded_files[0])

and finally we can see what the map looks like with::

    # and use sunpy's built-in features to plot it
    m.peek()
    plt.show()

which gives us the following awesome plot!

.. plot::

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


Simulating a transit
--------------------

Now that we have an image of the Sun to work with, let's simulate a transit of
an Earth-like planet transiting the Sun.

As with before, we'll have a bunch of packages to import things from::

    import astropy.units as u
    from astropy.constants import R_earth, R_sun
    from sunpy.map import Map
    from stash import simulate_lightcurve
    import matplotlib.pyplot as plt

The first thing we'll do is open up the image that we downloaded in the previous
example::

    # Load image from the `fetch_and_plot_example.py` script
    image = Map('data/hmi_ic_45s_2013_05_13_15_56_15_tai_continuum.fits').data

That's the image of the Sun that we'll simulate a transit on top of. Now let's
define some characteristics of the planet that will be doing the transiting::

    # Set up planet parameters
    orbital_period = 365 * u.day
    semimajor_axis = 1 * u.AU
    impact_parameter = 0.2
    R_planet = R_earth
    R_star = R_sun

The impact parameter ``b`` is defined on [-1, 1], and describes how far from the
center of the solar disk, in units of solar radii, the planet appears to cross
over. A planet that transits the solar equator has ``b=0``. A planet that just
barely grazes the tip of the Sun has ``b=1``. (For well-aligned planets, you
can convert between the solar latitude being occulted and the impact parameter
by noting that :math:`b = \sin \theta`, where :math:`\theta` is the latitude.)

Now we call the `~stash.simulate_lightcurve` function to simulate a light curve
of a transit of this planet on our image, ``image``::

    # Simulate a light curve for that system, return a `LightCurve` object
    lc = simulate_lightcurve(image, orbital_period, semimajor_axis,
                             impact_parameter, R_planet, R_star)

We can plot the light curve using the built-in plot function::

    # Plot the resulting light curve
    lc.plot()

    # Show me the plot!
    plt.show()

and we'll see something like this:

.. plot::

    import astropy.units as u
    from astropy.constants import R_earth, R_sun
    from sunpy.map import Map
    from stash import simulate_lightcurve
    import matplotlib.pyplot as plt

    # Load image from the `fetch_and_plot_example.py` script
    image = Map('data/hmi_ic_45s_2013_05_13_15_56_15_tai_continuum.fits').data

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


Would you look at that –– the planet occulted a starspot, causing the apparent
brightness to temporarily increase during the transit, because less flux was
missing when the planet was over the spot, compared to when the planet is over
the typically bright photosphere. Now we're cooking!

Simulating a bunch of transits
------------------------------

Now this time, let's iterate over impact parameter and see all of the different
transit light curves we could get as we vary :math:`b \in [-1, 0]`:

.. plot::

    import astropy.units as u
    from astropy.constants import R_earth, R_sun
    from sunpy.map import Map
    from stash import simulate_lightcurve
    import matplotlib.pyplot as plt
    import numpy as np

    # Load image from the `fetch_and_plot_example.py` script
    image = Map('data/hmi_ic_45s_2013_05_13_15_56_15_tai_continuum.fits').data

    # Set up planet parameters
    orbital_period = 365 * u.day
    semimajor_axis = 1 * u.AU
    R_planet = R_earth
    R_star = R_sun

    # Iterate over a range of impact parameters:
    for impact_parameter in np.arange(-0.8, 0, 0.05):
        # Simulate a light curve for that system, return a `LightCurve` object
        lc = simulate_lightcurve(image, orbital_period, semimajor_axis,
                                 impact_parameter, R_planet, R_star)

        # Plot the resulting light curve
        lc.plot()

    plt.xlabel('Time [d]')
    plt.ylabel('Flux')

    # Show me the plot!
    plt.show()

You can see that the planet occulted the big starspot at one of the impact
parameters that we swept through in the ``for`` loop.


