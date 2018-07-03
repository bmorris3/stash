
****************
Installing stash
****************

The easy way
------------

You can install stash via pip like so::

    pip install git+git://github.com/bmorris3/stash.git


Manually building stash
-----------------------

Stash requires some python packages to be installed first, which you can get
via pip like so::

    pip install numpy batman-package matplotlib sunpy astropy


Once you have those installed, you can install the latest version of ``stash``
by building from the source::

    git clone git@github.com:bmorris3/stash.git
    cd stash
    python setup.py install

