import numpy as np


def test_simple():

    from ..lightcurve import LightCurve

    lc = LightCurve(np.arange(10), np.ones(10))

    assert len(lc.times) == 10