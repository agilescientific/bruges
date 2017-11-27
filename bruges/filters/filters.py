# -*- coding: utf 8 -*-
"""
Smoothers.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
import scipy.ndimage

from bruges.util import nearest


def snn(horizon, size=5):

    # TODO: See how it handles Nans,
    # Consider removing them, interpolate over them,
    # and put them back at end.

    def func(this, pairs):
        centre = this.flat[this.size // 2]
        select = [nearest(this.flat[p], centre) for p in pairs]
        return np.mean(select)

    pairs = [[i, size**2-1 - i] for i in range(size**2 // 2)]
    return scipy.ndimage.generic_filter(horizon,
                                        func,
                                        size=size,
                                        extra_keywords={'pairs': pairs}
                                        )
