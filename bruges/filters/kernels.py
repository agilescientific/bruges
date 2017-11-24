# -*- coding: utf 8 -*-
"""
2D kernels for image processing.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np


def gaussian_kernel(size, size_y=None):
    size = int(size)
    if not size_y:
        size_y = size
    else:
        size_y = int(size_y)
    x, y = np.mgrid[-size:size+1, -size_y:size_y+1]
    g = np.exp(-(x**2/float(size)+y**2/float(size_y)))
    return g / g.sum()
