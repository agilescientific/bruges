# -*- coding: utf 8 -*-
"""
Smoothers.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np


def snn(horizon, kernel_size=5):
    """
    Slow but sure Symmetric Nearest Neighbours.

    http://subsurfwiki.org/wiki/Symmetric_nearest_neighbour_filter
    """
    # make an empty array the same shape as the input
    output = np.empty_like(horizon)
    inlines, xlines = horizon.shape
    offset = np.floor(kernel_size/2.)

    for x in range(inlines - 2*offset):
        x += offset  # Correct for offset
        for y in range(xlines - 2*offset):
            y += offset
            centre = horizon[(x), (y)]
            nearest_neighbours = [centre]
            for a in range(kernel_size**2. / 2.):
                i, j = divmod(a, kernel_size)
                i -= offset  # transform to relative coordinates in kernel
                j -= offset
                value = horizon[x+i, y+j]
                opposite = horizon[x-i, y-j]
                closest = min(value-centre, opposite-centre)
                nearest_neighbours.append(closest+centre)
            output[x, y] = np.mean(nearest_neighbours)

    return output
