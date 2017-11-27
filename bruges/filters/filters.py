# -*- coding: utf 8 -*-
"""
Smoothers.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
import scipy.ndimage

from bruges.util import nearest


def snn(arr, size=5, include=True):
    """
    Symmetric nearest neighbour, a nonlinear smoothing filter.
    http://subsurfwiki.org/wiki/Symmetric_nearest_neighbour_filter

    Args:
        arr (ndarray): a 2D array, such as a seismic horizon.
        size (int): the kernel size, e.g. 5 for 5x5. Should be odd,
            rounded up if not.
        include (bool): whether to include the central pixel itself.

    Returns:
        ndarray: the resulting smoothed array.

    TODO:
        See how it handles Nans, consider removing, interpolating, replacing.
    """
    def func(this, pairs):
        centre = this[this.size // 2]
        select = [nearest(this[p], centre) for p in pairs]
        if include:
            select += [centre]
        return np.mean(select)

    if not size // 2:
        size += 1

    pairs = [[i, size**2-1 - i] for i in range(size**2 // 2)]
    return scipy.ndimage.generic_filter(arr,
                                        func,
                                        size=size,
                                        extra_keywords={'pairs': pairs}
                                        )


def kuwahara(arr, size=5):
    """
    Kuwahara, a nonlinear smoothing filter.
    http://subsurfwiki.org/wiki/Kuwahara_filter

    Args:
        arr (ndarray): a 2D array, such as a seismic horizon.
        size (int): the kernel size, e.g. 5 for 5x5. Should be odd,
            rounded up if not.

    Returns:
        ndarray: the resulting smoothed array.

    TODO:
        See how it handles Nans, consider removing, interpolating, replacing.
    """

    def func(this):
        k = np.ceil(this.shape[0] / 2)
        t = this.reshape((k, k))
        sub = np.array([t[:k, :k].flatten(),
                        t[:k, k-1:].flatten(),
                        t[k-1:, :k].flatten(),
                        t[k-1:, k-1:].flatten()]
                       )
        select = sub[np.argmin(np.var(sub, axis=1))]
        return np.mean(select)

    if not size // 2:
        size += 1

    return scipy.ndimage.generic_filter(arr,
                                        func,
                                        size=size,
                                        )
