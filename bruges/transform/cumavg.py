# -*- coding: utf-8 -*-
"""
Average velocity equations.

All from Chris Liner, Elements of
3D Seismology, PennWell Press, 2004.
"""
import numpy as np


def v_rms(v, depth=None, time=None):
    """
    Cumulative RMS mean of a velocity log. You must provide either
    a depth or a time basis for the log.

    Args:
        v (ndarray): The velocity log.
        depth (ndarray): The depth values corresponding to the log.
        time (ndarray): The time values corresponding to the log.

    Returns:
        ndarray: The V_rms log.
    """
    if (depth is None) and (time is None):
        raise TypeError("You must provide a depth or time array")

    if depth is None:
        return np.sqrt(np.cumsum(v**2 / time) / np.cumsum(time))
    else:
        return np.sqrt(np.cumsum(depth * v) / np.cumsum(depth / v))


def v_avg(v, depth=None, time=None):
    """
    Cumulative average of a velocity log. You must provide either
    a depth or a time basis for the log.

    Args:
        v (ndarray): The velocity log.
        depth (ndarray): The depth values corresponding to the log.
        time (ndarray): The time values corresponding to the log.

    Returns:
        ndarray: The V_avg log.
    """
    if (depth is None) and (time is None):
        raise TypeError("You must provide a depth or time array")

    if depth is None:
        return np.cumsum(v * time) / np.cumsum(time)
    else:
        return np.cumsum(depth) / np.cumsum(depth / v)


def v_bac(v, rho, depth):
    """
    Cumulative Backus average of a velocity log. You must provide
    either a depth or a time basis for the log.

    For a non-cumulative version that can also provide sclaing for the
    V_s log, as well as quality factor, see bruges.anisotropy.backus.

    Args:
        v (ndarray): The velocity log.
        rho (ndarray): The density log.
        depth (ndarray): The depth values corresponding to the logs.

    Returns:
        ndarray: The V_bac log.
    """
    num = np.cumsum(depth**2)
    den = np.cumsum(rho * depth) * np.cumsum(depth/(v**2 * rho))
    return np.sqrt(num / den)
