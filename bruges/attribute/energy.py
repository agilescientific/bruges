# -*- coding: utf-8 -*-
"""
Mean-squared energy measurement.

:copyright: 2019 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
from bruges.filters import convolve


def energy(traces, duration, dt=1):
    """
    Compute an mean-squared energy measurement on seismic data.

    The attribute is computed over the last dimension. That is, time should
    be in the last dimension, so a 100 inline, 100 crossline seismic volume
    with 250 time slices would have shape (100, 100, 250).

    Args:
        traces (ndarray): The data array to use for calculating energy.
        duration (float): the time duration of the window (in seconds), or
            samples if dt=1.
        dt (float): the sample interval of the data (in seconds). Defaults
            to 1 so duration can be in samples.
    Returns:
        ndarray: An array the same dimensions as the input array.
    """
    data = traces.astype(np.float).reshape(-1, traces.shape[-1])
    n_samples = int(duration / dt)
    window = np.ones(n_samples) / n_samples
    energy = convolve(data**2, window)
    return energy.reshape(traces.shape)
  