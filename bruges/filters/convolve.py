"""
Convolution in n-dimensions.

:copyright: 2022 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
from bruges.util import apply_along_axis


def convolve(reflectivity, wavelet, axis=-1):
    """
    Convolve n-dimensional reflectivity with a 1D wavelet or 2D wavelet bank.
    
    Args:
        reflectivity (ndarray): The reflectivity trace, or 2D section, or volume.
        wavelet (ndarray): The wavelet, must be 1D function or a 2D wavelet 'bank'.
            If a wavelet bank, time should be on the last axis.
        axis (int): The time axis of the reflectivity data. In other words,
            the axis corresponding to a single 'trace'. If you index into this
            axis, you will get a single 'trace'.

    Returns:
        ndarray: Discrete, linear convolution of `reflectivity` and `wavelet`.
    """
    if wavelet.shape[-1] > reflectivity.shape[axis]:
        raise TypeError("Wavelet must shorter in time than the reflectivity.")

    reflectivity_ = np.moveaxis(np.asanyarray(reflectivity), axis, 0)

    # Compute the target shape of the final synthetic.
    outshape = wavelet.shape[:-1] + reflectivity_.shape

    # Force wavelet and reflectivity to both be 2D.
    bank = np.atleast_2d(wavelet)   
    reflectivity_2d = reflectivity_.reshape((-1, reflectivity_.shape[-1]))

    # Compute synthetic.
    syn = np.array([apply_along_axis(np.convolve, reflectivity_2d, w, mode='same') for w in bank])

    pos = wavelet.ndim - 1

    return np.moveaxis(syn.reshape(outshape), pos, axis)
