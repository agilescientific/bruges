"""
Convolution in n-dimensions.

:copyright: 2019 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
from bruges.util import apply_along_axis


def convolve(reflectivity, wavelet):
    """
    Convolve n-dimensional reflectivity with a 1D wavelet or 2D wavelet bank.
    
    Args
    reflectivity (ndarray): The reflectivity trace, or 2D section, or volume.
    wavelet (ndarray): The wavelet, must be 1D function or a 2D wavelet 'bank'.
        If a wavelet bank, time should be on the last axis.
    """
    # Compute the target shape of the final synthetic.
    outshape = wavelet.shape[:-1] + reflectivity.shape

    # Force wavelet and reflectivity to both be 2D.
    bank = np.atleast_2d(wavelet)   
    reflectivity_2d = reflectivity.reshape((-1, reflectivity.shape[-1]))

    # Compute synthetic, which will always be 3D.
    syn = np.array([apply_along_axis(np.convolve, reflectivity_2d, w, mode='same') for w in bank])

    return syn.reshape(outshape)


def apply_along_axis(func_1d, arr, *args, **kwargs):
    mapobj = map(lambda tr: func_1d(tr, *args, **kwargs), arr)
    return np.array(list(mapobj))
