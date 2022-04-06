"""
Convolution in n-dimensions.

:copyright: 2022 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
from bruges.util import apply_along_axis


def convolve(arr, v, axis=-1, verbose=False):
    """
    Convolve n-dimensional arr with a 1D wavelet or 2D wavelet bank.
    
    Args:
        arr (ndarray): The trace, or 2D section, or volume.
        v (ndarray): The wavelet, must be 1D function or a 2D wavelet 'bank'.
            If a wavelet bank, time should be on the last axis.
        axis (int): The time axis of the arr data. In other words, the axis
            corresponding to a single 'trace'. If you index into this axis,
            you will get a single 'trace'.
        verbose (bool): If True, print out the shapes of the inputs and output.

    Returns:
        ndarray: Discrete, linear convolution of `arr` and `v`.
    """
    if v.shape[-1] > arr.shape[axis]:
        raise TypeError("v (e.g. wavelet) must be shorter in time than arr (e.g. reflectivity).")

    arr_ = np.moveaxis(np.asanyarray(arr), axis, -1)

    # Compute the target shape of the final synthetic.
    outshape = v.shape[:-1] + arr_.shape

    # Force v and arr to both be 2D.
    bank = np.atleast_2d(v)   
    arr_2d = arr_.reshape((-1, arr_.shape[-1]))

    # Compute synthetic.
    syn = np.array([apply_along_axis(np.convolve, arr_2d, w, mode='same') for w in bank])

    pos = axis + v.ndim - 1

    out = np.moveaxis(syn.reshape(outshape), -1, pos)

    # Show the shapes of the data we're handling.
    if verbose:
        print(arr.shape, ' * ', v.shape, ' -> ', out.shape)

    return out
