#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utility functions.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import scipy.signal
import numpy as np


def rms(a):
    """
    Calculates the RMS of an array.

    :param a: An array.

    :returns: The RMS of the array.

    """

    return np.sqrt(np.sum(a**2.0)/a.size)


def moving_average(a, length, mode='valid'):
    """
    Computes the mean in a moving window. Naive implementation.

    Example:
        >>> test = np.array([1,9,9,9,9,9,9,2,3,9,2,2,3,1,1,1,1,3,4,9,9,9,8,3])
        >>> moving_average(test, 7, mode='same')
        [ 4.4285714  5.571428  6.7142857  7.8571428  8.          7.1428571
          7.1428571  6.142857  5.1428571  4.2857142  3.1428571  3.
          2.7142857  1.571428  1.7142857  2.          2.857142  4.
          5.1428571  6.142857  6.4285714  6.1428571  5.7142857  4.5714285 ]

    TODO:
        Other types of average.

    """
    pad = np.floor(length/2)

    if mode == 'full':
        pad *= 2

    # Make a padded version, paddding with first and last values
    r = np.empty(a.shape[0] + 2*pad)
    r[:pad] = a[0]
    r[pad:-pad] = a
    r[-pad:] = a[-1]

    # Cumsum with shifting trick
    s = np.cumsum(r, dtype=float)
    s[length:] = s[length:] - s[:-length]
    out = s[length-1:]/length

    # Decide what to return
    if mode == 'same':
        if out.shape[0] != a.shape[0]:
            # If size doesn't match, then interpolate.
            out = (out[:-1, ...] + out[1:, ...]) / 2
        return out
    elif mode == 'valid':
        return out[pad:-pad]
    else:  # mode=='full' and we used a double pad
        return out


def moving_avg_conv(a, length):
    """
    Moving average via convolution. Seems slower than naive.

    """
    boxcar = np.ones(length)/length
    return np.convolve(a, boxcar, mode="same")


def moving_avg_fft(a, length):
    """
    Moving average via FFT convolution. Seems slower than naive.

    """
    boxcar = np.ones(length)/length
    return scipy.signal.fftconvolve(a, boxcar, mode="same")


def normalize(a, new_min=0.0, new_max=1.0):
    """
    Normalize an array to [0,1] or to
    arbitrary new min and max.

    :param a: An array.
    :param new_min: A float to be the new min, default 0.
    :param new_max: A float to be the new max, default 1.

    :returns: The normalized array.
    """

    n = (a - np.amin(a)) / np.amax(a - np.amin(a))
    return n * (new_max - new_min) + new_min


def next_pow2(num):
    """
    Calculates the next nearest power of 2 to the input. Uses
      2**ceil( log2( num ) ).

    :param num: The number to round to the next power if two.

    :returns: the next power of 2 closest to num.
    """

    return int(2**np.ceil(np.log2(num)))


def top_and_tail(a, b=np.array([]), c=np.array([])):
    """
    Top and tail up to 3 arrays to the non-NaN extent of the first array.

    E.g. crop the NaNs from the top and tail of a well log.

    TODO:
        Make this work for an arbitrary number of arrays.

    """
    nans = np.where(~np.isnan(a))[0]
    first, last = nans[0], nans[-1]
    a = a[first:last]
    if b.any():
        b = b[first:last]
        if c.any():
            c = c[first:last]
            return a, b, c
        return a, b
    return a


def extrapolate(a):
    """
    Extrapolate up and down an array from the first and last non-NaN samples.

    E.g. Continue the first and last non-NaN values of a log up and down.

    """
    nans = np.where(~np.isnan(a))[0]
    first, last = nans[0], nans[-1]
    a[:first] = a[first]
    a[last + 1:] = a[last]
    return a
