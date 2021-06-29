"""
Smoothers.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
import scipy.ndimage
import scipy.signal

from bruges.bruges import BrugesError
from bruges.filters import apply_along_axis
from bruges.util import nearest
from bruges.util import rms as rms_


def mean(arr, size=5):
    """
    A linear n-D smoothing filter. Can be used as a moving average on 1D data.

    Args:
        arr (ndarray): an n-dimensional array, such as a seismic horizon.
        size (int): the kernel size, e.g. 5 for 5x5. Should be odd,
            rounded up if not.

    Returns:
        ndarray: the resulting smoothed array.
    """
    arr = np.array(arr, dtype=np.float)

    if not size // 2:
        size += 1

    return scipy.ndimage.generic_filter(arr, np.mean, size=size)


def rms(arr, size=5):
    """
    A linear n-D smoothing filter. Can be used as a moving average on 1D data.

    Args:
        arr (ndarray): an n-dimensional array, such as a seismic horizon.
        size (int): the kernel size, e.g. 5 for 5x5. Should be odd,
            rounded up if not.

    Returns:
        ndarray: the resulting smoothed array.
    """
    arr = np.array(arr, dtype=np.float)

    if not size // 2:
        size += 1

    return scipy.ndimage.generic_filter(arr, rms_, size=size)


def median(arr, size=5):
    """
    A nonlinear n-D edge-preserving smoothing filter.

    Args:
        arr (ndarray): an n-dimensional array, such as a seismic horizon.
        size (int): the kernel size, e.g. 5 for 5x5. Should be odd,
            rounded up if not.

    Returns:
        ndarray: the resulting smoothed array.
    """
    arr = np.array(arr, dtype=np.float)

    if not size // 2:
        size += 1

    return scipy.ndimage.generic_filter(arr, np.median, size=size)


def mode(arr, size=5, tie='smallest'):
    """
    A nonlinear n-D categorical smoothing filter. Use this to filter non-
    continuous variables, such as categorical integers, e.g. to label facies.

    Args:
        arr (ndarray): an n-dimensional array, such as a seismic horizon.
        size (int): the kernel size, e.g. 5 for 5x5. Should be odd,
            rounded up if not.
        tie (str): `'smallest'` or `'largest`'. In the event of a tie (i.e. two
            or more values having the same count in the kernel), whether to
            give back the smallest of the tying values, or the largest.

    Returns:
        ndarray: the resulting smoothed array.
    """
    def func(this, tie):
        if tie == 'smallest':
            m, _ = scipy.stats.mode(this)
        else:
            m, _ = -scipy.stats.mode(-this)
        return np.squeeze(m)

    arr = np.array(arr, dtype=np.float)

    if not size // 2:
        size += 1

    return scipy.ndimage.generic_filter(arr, func, size=size,
                                        extra_keywords={'tie': tie}
                                       )


def snn(arr, size=5, include=True):
    """
    Symmetric nearest neighbour, a nonlinear 2D smoothing filter.
    http://subsurfwiki.org/wiki/Symmetric_nearest_neighbour_filter

    Args:
        arr (ndarray): a 2D array, such as a seismic horizon.
        size (int): the kernel size, e.g. 5 for 5x5. Should be odd,
            rounded up if not.
        include (bool): whether to include the central pixel itself.

    Returns:
        ndarray: the resulting smoothed array.
    """
    def func(this, pairs, include):
        """
        Deal with this patch.
        """
        centre = this[this.size // 2]
        select = [nearest(this[p], centre) for p in pairs]
        if include:
            select += [centre]
        return np.mean(select)

    arr = np.array(arr, dtype=np.float)
    if arr.ndim != 2:
        raise BrugesError("arr must have 2-dimensions")

    if not size // 2:
        size += 1

    pairs = [[i, size**2-1 - i] for i in range(size**2 // 2)]
    return scipy.ndimage.generic_filter(arr,
                                        func,
                                        size=size,
                                        extra_keywords={'pairs': pairs,
                                                        'include': include}
                                       )


def kuwahara(arr, size=5):
    """
    Kuwahara, a nonlinear 2D smoothing filter.
    http://subsurfwiki.org/wiki/Kuwahara_filter

    Args:
        arr (ndarray): a 2D array, such as a seismic horizon.
        size (int): the kernel size, e.g. 5 for 5x5. Should be odd,
            rounded up if not.

    Returns:
        ndarray: the resulting smoothed array.
    """
    def func(this, s, k):
        """
        Deal with this patch.
        """
        t = this.reshape((s, s))
        sub = np.array([t[:k, :k].flatten(),
                        t[:k, k-1:].flatten(),
                        t[k-1:, :k].flatten(),
                        t[k-1:, k-1:].flatten()]
                      )
        select = sub[np.argmin(np.var(sub, axis=1))]
        return np.mean(select)

    arr = np.array(arr, dtype=np.float)
    if arr.ndim != 2:
        raise BrugesError("arr must have 2-dimensions")

    if not size // 2:
        size += 1

    k = int(np.ceil(size / 2))

    return scipy.ndimage.generic_filter(arr,
                                        func,
                                        size=size,
                                        extra_keywords={'s': size,
                                                        'k': k,
                                                       }
                                       )


def conservative(arr, size=5, supercon=False):
    """
    Conservative, a nonlinear n-D despiking filter. Very conservative! Only
    changes centre value if it is outside the range of all the other values
    in the kernel. Read http://subsurfwiki.org/wiki/Conservative_filter

    Args:
        arr (ndarray): an n-dimensional array, such as a seismic horizon.
        size (int): the kernel size, e.g. 5 for 5x5 (in a 2D arr). Should be
            odd, rounded up if not.
        supercon (bool): whether to be superconservative. If True, replaces
            pixel with min or max of kernel. If False (default), replaces pixel
            with mean of kernel.

    Returns:
        ndarray: the resulting smoothed array.
    """
    def func(this, k, supercon):
        this = this.flatten()
        centre = this[k]
        rest = [this[:k], this[-k:]]
        mi, ma = np.nanmin(rest), np.nanmax(rest)
        if centre < mi:
            return mi if supercon else np.mean(rest)
        elif centre > ma:
            return ma if supercon else np.mean(rest)
        else:
            return centre

    arr = np.array(arr, dtype=np.float)

    if not size // 2:
        size += 1

    k = int(np.floor(size**arr.ndim / 2))

    return scipy.ndimage.generic_filter(arr,
                                        func,
                                        size=size,
                                        extra_keywords={'k': k,
                                                        'supercon': supercon,
                                                       }
                                       )


def rotate_phase(s, phi, degrees=False):
    """
    Performs a phase rotation of wavelet or wavelet bank using:

    .. math::

        A = w(t)\cos(\phi) - h(t)\sin(\phi)

    where w(t) is the wavelet and h(t) is its Hilbert transform.

    The analytic signal can be written in the form S(t) = A(t)exp(j*theta(t))
    where A(t) = magnitude(hilbert(w(t))) and theta(t) = angle(hilbert(w(t))
    then a constant phase rotation phi would produce the analytic signal
    S(t) = A(t)exp(j*(theta(t) + phi)). To get the non analytic signal
    we take real(S(t)) == A(t)cos(theta(t) + phi)
    == A(t)(cos(theta(t))cos(phi) - sin(theta(t))sin(phi)) <= trig identity
    == w(t)cos(phi) - h(t)sin(phi)

    Args:
        w (ndarray): The wavelet vector, can be a 2D wavelet bank.
        phi (float): The phase rotation angle (in radians) to apply.
        degrees (bool): If phi is in degrees not radians.

    Returns:
        The phase rotated signal (or bank of signals).
    """
    # Make sure the data is at least 2D to apply_along
    data = np.atleast_2d(s)

    # Get Hilbert transform. This will be 2D.
    a = apply_along_axis(scipy.signal.hilbert, data, axis=0)

    # Transform angles into what we need.
    phi = np.asanyarray(phi).reshape(-1, 1, 1)
    if degrees:
        phi = np.radians(phi)
        
    rotated = np.real(a) * np.cos(phi)  -  np.imag(a) * np.sin(phi)
    return np.squeeze(rotated)
