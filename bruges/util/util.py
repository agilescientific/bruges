"""
Utility functions.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import functools
import inspect
import warnings

import scipy.signal
import numpy as np


def deprecated(instructions):
    """
    Flags a method as deprecated. This decorator can be used to mark functions
    as deprecated. It will result in a warning being emitted when the function
    is used.
    Args:
        instructions (str): A human-friendly string of instructions, such
            as: 'Please migrate to add_proxy() ASAP.'
    Returns:
        The decorated function.
    """
    def decorator(func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            message = 'Call to deprecated function {}. {}'.format(
                func.__name__,
                instructions)

            frame = inspect.currentframe().f_back

            warnings.warn_explicit(message,
                                   category=DeprecationWarning,
                                   filename=inspect.getfile(frame.f_code),
                                   lineno=frame.f_lineno)

            return func(*args, **kwargs)

        return wrapper

    return decorator


greek = {
    'Alpha': 'Α',
    'Beta': 'Β',
    'Gamma': 'Γ',
    'Delta': 'Δ',
    'Epsilon': 'Ε',
    'Zeta': 'Ζ',
    'Eta': 'Η',
    'Kappa': 'Κ',
    'Lambda': 'Λ',
    'Mu': 'Μ',
    'Nu': 'Ν',
    'Phi': 'Φ',
    'Pi': 'Π',
    'Rho': 'Ρ',
    'Sigma': 'Σ',
    'Tau': 'Τ',
    'Upsilon': 'Υ',
    'Theta': 'Θ',
    'Chi': 'Χ',
    'Psi': 'Ψ',
    'Omega': 'Ω',
    'alpha': 'α',
    'beta': 'β',
    'gamma': 'γ',
    'delta': 'δ',
    'epsilon': 'ε',
    'zeta': 'ζ',
    'eta': 'η',
    'theta': 'θ',
    'kappa': 'κ',
    'lambda': 'λ',
    'mu': 'μ',
    'nu': 'ν',
    'pi': 'π',
    'rho': 'ρ',
    'sigma': 'σ',
    'tau': 'τ',
    'upsilon': 'υ',
    'phi': 'φ',
    'chi': 'χ',
    'psi': 'ψ',
    'omega': 'ω',
}


def rms(a, axis=None):
    """
    Calculates the RMS of an array.

    Args:
        a (ndarray). A sequence of numbers to apply the RMS to.
        axis (int). The axis along which to compute. If not given or None,
            the RMS for the whole array is computed.

    Returns:
        ndarray: The RMS of the array along the desired axis or axes.
    """
    a = np.array(a)
    if axis is None:
        div = a.size
    else:
        div = a.shape[axis]
    ms = np.sum(a**2.0, axis=axis) / div
    return np.sqrt(ms)


def moving_average(a, length, mode='same'):
    """
    Computes the mean in a moving window using convolution. For an alternative,
    as well as other kinds of average (median, mode, etc.), see bruges.filters.

    Example:
        >>> test = np.array([1,1,9,9,9,9,9,2,3,9,2,2,np.nan,1,1,1,1])
        >>> moving_average(test, 5, mode='same')
        array([ 2.2,  4. ,  5.8,  7.4,  9. ,  7.6,  6.4,  6.4,  5. ,  3.6,  nan,
                nan,  nan,  nan,  nan,  0.8,  0.6])
    """
    padded = np.pad(a, int(length/2), mode='edge')
    boxcar = np.ones(int(length))/length
    smooth = np.convolve(padded, boxcar, mode='same')
    return smooth[int(length/2):-int(length/2)]


@deprecated("Use bruges.filters() for moving linear and nonlinear statistics")
def moving_avg_conv(a, length, mode='same'):
    """
    Moving average via convolution. Keeping it for now for compatibility.
    """
    boxcar = np.ones(length)/length
    return np.convolve(a, boxcar, mode=mode)


@deprecated("Use bruges.filters() for moving linear and nonlinear statistics")
def moving_avg_fft(a, length, mode='same'):
    """
    Moving average via FFT convolution. Keeping it for now for compatibility.

    """
    boxcar = np.ones(length)/length
    return scipy.signal.fftconvolve(a, boxcar, mode=mode)


def normalize(a, new_min=0.0, new_max=1.0):
    """
    Normalize an array to [0,1] or to arbitrary new min and max.

    Args:
        a (ndarray): An array.
        new_min (float): The new min to scale to, default 0.
        new_max (float): The new max to scale to, default 1.

    Returns:
        ndarray. The normalized array.
    """
    a = np.array(a, dtype=np.float)
    n = (a - np.nanmin(a)) / np.nanmax(a - np.nanmin(a))
    return n * (new_max - new_min) + new_min


def nearest(a, num):
    """
    Finds the array's nearest value to a given num.

    Args:
        a (ndarray): An array.
        num (float): The value to find the nearest to.

    Returns:
        float. The normalized array.
    """
    a = np.array(a, dtype=float)
    return a.flat[np.abs(a - num).argmin()]


def next_pow2(num):
    """
    Calculates the next nearest power of 2 to the input. Uses
      2**ceil( log2( num ) ).

    Args:
        num (number): The number to round to the next power if two.

    Returns:
        number. The next power of 2 closest to num.
    """

    return int(2**np.ceil(np.log2(num)))


def top_and_tail(*arrays):
    """
    Top and tail all arrays to the non-NaN extent of the first array.

    E.g. crop the NaNs from the top and tail of a well log.

    Args:
        arrays (list): A list of arrays to treat.

    Returns:
        list: A list of treated arrays.
    """
    if len(arrays) > 1:
        for arr in arrays[1:]:
            assert len(arr) == len(arrays[0])
    nans = np.where(~np.isnan(arrays[0]))[0]
    first, last = nans[0], nans[-1]
    return [array[first:last+1] for array in arrays]


def extrapolate(a):
    """
    Extrapolate up and down an array from the first and last non-NaN samples.

    E.g. Continue the first and last non-NaN values of a log up and down.

    Args:
        a (ndarray): The array to treat.

    Returns:
        ndarray: The treated array.
    """
    a = np.array(a)
    nans = np.where(~np.isnan(a))[0]
    first, last = nans[0], nans[-1]
    a[:first] = a[first]
    a[last + 1:] = a[last]
    return a


def error_flag(pred, actual, dev = 1.0, method = 1):
    """Calculate the difference between a predicted and an actual curve 
    and return a log flagging large differences based on a user-defined distance 
    (in standard deviation units) from the mean difference

    Matteo Niccoli, October 2018
    
    Args:
        predicted (ndarray) = predicted log
        actual (ndarray) =  original log  
        dev  (float) = standard deviations to use, default 1
        error calcluation method (int), default 1
            1: difference between curves larger than mean difference plus dev
            
            2: curve slopes have opposite sign. Will require depth log for .diff method
            3: curve slopes of opposite sign OR difference larger than mean plus dev
    

    Returns:
    flag (ndarray) =  error flag curve"""
    
    flag = np.zeros(len(pred))
    err = np.abs(pred-actual)
    err_mean = np.mean(err)
    err_std = np.std(err)

    if method == 1:
        flag[np.where(err>(err_mean + (dev*err_std)))] = 1
     
     ###
     # add methods 2 and 3
     ###
    return flag


def apply_along_axis(func_1d, arr, kernel, **kwargs):
    """
    Apply 1D function across 2D slice as efficiently as possible.

    Although `np.apply_along_axis` seems to do well enough, map usually
    seems to end up beig a bit faster.

    Args:
        func_1d (function): the 1D function to apply, e.g. np.convolve. Should
            take 2 or more arguments: the

    Example
    >>> apply_along_axes(np.convolve, reflectivity_2d, wavelet, mode='same') 
    """
    mapobj = map(lambda tr: func_1d(tr, kernel, **kwargs), arr)
    return np.array(list(mapobj))


def sigmoid(start, stop, num):
    """
    Nonlinear space following a logistic function.
    
    The function is asymptotic; the parameters used in the sigmoid
    gets within 0.5% of the target thickness in a wedge increasing
    from 0 to 2x the original thickness.
    """
    x = np.linspace(-5.293305, 5.293305, num)
    return start + (stop-start) / (1 + np.exp(-x))
