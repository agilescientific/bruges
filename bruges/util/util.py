"""
Utility functions.

:copyright: 2021 Agile Scientific
:license: Apache 2.0
"""
import functools
import inspect
import warnings

import scipy.signal
import numpy as np

from skimage.util.shape import view_as_windows
from typing import Tuple, Union, cast

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


@deprecated("Use bruges.filters for moving linear and nonlinear statistics")
def moving_average(a, length, mode='same'):
    """
    Computes the mean in a moving window using convolution. For an alternative,
    as well as other kinds of average (median, mode, etc.), see bruges.filters.

    Args:
        a (ndarray): The array to average.
        length (int): The length of the moving window.
        mode (str): The mode of the convolution, see np.convolve.

    Returns:
        ndarray: The moving average.

    Example:
        >>> test = np.array([1,1,9,9,9,9,9,2,3,9,2,2,np.nan,1,1,1,1])
        >>> moving_average(test, 5, mode='same')
        array([ 2.2,  4. ,  5.8,  7.4,  9. ,  7.6,  6.4,  6.4,  5. ,  3.6,  nan,
                nan,  nan,  nan,  nan,  0.8,  0.6])
    """
    warnings.warn("Use bruges.filters for moving linear and nonlinear statistics", DeprecationWarning)
    warnings.warn("The mode argument is ignored in this function. Use bruges.filters for moving linear and nonlinear statistics", stacklevel=2)
    padded = np.pad(a, int(length/2), mode='edge')
    boxcar = np.ones(int(length))/int(length)
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


def error_flag(pred, actual, dev=1.0, method=1):
    """
    Calculate the difference between a predicted and an actual curve
    and return a log flagging large differences based on a user-defined
    distance (in standard deviation units) from the mean difference.

    Author:
        Matteo Niccoli, 2018

    Args:
        predicted (ndarray): predicted log.
        actual (ndarray):  original log.
        dev (float): standard deviations to use, default 1
        error calcluation method (int): default 1
            1: difference between curves larger than mean difference plus dev
            2: curve slopes have opposite sign. Will require depth log for .diff method
            3: curve slopes of opposite sign OR difference larger than mean plus dev
    Returns:
        flag (ndarray) =  error flag curve
    """
    flag = np.zeros(len(pred))
    err = np.abs(pred-actual)
    err_mean = np.mean(err)
    err_std = np.std(err)

    if method == 1:
        flag[np.where(err > (err_mean + (dev * err_std)))] = 1

    return flag


def apply_along_axis(func_1d, arr, *args, **kwargs):
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
    mapobj = map(lambda tr: func_1d(tr, *args, **kwargs), arr)
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


def root(start, stop, num):
    """
    Nonlinear space following a sqrt function.
    """
    x = np.linspace(0, 1, num)
    y = np.sqrt(x)
    return min(start, stop) + abs(stop-start) * y


def power(start, stop, num):
    """
    Nonlinear space following a power function.
    """
    x = np.linspace(0, 8, num)
    y = 1 - 2**-x
    return min(start, stop) + abs(stop-start) * y


Imsize = Union[Tuple[int, int], Tuple[int, int, int]]

def patchify(image: np.ndarray, patch_size: Imsize, step: int = 1) -> np.ndarray:
    """
    Split a 2D or 3D image into small patches given the patch size.
    Parameters
    ----------
    image: the image to be split. It can be 2d (m, n) or 3d (k, m, n)
    patch_size: the size of a single patch
    step: the step size between patches
    Examples
    --------
    >>> image = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    >>> patches = patchify(image, (2, 2), step=1)  # split image into 2*3 small 2*2 patches.
    >>> assert patches.shape == (2, 3, 2, 2)
    >>> reconstructed_image = unpatchify(patches, image.shape)
    >>> assert (reconstructed_image == image).all()
    """
    return view_as_windows(image, patch_size, step)


def unpatchify(patches: np.ndarray, imsize: Imsize) -> np.ndarray:
    """
    Merge patches into the orignal image
    Parameters
    ----------
    patches: the patches to merge. It can be patches for a 2d image (k, l, m, n)
             or 3d volume (i, j, k, l, m, n)
    imsize: the size of the original image or volume
    Examples
    --------
    >>> image = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    >>> patches = patchify(image, (2, 2), step=1)  # split image into 2*3 small 2*2 patches.
    >>> assert patches.shape == (2, 3, 2, 2)
    >>> reconstructed_image = unpatchify(patches, image.shape)
    >>> assert (reconstructed_image == image).all()
    """

    assert len(patches.shape) / 2 == len(
        imsize
    ), "The patches dimension is not equal to the original image size"

    if len(patches.shape) == 4:
        return _unpatchify2d(patches, cast(Tuple[int, int], imsize))
    elif len(patches.shape) == 6:
        return _unpatchify3d(patches, cast(Tuple[int, int, int], imsize))
    else:
        raise NotImplementedError(
            "Unpatchify only supports a matrix of 2D patches (k, l, m, n)"
            f"or 3D volumes (i, j, k, l, m, n), but got: {patches.shape}"
        )


def _unpatchify2d(  # pylint: disable=too-many-locals
    patches: np.ndarray, imsize: Tuple[int, int]
) -> np.ndarray:

    assert len(patches.shape) == 4

    i_h, i_w = imsize
    image = np.zeros(imsize, dtype=patches.dtype)

    n_h, n_w, p_h, p_w = patches.shape

    s_w = 0 if n_w <= 1 else (i_w - p_w) / (n_w - 1)
    s_h = 0 if n_h <= 1 else (i_h - p_h) / (n_h - 1)

    # The step size should be same for all patches, otherwise the patches are unable
    # to reconstruct into a image
    if int(s_w) != s_w:
        raise NonUniformStepSizeError(i_w, n_w, p_w, s_w)
    if int(s_h) != s_h:
        raise NonUniformStepSizeError(i_h, n_h, p_h, s_h)
    s_w = int(s_w)
    s_h = int(s_h)

    i, j = 0, 0

    while True:
        i_o, j_o = i * s_h, j * s_w

        image[i_o : i_o + p_h, j_o : j_o + p_w] = patches[i, j]

        if j < n_w - 1:
            j = min((j_o + p_w) // s_w, n_w - 1)
        elif i < n_h - 1 and j >= n_w - 1:
            # Go to next row
            i = min((i_o + p_h) // s_h, n_h - 1)
            j = 0
        elif i >= n_h - 1 and j >= n_w - 1:
            # Finished
            break
        else:
            raise RuntimeError("Unreachable")

    return image


def _unpatchify3d(  # pylint: disable=too-many-locals
    patches: np.ndarray, imsize: Tuple[int, int, int]
) -> np.ndarray:

    assert len(patches.shape) == 6

    i_h, i_w, i_c = imsize
    image = np.zeros(imsize, dtype=patches.dtype)

    n_h, n_w, n_c, p_h, p_w, p_c = patches.shape

    s_w = 0 if n_w <= 1 else (i_w - p_w) / (n_w - 1)
    s_h = 0 if n_h <= 1 else (i_h - p_h) / (n_h - 1)
    s_c = 0 if n_c <= 1 else (i_c - p_c) / (n_c - 1)

    # The step size should be same for all patches, otherwise the patches are unable
    # to reconstruct into a image
    if int(s_w) != s_w:
        raise NonUniformStepSizeError(i_w, n_w, p_w, s_w)
    if int(s_h) != s_h:
        raise NonUniformStepSizeError(i_h, n_h, p_h, s_h)
    if int(s_c) != s_c:
        raise NonUniformStepSizeError(i_c, n_c, p_c, s_c)

    s_w = int(s_w)
    s_h = int(s_h)
    s_c = int(s_c)

    i, j, k = 0, 0, 0

    while True:

        i_o, j_o, k_o = i * s_h, j * s_w, k * s_c

        image[i_o : i_o + p_h, j_o : j_o + p_w, k_o : k_o + p_c] = patches[i, j, k]

        if k < n_c - 1:
            k = min((k_o + p_c) // s_c, n_c - 1)
        elif j < n_w - 1 and k >= n_c - 1:
            j = min((j_o + p_w) // s_w, n_w - 1)
            k = 0
        elif i < n_h - 1 and j >= n_w - 1 and k >= n_c - 1:
            i = min((i_o + p_h) // s_h, n_h - 1)
            j = 0
            k = 0
        elif i >= n_h - 1 and j >= n_w - 1 and k >= n_c - 1:
            # Finished
            break
        else:
            raise RuntimeError("Unreachable")

    return image


class NonUniformStepSizeError(RuntimeError):
    def __init__(
        self, imsize: int, n_patches: int, patch_size: int, step_size: float
    ) -> None:
        super().__init__(imsize, n_patches, patch_size, step_size)
        self.n_patches = n_patches
        self.patch_size = patch_size
        self.imsize = imsize
        self.step_size = step_size

    def __repr__(self) -> str:
        return f"Unpatchify only supports reconstructing image with a uniform step size for all patches. \
However, reconstructing {self.n_patches} x {self.patch_size}px patches to an {self.imsize} image requires {self.step_size} as step size, which is not an integer."

    def __str__(self) -> str:
        return self.__repr__()
