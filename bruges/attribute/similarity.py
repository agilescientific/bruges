# -*- coding: utf-8 -*-
"""
Various methods to compute seismic similarity. 

:copyright: 2015 Joe Kington, SEG Tutorial https://library.seg.org/doi/pdf/10.1190/tle34121510.1
:license: Apache 2.0
"""
import numpy as np
from scipy.ndimage import generic_filter
from scipy.ndimage import gaussian_filter1d


def moving_window(traces, func, window):
    """
    Helper function for multi-trace attribute generation.
    This function applies a 3D function func to process a
    region of shape `window` over a dataset `data`.
    """
    wrapped = lambda x: func(x.reshape(window))
    return generic_filter(traces, wrapped, window)


def marfurt(traces):
    """
    Marfurt, K., V. Sudhaker, A. Gersztenkorn, K. D. Crawford, and S. E. Nissen, 1999,
    Coherency calculations in the presence of structural dip: GEOPHYSICS, 64, 104-111.
    doi:10.1190/1.1444508
    """
    i, x, t = traces.shape
    traces = traces.reshape(-1, t)
    square_sums = np.sum(traces, axis=0)**2
    sum_squares = np.sum(traces**2, axis=0)
    c = square_sums.sum() / (sum_squares.sum() + 1e-12)
    return c / (i * x)


def gersztenkorn(traces):
    """
    Gersztenkorn, A., and K. J. Marfurt, 1999, Eigenstructure‐based coherence
    computations as an aid to 3-D structural and stratigraphic mapping:
    GEOPHYSICS, 64, 1468-1479. doi:10.1190/1.1444651
    """
    # Stack traces in 3D traces into 2D array.
    traces = traces.reshape(-1, traces.shape[-1])

    # Calculate eigenvalues of covariance matrix.
    cov = traces.dot(traces.T)
    vals = np.linalg.eigvalsh(cov)
    return vals.max() / (vals.sum() + 1e-12)


def gradients(traces, sigma):
    grads = []
    for axis in range(3):
        grad = gaussian_filter1d(traces, sigma, axis=axis, order=1)
        grads.append(grad[..., np.newaxis])
    return np.concatenate(grads, axis=3)


def moving_window4d(grad, window, func):
    """Applies the given function *func* over a moving *window*, reducing 
    the input *grad* array from 4D to 3D."""
    # Pad in the spatial dimensions, but leave the gradient dimension unpadded.
    half_window = [(x // 2, x // 2) for x in window] + [(0, 0)]
    padded = np.pad(grad, half_window, mode='reflect')
    
    out = np.empty(grad.shape[:3], dtype=float)
    for i, j, k in np.ndindex(out.shape):
        region = padded[i:i+window[0], j:j+window[1], k:k+window[2], :]
        out[i,j,k] = func(region)
    return out


def gst_coherence_calc(region):
    region = region.reshape(-1, 3)
    gst = region.T.dot(region)
    eigs = np.sort(np.linalg.eigvalsh(gst))[::-1]
    return (eigs[0]-eigs[1]) / (eigs[0]+eigs[1])


def gradient_structure_tensor(seismic, window, sigma=1):
    """
    Randen, T., E. Monsen, C. Singe, A. Abrahamsen, J. Hansen, T. Saeter, and J. Schlaf, 2000,
    Three-dimensional texture attributes for seismic data analysis, 70th Annual International Meeting,
    SEG, Expanded Abstracts, 668-671.
    """
    grad = gradients(seismic, sigma)
    return moving_window4d(grad, window, gst_coherence_calc)


def similarity(traces, duration, dt, step_out=1, kind='gst', sigma=1):
    """
    Compute similarity for a seismic section using one of various methods.

    Expects time or depth to be in the last axis of a 2D or 3D input.

    :param traces: A 2D or 3D NumPy array arranged as (cdp, twt) or
        (iline, xline, twt).
    :param duration: The length in seconds of the window trace kernel
        used to calculate the similarity.
    :keyword dt (default=1): The sample interval of the traces in sec.
        (eg. 0.001, 0.002, ...). Will default to one, allowing
        duration to be given in samples.
    :keyword step_out (default=1):
        The number of adjacent traces to the kernel to compute similarity over.
    :keyword kind (default='gst'):
        The method to use for the computation. Can be "marfurt", "gersztenkorn"
        or "gst" (gradient structure tensor).
    :keyword sigma (default=1):
        The width of the Gaussian function used to compute gradients.
    """
    if traces.ndim == 2:
        window = 2*step_out+1, int(duration / dt)
    elif traces.ndim == 3:
        window = 2*step_out+1, int(duration / dt), 2*step_out+1
    else:
        raise NotImplementedError("Expected 2D or 3D seismic data.")

    methods = {
        "marfurt": moving_window(traces, marfurt, window),
        "gersztenkorn": moving_window(traces, gersztenkorn, window),
        "gst": gradient_structure_tensor(traces, window, sigma)
    }

    return methods[kind]
