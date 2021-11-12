# -*- coding: utf-8 -*-
"""
Complex trace attributes.

:copyright: 2021 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
from scipy.signal import hilbert


def instantaneous_amplitude(traces):
    """
    Compute instantaneous amplitude, also known as the envelope or
    reflection strength.

    The attribute is computed over the last dimension. That is, time should
    be in the last dimension, so a 100 inline, 100 crossline seismic volume
    with 250 time slices would have shape (100, 100, 250).

    Args:
        traces (ndarray): The data array to use for calculating energy.
    Returns:
        ndarray: An array the same dimensions as the input array.
    """
    return np.abs(hilbert(traces))

envelope = instantaneous_amplitude
reflection_strength = instantaneous_amplitude


def quadrature(traces):
    """
    Compute the quadrature trace.

    See https://wiki.seg.org/wiki/Instantaneous_attributes.

    Args:
        traces (ndarray): The data array to use for calculating energy.

    Returns:
        ndarray: An array the same dimensions as the input array.
    """
    h = hilbert(traces)
    return np.abs(h) * np.sin(np.log(h).imag)


def instantaneous_phase(traces):
    """
    Compute the instantaneous phase of the data.

    .. math::

        \phi(t) = {\rm Im}[\ln h(t)]

    See https://wiki.seg.org/wiki/Instantaneous_attributes.

    Args:
        traces (ndarray): The data array to use for calculating energy.

    Returns:
        ndarray: An array the same dimensions as the input array.
    """
    return np.angle(hilbert(traces))


def _inst_freq_claerbout(traces, dt):
    """
    Compute the instantaneous frequency using Claerbout's (1985) approximation.
    This is also the formulation given in Yilmaz (2001).

    Formulation from Barnes, A, 2016, Handbook of Poststack Seismic Attributes,
    SEG Books.

    Args:
        traces (ndarray): The data array to use for calculating energy.
        dt (float): The sample interval in seconds, e.g. 0.004 for 4 ms sample
            interval (250 Hz sample frequency).
    Returns:
        ndarray: An array the same dimensions as the input array.
    """
    h = hilbert(traces)
    term = (h[1:] - h[:-1]) / (h[1:] + h[:-1])
    return (1 / (np.pi * dt)) * np.imag(term)


def _inst_freq_scheuer_oldenburg(traces, dt):
    """Instantaneous frequency after Scheuer & Oldenburg (1988).
    
    Scheuer, TE and DW Oldenburg (1988). Local phase velocity from complex seismic data.
    Geophysics 53 (12), p1503. DOI: http://dx.doi.org/10.1190/1.1442431.
    
    Formulation from Barnes, A, 2016, Handbook of Poststack Seismic Attributes,
    SEG Books:
    
    .. math::

        f_i(t) = \frac{1}{2\pi} \ \mathrm{Im} \left[\frac{h'(t)}{h(t)} \right]
               \approx \frac{1}{\pi T} \ \mathrm{Im} \left[\frac{h(t+T) - h(t)}{h(t+T) + h(t)} \right] 

    Args:
        traces (ndarray): The data array to use for calculating energy.
        dt (float): The sample interval in seconds, e.g. 0.004 for 4 ms sample
            interval (250 Hz sample frequency).
    Returns:
        ndarray: An array the same dimensions as the input array.
    """
    y = quadrature(traces)
    expr = (traces[:-1] * y[1:] - traces[1:] * y[:-1]) / (traces[:-1] * traces[1:] + y[1:] * y[:-1])
    return (1 / (2 * np.pi * dt)) * np.arctan(expr)


def instantaneous_frequency(traces, dt, kind='so', percentile_clip=99):
    """
    Compute instantaneous frequency with a discrete approximation.

    The attribute is computed over the last dimension. That is, time should
    be in the last dimension, so a 100 inline, 100 crossline seismic volume
    with 250 time slices would have shape (100, 100, 250).

    These attributes can be noisy so a percentile clips is applied.

    Args:
        traces (ndarray): The data array to use for calculating energy.
        dt (float): The sample interval in seconds, e.g. 0.004 for 4 ms sample
            interval (250 Hz sample frequency).
        kind (str): "scipy", "claerbout" or "so" to denote a naive method from
            the SciPy docs, Claerbout's (1985) method or that of Scheuer & Oldenburg
            (1988). Claerbout's approximation is not stable above about half the
            Nyquist frequency (i.e. one quarter of the sampling frequency). The
            SciPy implementation is not recommended for seismic data.
        percentile_clip (float): Percentile at which to clip the data.
            Computed from the absolute values, clipped symmetrically
            at -p and +p, where p is the value at the 98th percentile.
    Returns:
        ndarray: An array the same dimensions as the input array.
    """
    methods = {'claerbout': _inst_freq_claerbout,
               'so': _inst_freq_scheuer_oldenburg,
               'scipy': lambda traces, dt: np.diff(instantaneous_phase(traces)) / (2.0 * np.pi * dt),
               }
    func = methods.get(kind)
    if func is None:
        m = f'{kind} is not supported, use "so" (Scheuer-Oldenburg, recommended), "claerbout" or "scipy".'
        raise NotImplementedError(m)
    f = func(traces, dt)
    p = np.percentile(np.abs(f), percentile_clip)
    return np.clip(f, a_min=-p, a_max=p)
