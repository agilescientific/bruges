# -*- coding: utf-8 -*-
"""
Seismic wavelets.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
from collections import namedtuple

import numpy as np
from scipy.signal import hilbert
from scipy.signal import chirp


def sinc(duration, dt, f, return_t=False, taper='blackman'):
    """
    sinc function centered on t=0, with a dominant frequency of f Hz.

    If you pass a 1D array of frequencies, you get a wavelet bank in return.
    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (ndarray): Dominant frequency of the wavelet in Hz. If a sequence is
            passed, you will get a 2D array in return, one row per frequency.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis, where time is the range from -duration/2 to
            duration/2 in steps of dt.
        taper (str or function): The window or tapering function to apply.
            To use one of NumPy's functions, pass 'bartlett', 'blackman' (the
            default), 'hamming', or 'hanning'; to apply no tapering, pass
            'none'. To apply your own function, pass a function taking only
            the length of the window and returning the window function.

    Returns:
        ndarray. sinc wavelet(s) with centre frequency f sampled on t.
    """
    f = np.asanyarray(f).reshape(-1, 1)
    t = np.arange(-duration/2., duration/2., dt)
    t[t == 0] = 1e-12  # Avoid division by zero.
    f[f == 0] = 1e-12  # Avoid division by zero.
    w = np.squeeze(np.sin(2*np.pi*f*t) / (2*np.pi*f*t))

    if taper:
        funcs = {
            'bartlett': np.bartlett,
            'blackman': np.blackman,
            'hamming': np.hamming,
            'hanning': np.hanning,
            'none': lambda x: x,
        }
        func = funcs.get(taper, taper)
        w *= func(t.size)

    if return_t:
        RickerWavelet = namedtuple('RickerWavelet', ['amplitude', 'time'])
        return RickerWavelet(w, t)
    else:
        return w


def ricker(duration, dt, f, return_t=False):
    """
    Also known as the mexican hat wavelet, models the function:
    A =  (1-2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}

    If you pass a 1D array of frequencies, you get a wavelet bank in return.

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (ndarray): Centre frequency of the wavelet in Hz. If a sequence is
            passed, you will get a 2D array in return, one row per frequency.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis, where time is the range from -duration/2 to
            duration/2 in steps of dt.

    Returns:
        ndarray. Ricker wavelet(s) with centre frequency f sampled on t.
    """
    f = np.asanyarray(f).reshape(-1, 1)
    t = np.arange(-duration/2, duration/2, dt)
    pft2 = (np.pi * f * t)**2
    w = np.squeeze((1 - (2 * pft2)) * np.exp(-pft2))

    if return_t:
        RickerWavelet = namedtuple('RickerWavelet', ['amplitude', 'time'])
        return RickerWavelet(w, t)
    else:
        return w


def sweep(duration, dt, f,
          autocorrelate=True,
          return_t=False,
          taper='blackman',
          **kwargs):
    """
    Generates a linear frequency modulated wavelet (sweep). Wraps
    scipy.signal.chirp, adding dimensions as necessary.

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): is the sample interval in seconds (usually 0.001, 0.002,
            or 0.004)
        f (ndarray): Any sequence like (f1, f2). A list of lists will create a
            wavelet bank.
        autocorrelate (bool): Whether to autocorrelate the sweep(s) to create
            a wavelet. Default is `True`.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis, where time is the range from -duration/2 to
            duration/2 in steps of dt.
        taper (str or function): The window or tapering function to apply.
            To use one of NumPy's functions, pass 'bartlett', 'blackman' (the
            default), 'hamming', or 'hanning'; to apply no tapering, pass
            'none'. To apply your own function, pass a function taking only
            the length of the window and returning the window function.
        **kwargs: Further arguments are passed to scipy.signal.chirp. They are
            `method` ('linear','quadratic','logarithmic'), `phi` (phase offset
            in degrees), and `vertex_zero`.

    Returns:
        ndarray: The waveform.
    """
    t0, t1 = -duration/2, duration/2
    t = np.arange(t0, t1, dt)

    f = np.asanyarray(f).reshape(-1, 1)
    f1, f2 = f

    c = [chirp(t, f1_+(f2_-f1_)/2., t1, f2_, **kwargs)
         for f1_, f2_
         in zip(f1, f2)]

    if autocorrelate:
        w = [np.correlate(c_, c_, mode='same') for c_ in c]

    w = np.squeeze(w) / np.amax(w)

    if taper:
        funcs = {
            'bartlett': np.bartlett,
            'blackman': np.blackman,
            'hamming': np.hamming,
            'hanning': np.hanning,
            'none': lambda x: x,
        }
        func = funcs.get(taper, taper)
        w *= func(t.size)

    if return_t:
        Sweep = namedtuple('Sweep', ['amplitude', 'time'])
        return Sweep(w, t)
    else:
        return w


def ormsby(duration, dt, f, return_t=False):
    """
    The Ormsby wavelet requires four frequencies which together define a
    trapezoid shape in the spectrum. The Ormsby wavelet has several sidelobes,
    unlike Ricker wavelets.

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (usually 0.001, 0.002,
            or 0.004).
        f (ndarray): Sequence of form (f1, f2, f3, f4), or list of lists of
            frequencies, which will return a 2D wavelet bank.

    Returns:
        ndarray: A vector containing the Ormsby wavelet, or a bank of them.
    """
    f = np.asanyarray(f).reshape(-1, 1)

    try:
        f1, f2, f3, f4 = f
    except ValueError:
        raise ValueError("The last dimension must be 4")

    def numerator(f, t):
        return (np.sinc(f * t)**2) * ((np.pi * f) ** 2)

    pf43 = (np.pi * f4) - (np.pi * f3)
    pf21 = (np.pi * f2) - (np.pi * f1)

    t = np.arange(-duration/2, duration/2, dt)

    w = ((numerator(f4, t)/pf43) - (numerator(f3, t)/pf43) -
         (numerator(f2, t)/pf21) + (numerator(f1, t)/pf21))

    w = np.squeeze(w) / np.amax(w)

    if return_t:
        OrmsbyWavelet = namedtuple('OrmsbyWavelet', ['amplitude', 'time'])
        return OrmsbyWavelet(w, t)
    else:
        return w


def rotate_phase(w, phi, degrees=False):
    """
    Performs a phase rotation of wavelet or wavelet bank using:

    The analytic signal can be written in the form S(t) = A(t)exp(j*theta(t))
    where A(t) = magnitude(hilbert(w(t))) and theta(t) = angle(hilbert(w(t))
    then a constant phase rotation phi would produce the analytic signal
    S(t) = A(t)exp(j*(theta(t) + phi)). To get the non analytic signal
    we take real(S(t)) == A(t)cos(theta(t) + phi)
    == A(t)(cos(theta(t))cos(phi) - sin(theta(t))sin(phi)) <= trig idenity
    == w(t)cos(phi) - h(t)sin(phi)

    A = w(t)Cos(phi) - h(t)Sin(phi)

    Where w(t) is the wavelet and h(t) is its Hilbert transform.

    Args:
        w (ndarray): The wavelet vector, can be a 2D wavelet bank.
        phi (float): The phase rotation angle (in radians) to apply.
        degrees (bool): If phi is in degrees not radians.

    Returns:
        The phase rotated signal (or bank of signals).
    """
    if degrees:
        phi = phi * np.pi / 180.0
    a = hilbert(w, axis=0)
    w = (np.real(a) * np.cos(phi) - np.imag(a) * np.sin(phi))
    return w
