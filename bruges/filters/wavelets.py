# -*- coding: utf-8 -*-
"""
Seismic wavelets.

:copyright: 2019 Agile Geoscience
:license: Apache 2.0
"""
from collections import namedtuple
import warnings

import numpy as np
import scipy.signal


def generic(func, duration, dt, f, return_t=False, taper='blackman'):
    """
    Generic wavelet generator: applies a window to a continuous function.

    Args:
        func (function): The continuous function, taking t, f as arguments.
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (array-like): Dominant frequency of the wavelet in Hz. If a sequence is
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
        ndarray. wavelet(s) with centre frequency f sampled on t. If you
            passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    if not return_t:
        warnings.warn("In future releases, return_t will be True by default.", FutureWarning)

    f = np.asanyarray(f).reshape(-1, 1)
    t = np.arange(-duration/2., duration/2., dt)
    t[t == 0] = 1e-12  # Avoid division by zero.
    f[f == 0] = 1e-12  # Avoid division by zero.

    w = np.squeeze(func(t, f))

    if taper:
        tapers = {
            'bartlett': np.bartlett,
            'blackman': np.blackman,
            'hamming': np.hamming,
            'hanning': np.hanning,
            'none': lambda _: 1,
        }
        taper = tapers.get(taper, taper)
        w *= taper(t.size)

    if return_t:
        Wavelet = namedtuple('Wavelet', ['amplitude', 'time'])
        return Wavelet(w, t)
    else:
        return w


def sinc(duration, dt, f, return_t=False, taper='blackman'):
    """
    sinc function centered on t=0, with a dominant frequency of f Hz.

    If you pass a 1D array of frequencies, you get a wavelet bank in return.

    .. plot::
        plt.plot(bruges.filters.sinc(.5, 0.002, 40))

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (array-like): Dominant frequency of the wavelet in Hz. If a sequence is
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
        ndarray. sinc wavelet(s) with centre frequency f sampled on t. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    def func(t_, f_):
        return np.sin(2*np.pi*f_*t_) / (2*np.pi*f_*t_)

    return generic(func, duration, dt, f, return_t, taper)


def cosine(duration, dt, f, return_t=False, taper='gaussian', sigma=None):
    """
    With the default Gaussian window, equivalent to a 'modified Morlet'
    also sometimes called a 'Gabor' wavelet. The `bruges.filters.gabor`
    function returns a similar shape, but with a higher mean frequancy,
    somewhere between a Ricker and a cosine (pure tone).

    If you pass a 1D array of frequencies, you get a wavelet bank in return.

    .. plot::
        plt.plot(bruges.filters.cosine(.5, 0.002, 40))

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (array-like): Dominant frequency of the wavelet in Hz. If a sequence is
            passed, you will get a 2D array in return, one row per frequency.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis, where time is the range from -duration/2 to
            duration/2 in steps of dt.
        taper (str or function): The window or tapering function to apply.
            To use one of NumPy's functions, pass 'bartlett', 'blackman' (the
            default), 'hamming', or 'hanning'; to apply no tapering, pass
            'none'. To apply your own function, pass a function taking only
            the length of the window and returning the window function.
        sigma (float): Width of the default Gaussian window, in seconds.
            Defaults to 1/8 of the duration.

    Returns:
        ndarray. sinc wavelet(s) with centre frequency f sampled on t. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    if sigma is None:
        sigma = duration / 8

    def func(t_, f_):
        return np.cos(2 * np.pi * f_ * t_)

    def taper(length):
        return scipy.signal.gaussian(length, sigma/dt)

    return generic(func, duration, dt, f, return_t, taper)


def gabor(duration, dt, f, return_t=False):
    """
    Generates a Gabor wavelet with a peak frequency f0 at time t.

    https://en.wikipedia.org/wiki/Gabor_wavelet

    If you pass a 1D array of frequencies, you get a wavelet bank in return.

    .. plot::
        plt.plot(bruges.filters.gabor(.5, 0.002, 40))

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (array-like): Centre frequency of the wavelet in Hz. If a sequence is
            passed, you will get a 2D array in return, one row per frequency.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis, where time is the range from -duration/2 to
            duration/2 in steps of dt.

    Returns:
        ndarray. Gabor wavelet(s) with centre frequency f sampled on t. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    def func(t_, f_):
        return np.exp(-2 * f_**2 * t_**2) * np.cos(2 * np.pi * f_ * t_)

    return generic(func, duration, dt, f, return_t)


def ricker(duration, dt, f, return_t=False):
    """
    Also known as the mexican hat wavelet, models the function:

    .. math::
        A =  (1 - 2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}

    If you pass a 1D array of frequencies, you get a wavelet bank in return.

    .. plot::
        plt.plot(bruges.filters.ricker(.5, 0.002, 40))

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (array-like): Centre frequency of the wavelet in Hz. If a sequence is
            passed, you will get a 2D array in return, one row per frequency.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis, where time is the range from -duration/2 to
            duration/2 in steps of dt.

    Returns:
        ndarray. Ricker wavelet(s) with centre frequency f sampled on t. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    if not return_t:
        warnings.warn("In future releases, return_t will be True by default.", FutureWarning)

    f = np.asanyarray(f).reshape(-1, 1)
    t = np.arange(-duration/2, duration/2, dt)
    pft2 = (np.pi * f * t)**2
    w = np.squeeze((1 - (2 * pft2)) * np.exp(-pft2))

    if return_t:
        RickerWavelet = namedtuple('RickerWavelet', ['amplitude', 'time'])
        return RickerWavelet(w, t)
    else:
        return w


def klauder(duration, dt, f,
            autocorrelate=True,
            return_t=False,
            taper='blackman',
            **kwargs):
    """
    By default, gives the autocorrelation of a linear frequency modulated
    wavelet (sweep). Uses scipy.signal.chirp, adding dimensions as necessary.

    .. plot::
        plt.plot(bruges.filters.klauder(.5, 0.002, [10, 80]))

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): is the sample interval in seconds (usually 0.001, 0.002,
            or 0.004)
        f (array-like): Upper and lower frequencies. Any sequence like (f1, f2).
            A list of lists will create a wavelet bank.
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
        ndarray: The waveform. If you passed `return_t=True` then a tuple of
            (wavelet, t) is returned.
    """
    if not return_t:
        warnings.warn("In future releases, return_t will be True by default.", FutureWarning)

    t0, t1 = -duration/2, duration/2
    t = np.arange(t0, t1, dt)

    f = np.asanyarray(f).reshape(-1, 1)
    f1, f2 = f

    c = [scipy.signal.chirp(t, f1_+(f2_-f1_)/2., t1, f2_, **kwargs)
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


sweep = klauder


def ormsby(duration, dt, f, return_t=False):
    """
    The Ormsby wavelet requires four frequencies which together define a
    trapezoid shape in the spectrum. The Ormsby wavelet has several sidelobes,
    unlike Ricker wavelets.

    .. plot::
        plt.plot(bruges.filters.ormsby(.5, 0.002, [5, 10, 40, 80]))

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (usually 0.001, 0.002,
            or 0.004).
        f (array-like): Sequence of form (f1, f2, f3, f4), or list of lists of
            frequencies, which will return a 2D wavelet bank.

    Returns:
        ndarray: A vector containing the Ormsby wavelet, or a bank of them. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.

    """
    if not return_t:
        warnings.warn("In future releases, return_t will be True by default.", FutureWarning)

    f = np.asanyarray(f).reshape(-1, 1)

    try:
        f1, f2, f3, f4 = f
    except ValueError:
        raise ValueError("The last dimension of the frequency array must be 4")

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


def berlage(duration, dt, f, n=2, alpha=180, phi=-np.pi/2, return_t=False):
    """
    Generates a Berlage wavelet with a peak frequency f. Implements

    .. math::

    w(t) = AH(t) t^n \mathrm{e}^{-\alpha t} \cos(2 \pi f_0 t + \phi_0)

    as described in Aldridge, DF (1990), The Berlage wavelet, GEOPHYSICS
    55 (11), p 1508-1511. Berlage wavelets are causal, minimum phase and
    useful for modeling marine airgun sources.

    If you pass a 1D array of frequencies, you get a wavelet bank in return.

    .. plot::
        plt.plot(bruges.filters.berlage(0.5, 0.002, 40))

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (array-like): Centre frequency of the wavelet in Hz. If a sequence is
            passed, you will get a 2D array in return, one row per frequency.
        n (float): The time exponent; non-negative and real.
        alpha(float): The exponential decay factor; non-negative and real.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis, where time is the range from -duration/2 to
            duration/2 in steps of dt.

    Returns:
        ndarray. Berlage wavelet(s) with centre frequency f sampled on t. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    if not return_t:
        warnings.warn("In future releases, return_t will be True by default.", FutureWarning)

    f = np.asanyarray(f).reshape(-1, 1)
    t = np.arange(-duration/2, duration/2, dt)

    H = np.heaviside(t, 0)
    w = H * t**n * np.exp(-alpha * t) * np.cos(2 * np.pi * f * t + phi)

    w = np.squeeze(w) / np.max(np.abs(w))

    if return_t:
        BerlageWavelet = namedtuple('BerlageWavelet', ['amplitude', 'time'])
        return BerlageWavelet(w, t)
    else:
        return w


def generalized(duration, dt, f, u=2, return_t=False, center=True, imag=False):
    """
    Wang's generalized wavelet, of which the Ricker is a special case where
    u = 2. The parameter u is the order of the time-domain derivative, which
    can be a fractional derivative.

    As given by Wang (2015), Generalized seismic wavelets. GJI 203, p 1172-78.
    DOI: https://doi.org/10.1093/gji/ggv346. I am using the (more accurate)
    frequency domain method (eq 4 in that paper).

    .. plot::
        plt.plot(bruges.filters.generalized(.5, 0.002, 40, u=1.0))

    Args:
        duration (float): The length of the wavelet, in s.
        dt (float): The time sample interval in s.
        f (float or array-like): The frequency or frequencies, in Hertz.
        u (float or array-like): The fractional derivative parameter u.
        return_t (bool): Whether to return the time basis array.
        center (bool): Whether to center the wavelet on time 0.
        imag (bool): Whether to return the imaginary component as well.

    Returns:
        ndarray. If f and u are floats, the resulting wavelet has duration/dt
            = A samples. If you give f as an array of length M and u as an
            array of length N, then the resulting wavelet bank will have shape
            (M, N, A). If f or u are floats, their size will be 1, and they
            will be squeezed out: the bank is always squeezed to its minimum
            number of dimensions. If you passed `return_t=True` then a tuple
            of (wavelet, t) is returned.
    """
    if not return_t:
        warnings.warn("In future releases, return_t will be True by default.", FutureWarning)

    # Make sure we can do banks.
    f = np.asanyarray(f).reshape(-1, 1)
    u = np.asanyarray(u).reshape(-1, 1, 1)

    # Basics.
    om0 = f * 2 * np.pi
    u2 = u / 2
    df = 1 / duration
    nyquist = (1 / dt) / 2
    nf = 1 + nyquist / df
    t0 = duration / 2
    om = 2 * np.pi * np.arange(0, nyquist, df)

    # Compute the spectrum from Wang's eq 4.
    exp1 = np.exp((-om**2 / om0**2) + u2)
    exp2 = np.exp(-1j*om*t0 + 1j*np.pi * (1 + u2))
    W = (u2**(-u2)) * (om**u / om0**u) * exp1 * exp2

    # Compute time domain response.
    if center:
        t = np.arange(-duration/2, duration/2, dt)
    else:
        t = np.arange(0, duration, dt)
    w = np.fft.ifft(W, t.size)
    if not imag:
        w = w.real

    # At this point the wavelet bank has the shape (u, f, a),
    # where u is the size of u, f is the size of f, and a is
    # the number of amplitude samples we generated.
    w_max = np.max(np.abs(w), axis=-1)[:, :, None]
    w = np.squeeze(w / w_max)

    if return_t:
        GeneralizedWavelet = namedtuple('GeneralizedWavelet', ['amplitude', 'time'])
        return GeneralizedWavelet(w, t)
    else:
        return w


def rotate_phase(w, phi, degrees=False):
    """
    Performs a phase rotation of wavelet or wavelet bank using:

    ..math::

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
    if degrees:
        phi = phi * np.pi / 180.0
    a = scipy.signal.hilbert(w, axis=0)
    w = np.real(a) * np.cos(phi)  -  np.imag(a) * np.sin(phi)
    return w
