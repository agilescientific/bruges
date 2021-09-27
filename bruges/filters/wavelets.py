"""
Seismic wavelets.

:copyright: 2021 Agile Geoscience
:license: Apache 2.0
"""
import warnings
from collections import namedtuple

import numpy as np
import scipy.signal
from bruges.util import deprecated


def _get_time(duration, dt, sym=None):
    """
    Make a time vector.

    If `sym` is `True`, the time vector will have an odd number of samples,
    and will be symmetric about 0. If it's False, and the number of samples
    is even (e.g. duration = 0.016, dt = 0.004), then 0 will bot be center.
    """
    if sym is None:
        m = "In future releases, the default legacy behaviour will be removed. "
        m += "We recommend setting sym=True. This will be the default in v0.5+."
        warnings.warn(m, category=FutureWarning, stacklevel=2)
        return np.arange(-duration/2, duration/2, dt)
    
    # This business is to avoid some of the issues with `np.arange`:
    # (1) unpredictable length and (2) floating point weirdness, like
    # 1.234e-17 instead of 0. Not using `linspace` because figuring out
    # the length and offset gave me even more of a headache than this.
    n = int(duration / dt)
    odd = n % 2
    k = int(10**-np.floor(np.log10(dt)))
    dti = int(k * dt)  # integer dt
        
    if (odd and sym):
        t = np.arange(n)
    if (not odd and sym):
        t = np.arange(n + 1)
    if (odd and not sym): 
        t = np.arange(n)
    if (not odd and not sym):
        t = np.arange(n) - 1
        
    t -= t[-1] // 2
    
    return dti * t / k


def _generic(func, duration, dt, f, t=None, return_t=False, taper='blackman', sym=None):
    """
    Generic wavelet generator: applies a window to a continuous function.

    Args:
        func (function): The continuous function, taking t, f as arguments.
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (often one of  0.001, 0.002,
            or 0.004).
        f (array-like): Dominant frequency of the wavelet in Hz. If a sequence is
            passed, you will get a 2D array in return, one row per frequency.
        t (array-like): The time series to evaluate at, if you don't want one
            to be computed. If you pass `t` then `duration` and `dt` will be
            ignored, so we recommend passing `None` for those arguments.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis.
        taper (str or function): The window or tapering function to apply.
            To use one of NumPy's functions, pass 'bartlett', 'blackman' (the
            default), 'hamming', or 'hanning'; to apply no tapering, pass
            'none'. To apply your own function, pass a function taking only
            the length of the window and returning the window function.
        sym (bool): If True (default behaviour before v0.5) then the wavelet
            is forced to have an odd number of samples and the central sample
            is at 0 time.

    Returns:
        ndarray. wavelet(s) with centre frequency f sampled on t. If you
            passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    if not return_t:
        m = "In future releases, return_t will be True by default."
        warnings.warn(m, FutureWarning, stacklevel=2)

    f = np.asanyarray(f).reshape(-1, 1)

    # Compute time domain response.
    if t is None:
        t = _get_time(duration, dt, sym=sym)
    else:
        if (duration is not None) or (dt is not None):
            m = "`duration` and `dt` are ignored when `t` is passed."
            warnings.warn(m, UserWarning, stacklevel=2)

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


def sinc(duration, dt, f, t=None, return_t=False, taper='blackman', sym=None):
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
        t (array-like): The time series to evaluate at, if you don't want one
            to be computed. If you pass `t` then `duration` and `dt` will be
            ignored, so we recommend passing `None` for those arguments.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis.
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

    return _generic(func, duration, dt, f, t, return_t, taper)


def cosine(duration, dt, f, t=None, return_t=False, taper='gaussian', sigma=None, sym=None):
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
        t (array-like): The time series to evaluate at, if you don't want one
            to be computed. If you pass `t` then `duration` and `dt` will be
            ignored, so we recommend passing `None` for those arguments.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis.
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

    return _generic(func, duration, dt, f, t, return_t, taper)


def gabor(duration, dt, f, t=None, return_t=False, sym=None):
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
        t (array-like): The time series to evaluate at, if you don't want one
            to be computed. If you pass `t` then `duration` and `dt` will be
            ignored, so we recommend passing `None` for those arguments.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis.

    Returns:
        ndarray. Gabor wavelet(s) with centre frequency f sampled on t. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    def func(t_, f_):
        return np.exp(-2 * f_**2 * t_**2) * np.cos(2 * np.pi * f_ * t_)

    return _generic(func, duration, dt, f, t, return_t)


def ricker(duration, dt, f, t=None, return_t=False, sym=None):
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
        t (array-like): The time series to evaluate at, if you don't want one
            to be computed. If you pass `t` then `duration` and `dt` will be
            ignored, so we recommend passing `None` for those arguments.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis.
        sym (bool): If True (default behaviour before v0.5) then the wavelet
            is forced to have an odd number of samples and the central sample
            is at 0 time.

    Returns:
        ndarray. Ricker wavelet(s) with centre frequency f sampled on t. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    if not return_t:
        m = "In future releases, return_t will be True by default."
        warnings.warn(m, FutureWarning, stacklevel=2)

    f = np.asanyarray(f).reshape(-1, 1)

    if t is None:
        t = _get_time(duration, dt, sym=sym)
    else:
        if (duration is not None) or (dt is not None):
            m = "`duration` and `dt` are ignored when `t` is passed."
            warnings.warn(m, UserWarning, stacklevel=2)

    pft2 = (np.pi * f * t)**2
    w = np.squeeze((1 - (2 * pft2)) * np.exp(-pft2))

    if return_t:
        RickerWavelet = namedtuple('RickerWavelet', ['amplitude', 'time'])
        return RickerWavelet(w, t)
    else:
        return w


def klauder(duration, dt, f,
            autocorrelate=True,
            t=None,
            return_t=False,
            taper='blackman',
            sym=None,
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
        t (array-like): The time series to evaluate at, if you don't want one
            to be computed. If you pass `t` then `duration` and `dt` will be
            ignored, so we recommend passing `None` for those arguments.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis.
        taper (str or function): The window or tapering function to apply.
            To use one of NumPy's functions, pass 'bartlett', 'blackman' (the
            default), 'hamming', or 'hanning'; to apply no tapering, pass
            'none'. To apply your own function, pass a function taking only
            the length of the window and returning the window function.
        sym (bool): If True (default behaviour before v0.5) then the wavelet
            is forced to have an odd number of samples and the central sample
            is at 0 time.
        **kwargs: Further arguments are passed to scipy.signal.chirp. They are
            `method` ('linear','quadratic','logarithmic'), `phi` (phase offset
            in degrees), and `vertex_zero`.

    Returns:
        ndarray: The waveform. If you passed `return_t=True` then a tuple of
            (wavelet, t) is returned.
    """
    if not return_t:
        m = "In future releases, return_t will be True by default."
        warnings.warn(m, FutureWarning, stacklevel=2)

    if t is None:
        t = _get_time(duration, dt, sym=sym)
    else:
        if (duration is not None) or (dt is not None):
            m = "`duration` and `dt` are ignored when `t` is passed. "
            m += "Pass None to suppress this warning."
            warnings.warn(m, UserWarning, stacklevel=2)

    t0, t1 = -duration/2, duration/2

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


def ormsby(duration, dt, f, t=None, return_t=False, sym=None):
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
        t (array-like): The time series to evaluate at, if you don't want one
            to be computed. If you pass `t` then `duration` and `dt` will be
            ignored, so we recommend passing `None` for those arguments.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis.
        sym (bool): If True (default behaviour before v0.5) then the wavelet
            is forced to have an odd number of samples and the central sample
            is at 0 time.

    Returns:
        ndarray: A vector containing the Ormsby wavelet, or a bank of them. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.

    """
    if not return_t:
        m = "In future releases, return_t will be True by default."
        warnings.warn(m, FutureWarning, stacklevel=2)

    f = np.asanyarray(f).reshape(-1, 1)

    try:
        f1, f2, f3, f4 = f
    except ValueError:
        raise ValueError("The last dimension of the frequency array must be of size 4.")

    def numerator(f, t):
        return (np.sinc(f * t)**2) * ((np.pi * f) ** 2)

    pf43 = (np.pi * f4) - (np.pi * f3)
    pf21 = (np.pi * f2) - (np.pi * f1)

    if t is None:
        t = _get_time(duration, dt, sym=sym)
    else:
        if (duration is not None) or (dt is not None):
            m = "`duration` and `dt` are ignored when `t` is passed."
            warnings.warn(m, UserWarning, stacklevel=2)

    w = ((numerator(f4, t)/pf43) - (numerator(f3, t)/pf43) -
         (numerator(f2, t)/pf21) + (numerator(f1, t)/pf21))

    w = np.squeeze(w) / np.amax(w)

    if return_t:
        OrmsbyWavelet = namedtuple('OrmsbyWavelet', ['amplitude', 'time'])
        return OrmsbyWavelet(w, t)
    else:
        return w


def ormsby_fft(duration, dt, f, P=(0, 0), return_t=True, sym=True):
    """
    Non-white Ormsby, with arbitary amplitudes.
    
    Can use as many points as you like. The power of f1 and f4 is assumed to be 0,
    so you only need to provide p2 and p3 (the corners). (You can actually provide
    as many f points as you like, as long as there are n - 2 matching p points.)

    .. plot::
        plt.plot(bruges.filters.ormsby(.5, 0.002, [5, 10, 40, 80]))

    Args:
        duration (float): The length in seconds of the wavelet.
        dt (float): The sample interval in seconds (usually 0.001, 0.002,
            or 0.004).
        f (array-like): Sequence of form (f1, f2, f3, f4), or list of lists of
            frequencies, which will return a 2D wavelet bank.
        P (tuple): The power of the f2 and f3 frequencies, in relative dB.
            (The magnitudes of f1 and f4 are assumed to be -∞ dB, i.e. a
            magnitude of 0.) The default power values of (0, 0) results in a
            trapezoidal spectrum and a conventional Ormsby wavelet. Pass, e.g.
            (0, -15) for a 'pink' wavelet, with more energy in the lower
            frequencies.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis.
        sym (bool): If True (default behaviour before v0.5) then the wavelet
            is forced to have an odd number of samples and the central sample
            is at 0 time.

    Returns:
        ndarray: A vector containing the Ormsby wavelet, or a bank of them. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    fs = 1 / dt
    fN = fs // 2
    n = int(duration / dt)
    a = map(lambda p: 10**(p/20), P)

    # Linear interpolation of points.
    x  = np.linspace(0, int(fN), int(10*n))
    xp = [  0.] + list(f) +  [fN]
    fp = [0., 0.] + list(a) + [0., 0.]
    W = np.interp(x, xp, fp)

    # Compute inverse FFT.
    w_ = np.fft.fftshift(np.fft.irfft(W))
    L = int(w_.size // 2)
    normalize = lambda d: d / np.max(abs(d))
    w = normalize(w_[L-n//2:L+n//2+sym])
    t = _get_time(duration, dt, sym=sym)

    if return_t:
        OrmsbyWavelet = namedtuple('OrmsbyWavelet', ['amplitude', 'time'])
        return OrmsbyWavelet(w, t)
    else:
        return w


def berlage(duration, dt, f, n=2, alpha=180, phi=-np.pi/2, t=None, return_t=False, sym=None):
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
        phi (float): The phase.
        t (array-like): The time series to evaluate at, if you don't want one
            to be computed. If you pass `t` then `duration` and `dt` will be
            ignored, so we recommend passing `None` for those arguments.
        return_t (bool): If True, then the function returns a tuple of
            wavelet, time-basis.
        sym (bool): If True (default behaviour before v0.5) then the wavelet
            is forced to have an odd number of samples and the central sample
            is at 0 time.

    Returns:
        ndarray. Berlage wavelet(s) with centre frequency f sampled on t. If
            you passed `return_t=True` then a tuple of (wavelet, t) is returned.
    """
    if not return_t:
        m = "In future releases, return_t will be True by default."
        warnings.warn(m, FutureWarning, stacklevel=2)

    f = np.asanyarray(f).reshape(-1, 1)
    if t is None:
        t = _get_time(duration, dt, sym=sym)
    else:
        if (duration is not None) or (dt is not None):
            m = "`duration` and `dt` are ignored when `t` is passed."
            warnings.warn(m, UserWarning, stacklevel=2)


    H = np.heaviside(t, 0)
    w = H * t**n * np.exp(-alpha * t) * np.cos(2 * np.pi * f * t + phi)

    w = np.squeeze(w) / np.max(np.abs(w))

    if return_t:
        BerlageWavelet = namedtuple('BerlageWavelet', ['amplitude', 'time'])
        return BerlageWavelet(w, t)
    else:
        return w


def generalized(duration, dt, f, u=2, t=None, return_t=False, imag=False, sym=None):
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
        t (array-like): The time series to evaluate at, if you don't want one
            to be computed. If you pass `t` then `duration` and `dt` will be
            ignored, so we recommend passing `None` for those arguments.
        return_t (bool): Whether to return the time basis array.
        center (bool): Whether to center the wavelet on time 0.
        imag (bool): Whether to return the imaginary component as well.
        sym (bool): If True (default behaviour before v0.5) then the wavelet
            is forced to have an odd number of samples and the central sample
            is at 0 time.

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
        m = "In future releases, return_t will be True by default."
        warnings.warn(m, FutureWarning, stacklevel=2)

    # Make sure we can do banks.
    f = np.asanyarray(f).reshape(-1, 1)
    u = np.asanyarray(u).reshape(-1, 1, 1)

    # Compute time domain response.
    if t is None:
        t = _get_time(duration, dt, sym=sym)
    else:
        if (duration is not None) or (dt is not None):
            m = "`duration` and `dt` are ignored when `t` is passed."
            warnings.warn(m, UserWarning, stacklevel=2)
        dt = t[1] - t[0]
        duration = len(t) * dt

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


@deprecated('bruges.filters.wavelets.rotate_phase() is deprecated. Please use bruges.filters.rotate_phase() instead.')
def rotate_phase(w, phi, degrees=False):
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
    if degrees:
        phi = phi * np.pi / 180.0
    a = scipy.signal.hilbert(w, axis=0)
    w = np.real(a) * np.cos(phi)  -  np.imag(a) * np.sin(phi)
    return w
