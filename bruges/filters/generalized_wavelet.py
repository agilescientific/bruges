# -*- coding: utf-8 -*-
import numpy as np


def generalized_wavelet(duration, dt, f, u=2, return_time=False):
    """
    Generalized wavelet, of which the Ricker is a special case where u = 2.
    The parameter u is the order of the time-domain derivative.

    Use the function get_derivative() to find the best value for f and u.

    As given by Wang (2015), Generalized seismic wavelets. GJI 203, p 1172-78.
    Using the more accurate frequency domain method (eq 4 in that paper).
    See in-bruges for paper and dev.
    """
    om0 = f * 2 * np.pi
    u2 = u / 2
    df = 1 / duration
    nyquist = (1 / dt) / 2
    nf = 1 + nyquist / df
    t0 = duration / 2
    om = np.linspace(0, nyquist, nf)

    # METHOD 1
    exp1 = np.exp((-om**2 / om0**2) + u2)
    exp2 = np.exp(-1j*om*t0 + 1j*np.pi * (1 + u2))
    W = (u2**(-u2)) * (om**u / om0**u) * exp1 * exp2
    W = np.abs(W)

    # METHOD 2
    G = -1 * np.exp(-om**2/om0**2) * np.exp(-1j*om*t0)
    G_ = G * (1j * om)**u
    PHI = G_ * om0**(-u) * (u2**(-u2)) * np.exp(u2)
    W = np.abs(PHI)

    # COMPUTE TIME DOMAIN RESPONSE
    # W_ = np.fft.ifftshift(W)
    w = np.fft.irfft(W)
    w = np.fft.ifftshift(w)
    t = np.linspace(-duration/2, duration/2, duration/dt)

    if return_time:
        return w, t
    return W, w
