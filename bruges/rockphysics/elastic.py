# -*- coding: utf-8 -*-
"""
Elastic impedance.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np


def elastic_impedance(vp, vs, rho, theta1,
                      k='auto',
                      normalize=False,
                      constants=None,
                      use_sin=False,
                      rho_term=False):
    """
    Returns the elastic impedance as defined by Connolly, 1999;
    we are using the formulation reported in Whitcombe et al. (2001).
    Inputs should be shape m x 1, angles should be n x 1. The
    result will be m x n.

    :param vp1: The P-wave velocity scalar or 1D array.
    :param vs1: The S-wave velocity scalar or 1D array.
    :param rho1: The bulk density scalar or 1D array.
    :param theta1: Incident angle(s), scalar or array [degrees].
    :param k: A constant, see Connolly (1999). Default is 'auto', which
        computes it from Vp and Vs
    :param normalize: if True, returns the normalized form of Whitcombe
        et. al (2001).
    :param constants: A sequence of 3 constants to use for normalization. If
        you don't provide this, then normalization just uses the means of the
        inputs. If these are scalars, the result will be the acoustic impedance
        (see Whitcombe, 2001).
    :param use_sin: Use sin^2 for the first term, instead of tan^2 (see
        Connolly 1999).
    :param rho_term: Alternative form, with Vp factored out; use in place of
        density in generating synthetics in other software (see Connolly 1999).
    """
    theta1 = np.array(np.radians(theta1)).reshape((-1, 1))
    if (np.nan_to_num(theta1) > np.pi/2.).any():
        raise ValueError("Incidence angle theta1 must be less than 90 deg.")

    alpha = np.array(vp, dtype=float)
    beta = np.array(vs, dtype=float)
    rho = np.array(rho, dtype=float)

    if use_sin:
        op = np.sin
    else:
        op = np.tan

    if k.lower() == 'auto':
        k = np.mean(beta**2.0 / alpha**2.0)
    # Otherwise, just use the k we were given.

    a = 1 + op(theta1)**2.0
    b = -8 * k * np.sin(theta1)**2.0
    c = 1 - 4 * k * np.sin(theta1)**2.0

    ei = alpha**a * beta**b * rho**c

    if normalize:
        if constants is None:
            # Use the means; this will return acoustic impedance for scalars.
            alpha_0 = np.nanmean(vp)
            beta_0 = np.nanmean(vs)
            rho_0 = np.nanmean(rho)
        else:
            try:
                alpha_0, beta_0, rho_0 = constants
            except ValueError as e:
                raise ValueError("You must provide a sequence of 3 constants.")
        ei *= alpha_0**(1 - a) * beta_0**(-b) * rho_0**(1 - c)

    if rho_term:
        ei /= alpha

    if ei.size == 1:
        return np.asscalar(ei)
    else:
        return np.squeeze(ei.T)
