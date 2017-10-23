#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Elastic impedance.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np


def elastic_impedance(vp, vs, rho, theta1, normalize=True):
    """
    Returns the elastic impedance (as defined by Connolly, 1999)
    for each incidence angle in theta1:

    :param vp1: The p-wave velocity or p-wave velocity array.
    :param vs1: The s-wave velocity of s-wave velocity array.
    :param rho1: The density (either scalar or array).
    :param theta1: An array of incident angles to use for reflectivity
        calculation [degrees].
    :param normalized: if True, returns the normalized form from
        Whitcombe et. al (2001).
    """
    theta1 = np.radians(theta1)
    # k = np.mean(vs)**2 / np.mean(vp)**2
    k = 0.25
    a = 1 + np.tan(theta1)**2
    b = -8 * k * np.sin(theta1)**2
    c = 1 - 4 * k * np.sin(theta1)**2

    ei = vp**a * vs**b * rho**c

    if normalize:
        n = vp**(1 - a) * vs**(-b) * rho**(1 - c)
        return n * ei
    else:
        return ei
