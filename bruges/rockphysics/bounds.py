#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Bounds on effective elastic moduli.
:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""

import numpy as np


def voigt_bound(f, m):
    """
    The upper bound on the effective elastic modulus, mv of a
    mixture of N material phases. This is defined at the arithmetic
    average of the constituents.

    Args:
        f: list or array of N volume fractions (must sum to 1 or 100).
        m: elastic modulus of N constituents (list or array).

     Returns:
        mv: Voigt upper bound.

    """

    f = np.array(f)

    if float(sum(f)) == 100.0:
        # fractions have been giving in percent: scale to 1.
        f /= 100.0

    m = np.array(m)
    mv = np.sum(f * m)

    return mv


def reuss_bound(f, m):
    """
    The lower bound on the effective elastic modulus, mv of a
    mixture of N material phases. This is defined at the harmonic
    average of the constituents.

    Args:
        f: list or array of N volume fractions (must sum to 1 or 100).
        m: elastic modulus of N constituents (list or array).

     Returns:
        mr: Reuss lower bound.

    """

    f = np.array(f)

    if float(sum(f)) == 100.0:
        # fractions have been giving in percent: scale to 1.
        f /= 100.0

    m = np.array(m)
    mr = 1.0 / np.sum(f / m)

    return mr


def hill_average(f, m):
    """
    The Hill average effective elastic modulus, mh of a
    mixture of N material phases. This is defined as the simple
    average of the Reuss (lower) and Voigt (upper) bounds.

    Args:
        f: list or array of N volume fractions (must sum to 1 or 100).
        m: elastic modulus of N constituents (list or array).

     Returns:
        mh: Hill average.

    """
    mv = voigt_bound(f, m)
    mr = reuss_bound(f, m)
    mh = (mv + mr) / 2.0

    return mh


def hashin_shtrikman(f, k, mu, bound='upper', moduli='bulk'):

    """
    Hashin-Shtrikman bounds for a mixture of two constituents.
    The best bounds for an isotropic elastic mixture, which give
    the narrowest possible range of elastic moduli without
    specifying anything about the geometries of the constituents.

    Args:
        f: list or array of volume fractions (<=1).
        k: bulk moduli of constituents (list or array).
        mu: shear moduli of constituents (list or array).
        moduli: A string specifying whether to return either the
            'bulk' or 'shear' HS bound.
        bound: A string specifying whether to return either
            the 'upper' or lower 'bound'.

    Returns:
        The Hashin Shtrikman modulus, mhs

    :source: Berryman, J.G., 1993, Mixture theories for rock properties
             Mavko, G., 1993, Rock Physics Formulas.

    : Written originally by Xingzhou 'Frank' Liu, in MATLAB
    : modified by Isao Takahashi, 4/27/99,
    : Translated into Python by Evan Bianco
    """

    f = np.array(f)

    if float(sum(f)) == 100.0:
        # fractions have been giving in percent: scale to 1.
        f /= 100.0

    k = np.array(k)
    mu = np.array(mu)

    c = 4.0 / 3

    kmx = max(k)
    kmn = min(k)
    umx = max(mu)
    umn = min(mu)

    if moduli == 'bulk':
        if bound == 'upper':
            mhs = 1.0 / np.sum(f / (k + c * umx)) - (c * umx)  # HS upper bound
        if bound == 'lower':
            mhs = 1.0 / np.sum(f / (k + c * umn)) - (c * umx)  # HS lower bound

    if moduli == 'shear':
        if bound == 'upper':
            etamx = umx * (9.0 * kmx + 8.0 * umx) / (kmx + 2.0 * umx) / 6.0
            mhs = 1.0 / sum(f / (mu + etamx)) - etamx    # HS upper bound

        if bound == 'lower':
            etamn = umn * (9.0 * kmn + 8.0 * umn) / (kmn + 2.0 * umn) / 6.0
            mhs = 1.0 / sum(f / (mu + etamn)) - etamn    # HS lower bound

    return mhs
