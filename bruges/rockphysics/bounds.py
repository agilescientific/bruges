#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Bounds on effective elastic moduli.
:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""

import numpy as np


def voight_reuss(f, k, mu):
    """
    Voight-Ruess bound of aggregate

    :param f: list or array of volume fractions (<=1).
    :param k: bulk moduli of constituents (list or array).
    :param mu: shear moduli of constituents (list or array).

    :returns: ku, kl, muu, mul - upper and lower bounds of the bulk
    and shear moduli, plus ka, mul - arithmetic average of the upper 
    and lower bounds  (equals the Hill average for Hashin-Shtrikman bounds)

    :source: Berryman, J.G., 1993, Mixture theories for rock properties
             Mavko, G., 1993, Rock Physics Formulas.

    : Written originally by Xingzhou 'Frank' Liu, in MATLAB 
    : modified by Isao Takahashi, 4/27/99,
    : Translated into Python by Evan Bianco
    """

    if len(f) != len(k) or len(f) != len(mu) or len(k) != len(mu):
        print ('Input f, k, and mu must have the same length')

    if sum(f) != 1.0:
        print ('F must sum up to 1')

    ku = np.sum(f * k)          # Voight bound
    kl = 1.0 / np.sum(f / k)    # Reuss bound

    muu = np.sum(f * mu)        # Voight bound
    mul = 1.0 / np.sum(f / mu)  # Reuss bound

    ka = (ku + kl) / 2.0        # Hill averages
    mua = (muu + mul) / 2.0 

    return ku, kl, muu, mul, ka, mua


def hashin_shtrikman(f, k, mu):

    """
    HashinShtrikman bounds of aggregate

    :param f: list or array of volume fractions (<=1).
    :param k: bulk moduli of constituents (list or array).
    :param mu: shear moduli of constituents (list or array).

    :returns: ku, kl, muu, mul - upper and lower bounds of the bulk
    and shear moduli, plus ka, mul - arithmetic average of the upper
    and lower bounds  (equals the Hill average for Hashin-Shtrikman bounds)

    :source: Berryman, J.G., 1993, Mixture theories for rock properties
             Mavko, G., 1993, Rock Physics Formulas.

    : Written originally by Xingzhou 'Frank' Liu, in MATLAB
    : modified by Isao Takahashi, 4/27/99,
    : Translated into Python by Evan Bianco
    """
    if len(f) != len(k) or len(f) != len(mu) or len(k) != len(mu):
        print ('Input f, k, and mu must have the same length')

    if sum(f) != 1.0:
        print ('F must sum up to 1')

    kmx = max(k)
    kmn = min(k)
    umx = max(k)
    umn = min(k)

    ku = 1.0 / np.sum(f / (k + c * umx)) - (c * umx)  # HS upper bound
    kl = 1.0 / np.sum(f / (k + c * umn)) - (c * umx)  # HS lower bound

    etamx = umx * (9 * kmx + 8 * umx) / (kmx + 2 * umx) / 6
    etamn = umn * (9 * kmn + 8 * umn) / (kmn + 2 * umn) / 6

    muu = 1.0 / sum(f / (mu + etamx)) - etamx    # HS upper bound
    mul = 1.0 / sum(f / (mu + etamn)) - etamn    # HS lower bound

    ka = (ku + kl) / 2.0        # Hill averages
    mua = (muu + mul) / 2.0

    return ku, kl, muu, mul, ka, mua
