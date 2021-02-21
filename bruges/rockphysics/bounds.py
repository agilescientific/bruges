"""
Bounds on effective elastic modulus.
:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
from collections import namedtuple

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
    f = np.array(f).astype(float)

    if float(sum(f)) == 100.0:
        # fractions have been given in percent: scale to 1.
        f /= 100.0

    m = np.array(m)
    mv = np.sum(f * m)

    return mv


def reuss_bound(f, m):
    """
    The lower bound on the effective elastic modulus of a
    mixture of N material phases. This is defined at the harmonic
    average of the constituents. Same as Wood's equation for homogeneous mixed fluids.

    Args:
        f: list or array of N volume fractions (must sum to 1 or 100).
        m: elastic modulus of N constituents (list or array).

     Returns:
        mr: Reuss lower bound.
    """
    f = np.array(f).astype(float)

    if float(sum(f)) == 100.0:
        # fractions have been given in percent: scale to 1.
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


def hashin_shtrikman(f, k, mu, modulus='bulk'):
    """
    Hashin-Shtrikman bounds for a mixture of two constituents.
    The best bounds for an isotropic elastic mixture, which give
    the narrowest possible range of elastic modulus without
    specifying anything about the geometries of the constituents.

    Args:
        f: list or array of volume fractions (must sum to 1.00 or 100%).
        k: bulk modulus of constituents (list or array).
        mu: shear modulus of constituents (list or array).
        modulus: A string specifying whether to return either the
            'bulk' or 'shear' HS bound.

    Returns:
        namedtuple: The Hashin Shtrikman (lower, upper) bounds.

    :source: Berryman, J.G., 1993, Mixture theories for rock properties
             Mavko, G., 1993, Rock Physics Formulas.

    : Written originally by Xingzhou 'Frank' Liu, in MATLAB
    : modified by Isao Takahashi, 4/27/99,
    : Translated into Python by Evan Bianco
    """
    def z_bulk(k, mu):
        return (4/3.) * mu

    def z_shear(k, mu):
        return mu * (9 * k + 8 * mu) / (k + 2 * mu) / 6

    def bound(f, k, z):
        return 1 / sum(f / (k + z)) - z

    f = np.array(f)
    if sum(f) == 100:
        f /= 100.0

    func = {'shear': z_shear,
            'bulk': z_bulk}

    k, mu = np.array(k), np.array(mu)
    z_min = func[modulus](np.amin(k), np.amin(mu))
    z_max = func[modulus](np.amax(k), np.amax(mu))

    fields = ['lower_bound', 'upper_bound']
    HashinShtrikman = namedtuple('HashinShtrikman', fields)
    return HashinShtrikman(bound(f, k, z_min), bound(f, k, z_max))
