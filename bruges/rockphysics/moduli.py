# -*- coding: utf-8 -*-
'''
===================
moduli.py
===================

Converts between various acoustic/eslatic parameters,
and provides a way to calculate all the elastic moduli
from Vp, Vs, and rho

Created June 2014

@author: Matt Hall

Using equations http://www.subsurfwiki.org/wiki/Elastic_modulus
from Mavko, G, T Mukerji and J Dvorkin (2003), The Rock Physics Handbook,
Cambridge University Press

'''
import numpy as np


def youngs(vp=None, vs=None, rho=None, mu=None, lam=None, bulk=None, pr=None,
           pmod=None):
    '''
    Computes Young's modulus given either Vp, Vs, and rho, or
    any two elastic moduli (e.g. lambda and mu, or bulk and P
    moduli).

    SI units only.

    Args:
        vp, vs, and rho
        or any 2 from lam, mu, bulk, pr, and pmod

    Returns:
        Young's modulus in pascals, Pa

    '''
    if None not in (vp, vs, rho):
        return rho * vs**2 * (3.*vp**2 - 4.*vs**2) / (vp**2 - vs**2)

    elif None not in (mu, lam):
        return mu * (3.*lam + 2*mu) / (lam + mu)

    elif None not in (bulk, lam):
        return 9.*bulk * (bulk - lam) / (3.*bulk - lam)

    elif None not in (bulk, mu):
        return 9.*bulk*mu / (3.*bulk + mu)

    elif None not in (lam, pr):
        return lam * (1+pr) * (1 - 2*pr) / pr

    elif None not in (pr, mu):
        return 2. * mu * (1+pr)

    elif None not in (pr, bulk):
        return 3. * bulk * (1 - 2*pr)

    else:
        return None


def bulk(vp=None, vs=None, rho=None, mu=None, lam=None, youngs=None, pr=None,
         pmod=None):
    '''
    Computes bulk modulus given either Vp, Vs, and rho, or
    any two elastic moduli (e.g. lambda and mu, or Young's
    and P moduli).

    SI units only.

    Args:
        vp, vs, and rho
        or any 2 from lam, mu, youngs, pr, and pmod

    Returns:
        Bulk modulus in pascals, Pa

    '''

    if None not in (vp, vs, rho):
        return rho * (vp**2 - (4./3.)*(vs**2))

    elif None not in (mu, lam):
        return lam + 2*mu/3.

    elif None not in (mu, youngs):
        return youngs * mu / (9.*mu - 3.*youngs)

    elif None not in (lam, pr):
        return lam * (1+pr) / 3.*pr

    elif None not in (pr, mu):
        return 2. * mu * (1+pr) / (3. - 6.*pr)

    elif None not in (pr, youngs):
        return youngs / (3. - 6.*pr)

    else:
        return None


def pr(vp=None, vs=None, rho=None, mu=None, lam=None, youngs=None, bulk=None,
       pmod=None):
    '''
    Computes Poisson ratio given either Vp, Vs, and rho, or
    any two elastic moduli (e.g. lambda and mu, or Young's
    and P moduli).

    SI units only.

    Args:
        vp, vs, and rho
        or any 2 from lam, mu, youngs, bulk, and pmod

    Returns:
        Poisson's ratio, dimensionless

    '''

    if None not in (vp, vs):
        return (vp**2. - 2.*vs**2) / (2. * (vp**2 - vs**2))

    elif None not in (mu, lam):
        return lam / (2. * (lam+mu))

    elif None not in (mu, youngs):
        return (youngs / (2.*mu)) - 1

    elif None not in (lam, bulk):
        return lam / (3.*bulk - lam)

    elif None not in (bulk, mu):
        return (3.*bulk - 2*mu) / (6.*bulk + 2*mu)

    elif None not in (bulk, youngs):
        return (3.*bulk - youngs) / (6.*bulk)

    else:
        return None


def mu(vp=None, vs=None, rho=None, pr=None, lam=None, youngs=None, bulk=None,
       pmod=None):
    '''
    Computes shear modulus given either Vp, Vs, and rho, or
    any two elastic moduli (e.g. lambda and bulk, or Young's
    and P moduli).

    SI units only.

    Args:
        vp, vs, and rho
        or any 2 from lam, bulk, youngs, pr, and pmod

    Returns:
        Shear modulus in pascals, Pa

    '''

    if None not in (vs, rho):
        return rho * vs**2

    elif None not in (bulk, lam):
        return 3. * (bulk - lam) / 2.

    elif None not in (bulk, youngs):
        return 3. * bulk * youngs / (9.*bulk - youngs)

    elif None not in (lam, pr):
        return lam * (1 - 2.*pr) / (2.*pr)

    elif None not in (pr, youngs):
        return youngs / (2. * (1 + pr))

    elif None not in (pr, bulk):
        return 3. * bulk * (1 - 2*pr) / (2. * (1 + pr))

    else:
        return None


def lam(vp=None, vs=None, rho=None, pr=None,  mu=None, youngs=None, bulk=None,
        pmod=None):
    '''
    Computes lambda given either Vp, Vs, and rho, or
    any two elastic moduli (e.g. bulk and mu, or Young's
    and P moduli).

    SI units only.

    Args:
        vp, vs, and rho
        or any 2 from bulk, mu, youngs, pr, and pmod

    Returns:
        Lambda in pascals, Pa

    '''
    if None not in (vp, vs, rho):
        return rho * (vp**2 - 2.*vs**2.)

    elif None not in (youngs, mu):
        return mu * (youngs - 2.*mu) / (3.*mu - youngs)

    elif None not in (bulk, mu):
        return bulk - (2.*mu/3.)

    elif None not in (bulk, youngs):
        return 3. * bulk * (3*bulk - youngs) / (9*bulk - youngs)

    elif None not in (pr, mu):
        return 2. * pr * mu / (1 - 2.*pr)

    elif None not in (pr, youngs):
        return pr * youngs / ((1+pr) * (1-2*pr))

    elif None not in (pr, bulk):
        return 3. * bulk * pr / (1+pr)

    else:
        return None


def pmod(vp=None, vs=None, rho=None, pr=None, mu=None, lam=None, youngs=None,
         bulk=None):
    '''
    Computes P-wave modulus given either Vp, Vs, and rho, or
    any two elastic moduli (e.g. lambda and mu, or Young's
    and bulk moduli).

    SI units only.

    Args:
        vp, vs, and rho
        or any 2 from lam, mu, youngs, pr, and bulk

    Returns:
        P-wave modulus in pascals, Pa

    '''

    if None not in (vp, rho):
        return rho * vp**2

    elif None not in (lam, mu):
        return lam + 2*mu

    elif None not in (youngs, mu):
        return mu * (4.*mu - youngs) / (3.*mu - youngs)

    elif None not in (bulk, lam):
        return 3*bulk - 2.*lam

    elif None not in (bulk, mu):
        return bulk + (4.*mu/3.)

    elif None not in (bulk, youngs):
        return 3. * bulk * (3*bulk + youngs) / (9*bulk - youngs)

    elif None not in (lam, pr):
        return lam * (1 - pr) / pr

    elif None not in (pr, mu):
        return 2. * pr * mu * (1-pr) / (1 - 2.*pr)

    elif None not in (pr, youngs):
        return (1-pr) * youngs / ((1+pr) * (1 - 2.*pr))

    elif None not in (pr, bulk):
        return 3. * bulk * (1-pr) / (1+pr)

    else:
        return None


def vp(youngs=None, vs=None, rho=None, mu=None, lam=None, bulk=None, pr=None,
       pmod=None):
    '''
    Computes Vp given bulk density and any two elastic moduli
    (e.g. lambda and mu, or Young's and P moduli).

    SI units only.

    Args:
        Any 2 from lam, mu, youngs, pr, pmod, bulk
        Rho

    Returns:
        Vp in m/s

    '''

    if None not in (mu, lam, rho):
        return np.sqrt((lam + 2.*mu) / rho)

    elif None not in (youngs, mu and rho):
        return np.sqrt(mu * (youngs - 4.*mu) / (rho * (youngs - 3.*mu)))

    elif None not in (youngs, pr and rho):
        return np.sqrt(youngs * (1 - pr) / (rho * (1+pr) * (1 - 2.*pr)))

    elif None not in (bulk, lam and rho):
        return np.sqrt((9.*bulk - 2.*lam) / rho)

    elif None not in (bulk, mu and rho):
        return np.sqrt((bulk + 4.*mu/3.) / rho)

    elif None not in (lam, pr and rho):
        return np.sqrt(lam * (1. - pr) / (pr*rho))

    else:
        return None


def vs(youngs=None, vp=None, rho=None, mu=None, lam=None, bulk=None, pr=None,
       pmod=None):
    '''
    Computes Vs given bulk density and shear modulus.

    SI units only.

    Args:
        Mu
        Rho

    Returns:
        Vs in m/s

    '''

    if None not in (mu, rho):
        return np.sqrt(mu / rho)

    else:
        return None


def moduli_dict(vp, vs, rho):
    '''
    Computes elastic moduli given Vp, Vs, and rho.

    SI units only.

    Args:
        Vp, Vs, and rho

    Returns:
        A dict of elastic moduli, plus P-wave impedance.

    '''

    mod = {}

    mod['imp'] = vp * rho

    mod['mu'] = mu(vs=vs, rho=rho)
    mod['pr'] = pr(vp=vp, vs=vs, rho=rho)
    mod['lam'] = lam(vp=vp, vs=vs, rho=rho)
    mod['bulk'] = bulk(vp=vp, vs=vs, rho=rho)
    mod['pmod'] = pmod(vp=vp, rho=rho)
    mod['youngs'] = youngs(vp=vp, vs=vs, rho=rho)

    return mod
