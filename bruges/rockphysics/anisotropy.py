#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Anisotropy effects.

Backus anisotropy is from thin layers.

Hudson anisotropy is from crack defects.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np

from bruges.rockphysics import moduli
from bruges.util import moving_average


def backus_parameters(vp, vs, rho, lb, dz):
    """
    Intermediate parameters for Backus averaging. This is expected to be a
    private function. You probably want backus() and not this.

    Args:
        vp (ndarray): P-wave interval velocity.
        vs (ndarray): S-wave interval velocity.
        rho (ndarray): Bulk density.
        lb (float): The Backus averaging length in m.
        dz (float): The depth sample interval in m.

    Returns:
        tuple: Liner's 5 intermediate parameters: A, C, F, L and M.

    Notes:
        Liner, C (2014), Long-wave elastic attenuation produced by horizontal
        layering. The Leading Edge, June 2014, p 634-638.

    """
    lam = moduli.lam(vp, vs, rho)
    mu = moduli.mu(vp, vs, rho)

    # Compute the layer parameters from Liner (2014) equation 2:
    a = rho * np.power(vp, 2.0)  # Acoustic impedance

    # Compute the Backus parameters from Liner (2014) equation 4:
    A1 = 4 * moving_average(mu*(lam+mu)/a, lb/dz, mode='same')
    A = A1 + np.power(moving_average(lam/a, lb/dz, mode='same'), 2.0)\
        / moving_average(1.0/a, lb/dz, mode='same')
    C = 1.0 / moving_average(1.0/a, lb/dz, mode='same')
    F = moving_average(lam/a, lb/dz, mode='same')\
        / moving_average(1.0/a, lb/dz, mode='same')
    L = 1.0 / moving_average(1.0/mu, lb/dz, mode='same')
    M = moving_average(mu, lb/dz, mode='same')

    return A, C, F, L, M


def backus(vp, vs, rho, lb, dz, return_rho=False):
    """
    Backus averaging. Using Liner's algorithm (2014; see Notes).

    Args:
        vp (ndarray): P-wave interval velocity.
        vs (ndarray): S-wave interval velocity.
        rho (ndarray): Bulk density.
        lb (float): The Backus averaging length in m.
        dz (float): The depth sample interval in m.
        return_rho (bool): Whether to return the smoothed rho log too.

    Returns:
        tuple: the smoothed logs. Default: vp and vs. Optional (see flag):
            vp, vs, plus rho. Useful for computing other elastic parameters
            at a seismic scale.

    Notes:
        Liner, C (2014), Long-wave elastic attenuation produced by horizontal
        layering. The Leading Edge, June 2014, p 634-638.

    """
    # Compute the Backus parameters:
    A, C, F, L, M = backus_parameters(vp, vs, rho, lb, dz)

    # Compute the vertical velocities from Liner (2014) equation 5:
    R = moving_average(rho, lb/dz, mode='same')
    vp0 = np.sqrt(C / R)
    vs0 = np.sqrt(L / R)

    if return_rho:
        return vp0, vs0, R
    else:
        return vp0, vs0


def backus_quality_factor(vp, vs, rho, lb, dz):
    """
    Compute Qp and Qs from Liner (2014) equation 10.

    """
    vp0, vs0 = backus(vp, vs, rho, lb, dz)

    ptemp = np.pi * np.log(vp0 / vp) / (np.log(vp0 / vp) + np.log(lb/dz))
    Qp = 1.0 / np.tan(ptemp)

    stemp = np.pi * np.log(vs0 / vs) / (np.log(vs0 / vs) + np.log(lb/dz))
    Qs = 1.0 / np.tan(stemp)

    return Qp, Qs


def thomsen_parameters(vp, vs, rho, lb, dz):
    """
    Liner, C, and T Fei (2006). Layer-induced seismic anisotropy from
    full-wave sonic logs: Theory, application, and validation.
    Geophysics 71 (6), p D183–D190. DOI:10.1190/1.2356997

    """
    A, C, F, L, M = backus_parameters(vp, vs, rho, lb, dz)

    delta = ((F + L)**2.0 - (C - L)**2.0) / (2.0 * C * (C - L))
    epsilon = (A - C) / (2.0 * C)
    gamma = (M - L) / (2.0 * L)

    return delta, epsilon, gamma


def dispersion_parameter(qp):
    """
    Kjartansson (1979). Journal of Geophysical Research, 84 (B9),
    4737-4748. DOI: 10.1029/JB084iB09p04737.
    """
    return np.arctan(1/qp) / np.pi


def blangy(vp1, vs1, rho1, d1, e1, vp0, vs0, rho0, d0, e0, theta):
    """
    Blangy, JP, 1994, AVO in transversely isotropic media-An overview.
    Geophysics 59 (5), 775-781. DOI: 10.1190/1.1443635

    Provide Vp, Vs, rho, delta, epsilon for the upper and lower intervals,
    and theta, the incidence angle.

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.
    :param d1: Thomsen's delta of the upper medium.
    :param e1: Thomsen's epsilon of the upper medium.

    :param vp0: The p-wave velocity of the lower medium.
    :param vs0: The s-wave velocity of the lower medium.
    :param rho0: The density of the lower medium.
    :param d0: Thomsen's delta of the lower medium.
    :param e0: Thomsen's epsilon of the lower medium.

    :param theta: A scalar [degrees].

    :returns: the isotropic and anisotropic reflectivities in a tuple. The
        isotropic result is equivalent to Aki-Rochards.


    TODO
        Use rocks.
    """
    # I already had code set up like this, so build some convenience dicts.
    # These should be rock objects.
    lower = {'vp': vp0,
             'vs': vs0,
             'rho': rho0,
             'd': d0,
             'e': e0,
             }

    upper = {'vp': vp1,
             'vs': vs1,
             'rho': rho1,
             'd': d1,
             'e': e1,
             }

    # Redefine theta
    inc_angle = np.radians(theta)
    trans_angle = np.arcsin(np.sin(inc_angle) * lower['vp']/upper['vp'])
    theta = 0.5 * (inc_angle + trans_angle)

    vp = (upper['vp'] + lower['vp'])/2.0
    vs = (upper['vs'] + lower['vs'])/2.0
    rho = (upper['rho'] + lower['rho'])/2.0

    dvp = lower['vp'] - upper['vp']
    dvs = lower['vs'] - upper['vs']
    drho = lower['rho'] - upper['rho']
    dd = lower['d'] - upper['d']
    de = lower['e'] - upper['e']

    A = 0.5 * (drho/rho + dvp/vp)
    B = 2.0 * (vs**2 / vp**2) * ((drho/rho + 2 * dvs/vs)) * np.sin(theta)**2
    C = 0.5 * (dvp/vp) * np.tan(theta)**2
    D = 0.5 * dd * np.sin(theta)**2
    E = 0.5 * (dd - de) * np.sin(theta)**2 * np.tan(theta)**2

    isotropic = A - B + C
    anisotropic = A - B + C + D - E

    return isotropic, anisotropic


def ruger(vp1, vs1, rho1, d1, e1, vp2, vs2, rho2, d2, e2, theta):
    """
    Coded by Alessandro Amato del Monte and (c) 2016 by him
    https://github.com/aadm/avo_explorer/blob/master/avo_explorer_v2.ipynb

    Rüger, A., 1997, P -wave reflection coefficients for transversely
    isotropic models with vertical and horizontal axis of symmetry:
    Geophysics, v. 62, no. 3, p. 713–722.

    Provide Vp, Vs, rho, delta, epsilon for the upper and lower intervals,
    and theta, the incidence angle.

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.
    :param d1: Thomsen's delta of the upper medium.
    :param e1: Thomsen's epsilon of the upper medium.

    :param vp0: The p-wave velocity of the lower medium.
    :param vs0: The s-wave velocity of the lower medium.
    :param rho0: The density of the lower medium.
    :param d0: Thomsen's delta of the lower medium.
    :param e0: Thomsen's epsilon of the lower medium.

    :param theta: A scalar [degrees].

    :returns: anisotropic reflectivity.

    """

    a = np.radians(theta)
    vp = np.mean([vp1, vp2])
    vs = np.mean([vs1, vs2])
#    rho = np.mean([rho1, rho2])
    z = np.mean([vp1*rho1, vp2*rho2])
    g = np.mean([rho1*vs1**2, rho2*vs2**2])
    dvp = vp2-vp1
#    dvs = vs2-vs1
#    drho = rho2-rho1
    z2, z1 = vp2*rho2, vp1*rho1
    dz = z2-z1
    dg = rho2*vs2**2 - rho1*vs1**2
    dd = d2-d1
    de = e2-e1
    A = 0.5*(dz/z)
    B = 0.5*(dvp/vp - (2*vs/vp)**2 * (dg/g) + dd) * np.sin(a)**2
    C = 0.5*(dvp/vp + de) * np.sin(a)**2 * np.tan(a)**2
    R = A+B+C

    return R


def crack_density(porosity, aspect):
    """
    Returns crack density from porosity and aspect ratio, phi and alpha
    respectively in the unnumbered equation between 15.40 and 15.41 in
    Dvorkin et al. 2014.

    Args:
        porosity (float): Fractional porosity.
        aspect (float): Aspect ratio.

    Returns:
        float: Crack density.
    """
    if porosity >= 1:
        porosity /= 100.

    return 3 * porosity / (4 * np.pi * aspect)


def hudson_delta_M(porosity, aspect, mu, lam=None, pmod=None):
    """
    The approximate reduction in compressional modulus M in the direction
    normal to a set of aligned cracks. Eqn 15.40 in Dvorkin et al (2014).

    Args:
        porosity (float): Fractional porosity, phi.
        aspect (float): Aspect ratio, alpha.
        mu (float): Shear modulus, sometimes called G.
        lam (float): Lame's first parameter.
        pmod (float): Compressional modulus, M.

    Returns:
        float: M_inf - M_0 = \Delta c_11.
    """
    epsilon = crack_density(porosity, aspect)
    if lam:
        return epsilon * (lam**2 / mu) * 4*(lam + 2*mu)/(3*lam + 3*mu)
    else:
        return (4*epsilon/3) * ((pmod - 2*mu)**2 / mu) * (pmod/(pmod-mu))


def hudson_delta_G(porosity, aspect, mu, lam=None, pmod=None):
    """
    The approximate reduction in shear modulus G (or mu) in the direction
    normal to a set of aligned cracks. Eqn 15.42 in Dvorkin et al (2014).

    Args:
        porosity (float): Fractional porosity, phi.
        aspect (float): Aspect ratio, alpha.
        mu (float): Shear modulus, sometimes called G.
        lam (float): Lame's first parameter, lambda.
        pmod (float): Compressional modulus, M.

    Returns:
        float: M_inf - M_0 = \Delta c_11.
    """
    epsilon = crack_density(porosity, aspect)
    if lam:
        return epsilon * mu * 16*(lam + 2*mu)/(9*lam + 12*mu)
    else:
        return (16*mu*epsilon/3) * pmod / (3*pmod - 2*mu)


def hudson_quality_factor(porosity, aspect, mu, lam=None, pmod=None):
    """
    Returns Q_p and Q_s for cracked media. Equations 15.41 and 15.43 in
    Dvorkin et al. (2014).
    """
    Qp = 2*mu / hudson_delta_M(porosity, aspect, mu, lam, pmod)
    Qs = 2*mu / hudson_delta_G(porosity, aspect, mu, lam, pmod)
    return Qp, Qs


def hudson_inverse_Q_ratio(mu=None,
                           pmod=None,
                           pr=None,
                           vp=None,
                           vs=None,
                           aligned=True):
    """
    Dvorkin et al. (2014), Eq 15.44 (aligned) and 15.48 (not aligned).
    """
    if pr:
        x = (2 - 2*pr) / (1 - 2*pr)
    elif vp and vs:
        x = vp**2 / vs**2
    elif mu and pmod:
        x = pmod / mu
    else:
        raise Exception

    if aligned:
        return 0.25 * (x - 2)**2 * (3*x - 2) / (x**2 - x)
    else:
        a = 2*x / (3*x - 2)
        b = x / 3*(x - 1)
        return 1.25 * ((x - 2)**2 / (x - 1)) / (a + b)
