#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Thin-layer anisotropy.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np

from bruges.rockphysics import moduli
from bruges.util import moving_average


def backus_parameters(vp, vs, rho, lb, dz):
    """
    Backus parameters.

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


def backus(vp, vs, rho, lb, dz):
    """
    Backus averaging.

    Liner, C (2014), Long-wave elastic attenuation produced by horizontal
        layering. The Leading Edge, June 2014, p 634-638.

    """
    # Compute the Backus parameters:
    A, C, F, L, M = backus_parameters(vp, vs, rho, lb, dz)

    # Compute the vertical velocities from Liner (2014) equation 5:
    R = moving_average(rho, lb/dz, mode='same')
    vp0 = np.sqrt(C / R)
    vs0 = np.sqrt(L / R)

    return vp0, vs0


def quality_factor(vp, vs, rho, lb, dz):
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


def blangy(vp0, vs0, rho0, d0, e0, vp1, vs1, rho1, d1, e1, theta):
    """
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
