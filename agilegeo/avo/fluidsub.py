# -*- coding: utf-8 -*-
'''
===================
fluidsub.py
===================

Calculates various parameters for fluid substitution
from Vp, Vs, and rho

Created July 2014

@author: Matt Hall, Evan Bianco

Using http://www.subsurfwiki.org/wiki/Gassmann_equation

The algorithm is from Avseth et al (2006), per the wiki page.

Informed by Smith et al, Geophysics 68(2), 2003. 

At some point we should do Biot too, per Russell...
http://cseg.ca/symposium/archives/2012/presentations/Biot_Gassmann_and_me.pdf


'''

import numpy as np
import moduli

def gassmann(ksat1, kf1, kf2, kmin, phi):
    """
    Applies the Gassmann_equation. Takes Ksat1,
    Kfluid1, Kfluid2, Kmineral, and phi.

    Returns Ksat2.
    """

    s  = ksat1 / (kmin - ksat1)
    f1 = kf1 / (phi * (kmin - kf1))
    f2 = kf2 / (phi * (kmin - kf2))

    ksat2 = kmin / ((1/(s - f1 + f2)) + 1)

    return ksat2


def vrh(kclay,  kqtz, vclay=None, vsh=None, vqtz=None):
    """
    Voigt-Reuss-Hill average to find Kmatrix from clay and qtz components.

    Must give one of vsh or vclay.

    From Smith et al, Geophysics 68(2), 2003.

    Works for any two components.

    Returns Kvrh, AKA Kmatrix.
    """
    if not vclay:
        vclay=0.7 * vsh

    if not vqtz:
        vqtz = 1 - vclay

    kreuss = 1. / (vclay/kclay + vqtz/kqtz)
    kvoigt = vclay*kclay + vqtz*kqtz
    kvrh  = 0.5 * (kreuss + kvoigt)

    return kvrh

def rhogas(gravity, temp, pressure):
    """ 
    From http://www.spgindia.org/geohorizon/jan_2006/dhananjay_paper.pdf
    """

    # Compute pseudo-reduced temp and pressure:
    tpr = (temp + 273.15) / (gravity * (94.72 + 170.75))
    ppr = pressure / (4.892 - 0.4048*gravity)

    exponent = -1 * (0.45 + 8 * (0.56 - (1/tpr))**2.0) * ppr**1.2 / tpr
    bige = 0.109 * (3.85 - tpr)**2.0 * np.exp(exponent)

    z = ppr * (0.03 + 0.00527 * (3.5 - tpr)) + (0.642*tpr - 0.007*tpr**4.0 - 0.52) + bige

    rhogas = 28.8 * gravity * pressure / (z * r * (temp + 273.15))

    return rhogas

def rhosat(phi, sw, rhomin, rhow, rhohc):
    """ 
    Density of partially saturated rock.

    """
    a = rhomin * (1 - phi)        # grains
    b = rhow * sw * phi           # brine
    c = rhohc * (1 - sw) * phi    # hydrocarbon

    return a + b + c


def fluidsub(vp, vs, rho, phi, rhof1, rhof2, kmin, kf1, kf2):
    """
    Naive fluid substitution.
    No pressure/temperature correction.

    :param vp: P-wave velocity 
    :param vs: S-wave velocity
    :param rho: bulk density
    :param phi: porosity (i.e. 0.20)
    :param rhof1: bulk density of original fluid (base case)
    :param rhof2: bulk density of substitute fluid (subbed case)
    :param kmin: bulk modulus of solid mineral(s)
    :param kf1: bulk modulus of original fluid
    :param kf2: bulk modulus of substitue fluid
    
    (only works for SI units right now)
    
    e.g.: 
    phi = 0.20
    kmin = 36000000000 Pa
    vp = 3800 m/s
    vs = 2400 m/s
    rho = 2450 kg / m^3 
    rhof1 = 1100 kg / m^3 
    rhof2 = 25 kg / m^3 

    returns Vp, Vs, and rho for the substituted case
    """
    
    # Step 1: Extract the dynamic bulk and shear moduli
    ksat1 = moduli.bulk(vp=vp, vs=vs, rho=rho)
    musat1 = moduli.mu(vp=vp, vs=vs, rho=rho)
        
    # Step 2: Apply Gassmann's relation
    ksat2 = gassmann(ksat1=ksat1, kf1=kf1, kf2=kf2, kmin=kmin, phi=phi)

    # Step 3: Leave the shear modulus unchanged
    musat2 = musat1

    # Step 4: Correct the bulk density for the change in fluid
    rho2 = rho + phi * (rhof2 - rhof1)

    # Step 5: recompute the fluid substituted velocities
    vp2 = moduli.vp(bulk=ksat2, mu=musat2, rho=rho2)
    vs2 = moduli.vs(mu=musat2, rho=rho2)

    return vp2, vs2, rho2

