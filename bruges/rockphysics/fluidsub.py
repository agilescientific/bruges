"""
:copyright: 2015 Agile Geoscience
:license: Apache 2.0

===================
fluidsub.py
===================

Calculates various parameters for fluid substitution
from Vp, Vs, and rho

Matt Hall, Evan Bianco, July 2014

Using http://www.subsurfwiki.org/wiki/Gassmann_equation

The algorithm is from Avseth et al (2006), per the wiki page.

Informed by Smith et al, Geophysics 68(2), 2003.

At some point we should do Biot too, per Russell...
http://cseg.ca/symposium/archives/2012/presentations/Biot_Gassmann_and_me.pdf
"""
from collections import namedtuple

import numpy as np
from . import moduli
from . import fluids


def avseth_gassmann(ksat1, kf1, kf2, k0, phi):
    """
    Applies the inverse Gassmann's equation to calculate the rock bulk modulus
    saturated with fluid with bulk modulus. Requires as inputs the insitu
    fluid bulk modulus, the insitu saturated rock bulk modulus, the mineral
    matrix bulk modulus and the porosity.

    Args:
        ksat1 (float): initial saturated rock bulk modulus.
        kf1 (float): initial fluid bulk modulus.
        kf2 (float): final fluid bulk modulus.
        k0 (float): mineral bulk modulus.
        phi (float): porosity.

    Returns:
        ksat2 (float): final saturated rock bulk modulus.
    """

    s = ksat1 / (k0 - ksat1)
    f1 = kf1 / (phi * (k0 - kf1))
    f2 = kf2 / (phi * (k0 - kf2))
    ksat2 = k0 / ((1/(s - f1 + f2)) + 1)

    return ksat2


def smith_gassmann(kdry, k0, kf, phi):
    """
    Applies the direct Gassmann's equation to calculate the saturated rock
    bulk modulus from porosity and the dry-rock, mineral and fluid bulk moduli.

    Args:
        kdry (float): dry-rock bulk modulus.
        k0 (float): mineral bulk modulus.
        kf (float): fluid bulk modulus.
        phi (float): Porosity.

    Returns:
        ksat (float): saturated rock bulk modulus.
    """

    a = (1 - kdry/k0)**2.0
    b = phi/kf + (1-phi)/k0 - (kdry/k0**2.0)
    ksat = kdry + (a/b)

    return ksat


def vrh(kclay, kqtz, vclay):
    """
    Voigt-Reuss-Hill average to find Kmatrix from clay and qtz components.
    Works for any two components, they don't have to be clay and quartz.

    From Smith et al, Geophysics 68(2), 2003.

    Args:
        kclay (float): K_clay.
        kqtz (float): K_quartz.
        vclay (float): V_clay.

    Returns:
        Kvrh, also known as Kmatrix.
    """

    vqtz = 1 - vclay

    kreuss = 1. / (vclay/kclay + vqtz/kqtz)
    kvoigt = vclay*kclay + vqtz*kqtz
    kvrh = 0.5 * (kreuss + kvoigt)

    return kvrh


def rhogas(gravity, temp, pressure):
    """
    From http://www.spgindia.org/geohorizon/jan_2006/dhananjay_paper.pdf
    """
    R = 8.3144621  # Gas constant in J.mol^-1.K^-1

    # Compute pseudo-reduced temp and pressure:
    tpr = (temp + 273.15) / (gravity * (94.72 + 170.75))
    ppr = pressure / (4.892 - 0.4048*gravity)

    exponent = -1 * (0.45 + 8 * (0.56 - (1/tpr))**2.0) * ppr**1.2 / tpr
    bige = 0.109 * (3.85 - tpr)**2.0 * np.exp(exponent)
    term2 = 0.642*tpr - 0.007*tpr**4.0 - 0.52

    Z = ppr * (0.03 + 0.00527 * (3.5 - tpr)) + term2 + bige

    rhogas = 28.8 * gravity * pressure / (Z * R * (temp + 273.15))

    return rhogas


def rhosat(phi, sw, rhomin, rhow, rhohc):
    """
    Density of partially saturated rock.
    """
    a = rhomin * (1 - phi)        # grains
    b = rhow * sw * phi           # brine
    c = rhohc * (1 - sw) * phi    # hydrocarbon

    return a + b + c


def avseth_fluidsub(vp, vs, rho, phi, rhof1, rhof2, kmin, kf1, kf2):
    """
    Naive fluid substitution from Avseth et al, section 1.3.1. Bulk modulus of
    fluid 1 and fluid 2 are already known, and the bulk modulus of the dry
    frame, Kmin, is known. Use SI units.

    Args:
        vp (float): P-wave velocity
        vs (float): S-wave velocity
        rho (float): bulk density
        phi (float): porosity (volume fraction, e.g. 0.20)
        rhof1 (float): bulk density of original fluid (base case)
        rhof2 (float): bulk density of substitute fluid (subbed case)
        kmin (float): bulk modulus of solid mineral (s)
        kf1 (float): bulk modulus of original fluid
        kf2 (float): bulk modulus of substitute fluid

    Returns:
        Tuple: Vp, Vs, and rho for the substituted case
    """

    # Step 1: Extract the dynamic bulk and shear moduli.
    ksat1 = moduli.bulk(vp=vp, vs=vs, rho=rho)
    musat1 = moduli.mu(vp=vp, vs=vs, rho=rho)

    # Step 2: Apply Gassmann's relation.
    ksat2 = avseth_gassmann(ksat1=ksat1, kf1=kf1, kf2=kf2, k0=kmin, phi=phi)

    # Step 3: Leave the shear modulus unchanged.
    musat2 = musat1

    # Step 4: Correct the bulk density for the change in fluid.
    rho2 = rho + phi * (rhof2 - rhof1)

    # Step 5: recompute the fluid substituted velocities.
    vp2 = moduli.vp(bulk=ksat2, mu=musat2, rho=rho2)
    vs2 = moduli.vs(mu=musat2, rho=rho2)

    FluidSubResult = namedtuple('FluidSubResult', ['Vp', 'Vs', 'rho'])
    return FluidSubResult(vp2, vs2, rho2)


def smith_fluidsub(vp, vs, rho, phi, rhow, rhohc,
                   sw, swnew, kw, khc, kclay, kqtz,
                   vclay,
                   rhownew=None, rhohcnew=None,
                   kwnew=None, khcnew=None
                   ):
    """
    Naive fluid substitution from Smith et al. 2003. No pressure/temperature
    correction. Only works for SI units right now.

    Args:
        vp (float): P-wave velocity
        vs (float): S-wave velocity
        rho (float): bulk density
        phi (float): porosity (v/v, fraction)
        rhow (float): density of water
        rhohc (float): density of HC
        sw (float): water saturation in base case
        swnew (float): water saturation in subbed case
        kw (float):  bulk modulus of water
        khc (float): bulk modulus of hydrocarbon
        kclay (float): bulk modulus of clay
        kqtz (float):  bulk modulus of quartz
        vclay (float): Vclay, v/v
        rhownew (float): density of water in subbed case (optional)
        rhohcnew (float): density of HC in subbed case (optional)
        kwnew (float):  bulk modulus of water in subbed case (optional)
        khcnew (float): bulk modulus of hydrocarbon in subbed case (optional)

    Returns:
        Tuple: Vp, Vs, and rho for the substituted case.
    """

    # Using the workflow in Smith et al., Table 2
    # Using Smith's notation, more or less (not the same
    # as Avseth's notation).
    #
    # Step 1: Log edits and interpretation.
    #
    # Step 2. Shear velocity estimation, if necessary.
    #
    # Step 3. Calculate K and G for the in-situ conditions.
    ksat = moduli.bulk(vp=vp, vs=vs, rho=rho)
    g = moduli.mu(vs=vs, rho=rho)

    # Step 4. Calculate K0 based on lithology estimates (VRH or HS mixing).
    k0 = vrh(kclay=kclay, kqtz=kqtz, vclay=vclay)

    # Step 5. Calculate fluid properties (K and ρ).
    # Step 6. Mix fluids for the in-situ case according to Sw.
    kfl = fluids.wood(kw, khc, sw)
    rhofl = sw * rhow + (1-sw)*rhohc

    # Step 7: Calculate K*.
    a = ksat * ((phi*k0/kfl) + 1 - phi) - k0
    b = (phi*k0/kfl) + (ksat/k0) - 1 - phi
    kstar = a / b

    # Step 8: Calculate new fluid properties (K and ρ) at the desired Sw
    # First set the new fluid properties, in case they are unchanged.
    if kwnew is None:
        kwnew = kw
    if rhownew is None:
        rhownew = rhow
    if khcnew is None:
        khcnew = khc
    if rhohcnew is None:
        rhohcnew = rhohc

    # Now calculate the new fluid properties.
    kfl2 = 1 / (swnew/kwnew + (1-swnew)/khcnew)
    rhofl2 = swnew * rhownew + (1-swnew)*rhohcnew

    # Step 9: Calculate the new saturated bulk modulus of the rock
    # using Gassmann.
    ksat2 = smith_gassmann(kdry=kstar, k0=k0, kf=kfl2, phi=phi)

    # Step 10: Calculate the new bulk density.
    # First we need the grain density...
    rhog = (rho - phi*rhofl) / (1-phi)
    # Now we can find the new bulk density
    rho2 = phi*rhofl2 + rhog*(1-phi)

    # Step 11: Calculate the new compressional velocity.
    # Remember, mu (G) is unchanged.
    vp2 = moduli.vp(bulk=ksat2, mu=g, rho=rho2)

    # Step 12: Calculate the new shear velocity.
    vs2 = moduli.vs(mu=g, rho=rho2)

    # Finish.
    FluidSubResult = namedtuple('FluidSubResult', ['Vp', 'Vs', 'rho'])
    return FluidSubResult(vp2, vs2, rho2)


def vels(Kdry, Gdry, K0, D0, Kf, Df, phi):
    '''
    Calculate velocities and densities of saturated rock via Gassmann equation.
    Provide all quantities in SI units.

    Parameters
    ----------   
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in Pa.
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in Pa.
    Kf, Df : float or array_like
        Fluid bulk modulus in Pa and density in kg/m^3.
    phi : float or array_like
        Porosity, v/v.

    Returns
    -------
    vp, vs : float or array_like
        Saturated rock P- and S-wave velocities in m/s.
    rho: float or array_like
        Saturated rock density in kg/m^3.
    K : float or array_like
        Saturated rock bulk modulus in Pa.
    '''
    rho = D0 * (1 - phi) + Df * phi
    with np.errstate(divide='ignore', invalid='ignore'):
        K = Kdry + (1 - Kdry / K0)**2 / ( (phi / Kf)
            + ((1 - phi) / K0) - (Kdry / K0**2) )
        vp = np.sqrt((K + 4/3 * Gdry) / rho)
        vs = np.sqrt(Gdry / rho)
    return vp, vs, rho, K
