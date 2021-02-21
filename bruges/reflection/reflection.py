# -*- coding: utf-8 -*-
"""
Various reflectivity algorithms.

:copyright: 2018 Agile Geoscience
:license: Apache 2.0
"""
from functools import wraps
from collections import namedtuple

import numpy as np

from bruges.rockphysics import moduli
from bruges.rockphysics import anisotropy
from bruges.util import deprecated


def critical_angles(vp1, vp2, vs2=None):
    """Compute critical angle at an interface

    Given the upper (vp1) and lower (vp2) velocities at an interface. If you want the PS-wave critical angle as well,
    pass vs2 as well.

    Args:
        vp1 (ndarray): Upper layer P-wave velocity.
        vp2 (ndarray): Lower layer P-wave velocity.
        vs2 (ndarray): Lower layer S-wave velocity (optional).

    Returns:
        tuple: The first and second critical angles at the interface, in
            degrees. If there isn't a critical angle, it returns np.nan.
    """
    ca1 = ca2 = np.nan

    if vp1 < vp2:
        ca1 = np.degrees(np.arcsin(vp1/vp2))

    if (vs2 is not None) and (vp1 < vs2):
        ca2 = np.degrees(np.arcsin(vp1/vs2))

    return ca1, ca2


def reflection_phase(reflectivity):
    """
    Compute the phase of the reflectivity. Returns an array (or float) of
    the phase, in positive multiples of 180 deg or pi rad. So 1 is opposite
    phase. A value of 1.1 would be +1.1 \times \pi rad.

    Args:
        reflectivity (ndarray): The reflectivity, eg from `zoeppritz()`.

    Returns:
        ndarray: The phase, strictly positive
    """
    ph = np.arctan2(np.imag(reflectivity), np.real(reflectivity)) / np.pi
    ph[ph == 1] = 0
    ph[ph < 0] = 2 + ph[ph < 0]
    return ph


def acoustic_reflectivity(vp, rho):
    """
    The acoustic reflectivity, given Vp and RHOB logs.

    Args:
        vp (ndarray): The P-wave velocity.
        rho (ndarray): The bulk density.

    Returns:
        ndarray: The reflectivity coefficient series.
    """
    upper = vp[:-1] * rho[:-1]
    lower = vp[1:] * rho[1:]
    return (lower - upper) / (lower + upper)

def reflectivity(vp, vs, rho, theta=0, method='zoeppritz_rpp'):
    """
    Offset reflectivity, given Vp, Vs, rho, and offset.

    Computes 'upper' and 'lower' intervals from the three provided arrays,
    then passes the result to the specified method to compute reflection
    coefficients.

    For acoustic reflectivity, either use the `acoustic_reflectivity()`
    function, or call `reflectivity()` passing any log as Vs, e.g. just give
    the Vp log twice (it won't be used anyway):

        reflectivity(vp, vp, rho)

    For anisotropic reflectivity, use either `anisotropy.blangy()` or
    `anisotropy.ruger()`.

    Args:
        vp (ndarray): The P-wave velocity; float or 1D array length m.
        vs (ndarray): The S-wave velocity; float or 1D array length m.
        rho (ndarray): The density; float or 1D array length m.
        theta (ndarray): The incidence angle; float or 1D array length n.
        method (str): The reflectivity equation to use; one of:

                - 'scattering_matrix': scattering_matrix
                - 'zoeppritz_element': zoeppritz_element
                - 'zoeppritz': zoeppritz
                - 'zoeppritz_rpp': zoeppritz_rpp
                - 'akirichards': akirichards
                - 'akirichards_alt': akirichards_alt
                - 'fatti': fatti
                - 'shuey': shuey
                - 'bortfeld': bortfeld
                - 'hilterman': hilterman

        Notes:

                - scattering_matrix gives the full solution
                - zoeppritz_element gives a single element which you specify
                - zoeppritz returns RPP element only; use zoeppritz_rpp instead
                - zoeppritz_rpp is faster than zoeppritz or zoeppritz_element

    Returns:
        ndarray. The result of running the specified method on the inputs.
            Will be a float (for float inputs and one angle), a 1 x n array
            (for float inputs and an array of angles), a 1 x m-1 array (for
            float inputs and one angle), or an m-1 x n array (for array inputs
            and an array of angles).
    """
    methods = {
        'scattering_matrix': scattering_matrix,
        'zoeppritz_element': zoeppritz_element,
        'zoeppritz': zoeppritz,
        'zoeppritz_rpp': zoeppritz_rpp,
        'akirichards': akirichards,
        'akirichards_alt': akirichards_alt,
        'fatti': fatti,
        'shuey': shuey,
        'bortfeld': bortfeld,
        'hilterman': hilterman,
    }
    func = methods[method.lower()]
    vp = np.asanyarray(vp, dtype=float)
    vs = np.asanyarray(vs, dtype=float)
    rho = np.asanyarray(rho, dtype=float)

    vp1, vp2 = vp[:-1], vp[1:]
    vs1, vs2 = vs[:-1], vs[1:]
    rho1, rho2 = rho[:-1], rho[1:]

    return func(vp1, vs1, rho1, vp2, vs2, rho2, theta)


def vectorize(func):
    """
    Decorator to make sure the inputs are arrays. We also add a dimension
    to theta to make the functions work in an 'outer product' way.

    Takes a reflectivity function requiring Vp, Vs, and RHOB for 2 rocks
    (upper and lower), plus incidence angle theta, plus kwargs. Returns
    that function with the arguments transformed to ndarrays.
    """
    @wraps(func)
    def wrapper(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, **kwargs):
        vp1 = np.asanyarray(vp1, dtype=float)
        vs1 = np.asanyarray(vs1, dtype=float) + 1e-12  # Prevent singular matrix.
        rho1 = np.asanyarray(rho1, dtype=float)
        vp2 = np.asanyarray(vp2, dtype=float)
        vs2 = np.asanyarray(vs2, dtype=float) + 1e-12  # Prevent singular matrix.
        rho2 = np.asanyarray(rho2, dtype=float)
        theta1 = np.asanyarray(theta1).reshape((-1, 1))
        return func(vp1, vs1, rho1, vp2, vs2, rho2, theta1, **kwargs)
    return wrapper


def preprocess(func):
    """
    Decorator to preprocess arguments for the reflectivity equations.

    Takes a reflectivity function requiring Vp, Vs, and RHOB for 2 rocks
    (upper and lower), plus incidence angle theta, plus kwargs. Returns
    that function with some arguments transformed.
    """
    @wraps(func)
    def wrapper(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, **kwargs):

        # Interpret tuple for theta1 as a linspace.
        if isinstance(theta1, tuple):
            if len(theta1) == 2:
                start, stop = theta1
                theta1 = np.linspace(start, stop, num=stop+1)
            elif len(theta1) == 3:
                start, stop, step = theta1
                steps = (stop / step) + 1
                theta1 = np.linspace(start, stop, num=steps)
            else:
                raise TypeError("Expected 2 or 3 parameters for theta1 expressed as range.")

        # Convert theta1 to radians and complex numbers.
        theta1 = np.radians(theta1).astype(complex)

        return func(vp1, vs1, rho1, vp2, vs2, rho2, theta1, **kwargs)
    return wrapper


@preprocess
@vectorize
def scattering_matrix(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0):
    """
    Full Zoeppritz solution, considered the definitive solution.
    Calculates the angle dependent p-wave reflectivity of an interface
    between two mediums.

    Originally written by: Wes Hamlyn, vectorized by Agile.

    Returns the complex reflectivity.

    Args:
        vp1 (float): The upper P-wave velocity.
        vs1 (float): The upper S-wave velocity.
        rho1 (float): The upper layer's density.
        vp2 (float): The lower P-wave velocity.
        vs2 (float): The lower S-wave velocity.
        rho2 (float): The lower layer's density.
        theta1 (ndarray): The incidence angle; float or 1D array length n.

    Returns:
        ndarray. The exact Zoeppritz solution for all modes at the interface.
            A 4x4 array representing the scattering matrix at the incident
            angle theta1.
    """
    theta1 *= np.ones_like(vp1)
    p = np.sin(theta1) / vp1  # Ray parameter.
    theta2 = np.arcsin(p * vp2)  # Trans. angle of P-wave.
    phi1 = np.arcsin(p * vs1)    # Refl. angle of converted S-wave.
    phi2 = np.arcsin(p * vs2)    # Trans. angle of converted S-wave.

    # Matrix form of Zoeppritz equations... M & N are matrices.
    M = np.array([[-np.sin(theta1), -np.cos(phi1), np.sin(theta2), np.cos(phi2)],
                  [np.cos(theta1), -np.sin(phi1), np.cos(theta2), -np.sin(phi2)],
                  [2 * rho1 * vs1 * np.sin(phi1) * np.cos(theta1),
                   rho1 * vs1 * (1 - 2 * np.sin(phi1) ** 2),
                   2 * rho2 * vs2 * np.sin(phi2) * np.cos(theta2),
                   rho2 * vs2 * (1 - 2 * np.sin(phi2) ** 2)],
                  [-rho1 * vp1 * (1 - 2 * np.sin(phi1) ** 2),
                   rho1 * vs1 * np.sin(2 * phi1),
                   rho2 * vp2 * (1 - 2 * np.sin(phi2) ** 2),
                   -rho2 * vs2 * np.sin(2 * phi2)]])

    N = np.array([[np.sin(theta1), np.cos(phi1), -np.sin(theta2), -np.cos(phi2)],
                  [np.cos(theta1), -np.sin(phi1), np.cos(theta2), -np.sin(phi2)],
                  [2 * rho1 * vs1 * np.sin(phi1) * np.cos(theta1),
                   rho1 * vs1 * (1 - 2 * np.sin(phi1) ** 2),
                   2 * rho2 * vs2 * np.sin(phi2) * np.cos(theta2),
                   rho2 * vs2 * (1 - 2 * np.sin(phi2) ** 2)],
                  [rho1 * vp1 * (1 - 2 * np.sin(phi1) ** 2),
                   -rho1 * vs1 * np.sin(2 * phi1),
                   - rho2 * vp2 * (1 - 2 * np.sin(phi2) ** 2),
                   rho2 * vs2 * np.sin(2 * phi2)]])

    M_ = np.moveaxis(np.squeeze(M), [0, 1], [-2, -1])
    A = np.linalg.inv(M_)
    N_ = np.moveaxis(np.squeeze(N), [0, 1], [-2, -1])
    Z_ = np.matmul(A, N_)

    return np.transpose(Z_, axes=list(range(Z_.ndim - 2)) + [-1, -2])


def zoeppritz_element(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, element='PdPu'):
    """
    Returns any mode reflection coefficients from the Zoeppritz
    scattering matrix. Pass in the mode as element, e.g. 'PdSu' for PS.

    Wraps scattering_matrix().

    Returns the complex reflectivity.

    Args:
        vp1 (float): The upper P-wave velocity.
        vs1 (float): The upper S-wave velocity.
        rho1 (float): The upper layer's density.
        vp2 (float): The lower P-wave velocity.
        vs2 (float): The lower S-wave velocity.
        rho2 (float): The lower layer's density.
        theta1 (ndarray): The incidence angle; float or 1D array length n.
        element (str): The name of the element to return, must be one of:
            'PdPu', 'SdPu', 'PuPu', 'SuPu', 'PdSu', 'SdSu', 'PuSu', 'SuSu',
            'PdPd', 'SdPd', 'PuPd', 'SuPd', 'PdSd', 'SdSd', 'PuSd', 'SuSd'.

    Returns:
        ndarray. Array length n of the exact Zoeppritz solution for the
            specified modes at the interface, at the incident angle theta1.
    """
    elements = np.array([['PdPu', 'SdPu', 'PuPu', 'SuPu'],
                         ['PdSu', 'SdSu', 'PuSu', 'SuSu'],
                         ['PdPd', 'SdPd', 'PuPd', 'SuPd'],
                         ['PdSd', 'SdSd', 'PuSd', 'SuSd']])

    Z = scattering_matrix(vp1, vs1, rho1, vp2, vs2, rho2, theta1).T

    return np.squeeze(Z[np.where(elements == element)].T)


def zoeppritz(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0):
    """
    Returns the PP reflection coefficients from the Zoeppritz
    scattering matrix. Wraps zoeppritz_element().

    Returns the complex reflectivity.

    Args:
        vp1 (float): The upper P-wave velocity.
        vs1 (float): The upper S-wave velocity.
        rho1 (float): The upper layer's density.
        vp2 (float): The lower P-wave velocity.
        vs2 (float): The lower S-wave velocity.
        rho2 (float): The lower layer's density.
        theta1 (ndarray): The incidence angle; float or 1D array length n.

    Returns:
        ndarray. Array length n of the exact Zoeppritz solution for the
            specified modes at the interface, at the incident angle theta1.
    """
    return zoeppritz_element(vp1, vs1, rho1, vp2, vs2, rho2, theta1, 'PdPu')


@preprocess
@vectorize
def zoeppritz_rpp(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0):
    """
    Exact Zoeppritz from expression.

    This is useful because we can pass arrays to it, which we can't do to
    scattering_matrix().

    Dvorkin et al. (2014). Seismic Reflections of Rock Properties. Cambridge.

    Returns the complex reflectivity.

    Args:
        vp1 (ndarray): The upper P-wave velocity; float or 1D array length m.
        vs1 (ndarray): The upper S-wave velocity; float or 1D array length m.
        rho1 (ndarray): The upper layer's density; float or 1D array length m.
        vp2 (ndarray): The lower P-wave velocity; float or 1D array length m.
        vs2 (ndarray): The lower S-wave velocity; float or 1D array length m.
        rho2 (ndarray): The lower layer's density; float or 1D array length m.
        theta1 (ndarray): The incidence angle; float or 1D array length n.

    Returns:
        ndarray. The exact Zoeppritz solution for P-P reflectivity at the
            interface. Will be a float (for float inputs and one angle), a
            1 x n array (for float inputs and an array of angles), a 1 x m
            array (for float inputs and one angle), or an n x m array (for
            array inputs and an array of angles).
    """
    p = np.sin(theta1) / vp1  # Ray parameter
    theta2 = np.arcsin(p * vp2)
    phi1 = np.arcsin(p * vs1)  # Reflected S
    phi2 = np.arcsin(p * vs2)  # Transmitted S

    a = rho2 * (1 - 2 * np.sin(phi2)**2.) - rho1 * (1 - 2 * np.sin(phi1)**2.)
    b = rho2 * (1 - 2 * np.sin(phi2)**2.) + 2 * rho1 * np.sin(phi1)**2.
    c = rho1 * (1 - 2 * np.sin(phi1)**2.) + 2 * rho2 * np.sin(phi2)**2.
    d = 2 * (rho2 * vs2**2 - rho1 * vs1**2)

    E = (b * np.cos(theta1) / vp1) + (c * np.cos(theta2) / vp2)
    F = (b * np.cos(phi1) / vs1) + (c * np.cos(phi2) / vs2)
    G = a - d * np.cos(theta1)/vp1 * np.cos(phi2)/vs2
    H = a - d * np.cos(theta2)/vp2 * np.cos(phi1)/vs1

    D = E*F + G*H*p**2

    rpp = (1/D) * (F*(b*(np.cos(theta1)/vp1) - c*(np.cos(theta2)/vp2)) \
                   - H*p**2 * (a + d*(np.cos(theta1)/vp1)*(np.cos(phi2)/vs2)))

    return np.squeeze(rpp)


@preprocess
@vectorize
def akirichards(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    The Aki-Richards approximation to the reflectivity.

    This is the formulation from Avseth et al., Quantitative seismic
    interpretation, Cambridge University Press, 2006. Adapted for a 4-term
    formula. See http://subsurfwiki.org/wiki/Aki-Richards_equation.

    Returns the complex reflectivity.

    Args:
        vp1 (ndarray): The upper P-wave velocity; float or 1D array length m.
        vs1 (ndarray): The upper S-wave velocity; float or 1D array length m.
        rho1 (ndarray): The upper layer's density; float or 1D array length m.
        vp2 (ndarray): The lower P-wave velocity; float or 1D array length m.
        vs2 (ndarray): The lower S-wave velocity; float or 1D array length m.
        rho2 (ndarray): The lower layer's density; float or 1D array length m.
        theta1 (ndarray): The incidence angle; float or 1D array length n.
        terms (bool): Whether or not to return a tuple of the terms of the
            equation. The first term is the acoustic impedance.

    Returns:
        ndarray. The Aki-Richards approximation for P-P reflectivity at the
            interface. Will be a float (for float inputs and one angle), a
            1 x n array (for float inputs and an array of angles), a 1 x m
            array (for float inputs and one angle), or an n x m array (for
            array inputs and an array of angles).
    """
    theta2 = np.arcsin(vp2/vp1*np.sin(theta1))
    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    meantheta = (theta1+theta2) / 2.0
    rho = (rho1+rho2) / 2.0
    vp = (vp1+vp2) / 2.0
    vs = (vs1+vs2) / 2.0

    # Compute the coefficients
    w = 0.5 * drho/rho
    x = 2 * (vs/vp1)**2 * drho/rho
    y = 0.5 * (dvp/vp)
    z = 4 * (vs/vp1)**2 * (dvs/vs)

    # Compute the terms
    term1 = w
    term2 = -1 * x * np.sin(theta1)**2
    term3 = y / np.cos(meantheta)**2
    term4 = -1 * z * np.sin(theta1)**2

    if terms:
        fields = ['term1', 'term2', 'term3', 'term4']
        AkiRichards = namedtuple('AkiRichards', fields)
        return AkiRichards(np.squeeze([term1 for _ in theta1]),
                           np.squeeze(term2),
                           np.squeeze(term3),
                           np.squeeze(term4)
                           )
    else:
        return np.squeeze(term1 + term2 + term3 + term4)


@preprocess
@vectorize
def akirichards_alt(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    This is another formulation of the Aki-Richards solution.
    See http://subsurfwiki.org/wiki/Aki-Richards_equation

    Returns the complex reflectivity.

    Args:
        vp1 (ndarray): The upper P-wave velocity; float or 1D array length m.
        vs1 (ndarray): The upper S-wave velocity; float or 1D array length m.
        rho1 (ndarray): The upper layer's density; float or 1D array length m.
        vp2 (ndarray): The lower P-wave velocity; float or 1D array length m.
        vs2 (ndarray): The lower S-wave velocity; float or 1D array length m.
        rho2 (ndarray): The lower layer's density; float or 1D array length m.
        theta1 (ndarray): The incidence angle; float or 1D array length n.
        terms (bool): Whether or not to return a tuple of the terms of the
            equation. The first term is the acoustic impedance.

    Returns:
        ndarray. The Aki-Richards approximation for P-P reflectivity at the
            interface. Will be a float (for float inputs and one angle), a
            1 x n array (for float inputs and an array of angles), a 1 x m
            array (for float inputs and one angle), or an n x m array (for
            array inputs and an array of angles).
    """
    theta2 = np.arcsin(vp2/vp1*np.sin(theta1))
    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    theta = (theta1+theta2)/2.0
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0

    # Compute the three terms
    term1 = 0.5 * (dvp/vp + drho/rho)
    term2 = (0.5*dvp/vp-2*(vs/vp)**2*(drho/rho+2*dvs/vs)) * np.sin(theta)**2
    term3 = 0.5 * dvp/vp * (np.tan(theta)**2 - np.sin(theta)**2)

    if terms:
        fields = ['term1', 'term2', 'term3']
        AkiRichards = namedtuple('AkiRichards', fields)
        return AkiRichards(np.squeeze([term1 for _ in theta1]),
                           np.squeeze(term2),
                           np.squeeze(term3)
                           )
    else:
        return np.squeeze(term1 + term2 + term3)


@preprocess
@vectorize
def fatti(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    Compute reflectivities with Fatti's formulation of the Aki-Richards
    equation, which does not account for the critical angle. See Fatti et al.
    (1994), Geophysics 59 (9). Real numbers only.

    Args:
        vp1 (ndarray): The upper P-wave velocity; float or 1D array length m.
        vs1 (ndarray): The upper S-wave velocity; float or 1D array length m.
        rho1 (ndarray): The upper layer's density; float or 1D array length m.
        vp2 (ndarray): The lower P-wave velocity; float or 1D array length m.
        vs2 (ndarray): The lower S-wave velocity; float or 1D array length m.
        rho2 (ndarray): The lower layer's density; float or 1D array length m.
        theta1 (ndarray): The incidence angle; float or 1D array length n.
        terms (bool): Whether or not to return a tuple of the terms of the
            equation. The first term is the acoustic impedance.

    Returns:
        ndarray. The Fatti approximation for P-P reflectivity at the
            interface. Will be a float (for float inputs and one angle), a
            1 x n array (for float inputs and an array of angles), a 1 x m
            array (for float inputs and one angle), or an n x m array (for
            array inputs and an array of angles).
    """
    theta1 = np.real(theta1)

    drho = rho2-rho1
    rho = (rho1+rho2) / 2.0
    vp = (vp1+vp2) / 2.0
    vs = (vs1+vs2) / 2.0
    dip = (vp2*rho2 - vp1*rho1)/(vp2*rho2 + vp1*rho1)
    dis = (vs2*rho2 - vs1*rho1)/(vs2*rho2 + vs1*rho1)
    d = drho/rho

    # Compute the three terms
    term1 = (1 + np.tan(theta1)**2) * dip
    term2 = -8 * (vs/vp)**2 * dis * np.sin(theta1)**2
    term3 = -1 * (0.5 * np.tan(theta1)**2 - 2 * (vs/vp)**2 * np.sin(theta1)**2) * d

    if terms:
        fields = ['term1', 'term2', 'term3']
        Fatti = namedtuple('Fatti', fields)
        return Fatti(np.squeeze(term1),
                     np.squeeze(term2),
                     np.squeeze(term3)
                     )
    else:
        return np.squeeze(term1 + term2 + term3)


@preprocess
@vectorize
def shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0,
          terms=False,
          return_gradient=False):
    """
    Compute Shuey approximation with 3 terms.
    http://subsurfwiki.org/wiki/Shuey_equation

    Args:
        vp1 (ndarray): The upper P-wave velocity; float or 1D array length m.
        vs1 (ndarray): The upper S-wave velocity; float or 1D array length m.
        rho1 (ndarray): The upper layer's density; float or 1D array length m.
        vp2 (ndarray): The lower P-wave velocity; float or 1D array length m.
        vs2 (ndarray): The lower S-wave velocity; float or 1D array length m.
        rho2 (ndarray): The lower layer's density; float or 1D array length m.
        theta1 (ndarray): The incidence angle; float or 1D array length n.
        terms (bool): Whether or not to return a tuple of the terms of the
            equation. The first term is the acoustic impedance.
        return_gradient (bool): Whether to return a tuple of the intercept
            and gradient (i.e. the second term divided by sin^2(theta)).

    Returns:
        ndarray. The Aki-Richards approximation for P-P reflectivity at the
            interface. Will be a float (for float inputs and one angle), a
            1 x n array (for float inputs and an array of angles), a 1 x m
            array (for float inputs and one angle), or an n x m array (for
            array inputs and an array of angles).
    """
    theta1 = np.real(theta1)

    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0

    # Compute three-term reflectivity
    r0 = 0.5 * (dvp/vp + drho/rho)
    g = 0.5 * dvp/vp - 2 * (vs**2/vp**2) * (drho/rho + 2 * dvs/vs)
    f = 0.5 * dvp/vp

    term1 = r0
    term2 = g * np.sin(theta1)**2
    term3 = f * (np.tan(theta1)**2 - np.sin(theta1)**2)

    if return_gradient:
        fields = ['intercept', 'gradient']
        Shuey = namedtuple('Shuey', fields)
        return Shuey(np.squeeze(r0), np.squeeze(g))
    elif terms:
        fields = ['R0', 'Rg', 'Rf']
        Shuey = namedtuple('Shuey', fields)
        return Shuey(np.squeeze([term1 for _ in theta1]),
                     np.squeeze(term2),
                     np.squeeze(term3)
                     )
    else:
        return np.squeeze(term1 + term2 + term3)


@deprecated('Please use shuey() instead.')
def shuey2(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0):
    """
    Compute Shuey approximation with 2 terms. Wraps `shuey()`. Deprecated,
    use `shuey()` instead.
    """
    r, g, _ = shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1=theta1, terms=True)
    return r + g


@deprecated('Please use shuey() instead.')
def shuey3(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    Compute Shuey approximation with 3 terms. Wraps `shuey()`. Deprecated,
    use `shuey()` instead.
    """
    return shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1=theta1)


@preprocess
@vectorize
def bortfeld(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    Compute Bortfeld approximation with three terms.
    http://sepwww.stanford.edu/public/docs/sep111/marie2/paper_html/node2.html
    Real numbers only.

    Args:
        vp1 (ndarray): The upper P-wave velocity; float or 1D array length m.
        vs1 (ndarray): The upper S-wave velocity; float or 1D array length m.
        rho1 (ndarray): The upper layer's density; float or 1D array length m.
        vp2 (ndarray): The lower P-wave velocity; float or 1D array length m.
        vs2 (ndarray): The lower S-wave velocity; float or 1D array length m.
        rho2 (ndarray): The lower layer's density; float or 1D array length m.
        theta1 (ndarray): The incidence angle; float or 1D array length n.
        terms (bool): Whether or not to return a tuple of the terms of the
            equation. The first term is the acoustic impedance.

    Returns:
        ndarray. The 3-term Bortfeld approximation for P-P reflectivity at the
            interface. Will be a float (for float inputs and one angle), a
            1 x n array (for float inputs and an array of angles), a 1 x m
            array (for float inputs and one angle), or an n x m array (for
            array inputs and an array of angles).
    """
    theta1 = np.real(theta1)

    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0
    k = (2 * vs/vp)**2
    rsh = 0.5 * (dvp/vp - k*drho/rho - 2*k*dvs/vs)

    # Compute three-term reflectivity
    term1 = 0.5 * (dvp/vp + drho/rho)
    term2 = rsh * np.sin(theta1)**2
    term3 = 0.5 * dvp/vp * np.tan(theta1)**2 * np.sin(theta1)**2

    if terms:
        fields = ['term1', 'term2', 'term3']
        Bortfeld = namedtuple('Bortfeld', fields)
        return Bortfeld(np.squeeze([term1 for _ in theta1]),
                        np.squeeze(term2),
                        np.squeeze(term3)
                        )
    else:
        return np.squeeze(term1 + term2 + term3)


@deprecated('Please use bortfeld() instead.')
def bortfeld2(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    The 2-term Bortfeld approximation for ava analysis. Wraps `shuey()`.
    Deprecated, use `bortfeld()` instead.

    Args:
        vp1 (ndarray): The upper P-wave velocity; float or 1D array length m.
        vs1 (ndarray): The upper S-wave velocity; float or 1D array length m.
        rho1 (ndarray): The upper layer's density; float or 1D array length m.
        vp2 (ndarray): The lower P-wave velocity; float or 1D array length m.
        vs2 (ndarray): The lower S-wave velocity; float or 1D array length m.
        rho2 (ndarray): The lower layer's density; float or 1D array length m.
        theta1 (ndarray): The incidence angle; float or 1D array length n.
        terms (bool): Whether or not to return a tuple of the terms of the
            equation. The first term is the acoustic impedance.

    Returns:
        ndarray. The 2-term Bortfeld approximation for P-P reflectivity at the
            interface. Will be a float (for float inputs and one angle), a
            1 x n array (for float inputs and an array of angles), a 1 x m
            array (for float inputs and one angle), or an n x m array (for
            array inputs and an array of angles).
    """
    theta1 = np.radians(theta1)
    theta2 = np.arcsin(vp2/vp1*np.sin(theta1))
    term1 = 0.5 * np.log((vp2*rho2*np.cos(theta1)) / (vp1*rho1*np.cos(theta2)))
    svp2 = (np.sin(theta1)/vp1)**2
    dvs2 = (vs1**2-vs2**2)
    term2 = svp2 * dvs2 * (2+np.log(rho2/rho1)/np.log(vs2/vs1))

    if terms:
        return term1, term2
    else:
        return (term1 + term2)


@deprecated('Please use bortfeld() instead.')
def bortfeld3(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    return bortfeld(vp1, vs1, rho1, vp2, vs2, rho2, theta1=theta1)


@preprocess
@vectorize
def hilterman(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    Not recommended, only seems to match Zoeppritz to about 10 deg.

    Hilterman (1989) approximation from Mavko et al. Rock Physics Handbook.
    According to Dvorkin: "arguably the simplest and a very convenient
    [approximation]." At least for small angles and small contrasts. Real
    numbers only.

    Args:
        vp1 (ndarray): The upper P-wave velocity; float or 1D array length m.
        vs1 (ndarray): The upper S-wave velocity; float or 1D array length m.
        rho1 (ndarray): The upper layer's density; float or 1D array length m.
        vp2 (ndarray): The lower P-wave velocity; float or 1D array length m.
        vs2 (ndarray): The lower S-wave velocity; float or 1D array length m.
        rho2 (ndarray): The lower layer's density; float or 1D array length m.
        theta1 (ndarray): The incidence angle; float or 1D array length n.
        terms (bool): Whether or not to return a tuple of the terms of the
            equation. The first term is the acoustic impedance.

    Returns:
        ndarray. The Hilterman approximation for P-P reflectivity at the
            interface. Will be a float (for float inputs and one angle), a
            1 x n array (for float inputs and an array of angles), a 1 x m
            array (for float inputs and one angle), or an n x m array (for
            array inputs and an array of angles).
    """
    theta1 = np.real(theta1)

    ip1 = vp1 * rho1
    ip2 = vp2 * rho2
    rp0 = (ip2 - ip1) / (ip2 + ip1)

    pr2, pr1 = moduli.pr(vp2, vs2), moduli.pr(vp1, vs1)
    pravg = (pr2 + pr1) / 2.
    pr = (pr2 - pr1) / (1 - pravg)**2.

    term1 = rp0 * np.cos(theta1)**2.
    term2 = pr * np.sin(theta1)**2.

    if terms:
        fields = ['term1', 'term2']
        Hilterman = namedtuple('Hilterman', fields)
        return Hilterman(np.squeeze(term1), np.squeeze(term2))
    else:
        return np.squeeze(term1 + term2)


def blangy(vp1, vs1, rho1, vp2, vs2, rho2,
           d1=0, e1=0, d2=0, e2=0,
           theta1=0):
    """Implements the Blangy equation with the same interface as the other
    reflectivity equations. Wraps bruges.anisotropy.blangy(), which you may
    prefer to use directly.

    Args:
        vp1 (ndarray): The upper P-wave velocity; float or 1D array length m.
        vs1 (ndarray): The upper S-wave velocity; float or 1D array length m.
        rho1 (ndarray): The upper layer's density; float or 1D array length m.
        vp2 (ndarray): The lower P-wave velocity; float or 1D array length m.
        vs2 (ndarray): The lower S-wave velocity; float or 1D array length m.
        rho2 (ndarray): The lower layer's density; float or 1D array length m.
        d1 (ndarray): The upper delta; float or 1D array length m.
        e1 (ndarray): The upper epsilon; float or 1D array length m.
        d2 (ndarray): The lower delta; float or 1D array length m.
        e2 (ndarray): The lower epsilon; float or 1D array length m.
        theta1 (ndarray): The incidence angle; float or 1D array length n.

    Returns:
        ndarray. The Blangy approximation for P-P reflectivity at the
            interface. Wraps `anisotropy.blangy()`.
    """
    _, anisotropic = anisotropy.blangy(vp1, vs1, rho1, d1, e1,  # UPPER
                                       vp2, vs2, rho2, d2, e2,  # LOWER
                                       theta1)
    return anisotropic
