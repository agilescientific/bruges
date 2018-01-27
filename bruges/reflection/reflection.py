# -*- coding: utf-8 -*-
"""
Various reflectivity algorithms.

:copyright: 2018 Agile Geoscience
:license: Apache 2.0
"""
from functools import wraps
from collections import namedtuple

import numpy as np
from numpy import tan, sin, cos

from bruges.rockphysics import moduli
from bruges.rockphysics import anisotropy
from bruges.util import deprecated


def reflectivity(vp, vs, rho, theta, method='zoeppritz_rpp'):
    """
    Compute 'upper' and 'lower' intervals from the three provided arrays,
    then pass the result to the specified method to compute reflection
    coefficients.

    :param vp: The P-wave velocity 1D array.
    :param vs: The S-wave velocity 1D array.
    :param rho: The density 1D array.
    :param theta: Incidence angle, a scalar or 1D array [degrees].
    :param method: Which method to use:

        - 'zoeppritz_rpp': zoeppritz_rpp,
        - 'akirichards': akirichards,
        - 'akirichards_alt': akirichards_alt,
        - 'fatti': fatti,
        - 'shuey': shuey,
        - 'bortfeld': bortfeld,
        - 'hilterman': hilterman,

    :returns: the result of running the specified method on the data.
    """
    methods = {
        'zoeppritz_rpp': zoeppritz_rpp,
        'akirichards': akirichards,
        'akirichards_alt': akirichards_alt,
        'fatti': fatti,
        'shuey': shuey,
        'bortfeld': bortfeld,
        'hilterman': hilterman,
    }
    func = methods[method]
    vp = np.array(vp).astype(float)
    vs = np.array(vs).astype(float)
    rho = np.array(rho).astype(float)
    vp1, vp2 = vp[:-1], vp[1:]
    vs1, vs2 = vs[:-1], vs[1:]
    rho1, rho2 = rhob[:-1], rhob[1:]
    return func(vp1, vs1, rho1, vp2, vs2, rho2, theta)


def vectorize(func):
    """
    Decorator to make sure the inputs are arrays. We also add a dimension
    to theta to make the functions work in an 'outer product' way.
    """
    @wraps(func)
    def wrapper(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, **kwargs):
        vp1 = np.array(vp1).astype(float)
        vp2 = np.array(vp2).astype(float)
        vs1 = np.array(vs1).astype(float)
        vs2 = np.array(vs2).astype(float)
        rho1 = np.array(rho1).astype(float)
        rho2 = np.array(rho2).astype(float)
        theta1 = np.array(theta1).reshape((-1, 1))
        return func(vp1, vs1, rho1, vp2, vs2, rho2, theta1, **kwargs)
    return wrapper


def scattering_matrix(vp1, vs1, rho1, vp0, vs0, rho0, theta1=0):
    '''
    Full Zoeppritz solution, considered the definitive solution.
    Calculates the angle dependent p-wave reflectivity of an interface
    between two mediums.

    Originally written by: Wes Hamlyn, vectorized by Agile.

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp0: The p-wave velocity of the lower medium.
    :param vs0: The s-wave velocity of the lower medium.
    :param rho0: The density of the lower medium.

    :param theta1: A scalar  [degrees].

    :returns: a 4x4 array representing the scattering matrix
                at the incident angle theta1.
    '''
    theta1 = np.radians(np.array(theta1))
    if theta1.size == 1:
        theta1 = np.expand_dims(theta1, axis=1)

    p = sin(theta1) / vp1  # Ray parameter.

    # Calculate reflection & transmission angles for Zoeppritz.
    theta2 = np.arcsin(p * vp0)  # Trans. angle of P-wave.
    phi1 = np.arcsin(p * vs1)    # Refl. angle of converted S-wave.
    phi2 = np.arcsin(p * vs0)    # Trans. angle of converted S-wave.

    # Matrix form of Zoeppritz Equations... M & N are matrices.
    M = np.array([[-sin(theta1), -cos(phi1), sin(theta2), cos(phi2)],
                  [cos(theta1), -sin(phi1), cos(theta2), -sin(phi2)],
                  [2 * rho1 * vs1 * sin(phi1) * cos(theta1),
                   rho1 * vs1 * (1 - 2 * sin(phi1) ** 2),
                   2 * rho0 * vs0 * sin(phi2) * cos(theta2),
                   rho0 * vs0 * (1 - 2 * sin(phi2) ** 2)],
                  [-rho1 * vp1 * (1 - 2 * sin(phi1) ** 2),
                   rho1 * vs1 * sin(2 * phi1),
                   rho0 * vp0 * (1 - 2 * sin(phi2) ** 2),
                   -rho0 * vs0 * sin(2 * phi2)]], dtype='float')

    N = np.array([[sin(theta1), cos(phi1), -sin(theta2), -cos(phi2)],
                  [cos(theta1), -sin(phi1), cos(theta2), -sin(phi2)],
                  [2 * rho1 * vs1 * sin(phi1) * cos(theta1),
                   rho1 * vs1 * (1 - 2 * sin(phi1) ** 2),
                   2 * rho0 * vs0 * sin(phi2) * cos(theta2),
                   rho0 * vs0 * (1 - 2 * sin(phi2) ** 2)],
                  [rho1 * vp1 * (1 - 2 * sin(phi1) ** 2),
                   -rho1 * vs1 * sin(2 * phi1),
                   - rho0 * vp0 * (1 - 2 * sin(phi2) ** 2),
                   rho0 * vs0 * sin(2 * phi2)]], dtype='float')

    A = np.linalg.inv(np.rollaxis(M, 2))
    Z = np.matmul(A, np.rollaxis(N, -1))

    return np.rollaxis(Z, 0, 3)


def zoeppritz_element(vp1, vs1, rho1, vp0, vs0, rho0, theta1=0, element='PdPu'):
    """
    Returns any mode reflection coefficients from the Zoeppritz
    scattering matrix. Pass in the mode as element, e.g. 'PdSu' for PS.

    Wraps scattering_matrix().

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp0: The p-wave velocity of the lower medium.
    :param vs0: The s-wave velocity of the lower medium.
    :param rho0: The density of the lower medium.

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    """
    elements = np.array([['PdPu', 'SdPu', 'PuPu', 'SuPu'],
                         ['PdSu', 'SdSu', 'PuSu', 'SuSu'],
                         ['PdPd', 'SdPd', 'PuPd', 'SuPd'],
                         ['PdSd', 'SdSd', 'PuSd', 'SuSd']])

    Z = scattering_matrix(vp1, vs1, rho1, vp0, vs0, rho0, theta1)

    return np.squeeze(Z[np.where(elements == element)])


def zoeppritz(vp1, vs1, rho1, vp0, vs0, rho0, theta1=0):
    '''
    Returns the PP reflection coefficients from the Zoeppritz
    scattering matrix. Wraps zoeppritz_element().

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp0: The p-wave velocity of the lower medium.
    :param vs0: The s-wave velocity of the lower medium.
    :param rho0: The density of the lower medium.

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    '''
    return zoeppritz_element(vp1, vs1, rho1, vp0, vs0, rho0, theta1, 'PdPu')


@vectorize
def zoeppritz_rpp(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    Exact Zoeppritz from expression.

    This is useful because we can pass arrays to it, which we can't do to
    scattering_matrix().

    Dvorkin et al. (2014). Seismic Reflections of Rock Properties. Cambridge.

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp0: The p-wave velocity of the lower medium.
    :param vs0: The s-wave velocity of the lower medium.
    :param rho0: The density of the lower medium.

    :returns: a matrix of length `len(vp1)` whose rows are of length `len(vp[RHOB])`.
    """
    theta1 = np.radians(theta1)

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

    return np.squeeze(rpp).T


@vectorize
def akirichards(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    This is the formulation from Avseth et al.,
    Quantitative seismic interpretation,
    Cambridge University Press, 2006. Adapted for a 4-term formula.
    See http://subsurfwiki.org/wiki/Aki-Richards_equation

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp2: The p-wave velocity of the lower medium.
    :param vs2: The s-wave velocity of the lower medium.
    :param rho2: The density of the lower medium.

    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle
    """
    theta1 = np.radians(theta1)
    if (np.nan_to_num(theta1) > np.pi/2.).any():
        raise ValueError("Incidence angle theta1 must be less than 90 deg.")

    # critical_angle = arcsin(vp1/vp2)
    theta2 = np.arcsin(vp2/vp1*sin(theta1))
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
    term1 = [w for _ in theta1]
    term2 = -1 * x * np.sin(theta1)**2
    term3 = y / np.cos(meantheta)**2
    term4 = -1 * z * np.sin(theta1)**2

    if terms:
        fields = ['term1', 'term2', 'term3', 'term4']
        AkiRichards = namedtuple('AkiRichards', fields)
        return AkiRichards(np.squeeze(term1),
                           np.squeeze(term2),
                           np.squeeze(term3),
                           np.squeeze(term4)
                          )
    else:
        return np.squeeze(term1 + term2 + term3 + term4).T


@vectorize
def akirichards_alt(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    This is another formulation of the Aki-Richards solution.
    See http://subsurfwiki.org/wiki/Aki-Richards_equation

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp2: The p-wave velocity of the lower medium.
    :param vs2: The s-wave velocity of the lower medium.
    :param rho2: The density of the lower medium.

    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    """
    theta1 = np.radians(theta1)

    # critical_angle = arcsin(vp1/vp2)
    theta2 = np.arcsin(vp2/vp1*sin(theta1))
    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    theta = (theta1+theta2)/2.0
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0

    # Compute the three terms
    term1 = 0.5 * (dvp/vp + drho/rho)
    term2 = (0.5*dvp/vp-2*(vs/vp)**2*(drho/rho+2*dvs/vs)) * sin(theta)**2
    term3 = 0.5 * dvp/vp * (tan(theta)**2 - sin(theta)**2)

    if terms:
        fields = ['term1', 'term2', 'term3']
        AkiRichards = namedtuple('AkiRichards', fields)
        return AkiRichards(np.squeeze(term1),
                           np.squeeze(term2),
                           np.squeeze(term3)
                           )
    else:
        return np.squeeze(term1 + term2 + term3).T


@vectorize
def fatti(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    Compute reflectivities with Fatti's formulation of the
    Aki-Richards equation, which does not account for the
    critical angle. Fatti et al (1994), Geophysics 59 (9).

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp2: The p-wave velocity of the lower medium.
    :param vs2: The s-wave velocity of the lower medium.
    :param rho2: The density of the lower medium.

    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    """
    theta1 = np.radians(theta1)

    drho = rho2-rho1
    rho = (rho1+rho2) / 2.0
    vp = (vp1+vp2) / 2.0
    vs = (vs1+vs2) / 2.0
    dip = (vp2*rho2 - vp1*rho1)/(vp2*rho2 + vp1*rho1)
    dis = (vs2*rho2 - vs1*rho1)/(vs2*rho2 + vs1*rho1)
    d = drho/rho

    # Compute the three terms
    term1 = (1 + tan(theta1)**2) * dip
    term2 = -8 * (vs/vp)**2 * dis * sin(theta1)**2
    term3 = -1 * (0.5 * tan(theta1)**2 - 2 * (vs/vp)**2 * sin(theta1)**2) * d

    if terms:
        fields = ['term1', 'term2', 'term3']
        Fatti = namedtuple('Fatti', fields)
        return Fatti(np.squeeze(term1),
                     np.squeeze(term2),
                     np.squeeze(term3)
                     )
    else:
        return np.squeeze(term1 + term2 + term3).T


@vectorize
def shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0,
          terms=False,
          return_gradient=False):
    """
    Compute Shuey approximation with 3 terms.
    http://subsurfwiki.org/wiki/Shuey_equation

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.
    :param vp2: The p-wave velocity of the lower medium.
    :param vs2: The s-wave velocity of the lower medium.
    :param rho2: The density of the lower medium.
    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].
    :param terms: bool. Whether to return a tuple of the 3 individual terms.
    :param return_gradient: bool. Whether to return a tuple of the intercept
                            and gradient (i.e. the second term divided by
                            sin^2(theta).

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    """
    theta1 = np.radians(theta1)

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
        return Shuey(np.squeeze(term1),
                     np.squeeze(term2),
                     np.squeeze(term3)
                     )
    else:
        return np.squeeze(term1 + term2 + term3).T


@deprecated('Please use shuey() instead.')
def shuey2(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0):
    """
    Compute Shuey approximation with 2 terms.
    """
    r, g, _ = shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1=theta1, terms=True)
    return r + g


@deprecated('Please use shuey() instead.')
def shuey3(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    Compute Shuey approximation with 3 terms.
    """
    return shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1=theta1)


@deprecated('Please use bortfeld() instead.')
def bortfeld2(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    The 2-term Bortfeld approximation for ava analysis.

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp2: The p-wave velocity of the lower medium.
    :param vs2: The s-wave velocity of the lower medium.
    :param rho2: The density of the lower medium.

    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    """
    theta1 = np.radians(theta1)
    # Bortfeld only needs one extra parameter
    theta2 = np.arcsin(vp2/vp1*np.sin(theta1))  # radians

    # This breaks if theta = 90 deg
    term1 = 0.5 * np.log((vp2*rho2*cos(theta1)) / (vp1*rho1*cos(theta2)))

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


@vectorize
def bortfeld(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    Compute Bortfeld approximation with three terms.
    http://sepwww.stanford.edu/public/docs/sep111/marie2/paper_html/node2.html

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp2: The p-wave velocity of the lower medium.
    :param vs2: The s-wave velocity of the lower medium.
    :param rho2: The density of the lower medium.

    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    """
    theta1 = np.radians(theta1)

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
        return Bortfeld(np.squeeze(term1),
                        np.squeeze(term2),
                        np.squeeze(term3)
                        )
    else:
        return np.squeeze(term1 + term2 + term3).T


@vectorize
def hilterman(vp1, vs1, rho1, vp2, vs2, rho2, theta1=0, terms=False):
    """
    Not recommended, only seems to match Zoeppritz to about 10 deg.

    Hilterman (1989) approximation from Mavko et al. Rock Physics Handbook.
    According to Dvorkin: "arguably the simplest and a very convenient
    [approximation]." At least for small angles and small contrasts.

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp2: The p-wave velocity of the lower medium.
    :param vs2: The s-wave velocity of the lower medium.
    :param rho2: The density of the lower medium.

    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    """
    theta1 = np.radians(theta1)

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
        return np.squeeze(term1 + term2).T


def blangy(vp1, vs1, rho1, vp2, vs2, rho2,
           d1=0, e1=0, d2=0, e2=0,
           theta1=0, terms=False):
    """Implements the Blangy equation with the same interface as the other
    reflectivity equations. Wraps bruges.anisotropy.blangy().

    Note that the anisotropic parameters come after the other rock properties,
    and all default to zero.

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.
    :param vp2: The p-wave velocity of the lower medium.
    :param vs2: The s-wave velocity of the lower medium.
    :param rho2: The density of the lower medium.
    :param d1: Thomsen's delta for the upper medium.
    :param e1: Thomsen's epsilon for the upper medium.
    :param d2: Thomsen's delta for the upper medium.
    :param e2: Thomsen's epsilon for the upper medium.
    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    """
    _, anisotropic = anisotropy.blangy(vp1, vs1, rho1, d1, e1,  # UPPER
                                       vp2, vs2, rho2, d2, e2,  # LOWER
                                       theta1)
    return anisotropic
