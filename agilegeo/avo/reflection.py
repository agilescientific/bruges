#!/usr/bin/env python
#
#   Python application for calculating angle-dependent p-wave reflectivity
#   using Zoeppritz equations & various approximations.
#
#   Originally written to give insight into limitations of the Zoeppritz
#   approximations and to get more familiar with GUI programming using wxPython
#   and Matplotlib.
#
#   Requires:   Python (2.6 or 2.7)
#               wxPython
#               Numpy & Matplotlib
#
#       Written by: Wes Hamlyn
#       Modified by: Sean Ross-Ross
#       Last Mod:   May 1, 2012
#
#   Use for whatever you like but at your own risk...
#

import numpy as np
from numpy import log, tan, sin, cos, arcsin, arccosh, radians, \
                  degrees

def zoeppritz(vp1, vs1, rho1, vp0, vs0, rho0, theta1):
    '''
    Full Zoeppritz solution, considered the definitive solution.
    Calculates the angle dependent p-wave reflectivity of an interface
    between two mediums.

    :param vp1: The p-wave velocity of the upper medium.
    :param vs1: The s-wave velocity of the upper medium.
    :param rho1: The density of the upper medium.

    :param vp0: The p-wave velocity of the lower medium.
    :param vs0: The s-wave velocity of the lower medium.
    :param rho0: The density of the lower medium.
    
    :param theta1: An array of incident angles to use for reflectivity
                   calculation [degrees].

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    '''
    
    # Set the ray paramter, p
    p = sin(radians(theta1)) / vp1 # ray parameter
    
    # Calculate reflection & transmission angles for Zoeppritz
    theta1 = radians(theta1)   # Convert theta1 to radians
    theta2 = arcsin(p * vp0);      # Trans. angle of P-wave
    phi1 = arcsin(p * vs1);      # Refl. angle of converted S-wave
    phi2 = arcsin(p * vs0);      # Trans. angle of converted S-wave
                
    # Matrix form of Zoeppritz Equations... M & N are matricies
    M = np.array([ \
        [-sin(theta1), -cos(phi1), sin(theta2), cos(phi2)], \
        [cos(theta1), -sin(phi1), cos(theta2), -sin(phi2)], \
        [2 * rho1 * vs1 * sin(phi1) * cos(theta1), rho1 * vs1 *\
            (1 - 2 * sin(phi1) ** 2), \
            2 * rho0 * vs0 * sin(phi2) * cos(theta2), \
            rho0 * vs0 * (1 - 2 * sin(phi2) ** 2)], \
        [-rho1 * vp1 * (1 - 2 * sin(phi1) ** 2), rho1 * vs1 * \
             sin(2 * phi1), \
            rho0 * vp0 * (1 - 2 * sin(phi2) ** 2), -rho0 * \
              vs0 * sin(2 * phi2)]
        ], dtype='float')
    
    N = np.array([ \
        [sin(theta1), cos(phi1), -sin(theta2), -cos(phi2)], \
        [cos(theta1), -sin(phi1), cos(theta2), -sin(phi2)], \
        [2 * rho1 * vs1 * sin(phi1) * cos(theta1), rho1 * vs1 * (1 - 2 * sin(phi1) ** 2), \
            2 * rho0 * vs0 * sin(phi2) * cos(theta2), rho0 * vs0 * (1 - 2 * sin(phi2) ** 2)], \
        [rho1 * vp1 * (1 - 2 * sin(phi1) ** 2), -rho1 * vs1 * sin(2 * phi1), \
            - rho0 * vp0 * (1 - 2 * sin(phi2) ** 2), rho0 * vs0 * sin(2 * phi2)]\
        ], dtype='float')
    
    # This is the important step, calculating
    # coefficients for all modes
    # and rays result is a 4x4 matrix, we want the R[0][0]
    # element for Rpp reflectivity only
    
    if M.ndim == 3:
        zoep = np.zeros([M.shape[-1]])
        for i in range(M.shape[-1]): 
            Mi = M[..., i]
            Ni = N[..., i]
            dt = np.dot(np.linalg.inv(Mi), Ni)
            zoep[i] = dt[0][0]
    else:
        dt = np.dot(np.linalg.inv(M), N)
        zoep = dt[0][0]
    
    return zoep

def akirichards(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
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
             value corresponding to each angle.
    
    """

    # We are not using this for anything, but will
    # critical_angle = arcsin(vp1/vp2)
    
    # Do we need to ensure that we get floats out before
    # computing sines?
    vp1 = float(vp1)

    theta2 = degrees(arcsin(vp2/vp1*sin(radians(theta1))))

    # Compute the various parameters
    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    meantheta = (theta1+theta2)/2.0
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0

    # Compute the coefficients 
    w = 0.5 * drho/rho
    x = 2 * (vs/vp1)**2 * drho/rho
    y = 0.5 * (dvp/vp)
    z = 4 * (vs/vp1)**2 * (dvs/vs)

    # Compute the terms
    term1 = w
    term2 = -1 * x * sin(radians(theta1))**2
    term3 = y / cos(radians(meantheta))**2
    term4 = -1 * z * sin(radians(theta1))**2
    
    return (term1 + term2 + term3 + term4)

def akirichards_alt(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
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

    # We are not using this for anything, but will
    # critical_angle = arcsin(vp1/vp2)
    
    # Do we need to ensure that we get floats out before
    # computing sines?
    vp1 = float(vp1)    

    theta2 = degrees(arcsin(vp2/vp1*sin(radians(theta1))))

    # Compute the various parameters
    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    theta = (theta1+theta2)/2.0
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0

    # Compute the three terms
    term1 = 0.5*(dvp/vp + drho/rho)
    term2 = (0.5*dvp/vp-2*(vs/vp)**2*(drho/rho+2*dvs/vs))*\
            sin(radians(theta))**2
    term3 = 0.5*dvp/vp*(tan(radians(theta))**2 - \
                        sin(radians(theta))**2)
    
    return (term1 + term2 + term3)

def fatti(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
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
    
    # Do we need to ensure that we get floats out before computing sines?
    vp1 = float(vp1)

    # Compute the various parameters
    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0
    dip = (vp2*rho2 - vp1*rho1)/(vp2*rho2 + vp1*rho1)
    dis = (vs2*rho2 - vs1*rho1)/(vs2*rho2 + vs1*rho1)

    # Compute the three terms
    term1 = (1 + tan(radians(theta1))**2) * dip
    term2 = -8 * (vs/vp)**2 * sin(radians(theta1))**2 * dis
    term3 = -1 * ( 0.5 * tan(radians(theta1))**2 - 2 * \
                   (vs/vp)**2 * sin(radians(theta1))**2 ) * drho/rho
    
    return (term1 + term2 + term3)

def shuey2(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    """
    Compute Shuey approximation with 2 terms.
    http://subsurfwiki.org/wiki/Shuey_equation

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

    # Compute some parameters
    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0

    # Compute two-term reflectivity

    r0 = 0.5 * (dvp/vp + drho/rho)
    g  = 0.5 * dvp/vp - 2 * (vs**2/vp**2) * ( drho/rho + 2 * dvs/vs)

    return r0 + g * sin(radians(theta1))**2

def shuey3(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
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

    :returns: a vector of len(theta1) containing the reflectivity
             value corresponding to each angle.
    
    """

    # Compute some parameters
    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0

    # Compute three-term reflectivity
    
    r0 = 0.5 * (dvp/vp + drho/rho)
    g  = 0.5 * dvp/vp - 2 * (vs**2/vp**2) * ( drho/rho + 2 * dvs/vs)
    f = 0.5 * dvp/vp

    return r0 + g * sin(radians(theta1))**2 + f * \
       (tan(radians(theta1))**2 - sin(radians(theta1))**2)

def bortfeld2(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
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
    
    # Bortfeld only needs one extra parameter
    theta2 = degrees(arcsin(vp2/vp1*sin(radians(theta1))))

    # This breaks if theta = 90 deg
    term1 = 0.5 * log( (vp2*rho2*cos(radians(theta1)))/ \
                       (vp1*rho1*cos(radians(theta2))) )
                       
    term2 = (sin(radians(theta1))/vp1)**2 * (vs1**2-vs2**2) * \
             (2+log(rho2/rho1)/log(vs2/vs1))

    return term1 + term2

def bortfeld3(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
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

    # Compute some parameters
    drho = rho2-rho1
    dvp = vp2-vp1
    dvs = vs2-vs1
    rho = (rho1+rho2)/2.0
    vp = (vp1+vp2)/2.0
    vs = (vs1+vs2)/2.0
    k = (2 * vs/vp)**2
    rsh = 0.5 * (dvp/vp - k*drho/rho -2*k*dvs/vs)

    # Compute three-term reflectivity

    term1 = 0.5 * (dvp/vp + drho/rho)
    term2 = rsh * sin(radians(theta1))**2
    term3 = 0.5 * dvp/vp * tan(radians(theta1))**2 * \
            sin(radians(theta1))**2

    return term1 + term2 + term3
