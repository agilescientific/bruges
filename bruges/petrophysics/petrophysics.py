# -*- coding: utf-8 -*-
'''
===================
petrophysics.py
===================

:copyright: 2016 Agile Geoscience
:license: Apache 2.0
'''


def gardner(vp, alpha=0.31, beta=0.25):
    '''
    Computes Gardner's density prediction.

    SI units only.

    Args:
        vp (this is actually slowness in ft / sec)

    Returns:
        rhob estimate.
    '''
    return alpha * vp ** beta


def convert_rhob_to_SI(rhob):
    return rhob*1000


def rhob(phi_rhob, Rho_matrix= 2650.0, Rho_fluid=1000.0):
    """
    Rho_matrix (sandstone) : 2650 g/cc
    Rho_matrix (Limestome): 2710 g/cc
    Rho_matrix (Dolomite): 2876 g/cc
    Rho_matrix (Anyhydrite): 2977 g/cc
    Rho_matrix (Salt): 20320 g/cc

    Rho_fluid (fresh water): 1000 g/cc (is this more mud-like?)
    Rho_fluid (salt water): 1100 g/cc
    see wiki.aapg.org/Density-neutron_log_porosity
    returns density porosity log """
    
    return Rho_matrix*(1 - phi_rhob) + Rho_fluid*phi_rhob

def slowness_to_velocity(sonic, imperial=True, nansub=-999.25):
    """
    This only works for sonic in us/ft, always converts to SI
    sonic :a numpy array of slowness (us/ft or us/m)
    """
    arr = sonic.as_matrix() 
    if imperial:
        arr[arr == nansub] = np.nan
        arr *= 3.2808  # units in us/m
        arr /= 1e6     # units in s/m
        velocity = 1/arr

    else:
        arr[arr == nansub] = np.nan
        arr /= 1e6     # units in s/m
        velocity = 1/sonic
    return velocity