# -*- coding: utf-8 -*-
"""
===================
petrophysics.py
===================

:copyright: 2018 Agile Geoscience
:license: Apache 2.0
"""


def gardner(vp, alpha=0.31, beta=0.25, fps=False):
    """
    Computes Gardner's density prediction from P-wave velocity.

    Args:
        vp (ndarray): P-wave velocity in m/s.
        alpha (float): The factor, 0.31 for m/s and 0.23 for fps.
        beta (float): The exponent, usually 0.25.
        fps (bool): Set to true for FPS and the equation will use the typical
            value for alpha.

    Returns:
        ndarray: RHOB estimate.
    """
    if fps:
        alpha = 0.23
    return alpha * vp ** beta


def porosity_to_density(phi_rhob, rho_matrix, rho_fluid):
    """
    Get density from a porosity log. Typical values:

        rho_matrix (sandstone) : 2650 g/cc
        rho_matrix (Limestome): 2710 g/cc
        rho_matrix (Dolomite): 2876 g/cc
        rho_matrix (Anyhydrite): 2977 g/cc
        rho_matrix (Salt): 20320 g/cc

        rho_fluid (fresh water): 1000 g/cc (is this more mud-like?)
        rho_fluid (salt water): 1100 g/cc

    See wiki.aapg.org/Density-neutron_log_porosity.

    Args:
        phi_rhob (ndarray): The density porosity log.
        rho_matrix (float)
        rho_fluid (float)

    Returns:
        Estimate of RHOB.
    """
    return rho_matrix*(1 - phi_rhob) + rho_fluid*phi_rhob


def slowness_to_velocity(slowness):
    """
    Convert a DT or DTS log in µs per unit depth, to velocity in unit depth
    per second.

    Args:
        slowness (ndarray): A DT value or sequence of values.

    Returns:
        ndarray: The velocity.
    """
    return 1e6 / np.array(slowness)
