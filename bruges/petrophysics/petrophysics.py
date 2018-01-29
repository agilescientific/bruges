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

        rho_matrix (sandstone) : 2650 kg/m^3
        rho_matrix (limestome): 2710 kg/m^3
        rho_matrix (dolomite): 2876 kg/m^3
        rho_matrix (anyhydrite): 2977 kg/m^3
        rho_matrix (salt): 20320 kg/m^3

        rho_fluid (fresh water): 1000 kg/m^3
        rho_fluid (salt water): 1100 kg/m^3

    See wiki.aapg.org/Density-neutron_log_porosity.

    Args:
        phi_rhob (ndarray): The density porosity log.
        rho_matrix (float)
        rho_fluid (float)

    Returns:
        Estimate of RHOB.
    """
    return rho_matrix * (1 - phi_rhob) + rho_fluid * phi_rhob


def density_to_porosity(density, rho_matrix, rho_fluid):
    """
    Get density from a porosity log. Typical values:

        rho_matrix (sandstone) : 2650 kg/m^3
        rho_matrix (limestome): 2710 kg/m^3
        rho_matrix (dolomite): 2876 kg/m^3
        rho_matrix (anyhydrite): 2977 kg/m^3
        rho_matrix (salt): 20320 kg/m^3

        rho_fluid (fresh water): 1000 kg/m^3
        rho_fluid (salt water): 1100 kg/m^3

    See wiki.aapg.org/Density-neutron_log_porosity.

    Args:
        density (ndarray): The density log or RHOB.
        rho_matrix (float)
        rho_fluid (float)

    Returns:
        Estimate of porosity as a volume fraction.
    """
    return (rho_matrix - density) / (rho_matrix - rho_fluid)


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
