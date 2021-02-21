"""
Elastic impedance.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np


def elastic_impedance(vp, vs, rho, theta1,
                      k=None,
                      normalize=False,
                      constants=None,
                      use_sin=False,
                      rho_term=False):
    """
    Returns the elastic impedance as defined by Connolly, 1999; we are using
    the formulation reported in Whitcombe et al. (2001). Inputs should be
    shape m x 1, angles should be n x 1. The result will be m x n.

    Args:
        vp (ndarray): The P-wave velocity scalar or 1D array.
        vs (ndarray): The S-wave velocity scalar or 1D array.
        rho (ndarray): The bulk density scalar or 1D array.
        theta1 (ndarray): Incident angle(s) in degrees, scalar or array.
        k (float): A constant, see Connolly (1999). Default is None, which
            computes it from Vp and Vs.
        normalize (bool): if True, returns the normalized form of Whitcombe.
        constants (tuple): A sequence of 3 constants to use for normalization.
            If you don't provide this, then normalization just uses the means
            of the inputs. If these are scalars, the result will be the
            acoustic impedance (see Whitcombe et al., 2001).
        use_sin (bool): If True, use sin^2 for the first term, instead of
            tan^2 (see Connolly).
        rho_term (bool): If True, provide alternative form, with Vp factored
            out; use in place of density in generating synthetics in other
            software (see Connolly). In other words, the result can be
            multipled with Vp to get the elastic impedance.

    Returns:
        ndarray: The elastic impedance log at the specficied angle or angles.
    """
    theta1 = np.radians(theta1).reshape(-1, 1)
    if (np.nan_to_num(theta1) > np.pi/2.).any():
        raise ValueError("Incidence angle theta1 must be less than 90 deg.")

    alpha = np.asanyarray(vp, dtype=float)
    beta = np.asanyarray(vs, dtype=float)
    rho = np.asanyarray(rho, dtype=float)
    op = np.sin if use_sin else np.tan
    k = np.mean(beta**2.0 / alpha**2.0) if k is None else k

    a = 1 + op(theta1)**2.0
    b = -8 * k * np.sin(theta1)**2.0
    c = 1 - 4 * k * np.sin(theta1)**2.0

    ei = alpha**a * beta**b * rho**c

    if normalize:
        if constants is None:
            # Use the means; this will return acoustic impedance for scalars.
            alpha_0 = np.nanmean(vp)
            beta_0 = np.nanmean(vs)
            rho_0 = np.nanmean(rho)
        else:
            try:
                alpha_0, beta_0, rho_0 = constants
            except ValueError as e:
                raise ValueError("You must provide a sequence of 3 constants.")
        ei *= alpha_0**(1 - a) * beta_0**(-b) * rho_0**(1 - c)

    if rho_term:
        ei /= alpha

    if ei.size == 1:
        return np.asscalar(ei)
    else:
        return np.squeeze(ei.T)
