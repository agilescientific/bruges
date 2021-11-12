# -*- coding: utf-8 -*-
"""
===================
petrophysics.py
===================

:copyright: 2018 Agile Geoscience
:license: Apache 2.0
"""


def gardner(vp, alpha=310, beta=0.25, fps=False):
    """
    Computes Gardner's density prediction from P-wave velocity.

    Args:
        vp (ndarray): P-wave velocity in m/s.
        alpha (float): The factor, 310 for m/s and 230 for ft/s.
        beta (float): The exponent, usually 0.25.
        fps (bool): Set to true for FPS and the equation will use the typical
            value for alpha. Overrides value for alpha, so if you want to use
            you own alpha, regardless of units, set this to False.

    Returns:
        ndarray: RHOB estimate in kg/m^3.
    """
    alpha = 230 if fps else alpha
    return alpha * vp ** beta


def inverse_gardner(rho, alpha=310, beta=0.25, fps=False):
    """
    Computes Gardner's density prediction from P-wave velocity.

    Args:
        rho (ndarray): Density in kg/m^3.
        alpha (float): The factor, 310 for m/s and 230 for fps.
        beta (float): The exponent, usually 0.25.
        fps (bool): Set to true for FPS and the equation will use the typical
            value for alpha. Overrides value for alpha, so if you want to use
            your own alpha, regardless of units, set this to False.

    Returns:
        ndarray: Vp estimate in m/s.
    """
    alpha = 230 if fps else alpha
    exponent = 1 / beta
    factor = 1 / alpha**exponent
    return factor * rho**exponent


velocity_to_density = gardner


density_to_velocity = inverse_gardner


def porosity_to_density(phi, rho_matrix, rho_fluid):
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
        phi (ndarray): The porosity log.
        rho_matrix (float)
        rho_fluid (float)

    Returns:
        Estimate of bulk density, rho.
    """
    return rho_matrix * (1 - phi) + rho_fluid * phi_rhob


def density_to_porosity(rho, rho_matrix, rho_fluid):
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
        rho (ndarray): The bulk density log or RHOB.
        rho_matrix (float)
        rho_fluid (float)

    Returns:
        Estimate of porosity as a volume fraction.
    """
    return (rho_matrix - rho) / (rho_matrix - rho_fluid)


def slowness_to_velocity(slowness):
    """
    Convert a slowness log in µs per unit depth, to velocity in unit depth
    per second.

    Args:
        slowness (ndarray): A value or sequence of values.

    Returns:
        ndarray: The velocity.
    """
    return 1e6 / np.array(slowness)


velocity_to_slowness = slowness_to_velocity


def gardner_param(vp,rhob):
     """
     Finds optimal alpha and beta parameters for the gardner

     Volodymyr Vragov, October 2018


     Args:
         rho(ndarray): Density.
         vp (ndarray): P-wave velocity.

     Returns:
         alpha (float): The factor.
         beta (float): The exponent, usually 0.25.
     """
     params, _ = curve_fit(optimizer_gardner, rhob, vp)
     return params[0], params[1]
     

def error_flag(pred, actual, dev = 1.0, method = 1):
    """Calculate the difference between a predicted and an actual curve 
    and return a log flagging large differences based on a user-defined distance 
    (in standard deviation units) from the mean difference

    Matteo Niccoli, October 2018
    
    Args:
        predicted (ndarray) = predicted log
        actual (ndarray) =  original log  
        dev  (float) = standard deviations to use, default 1
        error calcluation method (int), default 1
            1: difference between curves larger than mean difference plus dev
            2: curve slopes have opposite sign
            3: curve slopes of opposite sign OR difference larger than mean plus dev
    

    Returns:
    flag (ndarray) =  error flag curve
    """   
    flag = np.zeros(len(pred))
    err = np.abs(pred-actual)
    err_mean = np.mean(err)
    err_std = np.std(err)
    ss = np.sign(pred.diff().fillna(pred))
    ls = np.sign(actual.diff().fillna(actual))

    if method == 1:
        flag[np.where(err>(err_mean + (dev*err_std)))] = 1
    elif method == 2:      
        flag[np.where((ss + ls)==0 )]= 1
    elif method == 3:
        flag[np.where(np.logical_or(err>(err_mean + (dev*err_std)), (ss+ls)==0 ))]= 1
    return flag


def optimize_inverse_gardner(rho, alpha, beta):
     """ Wrapper function to pass inverse_gardner to scipy.curve_fit
     to get optimal alpha and beta parameters


     Matteo Niccoli and Volodymyr Vragov, October 2018

     Args:
         rho(ndarray): Density.
         alpha (float): The factor.
         beta (float): The exponent.

     Returns:
         inverse_gardner: this is passed to scipy.curve_fit as, for example:
         popt_synt, pcov = scipy.curve_fit(optimize_inverse_gardner, rho, vp)
         For a ull example, please read:
         mycarta.wordpress.com/2018/10/28/geophysics-python-sprint-2018-day-2-and-beyond-part-i
     """
     return inverse_gardner(rho, alpha=alpha, beta=beta)


