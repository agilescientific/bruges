"""
Fluid properties.

These functions implement equations from Batzle and Wang (1992), seismic
properties of pore fluids. GEOPHYSICS, VOL. 57, NO. 11; P. 1396-1408,

:copyright: 2018 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np

def wood(Kf1, Kf2, Sf1):
    """
    Wood's equation, per equation 35b in Batzle and Wang (1992).
    """
    return 1 / ((Sf1 / Kf1) + ((1 - Sf1) / Kf2))


def rho_water(temperature, pressure):
    """
    The density of pure water, as a function of temperature and pressure.
    Implements equation 27a from Batzle and Wang (1992).

    Use scalars or arrays; if you use arrays, they must be the same size.

    Args:
        temperature (array): The temperature in degrees Celsius.
        pressure (array): The pressure in pascals.

    Returns:
        array: The density in kg/m3.
    """
    # Align with the symbols and units in Batzle & Wang.
    T, P = np.asanyarray(temperature), np.asanyarray(pressure)*1e-6

    x = -80*T - 3.3*T**2 + 0.00175*T**3 + 489*P - 2*T*P \
        + 0.016*P*T**2 - 1.3e-5*P*T**3 - 0.333*P**2 - 0.002*T*P**2

    return 1000 + 1e-3 * x


def rho_brine(temperature, pressure, salinity):
    """
    The density of NaCl brine, given temperature, pressure, and salinity.
    The density of pure water is computed from rho_water(). Implements
    equation 27b from Batzle and Wang (1992).

    Use scalars or arrays; if you use arrays, they must be the same size.

    Args:
        temperature (array): The temperature in degrees Celsius.
        pressure (array): The pressure in pascals.
        salinity (array): The weight fraction of NaCl, e.g. 35e-3
            for 35 parts per thousand, or 3.5% (the salinity of
            seawater).
    Returns:
        array: The density in kg/m3.
    """
    # Align with the symbols and units in Batzle & Wang.
    T, P = np.asanyarray(temperature), np.asanyarray(pressure)*1e-6
    S = np.asanyarray(salinity)

    rho_w = rho_water(temperature, pressure) / 1000
    x = 300*P - 2400*P*S + T*(80 + 3*T - 3300*S - 13*P + 47*P*S)

    return rho_w + S*(0.668 + 0.44*S + 1e-6 * x)


def rho_gas(temperature, pressure, molecular_weight):
    """
    Gas density, given temperature (in deg C), pressure (in Pa), and molecular
    weight.

    Args:
        temperature (array): The temperature in degrees Celsius.
        pressure (array): The pressure in pascals.
        molecular_weight (array): The molecular weight.
    Returns:
        array: The density in kg/m3.
    """
    # Align with the symbols and units in Batzle & Wang.
    T, P = np.asanyarray(temperature), np.asanyarray(pressure)*1e-6
    M = np.asanyarray(molecular_weight)
    R = 8.3144598

    return M * P / (R * (T + 273.15))


def v_water(temperature, pressure):
    """
    The acoustic velocity of pure water, as a function of temperature
    and pressure. Implements equation 28 from Batzle and Wang (1992), using
    the coefficients in Table 1.

    Note that this function does not work at pressures above about 100 MPa.

    Use scalars or arrays; if you use arrays, they must be the same size.

    Args:
        temperature (array): The temperature in degrees Celsius.
        pressure (array): The pressure in pascals.

    Returns:
        array: The velocity in m/s.
    """
    w = np.array([[ 1.40285e+03,  1.52400e+00,  3.43700e-03, -1.19700e-05],
                  [ 4.87100e+00, -1.11000e-02,  1.73900e-04, -1.62800e-06],
                  [-4.78300e-02,  2.74700e-04, -2.13500e-06,  1.23700e-08],
                  [ 1.48700e-04, -6.50300e-07, -1.45500e-08,  1.32700e-10],
                  [-2.19700e-07,  7.98700e-10,  5.23000e-11, -4.61400e-13]])

    T, P = np.asanyarray(temperature), np.asanyarray(pressure) * 1e-6
    return sum(w[i, j] * T**i * P**j for i in range(5) for j in range(4))


def v_brine(temperature, pressure, salinity):
    """
    The acoustic velocity of brine, as a function of temperature (deg C),
    pressure (Pa), and salinity (weight fraction). Implements equation 29
    from Batzle and Wang (1992).

    Note that this function does not work at pressures above about 100 MPa.

    Use scalars or arrays; if you use arrays, they must be the same size.

    Args:
        temperature (array): The temperature in degrees Celsius.
        pressure (array): The pressure in pascals.
        salinity (array): The weight fraction of NaCl, e.g. 35e-3
            for 35 parts per thousand, or 3.5% (the salinity of
            seawater).

    Returns:
        array: The velocity in m/s.
    """
    # Align with the symbols and units in Batzle & Wang.
    T, P = np.asanyarray(temperature), np.asanyarray(pressure)*1e-6
    S = np.asanyarray(salinity)

    v_w = v_water(temperature, pressure)
    s1 = 1170 - 9.6*T + 0.055*T**2 - 8.5e-5*T**3 + 2.6*P - 0.0029*T*P - 0.0476*P**2
    s15 = 780 - 10*P + 0.16*P**2
    s2 = -820

    return v_w + s1 * S + s15 * S**1.5 + s2 * S**2
