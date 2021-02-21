"""
====================
rockphysicsmodels.py
====================

A bunch of rock physics models.
References are mentioned in docstrings of individual functions.
Docstrings follow numpy/scipy convention.

Alessandro Amato del Monte, March 2019
"""
import numpy as np


def critical_porosity(K0, G0, phi, phi_c=0.4):
    '''
    Critical porosity model.
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of dry sand for porosities below the critical porosity (phi_c).
    Above phi_c, the fluid phase is load-bearing, below phi_c the
    solid phase (mineral grains) is load-bearing.
    The equations here describe the variation of K and G
    for porosities below phi_c as straight lines.
    Critical porosity is usually 0.4 for sandstone,
    0.7 for chalk, 0.02-0.03 for granites.

    Parameters
    ----------   
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    phi_c : float, optional
        Critical porosity. Default: 0.4

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.

    References
    ----------
    Mavko et al. (2009), The Rock Physics Handbook, Cambridge University Press (p.370)
    '''
    Kdry = K0 * (1 - phi / phi_c)
    Gdry = G0 * (1 - phi / phi_c)
    return Kdry, Gdry


def hertz_mindlin(K0, G0, sigma, phi_c=0.4, Cn=8.6, f=1):
    '''
    Hertz-Mindlin model.
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of a dry pack of spheres subject to a hydrostatic confining pressure. 

    Parameters
    ----------   
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    sigma : float
        Effective stress in MPa.
    phi_c : float, optional
        Critical porosity. Default: 0.4
    Cn : float, optional
        Coordination number Default: 8.6.
    f : float, optional
        Shear modulus correction factor,
        f=1 for dry pack with perfect adhesion
        between particles and f=0 for dry frictionless pack.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.

    References
    ----------
    Mavko et al. (2009), The Rock Physics Handbook, Cambridge University Press (p.246)
    '''
    sigma0 = sigma / 1e3  # converts pressure in same units as solid moduli (GPa)
    pr0 =(3*K0-2*G0) / (6*K0+2*G0)  # poisson's ratio of mineral mixture
    Khm = (sigma0*(Cn**2*(1 - phi_c)**2*G0**2) / (18
        * np.pi**2 * (1 - pr0)**2))**(1/3)
    Ghm = ((2+3*f-pr0*(1+3*f)) / (5*(2-pr0))) * ((
        sigma0 * (3 * Cn**2 * (1 - phi_c)**2 * G0**2) / (
        2 * np.pi**2 * (1 - pr0)**2)))**(1/3)
    return Khm, Ghm


def soft_sand(K0, G0, phi, sigma, phi_c=0.4, Cn=8.6, f=1):
    '''
    Soft sand, or friable sand or uncemented sand model.
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of poorly sorted dry sand by interpolating with the lower Hashin-Shtrikman bound
    the two end members at zero porosity and critical porosity.
    The zero porosity end member has K and G equal to mineral.
    The end member at critical porosity has K and G given by Hertz-Mindlin model.

    Parameters
    ----------   
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    sigma : float
        Effective stress in MPa.
    phi_c : float, optional
        Critical porosity. Default: 0.4
    Cn : float, optional
        Coordination number Default: 8.6.
    f : float, optional
        Shear modulus correction factor,
        f=1 for dry pack with perfect adhesion
        between particles and f=0 for dry frictionless pack.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.

    References
    ----------
    Mavko et al. (2009), The Rock Physics Handbook, Cambridge University Press (p.258)
    '''
    Khm, Ghm = hertz_mindlin(K0, G0, sigma, phi_c, Cn, f)
    Kdry = -4/3 * Ghm + (((phi / phi_c) / (Khm + 4/3 * Ghm))
        + ((1 - phi / phi_c) / (K0 + 4/3 * Ghm)))**-1
    gxx = Ghm / 6 * ((9 * Khm + 8 * Ghm)  /  (Khm + 2 * Ghm))
    Gdry = -gxx + ((phi / phi_c) / (Ghm + gxx)
        + ((1 - phi / phi_c) / (G0 + gxx)))**-1
    return Kdry, Gdry


def stiff_sand(K0, G0, phi, sigma, phi_c=0.4, Cn=8.6, f=1):
    '''
    Stiff sand model.
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of stiff dry sands by interpolating with the upper Hashin-Shtrikman bound
    the two end members at zero porosity and critical porosity.
    The zero porosity end member has K and G equal to mineral.
    The end member at critical porosity has K and G given by Hertz-Mindlin model.

    Parameters
    ----------   
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    sigma : float
        Effective stress in MPa.
    phi_c : float, optional
        Critical porosity. Default: 0.4
    Cn : float, optional
        Coordination number Default: 8.6.
    f : float, optional
        Shear modulus correction factor,
        f=1 for dry pack with perfect adhesion
        between particles and f=0 for dry frictionless pack.  
    
    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa

    References
    ----------
    Mavko et al. (2009), The Rock Physics Handbook, Cambridge University Press (p.260)
    '''
    Khm, Ghm = hertz_mindlin(K0, G0, sigma, phi_c, Cn, f)
    Kdry = -4/3 * G0 + (((phi / phi_c) / (Khm + 4/3 * G0))
        + ((1 - phi / phi_c) / (K0 + 4/3 * G0)))**-1
    tmp = G0 / 6*((9 * K0 + 8 * G0)  /  (K0 + 2 * G0))
    Gdry = -tmp + ((phi / phi_c) / (Ghm + tmp)
        + ((1 - phi / phi_c) / (G0 + tmp)))**-1
    return Kdry, Gdry


def contact_cement(K0, G0, phi, phi_c=0.4, Cn=8.6, Kc=37, Gc=45, scheme=2):
    '''
    Contact cement or cemented sand model,.
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of dry sand where cement is deposited at grain contacts.
    The cement properties can be modified as well as the type of 
    cementation scheme.

    Parameters
    ----------
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    phi_c : float, optional
        Critical porosity. Default: 0.4
    Cn : float, optional
        Coordination number Default: 8.6.
    Kc, Gc : float, optional
        Cement bulk & shear modulus in GPa. Default: 37, 45.
    scheme : int, optional
        Cementation scheme, can be either 1 or 2:
        1: cement deposited at grain contacts
        2: cement as uniform layer around grains.
        Default: 2.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.

    References
    ----------
    Dvorkin-Nur (1996), Elasticity of High-Porosity Sandstones: Theory for Two North Sea Data Sets.
    Geophysics 61, no. 5 (1996).
    Mavko et al. (2009), The Rock Physics Handbook, Cambridge University Press (p.255)
    '''
    pr0 = (3 * K0 - 2 * G0) / (6 * K0 + 2 * G0)
    PRc = (3 * Kc - 2 * Gc) / (6 * Kc + 2 * Gc)
    if scheme == 1:  # scheme 1: cement deposited at grain contacts
        alpha = ((phi_c - phi) / (3 * Cn * (1 - phi_c)))**(1/4)
    else:  # scheme 2: cement evenly deposited on grain surface
        alpha = ((2 * (phi_c - phi)) / (3 * (1 - phi_c)))**(1/2)
    LambdaN = (2 * Gc * (1 - pr0) * (1 - PRc)) / (
        np.pi * G0 * (1 - 2 * PRc))
    N1 = -0.024153 * LambdaN**-1.3646
    N2 = 0.20405 * LambdaN**-0.89008
    N3 = 0.00024649 * LambdaN**-1.9864
    Sn = N1 * alpha**2 + N2 * alpha + N3
    LambdaT = Gc / (np.pi * G0)
    T1 = -10**-2 * (2.26 * pr0**2 + 2.07 * pr0
        + 2.3) * LambdaT**(0.079 * pr0**2 + 0.1754 * pr0 - 1.342)
    T2 = (0.0573 * pr0**2 + 0.0937 * pr0
        + 0.202) * LambdaT**(0.0274 * pr0**2 + 0.0529 * pr0 - 0.8765)
    T3 = 10**-4 * (9.654 * pr0**2 + 4.945 * pr0
        + 3.1) * LambdaT**(0.01867 * pr0**2 + 0.4011 * pr0 - 1.8186)
    St = T1 * alpha**2 + T2 * alpha + T3
    Kdry = 1 / 6 * Cn *(1 - phi_c ) * (Kc + 4/3 * Gc) * Sn
    Gdry = 3 / 5 * Kdry + 3 / 20 * Cn * (1 - phi_c) * Gc * St
    return Kdry, Gdry


def constant_cement(K0, G0, phi, phi_cem=0.38, phi_c=0.4, Cn=8.6, Kc=37, Gc=45, scheme=2):
    '''
    Constant cement model, Avseth et al. (2000).
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of high porosity dry sand with a certain initial cementation
    by interpolating with the lower Hashin-Shtrikman bound
    the two end members at zero porosity and critical porosity.
    The zero porosity end member has K and G equal to mineral.
    The high porosity end member has K and G given by the Contact Cement model.
    It is assumed that the porosity reduction is due to
    non-cementing material filling in the available pore space.

    Parameters
    ----------
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    phi_cem : float, optional
        Porosity at initial cementation. Default: 0.38.
    phi_c : float, optional
        Critical porosity. Default: 0.4.
    Cn : float, optional
        Coordination number. Default: 8.6.
    Kc, Gc : float, optional
        Cement bulk & shear modulus in GPa. Default: 37, 45.
    scheme : int, optional
        Cementation scheme, can be either 1 or 2:
        1: cement deposited at grain contacts
        2: cement as uniform layer around grains.
        Default: 2.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.
    
    References
    ----------
    Dvorkin et al. (2014), Seismic Reflections of Rock Properties, Cambridge University Press (p.30-31)
    '''
    # contact cement model
    Khi, Ghi = contact_cement(K0, G0, phi, phi_c=phi_c, Cn=Cn, Kc=Kc, Gc=Gc, scheme=scheme)
    # lower bound Hashin-Shtrikman starting from phi_cem
    Kcc, Gcc = contact_cement(K0, G0, phi_cem, phi_c=phi_c, Cn=Cn, Kc=Kc, Gc=Gc, scheme=scheme)
    Klo = -4/3 * Gcc + (((phi / phi_cem) / (Kcc + 4/3 * Gcc)) + (
        (1 - phi / phi_cem) / (K0 + 4/3 * Gcc)))**-1
    tmp = Gcc / 6* ((9 * Kcc + 8 * Gcc)  /  (Kcc + 2 * Gcc))
    Glo = -tmp + ((phi / phi_cem) / (Gcc + tmp) + (
        (1 - phi / phi_cem) / (G0 + tmp)))**-1
   
    # initialize empty vectors for K and G dry
    Kdry, Gdry = (np.full(phi.size, np.nan) for _ in range(2))
    # for porosities>phi_cem use [K,G]_HI = contact cement model
    # for porosities<=phi_cem use [K,G]_LO = constant cement model
    Kdry[phi > phi_cem] = Khi[phi > phi_cem]
    Kdry[phi <= phi_cem] = Klo[phi <= phi_cem]
    Gdry[phi > phi_cem] = Ghi[phi > phi_cem]
    Gdry[phi <= phi_cem] = Glo[phi <= phi_cem]
    return Kdry, Gdry


def increasing_cement(K0, G0, phi, phi_cem=0.38, phi_c=0.4, Cn=8.6, Kc=37, Gc=45, scheme=2):
    '''
    Increasing cement model (Modified Hashin-Shtrikman upper bound).
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of a dry sand with a certain initial cementation
    by interpolating with the upper Hashin-Shtrikman bound
    the two end members at zero porosity and critical porosity.
    The zero porosity end member has K and G equal to mineral.
    The high porosity end member has K and G given by the Contact Cement model.
    Probably best to avoid using if for porosities>phi_cem.
    Need to check references.

    Parameters
    ----------
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    phi_cem : float, optional
        Porosity at initial cementation. Default: 0.38.
    phi_c : float, optional
        Critical porosity. Default: 0.4.
    Cn : float, optional
        Coordination number. Default: 8.6.
    Kc, Gc : float, optional
        Cement bulk & shear modulus in GPa. Default: 37, 45.
    scheme : int, optional
        Cementation scheme, can be either 1 or 2:
        1: cement deposited at grain contacts
        2: cement as uniform layer around grains.
        Default: 2.

    Returns
    -------
    Kdry, Gdry : float or array_like
        dry rock bulk & shear modulus in GPa.
    '''
    Kcc, Gcc = contact_cement(K0, G0, phi_cem, phi_c=phi_c, Cn=Cn, Kc=Kc, Gc=Gc, scheme=scheme)
    Kdry = -4/3 * G0 + (((phi / phi_cem) / (Kcc + 4/3 * G0)) + (
        (1 - phi / phi_cem) / (K0 + 4/3 * G0)))**-1
    tmp = G0 / 6 * ((9 * K0 + 8 * G0)  /  (K0 + 2 * G0))
    Gdry = -tmp + ((phi / phi_cem) / (Gcc + tmp) + (
        (1 - phi / phi_cem) / (G0 + tmp)))**-1
    return Kdry, Gdry


def vernik_consol_sand(K0, G0, phi, sigma, b=10):
    '''
    Vernik & Kachanov Consolidated Sand Model.
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of consolidated dry sand subject to a hydrostatic confining pressure
    as a continuous solid containing pores and cracks.

    Parameters
    ----------
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    sigma : float
        Effective stress in MPa.
    b : float, optional
        Slope parameter in pore shape empirical equation, range: 8-12.
        Default: 10.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.
    
    References
    ----------
    Vernik & Kachanov (2010), Modeling elastic properties of siliciclastic rocks, Geophysics v.75 n.6
    '''
    # empirical pore shape factor:
    p = 3.6 + b * phi
    q = p # true if phi>0.03
    psf = phi / (1 - phi)  # psf = pore shape factor multiplier

    # matrix properties: assuming arenites with K=35.6 GPa, G=33 GPa, Poisson's ratio nu_m = 0.146
    nu_m = 0.146
    Avm = (16 * (1 - nu_m**2) ) / ( 9 * (1 - 2 * nu_m))       # nu_m=0.146 --> Avm=2.46
    Bvm = (32 * (1 - nu_m) * (5 - nu_m)) / (45 * (2 - nu_m))  # nu_m=0.146 --> Bvm=1.59

    # crack density: inversely correlated to effective stress
    eta0 = 0.3 + 1.6 * phi  # crack density at zero stress
    d = 0.07  # compaction coefficient
    d = 0.02 + 0.003 * sigma
    cd = (eta0 * np.exp(-d * sigma)) / (1 - phi)

    # note: the presence at denominator of the factor (1 - phi) in psf and cd is needed
    # to account for the interaction effects, i.e. the presence of pores raises the average stress
    # in the matrix increasing compliance contributions of pores and cracks
    # this correction is referred to as Mori-Tanaka's scheme.
    # in this way, the original model which is a NIA (non-interaction model)
    # is extended and becomes effectively a model which does take into account interactions.
    Kdry = K0 * (1 + p * psf + Avm * cd)**-1
    Gdry = G0 * (1 + q * psf + Bvm * cd)**-1
    return Kdry, Gdry


def vernik_soft_sand_1(K0, G0, phi, sigma, phi_c=0.36, phi_con=0.26, b=10, n=2.00, m=2.05):
    '''
    Vernik & Kachanov Soft Sand Model 1.
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of dry sand modeled as a granular material.
    Only applicable for porosities between the low-porosity end-member
    (at the consolidation porosity phi_con) and the high-porosity
    end-member (at the critical porosity phi_c).
    The low-porosity end member is calculated with Vernik's
    Consolidated Sand Model.

    Parameters
    ----------
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    sigma : float
        Effective stress in MPa.
    phi_c : float, optional
        Critical porosity, range 0.30-0.42. Default: 0.36.
    phi_con : float, optional
        Consolidation porosity, range 0.22-0.30. Default: 0.26.
    b : float, optional
        Slope parameter in pore shape empirical equation, range: 8-12.
        Default: 10.
    n, m : float
        Empirical factors. Default: 2.00, 2.05.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.
    
    References
    ----------
    Vernik & Kachanov (2010), Modeling elastic properties of siliciclastic rocks, Geophysics v.75 n.6
    '''
    if isinstance(phi, np.ndarray):
        phi_edit = phi.copy()
        phi_edit[(phi_edit < phi_con) | (phi_edit > phi_c)]=np.nan
    else:
        phi_edit = np.array(phi)
        if (phi_edit < phi_con) | (phi_edit > phi_c):
            return np.nan, np.nan
    M0 = K0 + 4/3 * G0
    K_con, G_con = vernik_consol_sand(K0, G0, phi_con, sigma, b)
    M_con = K_con + 4/3 * G_con
    T = (1 - (phi_edit - phi_con) / (phi_c - phi_con))
    Mdry = M_con * T**n
    Gdry = G_con * T**m
    Kdry = Mdry - 4/3 * Gdry
    return Kdry, Gdry


def vernik_soft_sand_2(K0, G0, phi, p=20, q=20):
    '''
    Vernik & Kachanov Soft Sand Model 2.
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of dry sand modeled as a granular material.
    Applicable in the entire porosity range.

    Parameters
    ----------
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    p, q : float, optional
        Pore shape factor for K and G, range: 10-45.
        Default: 20.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.
    
    References
    ----------
    Vernik & Kachanov (2010), Modeling elastic properties of siliciclastic rocks, Geophysics v.75 n.6
    '''
    M0 = K0 + 4/3 * G0
    Mdry = M0 * (1 + p * (phi / (1 - phi)))**-1
    Gdry = G0 * (1 + q * (phi / (1 - phi)))**-1
    Kdry = Mdry - 4/3 * Gdry
    return Kdry, Gdry


def vernik_sand_diagenesis(K0, G0, phi, sigma, phi_c=0.36, phi_con=0.26, b=10, n=2.00, m=2.05):
    '''
    Vernik & Kachanov Sandstone Diagenesis Model.
    This model describes the elastic behaviour (K=bulk and G=shear moduli)
    of dry sand modeled as a continuous solid containing pores and cracks
    for porosities below phi_con (consolidation porosity)
    using Vernik's Consolidated Sand Model,  and as a granular material
    for porosities above phi_con using Vernik's Soft Sand Model 1.

    Parameters
    ----------
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    sigma : float
        Effective stress in MPa.
    phi_c : float, optional
        Critical porosity, range 0.30-0.42. Default: 0.36.
    phi_con : float, optional
        Consolidation porosity, range 0.22-0.30. Default: 0.26.
    b : float, optional
        Slope parameter in pore shape empirical equation, range: 8-12.
        Default: 10.
    n, m : float
        Empirical factors. Default: 2.00, 2.05.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.

    References
    ----------
    Vernik & Kachanov (2010), Modeling elastic properties of siliciclastic rocks, Geophysics v.75 n.6
    '''
    Kdry, Gdry = vernik_consol_sand(K0, G0, phi, sigma, b)
    Kdry_soft, Gdry_soft = vernik_soft_sand_1(K0, G0, phi, sigma, phi_c, phi_con, b, n, m)
    if isinstance(phi, np.ndarray):
        uu = phi>=phi_con
        Kdry[uu] = Kdry_soft[uu]
        Gdry[uu] = Gdry_soft[uu]
        return Kdry, Gdry
    else:
        if phi <= phi_con:
            return Kdry, Gdry
        else:
            return Kdry_soft, Gdry_soft


def vernik_shale(vclay, phi, rhom=2.73, rhob=1, Mqz=96, c33_clay=33.4, A=0.00284):
    '''
    Vernik & Kachanov Shale Model.
    This model describes the elastic behaviour in terms of velocities
    and density of inorganic shales.

    Parameters
    ----------
    vclay : float or array_like
        Dry clay content volume fraction.
    phi : float or array_like
        Porosity, maximum 0.40.
    rhom : float, optional
        Shale matrix density in g/cc. Default: 2.73.
    rhob : float, optional
        Brine density in g/cc. Default: 1.
    Mqz : float, optional
        P-wave elastic modulus of remaining minerals in GPa
        Default: 96.
    c33_clay : float, optional
        Anisotropic clay constant in GPa. Default: 33.4.
    A : float, optional
        Empirical coefficient for Vs. Default is good for illite/smectite/chlorite,
        can be raised up to .006 for kaolinite-rich clays.
        Default: 0.00284.
        
    Returns
    -------
    vp, vs, density : float or array_like
        P- and S-wave velocities in m/s, density in g/cc.

    Notes
    -----
    Shale matrix density (rhom) averages 2.73 +/- 0.03 g/cc at porosities below 0.25.
    It gradually varies with compaction and smectite-to-illite transition.
    A more accurate estimate can be calculated with this equation:
    rhom = 2.76+0.001*((rho-2)-230*np.exp(-4*(rho-2)))

    References
    ----------
    Vernik & Kachanov (2010), Modeling elastic properties of siliciclastic rocks, Geophysics v.75 n.6
    '''
    rho_matrix = 2.65 * (1 - vclay) + rhom * vclay
    k = 5.2 - 1.3 * vclay
    B, C = 0.287, 0.79
    c33_min = (vclay / c33_clay + (1 - vclay) / Mqz)**-1
    c33 = c33_min * (1 - phi)**k
    vp = np.sqrt(c33 / (rhom * (1 - phi) + rhob * phi))
    vs = np.sqrt(A * vp**4 + B * vp**2 - C)
    rho = rho_matrix * (1 - phi) + rhob * phi
    return vp * 1e3,vs * 1e3, rho
