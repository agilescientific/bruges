import numpy as np
from scipy.signal import hilbert

def ricker( duration, dt, f ):
    """
    Also known as the mexican hat wavelet, models the function:
    A =  (1-2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}

    :param duration: The length in seconds of the wavelet.
    :param dt: is the sample interval in seconds (usually 0.001,
               0.002, 0.004)
    :params f: Center frequency of the wavelet (in Hz). If a list or tuple is
               passed, the first element will be used.

    :returns: The ricker wavelet with center frequency f sampled at t.
    """

    # Check size of f. API compatibility may allow for a list
    if isinstance(f, list) or isinstance(f,tuple): f = f[0]
     
    t = np.arange( -duration/2, duration/2 , dt)

    pi2 = (np.pi ** 2.0)
    fsqr = f ** 2.0
    tsqr = t ** 2.0
    pft = pi2 * fsqr * tsqr
    A = (1 - (2 * pft)) * np.exp(-pft)

    return A
    
    
def sweep( duration, dt, f, method = 'linear', phi = 0,
           vertex_zero = True, autocorrelate = True ):
    """
    Generates a linear frequency modulated wavelet (sweep)
    Does a wrapping of scipy.signal.chirp
    :param duration: The length in seconds of the wavelet.
    :param dt: is the sample interval in seconds (usually 0.001, 0.002, 0.004)
    :param f: Tuple of (f1, f2), or a similar list.
    :keyword method: {'linear','quadratic','logarithmic'}, optional
    :keyword phi: float, phase offset in degrees
    :keyword vertex_zero: bool, optional
        This parameter is only used when method is 'quadratic'. 
        It determines whether the vertex of the parabola that 
        is the graph of the frequency is at t=0 or t=t1.

    :returns: An LFM waveform.
    """
    
    t = np.arange( -duration/2, duration/2 , dt) 
    t0 = -duration/2
    t1 = duration/2
    
    f1 = f[0]
    f2 = f[1]
    
    from scipy.signal import chirp
    
    A = chirp( t , f1, t1, f2, method, phi, vertex_zero  )
    
    if autocorrelate:
        A = np.correlate(A, A, mode='full')
     
    return A / np.amax( A )


def ormsby(duration, dt, f):
    """
    The Ormsby wavelet requires four frequencies:
    f1 = low-cut frequency
    f2 = low-pass frequency
    f3 = high-pass frequency
    f4 = hi-cut frequency
    Together, the frequencies define a trapezoid shape in the spectrum.
    The Ormsby wavelet has several sidelobes, unlike Ricker wavelets
    which only have two, one either side.

    :param duration: The length in seconds of the wavelet.
    :param dt: is the sample interval in seconds (usually 0.001, 0.002, 0.004)
    :params f: Tuple of form (f1,f2,f3,f4), or a similar list. If fewer (or 
    more than 

    :returns: A vector containing the ormsby wavelet
    """
    # Try to handle some duck typing
    if not ( isinstance(f, list) or isinstance(f, tuple) ) : f = [f]
    
    # Deal with having fewer than 4 frequencies
    if len(f) == 4:
        f1 = f[0]
        f2 = f[1]
        f3 = f[2]
        f4 = f[3]
    else:
        # Cope with only having one frequency
        # This is an arbitrary hack, is this desirable?
        # Need a way to notify with warnings
        f1 = f[0]/4
        f2 = f[0]/2
        f3 = f[0]*2
        f4 = f[0]*2.5
    
    def numerator(f,t):
        return (np.sinc( f * t)**2) * ((np.pi * f) ** 2)

    pf43 = ( np.pi * f4 ) - ( np.pi * f3 )
    pf21 = ( np.pi * f2 ) - ( np.pi * f1 )
    
    t = np.arange( -duration/2, duration/2 , dt) 
    
    A = ( ( numerator(f4,t)/pf43 ) - ( numerator(f3,t)/pf43 ) -
         ( numerator(f2,t)/pf21 ) + ( numerator(f1,t)/pf21) )

    A /= np.amax( A )
    return A


def rotate_phase( w, phi ):
    """
    Performs a phase rotation of wavelet using:
    A = w(t)Cos(phi) + h(t)Sin(phi)
    Where w(t) is the wavelet and h(t) is it's hilbert transform.
    
    :params w: The wavelet vector.
    :params phi: The phase rotation angle (in Radians) to apply.

    :returns: The phase rotated signal.
    """

    # Get the analytic signal for the wavelet
    a = hilbert( w )

    A = (np.real( a ) * np.cos( phi ) +
         np.imag( a  ) * np.sin( phi ) )

    return A
    