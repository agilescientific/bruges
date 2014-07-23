import scipy as sp
import numpy as np

"""
These are generally useful function calls that don't necessarily fit
into a specific geophysics
"""
def rms(a):
    """
    Calculates the RMS of an array.

    :param a: An array.

    :returns: The RMS of the array.

    """

    return np.sqrt(np.sum(a**2.0)/a.size)


def normalize(a, new_min=0.0, new_max=1.0):
    """
    Normalize an array to [0,1] or to 
    arbitrary new min and max.

    :param a: An array.
    :param new_min: A float to be the new min, default 0.
    :param new_max: A float to be the new max, default 1.

    :returns: The normalized array.
    """
    
    n = (a - np.amin(a)) / np.amax(a - np.amin(a))
    return n * (new_max - new_min) + new_min


def next_pow2( num ):
    """
    Calculates the next nearest power of 2 to the input. Uses
      2**ceil( log2( num ) ).

    :param num: The number to round to the next power if two.

    :returns: the next power of 2 closest to num.
    """

    return( 2**np.ceil( np.log2( num ) ))


def noise_db(a, snr):
    """
    Takes an array of seismic amplitudes
    and SNR in dB.

    Returns an array of noise, the same 
    shape as the input.

    Note it does *not* return the input
    array with the noise added.

    """
    
    # Get the amplitude of the signal
    sigmean = rms(a)
    
    # Calculate the amp of the noise,
    # given the desired SNR
    noisemean = sigmean / 10.0**(snr/20.0)
    
    # Normal noise, centered on 0,
    # SD=sqrt(var), same shape as input
    noise = noisemean * np.random.normal(0.0, 1.0, a.shape)
    
    return noise
