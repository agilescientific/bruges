import scipy as sp
import numpy as np

"""
These are generally useful function calls that don't necessarily fit
into a specific geophysics
"""


def next_pow2( num ):
    """
    Calculates the next nearest power of 2 to the input. Uses
      2**ceil( log2( num ) ).

    :param num: The number to round to the next power if two.

    :returns: the next power of 2 closest to num.
    """

    return( 2**np.ceil( np.log2( num ) ))
