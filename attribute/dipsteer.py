import numpy as np

def dipsteer( data, window_length, stepout, maxlag ):
    """
    Calculates a dip field by finding the maximum correlation between
    adjacent traces.

    :param data: A 2D seismic section (samples,traces) used to
                 calculate dip.
    :param window_length: The length [samples] of the window to use
                          
    :param stepout: The number of traces on either side of each point
                    to average when calculating the dip.
    :param maxlag: The maximum amount time lag to use when correlating
                   the traces.

    :returns a dip field [samples/trace] of the same shape as the
             input data, and correlation coefficients corresponding
             to the data.
    """

    output = np.zeros( data.shape )
    
    # Define time windows

    # Loop over each time window

    # Normalize kernels
    
    # Correlate with adjacent traces

    # Find maximum sample

    # Check correlation threshold

    # Set dip in output

    return output


