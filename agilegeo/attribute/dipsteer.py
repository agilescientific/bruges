import numpy as np
from agilegeo.attribute import energy 

def dipsteer( data, window_length, stepout,
              maxlag, overlap=0., dt=1 ):
    """
    Calculates a dip field by finding the maximum correlation between
    adjacent traces.

    :param data: A 2D seismic section (samples,traces) used to
                 calculate dip.
    :param window_length: The length [seconds] of the window to use
    :param dt: The time sample interval of the traces.
                          
    :param stepout: The number of traces on either side of each point
                    to average when calculating the dip.
    :param maxlag: The maximum amount time lag to use when correlating
                   the traces.
    :keyword overlap: The fractional overlap for each window.
                      A value of 0 uses no redudant data, a value of 1
                      slides the dip correlator one sample at a time.
                      Defaults to 0.5
    

    :returns a dip field [samples/trace] of the same shape as the
             input data, and correlation coefficients corresponding
             to the data.
    """

    dip = np.zeros( data.shape )

    window_length = np.floor( window_length / dt )
    
    # Force the window length to be odd for index tracking
    if not ( window_length % 2 ): window_length+=1
    
    # Define time windows
    stride = window_length * ( 1 - overlap ) 
    n_windows = np.ceil( ( data.shape[0] - window_length ) /
                         stride ) + 1
 
    # Normalize each trace to the same RMS energy
    norm_factor = np.sqrt( energy( data, window_length ) )
    norm_data = data / norm_factor

    # Replace the 0/0 with 0
    norm_data = np.nan_to_num( norm_data )

    # Mid point in the data which corresponds to zero dip
    zero_dip = ( np.floor( window_length / 2.0) +  maxlag  )
    
    # Loop over each trace we can do a full calculation for
    for i in np.arange( stepout +1, data.shape[-1] - (stepout +1) ):
        
        # Loop over each time window
        for j in np.arange( 0, n_windows):

            start = (j * stride) + ( maxlag )
            end = start + window_length

            # Do not compute last samples if we don't
            # have a full window
            if ( end > (norm_data.shape[0]-maxlag) ): break
        
            # Get the  kernel
            kernel = norm_data[ start : end, i ]

            dips_j = 0
         
            # Correlate with adjacent traces
            for k in np.arange( 1, stepout + 1 ):

                # Do the trace on the right
                r_trace = norm_data[ start - (k*maxlag) : \
                                     end + (k*maxlag), i+k ]
               
                cor_r = np.correlate( kernel, r_trace, mode='same' )

                if ( np.amax( cor_r ) < .1 ): dip_r = 0
                else:
                    dip_r = (  np.argmax( cor_r ) - zero_dip ) / k

               
                # Do the left trace
                l_trace = norm_data[ start - (k*maxlag) : \
                                     end + (k*maxlag), i-k ]
                
                cor_l = np.correlate( kernel, l_trace, mode = 'same' )

                if ( np.amax( cor_l ) < .1 ): dip_l = 0
                else:
                    dip_l = -(  np.argmax( cor_l ) - zero_dip ) / k

          
                dips_j += dip_r
                dips_j += dip_l

              
            # Average the result
            dips_j /= ( 2. * stepout )

            # Update the output 
            dip[ start : start + stride, i ] = dips_j
        
    return dip


