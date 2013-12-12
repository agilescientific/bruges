import numpy
from scipy.signal import fftconvolve

def energy(traces, duration, dt=1 ):
    """
    Compute an mean-squared energy measurement for each point of a
    seismic section.
        
    :param traces: The data array to use for calculating MS energy.
                   Must be 1D or 2D numpy array.
    :param duration: the time duration of the window (in seconds), or
                     samples if dt=1.
    :param dt: the sample interval of the data (in seconds). Defaults
               to 1 so duration can be in samples.
    :returns: An array the same dimensions as the input array.
    """
    
    energy_data = numpy.zeros( traces.shape )
    signal = traces * traces 
    n_samples = int(duration / dt)
    
    window = numpy.ones( n_samples )
    
    if ( len( signal.shape ) ) == 1 : 
        ## Compute the sliding average using a convolution
        energy_data = fftconvolve( signal, window, mode='same' ) \
                     / n_samples
    
    elif ( len( signal.shape ) == 2 ):
        for trace in range(signal.shape[1]):
            energy_data[ :, trace ] = (
                    fftconvolve( signal[:,trace], window,
                                 mode='same' ) )
    
    else: raise ValueError( 'Array must be 1D or 2D' )
    
    return energy_data
    

    
    
