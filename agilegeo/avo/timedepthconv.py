from scipy.interpolate import griddata, interp1d
from numpy import  arange, amax, amin, asarray,zeros, cumsum, \
     transpose, gradient, mean

def time_to_depth( data,vmodel, dt, dz ):
    """
    Converts data from the time domain to the depth domain given a
    velocity model.

    :param data: The data to convert, will work with a 1 or 2D numpy
                 numpy array. array(samples,traces).
    :param vmodel: P-wave velocity model that corresponds to the data.
                   Must be the same shape as data.
    :param dt: The sample interval of the input data [s].
    :param dz: The sample interval of the output data [m].

    :returns: The data resampled in the depth domain.
    """

    # Do depth to time with inverted velocity profile
    return( depth_to_time( data, 1. / vmodel, dt, dz ) )


def depth_to_time( data,vmodel, dz, dt ):
    """
    Converts data from the time domain to the depth domain given a
    velocity model.

    :param data: The data to convert, will work with a 1 or 2D numpy
                 numpy array. array(samples,traces).
    :param vmodel: P-wave velocity model that corresponds to the data.
                   Must be the same shape as data.
    :param dz: The sample interval of the input data [m].
    :param dt: The sample interval of the output data [s].

    :returns: The data resampled in the time domain.
    """

    if( len( data.shape ) == 1 ):
        ntraces = 1
        nsamps = data.size
    else:
        ntraces = data.shape[-1]
        nsamps = data.shape[0]
        
    depths = transpose(asarray([ (arange( nsamps ) * dz) \
                       for i in range(ntraces) ] ) )

    v_avg = cumsum( vmodel, axis=0 ) / \
                    transpose( [arange( nsamps ) \
                                for i in range( ntraces )] )
    # convert depths to times
    times = depths / v_avg
    print( amax( depths ) )
    times_lin = arange( amin( times), amax( times ), dt )
    #for i in gradient( times[:,50] ) : print i
    if( ntraces == 1 ):
        inter = interp1d(times, data,
                         bounds_error=False)
        return( inter( times_lin ) )

    output = zeros( (times_lin.size, ntraces) )
    
    for i in range( ntraces ):
        
        inter = interp1d(times[:,i], data[:,i],
                         bounds_error=False,
                         fill_value = 0,
                         kind='nearest')
        output[:,i] += inter( times_lin )

    return( output )
         
        

                

    
    

    
