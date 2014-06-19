from scipy.interpolate import griddata, interp1d
from numpy import  arange, amax, amin, asarray,zeros, cumsum, \
     transpose, gradient, mean

def time_to_depth(data,vmodel, dt, dz, twt=True):
    """
    Converts data from the time domain to the depth domain given a
    velocity model.

    :param data: The data to convert, will work with a 1 or 2D numpy
                 numpy array. array(samples,traces).
    :param vmodel: P-wave velocity model that corresponds to the data.
                   Must be the same shape as data.
    :param dt: The sample interval of the input data [s].
    :param dz: The sample interval of the output data [m].

    :keyword twt: Use twt travel time, defaults to true

    :returns: The data resampled in the depth domain.
    """

    if twt:
        scale = 2.0
    else:
        scale = 1.0
    # Do depth to time with inverted velocity profile
    return convert(data, 1. / vmodel, dt, dz, scale)


def convert(data, vmodel, interval, interval_new,scale):
    """
    Generic function for converting between scales. Use either
    time to depth or depth to time
    """

    dz = interval
    dt = interval_new

    if( len( data.shape ) == 1 ):
        ntraces = 1
        nsamps = data.size
    else:
        ntraces = data.shape[-1]
        nsamps = data.shape[0]
        
    depths = transpose(asarray([(arange(nsamps) * dz) \
                       for i in range(ntraces)]))

    v_avg = cumsum( vmodel, axis=0 ) / \
                    transpose([arange( nsamps ) \
                            for i in range(ntraces)])

    

    
    # convert depths to times
    times = depths / v_avg

    times *= scale

    times_lin = arange(amin(times), amax(times ), dt)
    
    if( ntraces == 1 ):
        inter = interp1d(times, data,
                         bounds_error=False,
                         fill_value = data[-1],
                         kind='nearest')
        return(inter(times_lin))

    output = zeros((times_lin.size, ntraces))
    
    for i in range(ntraces):
        
        inter = interp1d(times[:,i], data[:,i],
                         bounds_error=False,
                         fill_value = data[-1,i],
                         kind='nearest')
        output[:,i] += inter(times_lin)


    return(output)
    
def depth_to_time(data,vmodel, dz, dt, twt=True):
    """
    Converts data from the depth domain to the time domain given a
    velocity model.

    :param data: The data to convert, will work with a 1 or 2D numpy
                 numpy array. array(samples,traces).
    :param vmodel: P-wave velocity model that corresponds to the data.
                   Must be the same shape as data.
    :param dz: The sample interval of the input data [m].
    :param dt: The sample interval of the output data [s].

    :keyword twt: Use twt travel time, defaults to true

    :returns: The data resampled in the time domain.
    """

    if twt:
        scale=0.5
    else: scale=1.0
    # Do depth to time with inverted velocity profile
    return convert(data, vmodel, dz, dt, scale)
         
        

                

    
    

    
