from numpy import zeros

def integrate_spectra( data, b1, b2, norm=True ):
    """
    Integrates across each slice of a spectrogram.

    :param data: A 2D numpy array (time,freq) of time-frequency data.
                 See spectrogram for details on how to create spectra.

    :param b1: The lower band limit to use for spectral integration.
    :param b2: The uppper band limit to use for spectral integration.

    :keyword norm :default=True: Normalize the spectra to unity before
               integrating the band.
    :returns a timeseries array of the integrated spectra.
    """

    if ( norm==True ):
        total = data.sum( axis=1 )
        
        output = data[ :, b1:b2].sum( axis = 1 ) / total
    else:
        output = data[ :, b1:b2].sum( axis = 1 ) 
        
        
    return output



        
