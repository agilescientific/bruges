from numpy import zeros, log2, ceil, arange, absolute
from scipy.fftpack import fft

from numpy import hanning, concatenate

def spectra( data, window=hanning(128), overlap=0.5 ):
    """
    Calculates a spectrogram using windowed STFTs.
     
    :param data: 1D numpy array to process into spectra. 
    :keyword window: Array of data defining the window to use for
                   the STFTs. (See scipy Tukey, Hann, Hamming, etc.).
                   Function will automatically pad to a power of 2.
                   Defaults to a 128 point Hann window.
                   
    :keyword overlap: The fractional overlap for each window.
                      A value of 0 uses no redudant data, a value of 1
                      slides the STFT window one sample at a time.
                      Defaults to 0.5
                      
    :returns A spectrogram of the data ([time, freq]).
            ( 2D array for 1D input )
    """

    
    # Calculate how many STFTs we need to do.
    stride = window.size * ( 1 - overlap ) + 1
    n_windows = ceil( ( data.size - window.size ) / stride ) + 1

    
    # Pad window to power of 2
    padded_window = zeros( 2**ceil( log2( window.size ) ) ) + window
    
    # Init the output array
    output = zeros( [n_windows, padded_window.size / 2.] )
    
    # Do the loop
    for i in arange( 0,n_windows, 1 ):
        
        start = i * stride
        end = start + padded_window.size 
        
        # Do not compute last samples if we don't have a full window
        if ( end > data.size-1 ): break        
 
        output[ i, : ] = ( 2.* \
                         absolute( fft( data[ start:end ] *\
                                        padded_window ) ) / \
                         window.size )[0 : padded_window.size /2. ]
    
    return output



    
        
        
        
        
