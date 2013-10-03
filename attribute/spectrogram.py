from numpy import zeros, log2, ceil, arange, absolute
from scipy.fftpack import fft
from traits.api import HasTraits, Array, Event, Range, Int
from numpy import hanning, concatenate

def compute_spectra( data, window, overlap=0.5 ):
    """
    Calculates a spectrogram using windowed STFTs. 
    :param data: Array of data to process into spectra.
    :param window: Array of data defining the window to use for
                   the STFTs. (See scipy Tukey, Hann, Hamming, etc.). Function
                   will automatically pad to a power of 2.
    :keyword overlap: The fractional overlap for each window. A value of 0
                      uses no redudant data, a value of 1 slides the STFT window
                      one sample at a time. Defaults to 0.5
    :returns A spectrogram of the data. ( 2D array for 1D input )
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
                         absolute( fft( data[ start:end ] * padded_window ) ) / \
                         window.size )[0 : padded_window.size /2. ]
    
    return output


class Spectrogram( HasTraits ):
    """
    Processing module that turns traces into spectra using hanning windowed 
    STFTs.
    """
    window = Array 
    updated = Event
    overlap = Range( 0.0,1.0, 0.5 )
    nsamps = Int
    
    def __init__( self, nsamps = 128, overlap=0.5 ):
        """
        Initializes the class
        :param nsamps: The number of samples to use for STFT window. 
        """
        
        super( Spectrogram, self ).__init__()
        
        self.nsamps = nsamps
        self.overlap = overlap
        
        self.window = hanning( nsamps )
    
    def _on_nsamps_changed( self ):
        
        self.window = hanning( self.nsamps )
        self.updated = True
    
    def compute( self, data ):
        
        
        for trace in range( data.shape[-1]):
            
            if trace == 0:
                first = compute_spectra( data[:,trace], self.window, 
                                         self.overlap ) 
                output = zeros( [ first.shape[0], data.shape[-1], \
                                  first.shape[1] ] )
                output[:,0, : ] = first        
            else: 
                spec = compute_spectra( data[:,trace], \
                                      self.window, 
                                      self.overlap )  
                output[ :, trace, : ] = spec
        return output
        
    
        
        
        
        