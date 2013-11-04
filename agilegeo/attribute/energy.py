import numpy
from scipy.signal import fftconvolve

from traits.api import HasTraits, Int, Event
from traitsui.api import View, Item

def compute_energy(traces, n_samples ):
    """
    Compute an mean-squared energy measurement for each point of a
    seismic section.
        
    :param traces: The data array to use for calculating MS energy.
                   Must be 1D or 2D numpy array.
    :param n_samples: The number of samples to use for the MS window.

    :returns An array the same dimensions as the input array.
    """
    
    energy_data = numpy.zeros( traces.shape )
    signal = traces * traces 
    window = numpy.zeros( n_samples ) + 1.0
    
    if ( len( signal.shape ) ) == 1 : 
        ## Compute the sliding average using a convolution
        energy_data = fftconvolve( signal, window, mode='same' ) \
                     / n_samples
    
    elif ( len( signal.shape ) == 2 ):
        for trace in range(signal.shape[1]):
            energy_data[ :, trace ] = (
                    fftconvolve( signal[:,trace], window, mode='same' ) )
    
    else: raise ValueError( 'Array must be 1D or 2D' )
    
    return energy_data
    
class Energy( HasTraits ):
    """
    Simple class for creating energy attribute calculator.
    """
    n_samples = Int
    updated = Event
    traits_view = View( Item( 'nsamples' ) )
    
    def __init__( self, n_samples ):
        
        super( Energy, self ).__init__()
        
        self.n_samples = n_samples
    
    def compute( self, traces ):
        
        data_out = compute_energy( traces,
                                   self.n_samples )
    
        return data_out
    
    def _n_samples_changed(self):
        self.updated = True
    
    
