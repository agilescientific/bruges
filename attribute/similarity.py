from traits.api import HasTraits, Int
from traitsui.api import View, Item
import numpy 

def compute_similarity(traces, n_samples):
    """ Compute an similarity for each point of a seismic cube.
    """

    half_n_samples = n_samples // 2

    similarity_cube = numpy.zeros_like(traces, dtype='float32')
    traces = numpy.nan_to_num(traces)

    for idx in xrange(similarity_cube.shape[0]):
        start_idx = max(0, idx-half_n_samples)
        stop_idx = min(similarity_cube.shape[0]-1, idx+half_n_samples)
        signal = traces[start_idx:stop_idx, :]

        squares = (signal*signal).sum(axis=0)
        squares_of_diff = ((signal[:,1:] - signal[:,:-1])**2.).sum(axis=0)

        sqrt=numpy.sqrt
        squares[squares==0.0] = 0.001
        squares_of_diff[squares_of_diff==0.0] = 0.001
        sim = 1.0 -  sqrt(squares_of_diff) / (sqrt(squares[1:]) + sqrt(squares[:-1]))

        similarity_cube[idx,1:] = sim

    return similarity_cube
    
class Similarity( HasTraits ):
    
    n_samples = Int
    traits_view = View( Item( 'n_samples' ))
    
    def __init__( self, n_samples ):
        
        self.n_samples = n_samples
        
    def compute( self, traces ):
        
        return compute_similarity( traces, self.n_samples ) 
    
    