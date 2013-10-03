from traits.api import HasTraits, Int, Event, Range
from traitsui.api import View, Item
import numpy 

def compute_similarity(traces, n_samples, step_out):
    """ Compute an similarity for each point of a seismic cube.
    """

    half_n_samples = n_samples // 2

    similarity_cube = numpy.zeros_like(traces, dtype='float32')
    traces = numpy.nan_to_num(traces)

    for i in range( step_out ):
        for idx in xrange(similarity_cube.shape[0]):
            start_idx = max(0, idx-half_n_samples)
            stop_idx = min(similarity_cube.shape[0]-1, idx+half_n_samples)
            signal = traces[start_idx:stop_idx, :]

            squares = (signal*signal).sum(axis=0)
            squares_of_diff = ((signal[:,1+i:] - signal[:,:-(1+i)])**2.).sum(axis=0)

            sqrt=numpy.sqrt
            squares[squares==0.0] = 0.001
            squares_of_diff[squares_of_diff==0.0] = 0.001
            sim = 1.0 -  sqrt(squares_of_diff) / (sqrt(squares[1+i:]) + sqrt(squares[:-(1+i)]))

            similarity_cube[idx,i+1:] += sim

    return (similarity_cube / step_out )
    
class Similarity( HasTraits ):
    
    n_samples = Range( 1,100,3)
    step_out = Range(1,10,3)
    updated = Event
    traits_view = View( Item( 'n_samples' ), Item( 'step_out'))
    
    def __init__( self, n_samples ):
        
        self.n_samples = n_samples
        
    def _n_samples_changed(self):
        self.updated=True
        
    def _step_out_changed(self):
        self.updated=True
    def compute( self, traces ):
        
        return compute_similarity( traces, self.n_samples, self.step_out) 
    
    