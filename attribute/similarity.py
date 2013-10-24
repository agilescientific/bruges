from traits.api import HasTraits, Int, Event, Range
from traitsui.api import View, Item
import numpy 

def compute_similarity(traces, n_samples, step_out=1, lag=0):
    """
    Compute similarity for each point of a seismic section using
    a variance method between traces.  

    For each point, a kernel of n_samples length is extracted from a
    trace. The similarity is calculated as a normalized variance
    between two adjacent trace sections, where a value of 1 is
    obtained by identical if the traces are identical. The step out
    will decide how many adjacent traces will be used for each kernel,
    and should be increased for poor quality data. The lag determines
    how much neighbouring traces can be shifted when calculating
    similiarity, which should be increased for dipping data.
    
    :param traces: A 2D numpy array arranged as [time, trace].
    :param n_samples: The length in samples of the trace kernel
                      used to calculate the similarity.
    :keyword step_out (default=1 ):
                       The number of adjacent traces to 
                       the kernel to check similarity. The maximum
                       similarity value will be chosen. 
    :keyword lag (default=0):
                  The maximum number of time samples adjacent traces
                  can be shifted by. The maximum similarity of
                  will be used.
                  
                    
    """

    half_n_samples = n_samples // 2

    similarity_cube = numpy.zeros_like(traces, dtype='float32')
    traces = numpy.nan_to_num(traces)

    for j in range( lag + 1 ):
        
        for i in range( step_out ):
            for idx in xrange(similarity_cube.shape[0]):

                start_idx = max(0, (idx+(j*i)-half_n_samples))
                stop_idx = min(similarity_cube.shape[0]-1, \
                               (idx-(i*j))+half_n_samples)
             
                signal = traces[start_idx:stop_idx, :]

                squares = (signal*signal).sum(axis=0)
                squares_of_diff = ((signal[:,1+i:] -
                                signal[:,:-(1+i)])**2.).sum(axis=0)

                sqrt=numpy.sqrt
                squares[squares==0.0] = 0.001
                squares_of_diff[squares_of_diff==0.0] = 0.001
                sim = 1.0 -  sqrt(squares_of_diff) / \
                                 (sqrt(squares[1+i:]) \
                                 + sqrt(squares[:-(1+i)]))

             
                similarity_cube[idx,(i+1):] = \
                      numpy.maximum(sim,similarity_cube[idx, (i+1):])

    return (similarity_cube)
    
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
    
    
