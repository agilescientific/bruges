import unittest
import numpy
from agilegeo.attribute import compute_similarity

class SimilarityTest( unittest.TestCase ):

    def test_same_data( self ):
        """
        Simple test to check if the algorithm works for the
        trivial case.
        """
        data = numpy.zeros( [1000, 1000] ) 
        check_data = data + 1.0
        
        data +=10.
        window_size = 20
        
        output = compute_similarity( data, window_size )
        same = numpy.allclose( check_data[1:,:], 
                               output[1:,:], .01 )
        
        self.assertTrue( same )


if __name__ == '__main__':
    unittest.main()
    