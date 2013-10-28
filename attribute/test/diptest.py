import unittest
from numpy import linspace, sin, pi, amax
from agilegeo.attribute import dipsteer

class DipTest( unittest.TestCase ):


    def test_spikes( self ):
        """
        Make a simple case with spikes offset on each trace.
        """

        test_data = np.zeros( (10,10 ) )

        index = np.arange( 0,10 )

        test_data[index,index] = 1.0

        dip, cor = dipsteer( test_data, 10,1,3 )


        self.assertTrue( dip[5,5] == 1.0 )
    

        
        
        
        
