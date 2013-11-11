import unittest
from agilegeo.util import next_pow2

class UtilityTest( unittest.TestCase ):

    def test_nextpow2( self ):

        num = 888
        ans = 1024

        self.assertEqual( ans, next_pow2( num ) )

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(UtilityTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
