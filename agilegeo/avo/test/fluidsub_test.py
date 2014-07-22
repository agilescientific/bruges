import unittest
import agilegeo.avo as avo
import numpy as np

class FluidsubTest( unittest.TestCase ):
    """
    Tests fluid sub calculations against Smith et al 2003.
    https://dl.dropboxusercontent.com/u/14965965/Smith_etal_2003.pdf
    """
    
    def test_fluidsub(self):

        vp_gas = 2429.0
        vs_gas = 1462.4
        rho_gas = 2080.

        vp_brine = 2850.5
        vs_shale = 1416.1
        rho_shale = 2210.

        # These are wrong
        phi = 0.2
        rhof1 = 800.
        rhof2 = 1100.
        kmin = 24000000000.
        kf1 = 30000000000.
        kf2 = 50000000000.

        self.assertAlmostEqual( avo.fluidsub(vp=vp_gas, vs=vs_gas, rho=rho_gas, phi=phi, rhof1=rhof1, rhof2=rhof2, kmin=kmin, kf1=kf1, kf2=kf2)[0], vp_brine,  places=0 )        
        self.assertAlmostEqual( avo.fluidsub(vp=vp_gas, vs=vs_gas, rho=rho_gas, phi=phi, rhof1=rhof1, rhof2=rhof2, kmin=kmin, kf1=kf1, kf2=kf2)[1], vs_brine,  places=0 )        
        self.assertAlmostEqual( avo.fluidsub(vp=vp_gas, vs=vs_gas, rho=rho_gas, phi=phi, rhof1=rhof1, rhof2=rhof2, kmin=kmin, kf1=kf1, kf2=kf2)[2], rho_brine, places=0 )        
        
if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(ModuliTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
