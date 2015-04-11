import unittest
import agilegeo.avo as avo
import numpy as np

class AvoTest( unittest.TestCase ):
    """
    Tests zoeppritz using a values from a spreadsheet, and also a
    qualitative comparison to plots made by the CREWES avo explorer
    web app. Other algorithms are then tested to be within 10% of the
    zoeppritz answer for angles < 45 degrees.
    """
    
    # 5 percent tolerance
    tolerance = 0.1
    
    def test_zoeppritz(self):

        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = 40.


        reflect = avo.zoeppritz( vp1,vs1,rho1,vp2,vs2,rho2,theta )

        # Number manually verified using
        # spreadsheet from http://tbberge.com/id63.html
        self.assertAlmostEquals( reflect, -0.112236, places=3 )
        
    
    def test_akirichards(self):
        
        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = np.arange( 45 )


        reflect = avo.akirichards( vp1,vs1,rho1,vp2,vs2,rho2,theta )
        reflect_zoep = avo.zoeppritz( vp1,vs1,rho1,vp2,
                                      vs2,rho2,theta )

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose( reflect, reflect_zoep,
                            rtol=self.tolerance )
        self.assertTrue(test)

        # Test it won't complain about arrays
        

        
        vp1 = np.ones(1000)*2
        vp2 = np.ones(1000)*3

        vs1 = np.ones(1000)*4
        vs2 = np.ones(1000)*5

        rho1 = np.ones(1000)*6
        rho2 = np.ones(1000)*7
        theta = 0
        reflect = avo.akirichards(vp1,vs1,rho1,vp2,vs2,rho2,theta)
        

    def test_akirichards_alt(self):

        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = np.arange( 45 )


        reflect = avo.akirichards_alt( vp1,vs1,rho1,vp2,
                                        vs2,rho2,theta )
        reflect_zoep = avo.zoeppritz( vp1,vs1,rho1,vp2,
                                      vs2,rho2,theta )

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose( reflect, reflect_zoep,
                            rtol=self.tolerance )
        self.assertTrue( test )

        vp1 = np.ones(1000)*2
        vp2 = np.ones(1000)*3

        vs1 = np.ones(1000)*4
        vs2 = np.ones(1000)*5

        rho1 = np.ones(1000)*6
        rho2 = np.ones(1000)*7
        theta = 0
        reflect = avo.akirichards_alt(vp1,vs1,rho1,vp2,vs2,rho2,theta)
        
    def test_fatti(self):

        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = np.arange( 45 )


        reflect = avo.fatti( vp1,vs1,rho1,vp2,vs2,rho2,theta )

        reflect_zoep = avo.zoeppritz( vp1,vs1,rho1,vp2,
                                      vs2,rho2,theta )

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose( reflect, reflect_zoep,
                            rtol=self.tolerance )
        self.assertTrue( test )

        vp1 = np.ones(1000)*2
        vp2 = np.ones(1000)*3

        vs1 = np.ones(1000)*4
        vs2 = np.ones(1000)*5

        rho1 = np.ones(1000)*6
        rho2 = np.ones(1000)*7

        theta = 0
        reflect = avo.fatti(vp1,vs1,rho1,vp2,vs2,rho2,theta)
        
    def test_shuey2(self):

        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = np.arange( 45 )


        reflect = avo.shuey2( vp1,vs1,rho1,vp2,vs2,rho2,theta )
        reflect_zoep = avo.zoeppritz( vp1,vs1,rho1,vp2,
                                      vs2,rho2,theta )

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose( reflect, reflect_zoep,
                            rtol=self.tolerance )
        self.assertTrue( test )
        

    def test_shuey3(self):
        
    
        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = np.arange(45)


        reflect = avo.shuey3( vp1,vs1,rho1,vp2,vs2,rho2,theta )

        reflect_zoep = avo.zoeppritz( vp1,vs1,rho1,vp2,
                                      vs2,rho2,theta )

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose( reflect, reflect_zoep,
                            rtol=self.tolerance )
        self.assertTrue( test )
        

    def test_bortfeld2(self):
        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = np.arange(45)


        reflect = avo.bortfeld2( vp1,vs1,rho1,vp2,vs2,rho2,theta )

        reflect_zoep = avo.zoeppritz( vp1,vs1,rho1,vp2,
                                      vs2,rho2,theta )

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose( reflect, reflect_zoep,
                            rtol=self.tolerance )
        self.assertTrue( test )

    def test_bortfeld3(self):
        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = np.arange(45)


        reflect = avo.bortfeld3( vp1,vs1,rho1,vp2,vs2,rho2,theta )

        reflect_zoep = avo.zoeppritz( vp1,vs1,rho1,vp2,
                                      vs2,rho2,theta )

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose( reflect, reflect_zoep,
                            rtol=self.tolerance )
        self.assertTrue( test )
        
if __name__ == '__main__':

    suite = \
      unittest.TestLoader().loadTestsFromTestCase(AvoTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
