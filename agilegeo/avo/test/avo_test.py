import unittest
import agilegeo.avo as avo

class AvoTest( unittest.TestCase ):

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

        theta = 40.


        reflect = avo.akirichards( vp1,vs1,rho1,vp2,vs2,rho2,theta )

        self.assertAlmostEquals( reflect, -0.112236, places=3 )

    def test_akirichards_alt(self):

        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = 40.


        reflect = avo.akirichards_alt( vp1,vs1,rho1,vp2,
                                        vs2,rho2,theta )

        self.assertAlmostEquals( reflect, -0.112236, places = 3 )

    def test_fatti(self):

        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = 40.


        reflect = avo.fatti( vp1,vs1,rho1,vp2,vs2,rho2,theta )

        self.assertAlmostEquals( reflect, -0.11155, places = 3 )

    def test_shuey2(self):

        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = 40.


        reflect = avo.shuey2( vp1,vs1,rho1,vp2,vs2,rho2,theta )

        self.assertAlmostEquals( reflect, -0.08955, places = 3 )

    def test_shuey3(self):
        
    
        vp1 = 12250.
        vp2 = 11600.

        vs1 = 6200.
        vs2 = 6650.

        rho1 = 2.66
        rho2 = 2.34

        theta = 40.


        reflect = avo.shuey3( vp1,vs1,rho1,vp2,vs2,rho2,theta )

        self.assertAlmostEquals( reflect, -0.0975, places = 3 )

   # def test_borfield2(self):

    #def test_borfield3(self):
        
if __name__ == '__main__':
    # sys.path.append( '~/agilegeo')
    unittest.main()
