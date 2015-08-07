import unittest
import numpy as np

from bruges.reflection import reflection as avo

vp1 = 12250.
vp2 = 11600.
vs1 = 6200.
vs2 = 6650.
rho1 = 2.66
rho2 = 2.34
theta = np.arange(45)

arr_vp1 = np.ones(1000)*2
arr_vp2 = np.ones(1000)*3
arr_vs1 = np.ones(1000)*4
arr_vs2 = np.ones(1000)*5
arr_rho1 = np.ones(1000)*6
arr_rho2 = np.ones(1000)*7
arr_theta = 0


class AvoTest(unittest.TestCase):
    """
    Tests zoeppritz using a values from a spreadsheet, and also a
    qualitative comparison to plots made by the CREWES avo explorer
    web app. Other algorithms are then tested to be within 10% of the
    zoeppritz answer for angles < 45 degrees.
    """

    tolerance = 0.1

    def test_zoeppritz(self):

        theta = 40.

        reflect = avo.zoeppritz(vp1, vs1, rho1, vp2, vs2, rho2, theta)

        # Number manually verified using
        # spreadsheet from http://tbberge.com/id63.html
        self.assertAlmostEquals(reflect, -0.112236, places=3)

    def test_zoeppritz_rpp(self):

        theta = 40.

        reflect = avo.zoeppritz(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_rpp = avo.zoeppritz_rpp(vp1, vs1, rho1, vp2, vs2, rho2, theta)

        # Should be the same as the exact solution.
        self.assertAlmostEquals(reflect_rpp, reflect, places=3)

    def test_akirichards(self):

        reflect = avo.akirichards(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_zoep = avo.zoeppritz(vp1, vs1, rho1, vp2,
                                     vs2, rho2, theta)

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose(reflect, reflect_zoep,
                           rtol=self.tolerance)
        self.assertTrue(test)

        # Test it won't complain about arrays
        reflect = avo.akirichards(arr_vp1, arr_vs1, arr_rho1,
                                  arr_vp2, arr_vs2, arr_rho2,
                                  arr_theta)

    def test_akirichards_alt(self):

        reflect = avo.akirichards_alt(vp1, vs1, rho1, vp2,
                                      vs2, rho2, theta)
        reflect_zoep = avo.zoeppritz(vp1, vs1, rho1, vp2,
                                     vs2, rho2, theta)

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose(reflect, reflect_zoep,
                           rtol=self.tolerance)
        self.assertTrue(test)

        reflect = avo.akirichards_alt(arr_vp1, arr_vs1, arr_rho1,
                                      arr_vp2, arr_vs2, arr_rho2,
                                      arr_theta)

    def test_fatti(self):
        reflect = avo.fatti(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_zoep = avo.zoeppritz(vp1, vs1, rho1, vp2,
                                     vs2, rho2, theta)

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose(reflect, reflect_zoep,
                           rtol=self.tolerance)
        self.assertTrue(test)

        reflect = avo.fatti(arr_vp1, arr_vs1, arr_rho1,
                            arr_vp2, arr_vs2, arr_rho2, arr_theta)

    def test_shuey2(self):
        reflect = avo.shuey2(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_zoep = avo.zoeppritz(vp1, vs1, rho1, vp2,
                                     vs2, rho2, theta)

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose(reflect, reflect_zoep,
                           rtol=self.tolerance)
        self.assertTrue(test)

    def test_shuey3(self):
        reflect = avo.shuey3(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_zoep = avo.zoeppritz(vp1, vs1, rho1, vp2,
                                     vs2, rho2, theta)

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose(reflect, reflect_zoep,
                           rtol=self.tolerance)
        self.assertTrue(test)

    def test_bortfeld2(self):
        reflect = avo.bortfeld2(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_zoep = avo.zoeppritz(vp1, vs1, rho1, vp2,
                                     vs2, rho2, theta)

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose(reflect, reflect_zoep,
                           rtol=self.tolerance)
        self.assertTrue(test)

    def test_bortfeld3(self):
        reflect = avo.bortfeld3(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_zoep = avo.zoeppritz(vp1, vs1, rho1, vp2,
                                     vs2, rho2, theta)

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose(reflect, reflect_zoep,
                           rtol=self.tolerance)
        self.assertTrue(test)

    def test_hilterman(self):
        reflect = avo.hilterman(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_zoep = avo.zoeppritz(vp1, vs1, rho1, vp2,
                                     vs2, rho2, theta)

        # See if it is within .1 of zoep for < 45 deg
        test = np.allclose(reflect, reflect_zoep,
                           rtol=self.tolerance)
        self.assertTrue(test)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(AvoTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
