# -*- coding: utf-8 -*-
import unittest
import numpy as np

from bruges.reflection import reflection as avo

vp1 = 12250.
vp2 = 11600.
vs1 = 6200.
vs2 = 6650.
rho1 = 2.66
rho2 = 2.34
theta = 25

arr_vp1 = np.ones(1000)*2
arr_vp2 = np.ones(1000)*3
arr_vs1 = np.ones(1000)*1
arr_vs2 = np.ones(1000)*1.5
arr_rho1 = np.ones(1000)*2.5
arr_rho2 = np.ones(1000)*3.5
arr_theta = np.arange(40)


class AvoTest(unittest.TestCase):
    """
    Tests zoeppritz using a values from a spreadsheet, and also a
    qualitative comparison to plots made by the CREWES avo explorer
    web app. Other algorithms are then tested to be within 1% of the
    zoeppritz answer for angles < 40 degrees.
    """

    tolerance = 0.01
    reflect_rpp = avo.zoeppritz_rpp(vp1, vs1, rho1, vp2, vs2, rho2, theta)

    def test_zoeppritz(self):
        theta = 40.
        reflect = avo.zoeppritz(vp1, vs1, rho1, vp2, vs2, rho2, theta)

        # Number manually verified using
        # spreadsheet from http://tbberge.com/id63.html
        self.assertAlmostEquals(reflect, -0.112236, places=5)

    def test_zoeppritz_rpp(self):
        theta = 40.
        reflect = avo.zoeppritz(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_rpp = avo.zoeppritz_rpp(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        self.assertAlmostEquals(reflect, reflect_rpp, places=5)

    def test_akirichards(self):
        reflect = avo.akirichards(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        test = np.allclose(reflect, self.reflect_rpp, rtol=self.tolerance)
        self.assertTrue(test)

        # Test it won't complain about arrays.
        reflect = avo.akirichards(arr_vp1, arr_vs1, arr_rho1,
                                  arr_vp2, arr_vs2, arr_rho2,
                                  arr_theta)
        self.assertTrue(reflect.shape == (40, 1000))

    def test_akirichards_alt(self):
        reflect = avo.akirichards_alt(vp1, vs1, rho1, vp2,
                                      vs2, rho2, theta)
        test = np.allclose(reflect, self.reflect_rpp, rtol=self.tolerance)
        self.assertTrue(test)

        reflect = avo.akirichards_alt(arr_vp1, arr_vs1, arr_rho1,
                                      arr_vp2, arr_vs2, arr_rho2,
                                      arr_theta)
        self.assertTrue(reflect.shape == (40, 1000))

    def test_fatti(self):
        reflect = avo.fatti(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        test = np.allclose(reflect, self.reflect_rpp, rtol=self.tolerance)
        self.assertTrue(test)
        reflect = avo.fatti(arr_vp1, arr_vs1, arr_rho1,
                            arr_vp2, arr_vs2, arr_rho2, arr_theta)

    def test_shuey(self):
        reflect = avo.shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        test = np.allclose(reflect, self.reflect_rpp, rtol=self.tolerance)
        self.assertTrue(test)

    def test_bortfeld(self):
        reflect = avo.bortfeld(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        test = np.allclose(reflect, self.reflect_rpp, rtol=self.tolerance)
        self.assertTrue(test)

    def test_hilterman(self):
        """
        Does not pass at theta = 25 deg.
        """
        theta = np.arange(10)
        reflect = avo.hilterman(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        reflect_rpp = avo.zoeppritz_rpp(vp1, vs1, rho1, vp2, vs2, rho2, theta)
        test = np.allclose(reflect, reflect_rpp, rtol=self.tolerance)
        self.assertTrue(test)

    def test_ca(self):
        """
        Test critical angles.
        """
        ca1, ca2 = avo.critical_angles(vp1=2300, vp2=2400, vs2=None)
        self.assertAlmostEquals(ca1, 73.40215786, places=5)
        self.assertTrue(np.isnan(ca2))

        ca1, ca2 = avo.critical_angles(vp1=1500, vp2=3200, vs2=1600)
        self.assertAlmostEquals(ca1, 27.9531869, places=5)
        self.assertAlmostEquals(ca2, 69.6358652, places=5)


if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(AvoTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
