# -*- coding: utf-8 -*-
"""
Tests.
"""
import unittest
import numpy as np

from bruges.rockphysics import elastic_impedance

vp = 2350.
vs = 1125.
rho = 2500.

theta2 = np.array([0.0, 10.0, 20.0, 30.0])

ones_vp = np.ones(1000) * vp
ones_vs = np.ones(1000) * vs
ones_rho = np.ones(1000) * rho

arr_vp = np.array([2350, 2690, 2430, 2400, 3700])
arr_vs = np.array([1125, 1300, 1240.0, 1180, 1900.0])
arr_rho = np.array([2500, 2550, 2480, 2420, 2690])

arr_vp_nan = arr_vp.astype(float)
arr_vp_nan[-2] = np.nan


class ElasticImpedanceTest(unittest.TestCase):
    """
    Tests moduli calculation against spreadsheet
    https://dl.dropboxusercontent.com/u/14965965/Elastic_moduli_formula_checking.xlsx
    """

    def test_ei_at_zero(self):
        """
        Checks that EI at zero degrees is equal to AI.
        """
        ei = elastic_impedance(vp, vs, rho, theta1=0)
        self.assertAlmostEqual(ei, vp * rho, places=-2)

    def test_ei_ninety(self):
        """
        Checks that EI at 90 degrees = (Vp/Vs)^2. This occurs when K = 0.25
        """
        ei = elastic_impedance(vp, vp/2.0, rho, theta1=90.0, use_sin=True)
        self.assertAlmostEqual(ei, (vp/(vp/2.0))**2, places=-2)

    def test_ei_log(self):
        """
        Checks that passing in 1-D arrays for vp, vs, and rho works for
        single values of theta1 as well as a 1-d vector.
        """
        x1 = elastic_impedance(ones_vp, ones_vs, ones_rho, theta1=np.array(30.0))
        x2 = elastic_impedance(ones_vp, ones_vs, ones_rho, theta1=theta2)
        x3 = elastic_impedance(arr_vp, arr_vs, arr_rho, theta1=theta2)

        self.assertEqual(x1.shape, (ones_vp.size,))
        self.assertEqual(x2.shape, (ones_vp.size, theta2.size))
        self.assertEqual(x3.shape, (arr_vp.size, theta2.size))

    def test_ei_nan(self):
        """
        Checks to see if we survuve NaNs.
        """
        ei = elastic_impedance(2350, np.nan, 2500, theta1=30.0)
        self.assertTrue(np.isnan(ei))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ElasticImpedanceTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
