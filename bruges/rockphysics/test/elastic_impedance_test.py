#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Tests.
"""
import unittest
import numpy as np

from bruges.rockphysics import elastic as elas

vp = 2350.
vs = 1125.
rho = 2500.

theta1 = np.arange(0,89)
theta2 = np.array([0.0, 10.0, 20.0, 30.0])
theta3 = theta2
theta3[-1] = np.nan

ones_vp = np.ones(1000) * vp 
ones_vs = np.ones(1000) * vs 
ones_rho = np.ones(1000) * rho

arr_vp = np.array([2350, 2690, 2430, 2400, 3700])
arr_vs = np.array([1125, 1300, 1240.0, 1180, 1900.0])
arr_rho =np.array([2500, 2550, 2480, 2420, 2690])

arr_vp_nan = arr_vp.astype(float)
arr_vp_nan[-2] = np.nan


class ElasticImpedanceTest(unittest.TestCase):
    """
    Tests moduli calculation against spreadsheet
    https://dl.dropboxusercontent.com/u/14965965/Elastic_moduli_formula_checking.xlsx
    """

    def test_ei_at_zero(self):
        """
        Checks to see if single elastic impedance at zero degrees is equal to acoustic impedance 
        for an integer of zero
        """
        self.assertAlmostEqual(elas.elastic_impedance(vp=vp, vs=vs, rho=rho, theta1=0, use_sin=True),
                               vp * rho, places=-2)
    
    def test_ei_ninety(self):
        """
        Checks to see if ei at 90 degrees = (Vp/Vs)^2. This occurs when K = 0.25
        """
        self.assertAlmostEqual(elas.elastic_impedance(vp=vp, vs=vp/2.0, rho=rho, theta1=90.0, use_sin=True),
                               (vp/(vp/2))**2, places=-2)

    def test_ei_log(self):
        """
        Checks that passing in 1-D arrays for vp, vs, and rho works for
        single values of theta1 as well as the a 1-d vector
        """
        x1 = elas.elastic_impedance(vp=ones_vp, vs=ones_vs, rho=ones_rho, theta1=np.array(30.0))
        x2 = elas.elastic_impedance(vp=ones_vp, vs=ones_vs, rho=ones_rho, theta1=theta2)
        x3 = elas.elastic_impedance(vp=arr_vp, vs=arr_vs, rho=arr_rho, theta1=theta2)

        self.assertEqual(x1.shape, (len(ones_vp),))
        self.assertEqual(x2.shape, (len(ones_vp), len(theta2)))
        self.assertEqual(x3.shape, (len(arr_vp), len(theta2)))

    def test_ei_nan(self):
        """
        Checks to see if ei at 90 degrees = (Vp/Vs)^2. This occurs when K = 0.25
        """
        self.assertIsNaN(elas.elastic_impedance(vp=2350, vs=np.nan, rho=2500, theta1=30.0))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ElasticImpedanceTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
