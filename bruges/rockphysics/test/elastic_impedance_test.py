#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Tests.
"""
import unittest

from bruges.rockphysics import elastic as elas

vp = 2350.
vs = 1125.
rho = 2500.

class ElasticImpedanceTest(unittest.TestCase):
    """
    Tests moduli calculation against spreadsheet
    https://dl.dropboxusercontent.com/u/14965965/Elastic_moduli_formula_checking.xlsx
    """

    def test_ei_at_zero(self):
        """
        Checks to see if single elastic impedance at zero degrees is equal to acoustic impedance 
        for a) an integer of zero, b) a float of zero, c) an array with one element (zero)
        and d) a list whose value is zero
        """
        self.assertAlmostEqual(elas.elastic_impedance(vp=vp, vs=vs, rho=rho, theta1=0),
                               vp * rho, places=-2)
        self.assertAlmostEqual(elas.elastic_impedance(vp=vp, vs=vs, rho=rho, theta1=np.array([0])),
                               vp * rho, places=-2)
        self.assertAlmostEqual(elas.elastic_impedance(vp=vp, vs=vs, rho=rho, theta1=np.array([0])),
                               vp * rho, places=-2)
        self.assertAlmostEqual(elas.elastic_impedance(vp=vp, vs=vs, rho=rho, theta1=[0]),
                               vp * rho, places=-2)
        self.assertAlmostEqual(elas.elastic_impedance(vp=vp, vs=vs, rho=rho, theta1=[0.0]),
                               vp * rho, places=-2)
    
    def test_ei_ninety(self):
        """
        Checks to see if ei at 90 degrees = (Vp/Vs)^2. This occurs when K = 0.25
        """
        self.assertAlmostEqual(elas.elastic_impedance(vp=vp, vs=vp/2.0, rho=rho, theta1=90.0),
                               (vp/vs)**2, places=-2)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ElasticImpedanceTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
