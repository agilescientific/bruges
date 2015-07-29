#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Tests.
"""
import unittest

from bruges.rockphysics import bounds as b

# Inputs.
f0 = [0.5, 0.5]
f1 = [50, 50]
f2 = [0.6, 0.4]    # suspension of quartz in water
k1 = [36.0, 4.0]   # Kquartz, Kclay
k2 = [36.0, 2.2]   # Kquartz, Kwater


# Expected outputs.
mv = 20.0    # Voigt average (hand-calculated)
mr = 7.199   # Reuss average (hand-calculated)
mr2 = 5.04   # Example from Mavko et al. RPH
mh = 13.599  # Voigt-Reuss Hill average (hand-calculated)


class BoundsTest(unittest.TestCase):
    """
    Tests calculations of elastic mixture bounds
    """

    def test_voigt(self):
        self.assertAlmostEqual(b.voigt_bound(f0, k1), mv, places=-2)
        self.assertAlmostEqual(b.voigt_bound(f1, k1), mv, places=-2)

    def test_reuss(self):
        self.assertAlmostEqual(b.reuss_bound(f0, k1), mr, places=-2)
        self.assertAlmostEqual(b.reuss_bound(f2, k2), mr2, places=-2)
        self.assertAlmostEqual(b.reuss_bound([0.6, 0.2, 0.2],
                                             [36.0, 2.2, 0.000131]), 0.00065, places=-2)

    def test_hill_average(self):
        self.assertAlmostEqual(b.hill_average(f0, k1), mh, places=-2)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(BoundsTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
