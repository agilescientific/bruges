#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Tests.
"""
import unittest

from bruges.rockphysics import bounds as b

# Inputs.
f1 = [0.5, 0.5]
f2 = [50, 50]
k1 = [36.0, 4.0]


# Expected outputs.
mv = 20.0    # Voigt average (hand-calculated)
mr = 7.199   # Reuss average (hand-calculated)
mh = 13.599  # Voigt-Reuss Hill average (hand-calculated)


class BoundsTest(unittest.TestCase):
    """
    Tests calculations of elastic mixture bounds
    """

    def test_voigt(self):
        self.assertAlmostEqual(b.voigt_bound(f1, k1), mv, places=-2)
        self.assertAlmostEqual(b.voigt_bound(f2, k1), mv, places=-2)

    def test_reuss(self):
        self.assertAlmostEqual(b.reuss_bound(f1, k1), mr, places=-2)

    def test_hill_average(self):
        self.assertAlmostEqual(b.hill_average(f1, k1), mh, places=-2)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(BoundsTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
