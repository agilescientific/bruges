#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Tests.
"""
import unittest

from bruges.rockphysics import bounds as b

# Inputs.
f = [0.5, 0.5]
k = [36.0, 4.0]
mu = [45, 5.0]

# Expected outputs.
mv = 20.0
mr = 7.199
mh = 13.599


class BoundsTest(unittest.TestCase):
    """
    Tests calculations of elastic mixture bounds
    """

    def test_voigt(self):
        self.assertAlmostEqual(b.voight(f, k), mv, places=-2)

    def test_reuss(self):
        self.assertAlmostEqual(b.reuss(f, k), mr, places=-2)

    def test_hill_average(self):
        self.assertAlmostEqual(b.hill_average(f, k), mh, places=-2)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(BoundsTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
