#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Testing util.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import unittest
import numpy as np
from bruges.util import next_pow2, rms


class UtilTest(unittest.TestCase):

    def test_nextpow2(self):
        num = 888
        ans = 1024

        self.assertEqual(ans, next_pow2(num))

    def test_rms(self):
        l = np.array([2, 3, 4.5])
        ans = 3.329164

        self.assertAlmostEqual(ans, rms(l), places=5)

    def noise_db(self):
        a = np.ones(1000000)
        ans10p = 3.32
        ans0 = 1.41
        ans10n = 1.05

        self.assertAlmostEqual(ans10p, noise_db(a, 10), places=2)
        self.assertAlmostEqual(ans0, noise_db(a, 10), places=2)
        self.assertAlmostEqual(ans10n, noise_db(a, -10), places=2)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(UtilityTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
