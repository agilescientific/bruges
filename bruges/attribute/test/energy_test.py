# -*- coding: utf-8 -*-
"""
Tests.
"""
import unittest
from numpy import linspace, sin, pi, amax

from bruges.attribute import energy


class EnergyTest(unittest.TestCase):

    def setUp(self):
        """
        Makes a simple sin wave with 1 amplitude to use as test data.
        """
        self.n_samples = 1001
        duration = 1.0
        freq = 10.0
        w = freq * 2 * pi
        t = linspace(0.0, duration, self.n_samples)

        self.data = sin(w * t)

    def test_amplitude(self):
        """
        Tests the basic algorithm returns the right amplitude
        location.
        """
        amplitude = energy(self.data, self.n_samples)
        max_amp = amax(amplitude)

        ms_sin = 0.5  # The MS energy of a sin wave
        self.assertAlmostEquals(ms_sin, max_amp, places=3)

        # Check that it is in the right location
        self.assertAlmostEqual(max_amp, amplitude[501], places=3)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(EnergyTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
