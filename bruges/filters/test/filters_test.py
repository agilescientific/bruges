# -*- coding: utf 8 -*-
import unittest
import numpy as np
from numpy import array

from bruges.filters import snn
from bruges.filters import kuwahara
from bruges.filters import conservative

horizon = np.arange(25).reshape((5, 5))


class FilterTest(unittest.TestCase):
    """
    Tests filters on a very small 2D horizon-like array.
    """
    def test_snn(self):

        x = array([[ 0,  0,  1,  2,  3],
                   [ 2,  2,  3,  4,  9],
                   [ 7,  7,  6, 10, 14],
                   [12, 14, 15, 17, 19],
                   [20, 20, 21, 22, 23]])
        horizon[2, 2] = -1
        y = snn(horizon, size=3)
        self.assertTrue((x == y).all())

    def test_kuwahara(self):

        x = array([[ 0,  0,  1,  2,  4],
                   [ 2,  3,  4,  5,  6],
                   [ 7,  8,  5, 11, 11],
                   [12, 13, 19, 16, 16],
                   [20, 20, 21, 22, 24]])
        horizon[2, 2] = -1
        y = kuwahara(horizon, size=3)
        self.assertTrue((x == y).all())

    def test_conservative(self):

        x = array([[ 0,  1,  2,  3,  4],
                   [ 5,  6,  7,  8,  9],
                   [10, 11,  6, 13, 14],
                   [15, 16, 17, 18, 19],
                   [20, 21, 22, 23, 24]])
        horizon[2, 2] = -1
        y = conservative(horizon, size=3, supercon=True)
        self.assertTrue((x == y).all())


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(FilterTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
