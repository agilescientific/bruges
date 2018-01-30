# -*- coding: utf-8 -*-
import unittest
import numpy as np
from numpy import array

from bruges.filters import snn
from bruges.filters import kuwahara
from bruges.filters import conservative
from bruges.filters import mean
from bruges.filters import rms
from bruges.filters import median
from bruges.filters import mode

horizon = np.arange(25).reshape((5, 5))
horizon[2, 2] = -1


class FilterTest(unittest.TestCase):
    """
    Tests filters on a very small 2D horizon-like array.
    """
    def test_snn(self):
        x = array([[ 0.2,  0.8,  1.8,  2.8,  3.8],
                   [ 2.2,  6.6,  8.6,  9.,   9.6],
                   [ 7.2,  7.6,  8.2,  8.8,  9.2],
                   [12.,  16.4, 19.8, 20.8, 21.8],
                   [20.2, 20.8, 21.8, 22.8, 23.8]])
        y = snn(horizon, size=3)
        self.assertTrue((x == y).all())

    def test_kuwahara(self):
        x = array( [[ 0.,    0.5,   1.5,   2.5,   4.  ],
                    [ 2.5,   3.,    4.,    5.,    6.5 ],
                    [ 7.5,   8.,    5.75, 11.,   11.5 ],
                    [12.5,  13.,   19.,   16.,   16.5 ],
                    [20.,   20.5,  21.5,  22.5,  24.  ]])
        y = kuwahara(horizon, size=3)
        self.assertTrue((x == y).all())

    def test_conservative(self):
        x = array([[ 0,  1,  2,  3,  4],
                   [ 5,  6,  7,  8,  9],
                   [10, 11,  6, 13, 14],
                   [15, 16, 17, 18, 19],
                   [20, 21, 22, 23, 24]])
        y = conservative(horizon, size=3, supercon=True)
        self.assertTrue((x == y).all())

    def test_mean(self):
        x = array([[ 4.28,  4.68,  5.48,  6.28,  6.68],
                   [ 6.28,  6.68,  7.48,  8.28,  8.68],
                   [10.28, 10.68, 11.48, 12.28, 12.68],
                   [14.28, 14.68, 15.48, 16.28, 16.68],
                   [16.28, 16.68, 17.48, 18.28, 18.68]])
        y = mean(horizon)
        self.assertTrue((x == y).all())

    def test_rms(self):
        x = array([[ 3.12694384,  3.65148372,  4.43471157,  5.29150262,  5.84997626],
                   [ 6.55743852,  6.79051626,  7.17247826,  7.62306442,  7.93725393],
                   [11.01514109, 11.45037602, 12.10142324, 12.77149604, 13.21615173],
                   [15.80084386, 16.3129124,  17.07174404, 17.83566963, 18.33939294],
                   [18.82079228, 19.49358869, 20.48576742, 21.47867159, 22.13092356]])
        y = rms(horizon, 3)
        test = np.allclose(x, y, rtol=1e-7)
        self.assertTrue(test)

    def test_median(self):
        x = array([[ 1.,  2.,  3.,  4.,  4.],
                   [ 5.,  5.,  6.,  7.,  9.],
                   [10., 10., 11., 13., 14.],
                   [15., 16., 17., 18., 19.],
                   [20., 20., 21., 22., 23.]])
        y = median(horizon, 3)
        self.assertTrue((x == y).all())

    def test_mode(self):
        x = array([[ 0.,  0.,  1.,  2.,  4.],
                   [10., 10., 10., 10., 10.],
                   [10., 10., 10., 10., 10.],
                   [10., 10., 10., 10., 10.],
                   [20., 20., 21., 22., 24.]])
        horizon[2] = np.array([10., 10., 10., 10., 10.])
        y = mode(horizon, 3)
        self.assertTrue((x == y).all())


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(FilterTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
