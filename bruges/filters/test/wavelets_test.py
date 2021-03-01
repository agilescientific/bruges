# -*- coding: utf-8 -*-
import unittest
import numpy as np
from numpy import array

from bruges.filters import sinc
from bruges.filters import ricker
from bruges.filters import sweep
from bruges.filters import berlage
from bruges.filters import ormsby


class WaveletTest(unittest.TestCase):
    """
    Tests wavelets.
    """
    def taper(self, l):
        w = np.zeros(l)
        w[int(l/4):-int(l/4)] = 1
        return w

    def test_sinc(self):
        x = array([-0.        ,  0.        ,  0.        ,  0.        , -0.15321277,
                   -0.1887466 ,  0.17251039,  0.72943393,  1.        ,  0.72943393,
                    0.17251039, -0.1887466 , -0.        ,  0.        ,  0.        ,  0.        ])
        y = sinc(0.064, 0.004, np.arange(50, 100), taper=self.taper)[3]
        self.assertTrue(np.allclose(x, y))

    def test_ricker(self):
        x = array([0.88273185,  0.90950743,  0.9330604 ,  0.95324475,  0.96993473,
                   0.9830259 ,  0.99243608,  0.99810603,  1.        ,  0.99810603,
                   0.99243608,  0.9830259 ,  0.96993473,  0.95324475,  0.9330604 ,
                   0.90950743])
        y = ricker(0.064, 0.004, np.arange(1, 10))[1]
        self.assertTrue(np.allclose(x, y))

    def test_sweep(self):
        x = array([-1.59650205e-18,   2.75409007e-03,  -3.01670223e-03,
                   -3.63692317e-02,  -1.63002095e-02,  -5.20943260e-02,
                   -2.83497662e-01,   2.76098252e-01,   9.82157437e-01,
                    2.38730442e-01,  -2.10312350e-01,  -3.25806534e-02,
                   -8.30581776e-03,  -1.39615614e-02,  -6.55916221e-04,
                   -2.28078034e-18])
        y = sweep(0.064, 0.004, [10, 100])
        self.assertTrue(np.allclose(x, y))

    def test_berlage(self):
        x = array([ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -0.00000000e+00,
                   -0.00000000e+00, -0.00000000e+00, -0.00000000e+00, -0.00000000e+00,
                    0.00000000e+00,  2.89837052e-01,  9.13081427e-01,  1.00000000e+00,
                    5.34807890e-01, -5.72248038e-16, -2.85099378e-01, -3.05622598e-01])
        y = berlage(0.064, 0.004, 25)
        self.assertTrue(np.allclose(x, y))

    def test_ormsby(self):
        x = array([ 0.03109604, -0.05796127, -0.15655105, -0.06436394, -0.06901495,
                   -0.39791603, -0.37263257,  0.43754754,  1.        ,  0.43754754,
                   -0.37263257, -0.39791603, -0.06901495, -0.06436394, -0.15655105,
                   -0.05796127])
        y = ormsby(0.064, 0.004, [10, 20, 60, 80])
        self.assertTrue(np.allclose(x, y))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(WaveletTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
