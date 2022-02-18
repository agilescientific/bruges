import unittest
import numpy as np

from bruges.filters import convolve


class ConvolveTest(unittest.TestCase):

    boxcar = np.array([1, 1, 1, 1, 1, 1, 1]) / 7
    sinc = np.array([-0.        ,  0.        ,  0.        ,  0.        , -0.15321277,
                  -0.1887466 ,  0.17251039,  0.72943393,  1.        ,  0.72943393,
                   0.17251039, -0.1887466 , -0.15321277,  0.        ,  0.        ,
                   0.        , -0.        ])

    def test_1d(self):
        rng = np.random.default_rng(42)
        data = rng.uniform(size=96)
        result = convolve(data, self.boxcar)
        self.assertTrue(result.shape == data.shape)

        data2 = np.zeros(100)
        data2[50] = 1
        result2 = convolve(data2, self.boxcar)
        self.assertTrue(result2.sum() - 7 < 1e-6)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ConvolveTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
