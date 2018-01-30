# -*- coding: utf-8 -*-
import unittest
import numpy as np

from bruges.attribute import dipsteer
from bruges.filters import ricker


class DipTest(unittest.TestCase):

    def test_spikes(self):
        """
        Make a simple case with spikes offset on each trace.
        """

        test_data = np.zeros((100, 100))
        index = np.arange(0, 100)

        test_data[index, index] = 1.0

        maxlag = 5
        dip = dipsteer(test_data, 1, 1, maxlag)

        # Output should be all -1 except for the edge effects
        same = np.allclose(dip[maxlag:-maxlag, maxlag:-maxlag],
                           -test_data[maxlag:-maxlag, maxlag:-maxlag],
                           .01)

        self.assertTrue(same)

        # Answer should not change if we increase the stepout
        stepout = 4
        maxlag = 1
        dip = dipsteer(test_data, 1, stepout, maxlag)

        diff = (dip[maxlag*stepout:100-(maxlag*stepout),
                stepout * maxlag:100-(stepout*maxlag)] +
                test_data[maxlag * stepout:100-(maxlag * stepout),
                          maxlag*stepout:100-(stepout*maxlag)])

        diff = diff[1:-1, 1:-1]

        # Output should be all -1 except for the edge effects
        same = np.allclose(diff, np.zeros(diff.shape), .01)

        self.assertTrue(same)

    def test_wavelet(self):

        # Make a dipping reflectivity field
        data = np.zeros((150, 100))

        data[0, :] = np.arange(100) / 10.0

        dips = np.arange(0, 100, 2)
        dips = np.concatenate([dips, dips[::-1]])

        window_size = 10
        freq = 20.0
        duration = 1
        wavelet = ricker(window_size, duration, freq)

        for i in range(dips.size):

            data[:, i] = np.convolve(wavelet,
                                     np.roll(data[:, i], dips[i]),
                                     mode='same')

        stepout = 1
        maxlag = 2
        dip = dipsteer(data, window_size, stepout=stepout,
                       maxlag=maxlag, overlap=0.5)

        undipped = np.zeros(dips.shape)
        for i in range(dips.size):
            undipped[i] = (dip[dips[i], i])

        check = undipped[maxlag*stepout:-maxlag*stepout]
        test = np.zeros(check.shape)
        test[0: test.size//2 - 1] = -2.0
        test[(test.size//2) - 1] = -1
        test[(test.size//2)] = 1
        test[(test.size//2)+1:] = 2

        self.assertTrue(np.allclose(test, check))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(DipTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
