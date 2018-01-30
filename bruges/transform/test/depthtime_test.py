# -*- coding: utf-8 -*-
"""
Tests.
"""
import unittest
import numpy as np

from bruges.transform import time_to_depth, depth_to_time


class TimeDepthTest(unittest.TestCase):
    """
    Tests the basic functionality of the time to depth and depth to
    time conversions in bruges.
    """

    def test_depthtotime(self):
        data = np.zeros((100, 100))
        vmodel = np.zeros((100, 100))

        data[0:50, :] += 100.0
        data[50:, :] += 200.0

        vmodel[0:50, :] += 1500.0
        vmodel[50:, :] += 3000.0

        dt = 0.001
        dz = 1.0

        face_change = int(np.floor(((49 * dz) / 1500.0) / dt))

        output = depth_to_time(data, vmodel, dz, dt, twt=False)

        self.assertTrue((output[face_change+1, 50] -
                         output[face_change, 50]) == 100)

    def test_timetodepth(self):
        data = np.zeros((100, 100))
        vmodel = np.zeros((100, 100))

        data[0:50, :] += 100.0
        data[50:, :] += 200.0

        vmodel[0:50, :] += 1500.0
        vmodel[50:, :] += 5000.0

        dt = 0.001
        dz = 1.0

        output = time_to_depth(data, vmodel, dt, dz, twt=False)

        face_change = int(np.floor(((49.5 * dt) * 1500.0) / dz))

        self.assertTrue((output[face_change+1, 50] -
                         output[face_change, 50]) == 100)

    def test_recip(self):
        data = np.zeros((100, 100))
        vmodel = np.zeros((100, 100))

        data[0:50, :] += 100.0
        data[50:, :] += 200.0

        vmodel[0:50, :] += 1500.0
        vmodel[50:, :] += 1800.0

        dt = 0.001
        dz = 1.

        out1 = depth_to_time(data, vmodel, dz, dt)
        v2 = depth_to_time(vmodel, vmodel, dz, dt)
        out2 = time_to_depth(out1, v2, dt, dz)

        # This test is not finished!


if (__name__ == '__main__'):

    suite = \
      unittest.TestLoader().loadTestsFromTestCase(TimeDepthTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
