# -*- coding: utf-8 -*-
import unittest
import numpy

from bruges.attribute import similarity


class SimilarityTest(unittest.TestCase):

    def test_same_data(self):
        """
        Simple test to check if the algorithm works for the trivial case.
        """
        # data2d = numpy.ones([3, 100])
        # output2d = similarity(data2d, duration=0.01, dt=0.001, kind='marfurt')
        # self.assertTrue(output2d.dims == 2)

        data3d = numpy.ones([3, 3, 100])
        output3d = similarity(data3d, duration=0.01, dt=0.001, kind='marfurt')
        self.assertTrue(output3d.ndim == 3)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(SimilarityTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
