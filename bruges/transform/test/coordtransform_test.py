# -*- coding: utf-8 -*-
"""
Tests.
"""
import unittest
import numpy as np

from bruges.transform import CoordTransform


# UTM coords of 3 unique inline, crossline locations.
corners_xy = np.array([[600938.125,  6073394.5],
                       [631226.3125, 6074241.0],
                       [630720.25,   6092358.5]])

# The inline, crossline locations you just provided.
corners_ix = np.array([[99,   104],
                       [99,  1316],
                       [824, 1316]])


class CoordTransformTest(unittest.TestCase):
    """
    Tests the basic functionality of coordinate transformation.
    """

    def test_coordtransform(self):
        transform = CoordTransform(corners_ix, corners_xy)
        ix = [440, 763]
        xy = np.array([617167, 6082379])
        result = transform(ix)
        self.assertAlmostEqual(result[0], xy[0], places=-1)
        self.assertAlmostEqual(result[1], xy[1], places=-1)

        result = transform.forward(ix)
        self.assertAlmostEqual(result[0], xy[0], places=-1)
        self.assertAlmostEqual(result[1], xy[1], places=-1)

        result = transform.reverse(xy)
        self.assertAlmostEqual(result[0], ix[0], places=0)
        self.assertAlmostEqual(result[1], ix[1], places=0)


if (__name__ == '__main__'):

    suite = \
      unittest.TestLoader().loadTestsFromTestCase(CoordTransformTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
