# -*- coding: utf-8 -*-
"""
Coordinate transformation. This module contains a class for
converting between seismic survey inline-xline coordinates
and real-world UTM coordinates.

:copyright: 2018 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
from ..util.transformations import affine_matrix_from_points


class CoordTransform(object):
    """
    A class for converting between seismic survey inline-xline
    coordinates and real-world UTM coordinates.

    Instantiate with a pair of at least 3 coordinates mapping
    one space to the other. Provide two array-likes of shape
    (3, 2). See below for example.

    After instantiation, the class is callable in the 'forward'
    (inline-xline to UTMx-UTMy) direction. You can also use
    ``CoordTransform.forwrd()``. The ``CoordTransform.reverse()``
    method converts in the other direction (UTMx-UTMy to
    inline-xline).

    Example
    >>> corner_ix = [[0,  0], [0, 950], [650, 950]]
    >>> corner_xy = [[605835.5, 6073556.5],
                     [629576.3, 6074220.0],
                     [629122.5, 6090463.2]]
    >>> transform = bruges.transform.CoordTransform(corner_ix, corner_xy)
    >>> transform([300, 400])
    array([ 615622.18016194, 6081332.72995951])
    >>> transform.forward([300, 400])
    array([ 615622.18016194, 6081332.72995951])
    >>> transform.reverse([ 615622.18016194, 6081332.72995951])
    ​array([300, 400])
    """
    def __init__(self, ix, xy):
        ix = np.array(ix)
        xy = np.array(xy)
        if ix.shape[1] == 2:
            ix = ix[:3].T
        if xy.shape[1] == 2:
            xy = xy[:3].T
        self.A = affine_matrix_from_points(ix, xy)

    def __call__(self, p):
        p = np.asanyarray(p)
        return np.dot(self.A, np.append(p, [1]))[:2]

    def forward(self, p):
        """
        Convert inline-xline to UTMx-UTM-y.

        Example
        >>> transform.forward([300, 400])
        array([ 615622.18016194, 6081332.72995951])
        """
        p = np.asanyarray(p)
        return self(p)

    def reverse(self, q):
        """
        Convert UTMx-UTM-y to inline-xline.

        Example
        >>> transform.reverse([ 615622.18016194, 6081332.72995951])
        ​array([300, 400])
        """
        p = np.dot(np.linalg.pinv(self.A), np.append(q, [1]))[:2]
        return np.rint(p).astype(np.int)
