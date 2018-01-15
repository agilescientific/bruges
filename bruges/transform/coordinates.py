# -*- coding: utf 8 -*-
"""
Seismic wavelets.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
from ..util.transformations import affine_matrix_from_points


class CoordTransform(object):
    def __init__(self, ix, xy):
        ix = np.array(ix)
        xy = np.array(xy)
        if ix.shape[1] == 2:
            ix = ix[:3].T
        if xy.shape[1] == 2:
            xy = xy[:3].T
        self.A = affine_matrix_from_points(ix, xy)

    def __call__(self, p):
        return (self.A @ np.append(p, [1]))[:2]

    def forward(self, p):
        return self(p)

    def reverse(self, q):
        p = (np.linalg.pinv(self.A) @ np.append(q, [1]))[:2]
        return np.rint(p).astype(np.int)
