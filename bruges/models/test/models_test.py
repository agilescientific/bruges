import unittest
import numpy as np
from numpy import array

from bruges.models import reconcile, interpolate, panel
from bruges.models import wedge


class ModelTest(unittest.TestCase):
    """
    Tests models.
    """
    def test_reconcile(self):
        a = np.array([2, 6, 7, 7, 3])
        b = np.array([3, 7, 3])
        A, B = reconcile(a, b, order=0)
        A_, B_ = array([2, 6, 7, 7, 3]), array([3, 7, 7, 3, 3])
        self.assertTrue(np.array_equal(A, A_))
        self.assertTrue(np.array_equal(B, B_))

    def test_interpolate(self):
        a = np.array([2, 6, 7, 7, 3])
        b = np.array([3, 7, 7, 3, 3])
        interp = interpolate(a, b, num=10)
        self.assertTrue(interp.shape == (5, 10))

    def test_panel(self):
        a = np.array([2, 6, 7, 7, 3])
        b = np.array([3, 7, 3])
        dists = (10,)
        out = panel(a, b, num=15, dists=dists)
        sample = out[:, 7]
        self.assertTrue(np.all(sample[:4] == array([2.5, 6.5, 5., 3.])))
        self.assertTrue(np.isnan(sample[-1]))

    def test_wedge(self):
        w, top, base, ref = wedge(depth=10, width=7, strat=(10, (20, 30), 40))
        col = array([10, 10, 10, 20, 20, 30, 40, 40, 40, 40])
        t = array([3., 3., 3., 3., 3., 3., 3.])
        b = array([3., 3., 3.6, 4.2, 4.8, 5.4, 6. ])
        self.assertTrue(np.all(w[:, -1] == col))
        self.assertTrue(w.sum() == 1990)
        self.assertTrue(np.allclose(top, t))
        self.assertTrue(np.allclose(base, b))
        self.assertTrue(ref == 6)

    def test_netgross(self):
        w, top, *_ = wedge(depth=10, width=7, breadth=3, strat=(10, (20, 30), 40))
        self.assertTrue(w.sum() == 6003)
        self.assertTrue(w.shape == (10, 7, 3))
        self.assertTrue(top.sum() == 63.0)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ModelTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
