# -*- coding: utf-8 -*-
"""
Tests.
"""
import unittest

from bruges.rockphysics import moduli as m

vp = 2350.
vs = 1125.
rho = 2500
lam = 7478125000.
mu = 3164062500.
youngs = 8551470048.4510355
pr = 0.3513434150638673
pmod = 13806250000.
bulk = 9587500000.


class ModuliTest(unittest.TestCase):
    """
    Tests moduli calculation against spreadsheet
    https://dl.dropboxusercontent.com/u/14965965/Elastic_moduli_formula_checking.xlsx
    """

    def test_vpvs(self):
        self.assertAlmostEqual(m.youngs(vp=vp, vs=vs, rho=rho),
                               youngs, places=-2)
        self.assertAlmostEqual(m.mu(vp=vp, vs=vs, rho=rho), mu, places=-2)
        self.assertAlmostEqual(m.pr(vp=vp, vs=vs, rho=rho), pr)
        self.assertAlmostEqual(m.lam(vp=vp, vs=vs, rho=rho), lam, places=-2)
        self.assertAlmostEqual(m.bulk(vp=vp, vs=vs, rho=rho), bulk, places=-2)
        self.assertAlmostEqual(m.pmod(vp=vp, vs=vs, rho=rho), pmod, places=-2)

    def test_lammu(self):
        self.assertAlmostEqual(m.vp(lam=lam, mu=mu, rho=rho), vp, places=2)
        self.assertAlmostEqual(m.youngs(lam=lam, mu=mu), youngs, places=-2)
        self.assertAlmostEqual(m.pr(lam=lam, mu=mu), pr)
        self.assertAlmostEqual(m.bulk(lam=lam, mu=mu), bulk, places=-2)
        self.assertAlmostEqual(m.pmod(lam=lam, mu=mu), pmod, places=-2)

    def test_youngspr(self):
        self.assertAlmostEqual(m.vp(youngs=youngs, pr=pr, rho=rho),
                               vp, places=2)
        self.assertAlmostEqual(m.mu(youngs=youngs, pr=pr), mu, places=-2)
        self.assertAlmostEqual(m.lam(youngs=youngs, pr=pr), lam, places=-2)
        self.assertAlmostEqual(m.bulk(youngs=youngs, pr=pr), bulk, places=-2)
        self.assertAlmostEqual(m.pmod(youngs=youngs, pr=pr), pmod, places=-2)

    def test_youngslam(self):
        self.assertAlmostEqual(m.mu(youngs=youngs, lam=lam), mu, places=-2)
        self.assertAlmostEqual(m.pr(youngs=youngs, lam=lam), pr, places=-2)
        self.assertAlmostEqual(m.bulk(youngs=youngs, lam=lam), bulk, places=-2)
        self.assertAlmostEqual(m.pmod(youngs=youngs, lam=lam), pmod, places=-2)

    def test_moduli(self):
        mod = {'imp': vp * rho}
        mod['mu'] = mu
        mod['pr'] = pr
        mod['lam'] = lam
        mod['bulk'] = bulk
        mod['pmod'] = pmod
        mod['youngs'] = youngs

        self.assertDictEqual(m.moduli_dict(vp=vp, vs=vs, rho=rho), mod)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ModuliTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
