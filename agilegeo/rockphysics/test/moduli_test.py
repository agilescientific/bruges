#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Tests.
"""
import unittest

from agilegeo.rockphysics import moduli

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
        self.assertAlmostEqual(moduli.youngs(vp=vp, vs=vs, rho=rho), youngs, places=-2)
        self.assertAlmostEqual(    moduli.mu(vp=vp, vs=vs, rho=rho), mu, places=-2)
        self.assertAlmostEqual(    moduli.pr(vp=vp, vs=vs, rho=rho), pr)
        self.assertAlmostEqual(   moduli.lam(vp=vp, vs=vs, rho=rho), lam, places=-2)
        self.assertAlmostEqual(  moduli.bulk(vp=vp, vs=vs, rho=rho), bulk, places=-2)
        self.assertAlmostEqual(  moduli.pmod(vp=vp, vs=vs, rho=rho), pmod, places=-2)

    def test_lammu(self):
        self.assertAlmostEqual(    moduli.vp(lam=lam, mu=mu, rho=rho), vp, places=2)
        self.assertAlmostEqual(moduli.youngs(lam=lam, mu=mu), youngs, places=-2)
        self.assertAlmostEqual(    moduli.pr(lam=lam, mu=mu), pr)
        self.assertAlmostEqual(  moduli.bulk(lam=lam, mu=mu), bulk, places=-2)
        self.assertAlmostEqual(  moduli.pmod(lam=lam, mu=mu), pmod, places=-2)

    def test_youngspr(self):
        self.assertAlmostEqual(  moduli.vp(youngs=youngs, pr=pr, rho=rho), vp, places=2)
        self.assertAlmostEqual(  moduli.mu(youngs=youngs, pr=pr), mu, places=-2)
        self.assertAlmostEqual( moduli.lam(youngs=youngs, pr=pr), lam, places=-2)
        self.assertAlmostEqual(moduli.bulk(youngs=youngs, pr=pr), bulk, places=-2)
        self.assertAlmostEqual(moduli.pmod(youngs=youngs, pr=pr), pmod, places=-2)

    def test_moduli(self):
        mod = {'imp': vp * rho}
        mod['mu'] = mu
        mod['pr'] = pr
        mod['lam'] = lam
        mod['bulk'] = bulk
        mod['pmod'] = pmod
        mod['youngs'] = youngs

        self.assertDictEqual(moduli.moduli(vp=vp, vs=vs, rho=rho), mod)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ModuliTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
