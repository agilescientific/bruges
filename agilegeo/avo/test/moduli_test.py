import unittest
import agilegeo.avo as avo
import numpy as np

vp = 2350.
vs = 1125.
rho = 2500
lam = 7478125000.
mu = 3164062500.
youngs = 8551470048.4510355
pr = 0.3513434150638673
pmod = 13806250000.
bulk = 9587500000.

class ModuliTest( unittest.TestCase ):
    """
    Tests moduli calculation against spreadsheet
    https://dl.dropboxusercontent.com/u/14965965/Elastic_moduli_formula_checking.xlsx
    """
    
    def test_vpvs(self):
        self.assertAlmostEqual( avo.youngs(vp=vp, vs=vs, rho=rho), youngs, places=-2 )
        self.assertAlmostEqual(     avo.mu(vp=vp, vs=vs, rho=rho), mu, places=-2 )
        self.assertAlmostEqual(     avo.pr(vp=vp, vs=vs, rho=rho), pr )
        self.assertAlmostEqual(    avo.lam(vp=vp, vs=vs, rho=rho), lam, places=-2 )
        self.assertAlmostEqual(   avo.bulk(vp=vp, vs=vs, rho=rho), bulk, places=-2 )
        self.assertAlmostEqual(   avo.pmod(vp=vp, vs=vs, rho=rho), pmod, places=-2 )

    def test_lammu(self):
        self.assertAlmostEqual(     avo.vp(lam=lam, mu=mu, rho=rho), vp, places=2 )
        self.assertAlmostEqual( avo.youngs(lam=lam, mu=mu), youngs, places=-2 )
        self.assertAlmostEqual(     avo.pr(lam=lam, mu=mu), pr )
        self.assertAlmostEqual(   avo.bulk(lam=lam, mu=mu), bulk, places=-2 )
        self.assertAlmostEqual(   avo.pmod(lam=lam, mu=mu), pmod, places=-2 )
        
    def test_youngspr(self):
        self.assertAlmostEqual(    avo.vp(youngs=youngs, pr=pr, rho=rho), vp, places=2 )
        self.assertAlmostEqual(    avo.mu(youngs=youngs, pr=pr), mu, places=-2 )
        self.assertAlmostEqual(   avo.lam(youngs=youngs, pr=pr), lam, places=-2 )
        self.assertAlmostEqual(  avo.bulk(youngs=youngs, pr=pr), bulk, places=-2 )
        self.assertAlmostEqual(  avo.pmod(youngs=youngs, pr=pr), pmod, places=-2 )

    def test_moduli(self):
        mod = {'imp': vp * rho}
        mod['mu'] = mu
        mod['pr'] = pr
        mod['lam'] = lam
        mod['bulk'] = bulk
        mod['pmod'] = pmod
        mod['youngs'] = youngs

        self.assertDictEqual( avo.moduli(vp=vp, vs=vs, rho=rho), mod )
        
if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(ModuliTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
