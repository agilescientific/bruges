import unittest
import agilegeo.avo as avo
import numpy as np

class ModuliTest( unittest.TestCase ):
    """
    Tests moduli calculation against spreadsheet
    https://dl.dropboxusercontent.com/u/14965965/Elastic_moduli_formula_checking.xlsx
    """
    
    def test_vpvs(self):

        vp = 2350.
        vs = 1125.
        rho = 2500

        self.assertAlmostEqual( avo.youngs(vp=vp, vs=vs, rho=rho), 8551470048., places=-2 )
        self.assertAlmostEqual(     avo.mu(vp=vp, vs=vs, rho=rho), 3164062500., places=-2 )
        self.assertAlmostEqual(     avo.pr(vp=vp, vs=vs, rho=rho),0.35134341506387 )
        self.assertAlmostEqual(    avo.lam(vp=vp, vs=vs, rho=rho), 7478125000., places=-2 )
        self.assertAlmostEqual(   avo.bulk(vp=vp, vs=vs, rho=rho), 9587500000., places=-2 )
        self.assertAlmostEqual(   avo.pmod(vp=vp, vs=vs, rho=rho),13806250000., places=-2 )
        

    def test_lammu(self):

        lam = 7478125000.
        mu = 3164062500.
        rho = 2500

        self.assertAlmostEqual(     avo.vp(lam=lam, mu=mu, rho=rho), 2350.0, places=2 )
        self.assertAlmostEqual( avo.youngs(lam=lam, mu=mu), 8551470048., places=-2 )
        self.assertAlmostEqual(     avo.pr(lam=lam, mu=mu),0.35134341506387 )
        self.assertAlmostEqual(   avo.bulk(lam=lam, mu=mu), 9587500000., places=-2 )
        self.assertAlmostEqual(   avo.pmod(lam=lam, mu=mu),13806250000., places=-2 )
        
    def test_youngspr(self):

        youngs = 8551470048.
        pr = 0.35134341506387
        rho = 2500

        self.assertAlmostEqual(    avo.vp(youngs=youngs, pr=pr, rho=rho), 2350.0, places=2 )
        self.assertAlmostEqual(    avo.mu(youngs=youngs, pr=pr), 3164062500., places=-2 )
        self.assertAlmostEqual(   avo.lam(youngs=youngs, pr=pr), 7478125000., places=-2 )
        self.assertAlmostEqual(  avo.bulk(youngs=youngs, pr=pr), 9587500000., places=-2 )
        self.assertAlmostEqual(  avo.pmod(youngs=youngs, pr=pr),13806250000., places=-2 )

    def test_moduli(self):

        vp = 2350.
        vs = 1125.
        rho = 2500

        mod = {'imp': vp * rho}
        mod['mu'] = 3164062500.
        mod['pr'] = 0.3513434150638673
        mod['lam'] = 7478125000.
        mod['bulk'] = 9587500000.
        mod['pmod'] = 13806250000.
        mod['youngs'] = 8551470048.4510355

        self.assertDictEqual( avo.moduli(vp=vp, vs=vs, rho=rho), mod )
        
        
if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(ModuliTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
