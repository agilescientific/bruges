# -*- coding: utf-8 -*-
"""
Tests.
"""
import unittest

from bruges.rockphysics import bounds as b

# Inputs.
f0 = [0.5, 0.5]
f1 = [50, 50]
f2 = [0.6, 0.4]     # suspension of quartz in water
k1 = [36.0, 4.0]    # Kquartz, Kclay [GPa]
k2 = [36.0, 2.2]    # Kquartz, Kwater [GPa]

# Hashin-Shtrikman tests
fh1 = [0.0, 1.0]     # end-member test case for HS bounds
fh2 = [1.0, 1.0]     # end-member test case for HS bounds
fh3 = [0.4, 0.6]    # middle-ish test case for HS bounds
khs = [21.0, 76.8]  # Bulk moduli [GPa] for HS tests, Clay, Calcite
ghs = [7.0, 32.0]   # Shear moduli [GPa] for HS tests, Clay, Calcite

# Expected outputs.
mv = 20.0       # Voigt average (hand-calculated)
mr = 7.199      # Reuss average (hand-calculated)
mr2 = 5.04      # Example from Mavko et al. RPH
mh = 13.599     # Voigt-Reuss Hill average (hand-calculated)
mhs1 = 76.8     # HS upper and lower bound at end-member 1
mhs2 = 21.0     # HS upper and lower bound at end-member 2
mhs3 = 35.5329  # HS lower bound at 50-50
mhs4 = 40.3989  # HS upper bound at 50-50


class BoundsTest(unittest.TestCase):
    """
    Tests calculations of elastic mixture bounds
    """

    def test_voigt(self):
        self.assertAlmostEqual(b.voigt_bound(f0, k1), mv, places=-2)
        self.assertAlmostEqual(b.voigt_bound(f1, k1), mv, places=-2)

    def test_reuss(self):
        self.assertAlmostEqual(b.reuss_bound(f0, k1), mr, places=-2)
        self.assertAlmostEqual(b.reuss_bound(f2, k2), mr2, places=-2)
        self.assertAlmostEqual(b.reuss_bound([0.6, 0.2, 0.2],
                                             [36.0, 2.2, 0.000131]),
                               0.00065, places=-2)

    def test_hill_average(self):
        self.assertAlmostEqual(b.hill_average(f0, k1), mh, places=-2)

    def test_hashin_shtrikman(self):
        self.assertAlmostEqual(
            b.hashin_shtrikman(fh1, khs, ghs, modulus='bulk')[1],
            mhs1, places=-2)
        self.assertAlmostEqual(
            b.hashin_shtrikman(fh1, khs, ghs, modulus='bulk')[0],
            mhs1, places=-2)
        self.assertAlmostEqual(
            b.hashin_shtrikman(fh2, khs, ghs, modulus='bulk')[1],
            mhs2, places=-2)
        self.assertAlmostEqual(
            b.hashin_shtrikman(fh2, khs, ghs, modulus='bulk')[0],
            mhs2, places=-2)
        self.assertAlmostEqual(
            b.hashin_shtrikman(fh3, khs, ghs, modulus='bulk')[1],
            mhs3, places=-2)
        self.assertAlmostEqual(
            b.hashin_shtrikman(fh3, khs, ghs, modulus='bulk')[0],
            mhs3, places=-2)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(BoundsTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
