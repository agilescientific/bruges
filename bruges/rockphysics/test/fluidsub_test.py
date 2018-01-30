# -*- coding: utf-8 -*-
"""
Tests.
"""
import unittest

from bruges.rockphysics import fluidsub

# Inputs... GAS case
vp_gas = 2429.0
vs_gas = 1462.4
rho_gas = 2080.

# Expected outputs... BRINE case
vp_brine = 2850.5
vs_brine = 1416.1
rho_brine = 2210.0

phi = 0.275           # Don't know this... reading from fig
rhohc = 250.0         # gas
rhow = 1040.0         # brine
sw = 0.3              # Don't know this... just guessing
swnew = 1.0           # Don't know this... just guessing
khc = 207000000.0     # gas
kw = 2950000000.0     # brine
kclay = 25000000000.0
kqtz = 37000000000.0
vclay = 0.05
kmin = 36266406250.0  # Don't know this... reading from fig


class FluidsubTest(unittest.TestCase):
    """
    Tests fluid sub calculations against Smith et al 2003.
    https://dl.dropboxusercontent.com/u/14965965/Smith_etal_2003.pdf
    """

    def test_avseth(self):
        # Base case: gas
        # Subbing with: brine

        sub = fluidsub.avseth_fluidsub(vp=vp_gas,
                                       vs=vs_gas,
                                       rho=rho_gas,
                                       phi=phi,
                                       rhof1=rhohc,
                                       rhof2=rhow,
                                       kmin=37000000000,
                                       kf1=khc,
                                       kf2=kw)

        self.assertAlmostEqual(sub[0], vp_brine, places=-1)  # Cannot match :(
        self.assertAlmostEqual(sub[1], vs_brine, places=-1)  # Cannot match :(
        self.assertAlmostEqual(sub[2], rho_brine, places=-1)  # Cannot match :(

    def test_smith(self):
        # Base case: gas
        # Subbing with: brine

        sub = fluidsub.smith_fluidsub(vp=vp_gas,
                                      vs=vs_gas,
                                      rho=rho_gas,
                                      phi=phi,
                                      rhohc=rhohc,
                                      rhow=rhow,
                                      sw=sw,
                                      swnew=swnew,
                                      khc=khc,
                                      kw=kw,
                                      kclay=kclay,
                                      kqtz=kqtz,
                                      vclay=vclay)

        self.assertAlmostEqual(sub[0], vp_brine,  places=-1)
        self.assertAlmostEqual(sub[1], vs_brine,  places=-1)
        self.assertAlmostEqual(sub[2], rho_brine, places=-1)  # Cannot match :(


if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(FluidsubTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
