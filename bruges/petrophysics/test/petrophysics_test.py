"""Regression tests for Gardner parameter estimation."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from bruges.petrophysics.petrophysics import gardner_param


class GardnerParamTest(unittest.TestCase):

    def check_coefficients(self, alpha, beta, low, high):
        vp = np.linspace(low, high, 40)
        # Construct observations directly, independently of the implementation.
        rhob = alpha * vp ** beta
        fitted = gardner_param(vp, rhob)
        self.assertIsInstance(fitted, tuple)
        self.assertEqual(len(fitted), 2)
        assert_allclose(fitted, (alpha, beta), rtol=1e-6)
        predicted = (rhob / fitted[0]) ** (1 / fitted[1])
        assert_allclose(predicted, vp, rtol=1e-7)

    def test_standard_coefficients(self):
        self.check_coefficients(310, 0.25, 1500, 6000)

    def test_nondefault_higher_exponent(self):
        self.check_coefficients(200, 0.30, 1800, 4500)

    def test_nondefault_lower_exponent(self):
        self.check_coefficients(450, 0.20, 2500, 6500)

    def test_noisy_fit_minimizes_velocity_residuals(self):
        rhob = np.linspace(1900, 2800, 12)
        vp = (rhob / 310) ** 4 + np.array(
            [20, -30, 45, -15, 35, -50, 10, 40, -25, 30, -40, 15]
        )
        alpha, beta = gardner_param(vp, rhob)
        predicted = (rhob / alpha) ** (1 / beta)
        residual = predicted - vp
        # Independent analytic derivatives of (rho / alpha)**(1 / beta).
        # Velocity least squares requires residuals orthogonal to both columns
        # of the Jacobian; density or log-space fitting has a different optimum.
        derivatives = (
            -predicted / (alpha * beta),
            -predicted * np.log(rhob / alpha) / beta ** 2,
        )
        for derivative in derivatives:
            correlation = np.dot(residual, derivative) / (
                np.linalg.norm(residual) * np.linalg.norm(derivative)
            )
            self.assertLess(abs(correlation), 1e-5)


if __name__ == '__main__':
    unittest.main()
