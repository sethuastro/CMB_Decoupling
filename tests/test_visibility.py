"""Optical depth and visibility function: numerics and analytic checks."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad, trapezoid

from cmb_decoupling.constants import C_LIGHT, SIGMA_T
from cmb_decoupling.visibility import last_scattering


def fully_ionised(z):
    return np.ones_like(np.asarray(z, dtype=float))


def test_tau_matches_quadrature(cosmo):
    ls = last_scattering(fully_ionised, cosmo, z_max=50.0, n=5001)

    def integrand(z):
        return cosmo.n_H(z) * SIGMA_T * C_LIGHT / ((1 + z) * cosmo.hubble(z))

    exact, _ = quad(integrand, 0.0, 50.0, epsabs=0, epsrel=1e-12)
    assert ls.tau[-1] == pytest.approx(exact, rel=1e-7)


def test_visibility_normalisation(cosmo, peebles):
    ls = last_scattering(peebles.x_e, cosmo)
    assert trapezoid(ls.g, ls.z) == pytest.approx(1 - np.exp(-ls.tau[-1]), rel=1e-6)


def test_grid_convergence(cosmo, peebles):
    coarse = last_scattering(peebles.x_e, cosmo, n=15001)
    fine = last_scattering(peebles.x_e, cosmo, n=60001)
    assert coarse.z_peak() == pytest.approx(fine.z_peak(), abs=0.05)
    assert coarse.z_tau_one() == pytest.approx(fine.z_tau_one(), abs=0.05)
    assert coarse.fwhm()[0] == pytest.approx(fine.fwhm()[0], abs=0.1)


def test_reproduces_camb_zstar_from_camb_xe(cosmo, camb_ref):
    # Feeding CAMB's own x_e through our tau integrator must recover CAMB's z_*.
    ls = last_scattering(camb_ref.x_e, cosmo)
    assert ls.z_tau_one() == pytest.approx(camb_ref.z_star, abs=0.1)


def test_conformal_peak_is_later_in_time(cosmo, peebles):
    ls = last_scattering(peebles.x_e, cosmo)
    assert ls.z_peak_conformal() > ls.z_peak()
