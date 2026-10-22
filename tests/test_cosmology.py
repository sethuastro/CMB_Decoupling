"""Background cosmology and constants."""

from __future__ import annotations

import numpy as np
import pytest

from cmb_decoupling.constants import E_ION_H, EV, LAMBDA_LYA


def test_hydrogen_ionisation_energy():
    # NIST: 13.598434 eV for 1H (reduced-mass corrected Rydberg).
    assert E_ION_H / EV == pytest.approx(13.5984, abs=1e-3)


def test_lyman_alpha_wavelength():
    # Lyman-alpha ~ 121.6 nm.
    assert LAMBDA_LYA * 1e9 == pytest.approx(121.57, abs=0.05)


def test_flatness(cosmo):
    assert cosmo.Omega_m + cosmo.Omega_r + cosmo.Omega_lambda == pytest.approx(1.0, abs=1e-14)


def test_radiation_density(cosmo):
    # Omega_gamma h^2 = 2.47e-5 for T_cmb = 2.7255 K (Dodelson & Schmidt 2020).
    assert cosmo.Omega_gamma * cosmo.h**2 == pytest.approx(2.47e-5, rel=5e-3)


def test_hubble_today(cosmo):
    assert cosmo.hubble(0.0) == pytest.approx(cosmo.H0, rel=1e-14)


def test_hydrogen_density_includes_helium(cosmo):
    # n_H0 is (1 - Y_p) of the all-protons estimate (up to m_H vs m_p).
    n_all = cosmo.Omega_b * cosmo.rho_crit0 / 1.67262192e-27
    assert cosmo.n_H0 / n_all == pytest.approx(1 - cosmo.Y_p, rel=1e-3)


def test_hubble_matches_camb(cosmo, camb_ref):
    z = np.array([0.0, 10.0, 500.0, 1100.0, 3000.0])
    np.testing.assert_allclose(cosmo.hubble(z), camb_ref.hubble(z), rtol=1e-6)
