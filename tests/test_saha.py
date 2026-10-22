"""Saha equilibrium: limiting cases and conservation."""

from __future__ import annotations

import numpy as np
import pytest

from cmb_decoupling.saha import (
    saha_redshift,
    saha_rhs_hydrogen,
    saha_species,
    x_e_saha_full,
    x_e_saha_hydrogen,
)


def test_hydrogen_limits(cosmo):
    assert x_e_saha_hydrogen(5000.0, cosmo) == pytest.approx(1.0, abs=1e-12)
    assert x_e_saha_hydrogen(500.0, cosmo) < 1e-10


def test_hydrogen_satisfies_saha_equation(cosmo):
    z = np.array([1200.0, 1350.0, 1500.0, 1700.0])
    x = x_e_saha_hydrogen(z, cosmo)
    np.testing.assert_allclose(x**2 / (1 - x), saha_rhs_hydrogen(z, cosmo), rtol=1e-10)


def test_hydrogen_no_cancellation_at_high_z(cosmo):
    # The naive quadratic root returned exactly 1 here (archive/Saha_py.py).
    assert 0.0 < 1.0 - x_e_saha_hydrogen(3000.0, cosmo) < 1e-8


def test_hydrogen_monotonic(cosmo):
    x = x_e_saha_hydrogen(np.linspace(800, 2000, 500), cosmo)
    assert np.all(np.diff(x) >= 0)


def test_saha_redshift_value(cosmo):
    # Regression: Saha x_e = 0.5 at z ~ 1369 with n_H (not n_b).
    assert saha_redshift(cosmo) == pytest.approx(1369.0, abs=1.0)


def test_full_saha_fully_ionised(cosmo):
    assert x_e_saha_full(1.0e4, cosmo) == pytest.approx(1 + 2 * cosmo.f_He, rel=1e-6)


@pytest.mark.parametrize("z", [1500.0, 2500.0, 5000.0, 8000.0])
def test_charge_conservation(cosmo, z):
    s = saha_species(z, cosmo)
    total = s["x_HII"] + cosmo.f_He * (s["y_HeII"] + 2 * s["y_HeIII"])
    assert s["x_e"] == pytest.approx(total, rel=1e-10)
    assert 0 <= s["y_HeII"] + s["y_HeIII"] <= 1


def test_full_reduces_to_hydrogen_when_helium_neutral(cosmo):
    z = np.array([200.0, 1000.0, 1400.0])
    np.testing.assert_allclose(x_e_saha_full(z, cosmo), x_e_saha_hydrogen(z, cosmo), rtol=1e-8)
