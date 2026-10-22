"""Peebles three-level atom: physical limits and numerical convergence."""

from __future__ import annotations

import numpy as np
import pytest

from cmb_decoupling import solve_peebles
from cmb_decoupling.peebles import alpha_B, peebles_C, peebles_rhs
from cmb_decoupling.saha import x_e_saha_hydrogen


def test_alpha_B_magnitude():
    # Case-B coefficient at 10^4 K is ~2.6e-19 m^3/s (Osterbrock & Ferland 2006).
    assert alpha_B(1.0e4) == pytest.approx(2.58e-19, rel=0.03)


def test_detailed_balance_at_saha(cosmo):
    # The net rate vanishes when x_e is the Saha value (detailed balance).
    z = 1800.0
    x = float(x_e_saha_hydrogen(z, cosmo))
    rate = peebles_rhs(z, np.array([x]), cosmo)[0]
    scale = peebles_rhs(z, np.array([0.99 * x]), cosmo)[0]
    assert abs(rate) < 1e-8 * abs(scale)


@pytest.mark.parametrize("z", [600.0, 1000.0, 1400.0])
def test_C_in_unit_interval(cosmo, peebles, z):
    C = peebles_C(z, float(peebles.x_e(z)), cosmo)
    assert 0.0 < C <= 1.0


def test_C_tends_to_one_late(cosmo, peebles):
    # Once photo-ionisation from n=2 is negligible, every n=2 atom reaches 1s.
    assert peebles_C(200.0, float(peebles.x_e(200.0)), cosmo) == pytest.approx(1.0, abs=1e-6)


def test_recombination_delayed_relative_to_saha(cosmo, peebles):
    z = np.linspace(400, 1500, 200)
    assert np.all(peebles.x_e(z) >= x_e_saha_hydrogen(z, cosmo))


def test_monotonic_decrease(peebles):
    z = np.linspace(0, 1799, 2000)
    assert np.all(np.diff(peebles.x_e(z)) >= -1e-15)


def test_continuous_at_switch(peebles):
    zs = peebles.z_switch
    assert peebles.x_e(zs - 1e-6) == pytest.approx(peebles.x_e(zs), rel=1e-4)


def test_insensitive_to_switch_redshift(cosmo, peebles):
    other = solve_peebles(cosmo, z_switch=1700.0)
    z = np.array([200.0, 800.0, 1100.0, 1400.0])
    np.testing.assert_allclose(other.x_e(z), peebles.x_e(z), rtol=1e-4)


def test_tolerance_convergence(cosmo, peebles):
    loose = solve_peebles(cosmo, rtol=1e-6, atol=1e-10)
    z = np.array([200.0, 800.0, 1100.0, 1400.0])
    np.testing.assert_allclose(loose.x_e(z), peebles.x_e(z), rtol=1e-4)


def test_fudge_speeds_recombination(cosmo, peebles):
    fast = solve_peebles(cosmo, fudge=1.14)
    assert fast.x_e(800.0) < peebles.x_e(800.0)
