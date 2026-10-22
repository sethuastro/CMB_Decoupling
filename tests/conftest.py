"""Shared fixtures: the expensive solutions are computed once per session."""

from __future__ import annotations

import pytest

from cmb_decoupling import PLANCK2018, solve_peebles


@pytest.fixture(scope="session")
def cosmo():
    """Planck 2018 background."""
    return PLANCK2018


@pytest.fixture(scope="session")
def peebles(cosmo):
    """Default (pure Peebles) ionisation history."""
    return solve_peebles(cosmo)


@pytest.fixture(scope="session")
def camb_ref(cosmo):
    """CAMB/RECFAST reference; skipped if camb is not installed."""
    pytest.importorskip("camb")
    from cmb_decoupling.reference import camb_reference

    return camb_reference(cosmo)
