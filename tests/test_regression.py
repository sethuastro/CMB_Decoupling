"""Regression against CAMB/RECFAST and against this code's own validated outputs.

The tolerances below describe the *measured* agreement of the simplified
model; they are not tuned targets. If physics changes move these numbers,
update them deliberately and explain why in the commit message.
"""

from __future__ import annotations

import pytest

from cmb_decoupling.results import key_results


@pytest.fixture(scope="module")
def ours(cosmo, peebles):
    return key_results(peebles.x_e, cosmo)


@pytest.fixture(scope="module")
def ref(cosmo, camb_ref):
    return key_results(camb_ref.x_e, cosmo)


def test_last_scattering_against_camb(ours, ref):
    assert ours.z_tau_one == pytest.approx(ref.z_tau_one, abs=2.0)
    assert ours.z_peak == pytest.approx(ref.z_peak, abs=2.0)
    assert ours.fwhm == pytest.approx(ref.fwhm, abs=3.0)
    assert ours.z_decoupling == pytest.approx(ref.z_decoupling, abs=2.0)


def test_freeze_out_against_camb(ours, ref):
    # Pure Peebles overestimates the residual ionisation (~14% at z=200).
    assert ours.x_e_freeze_out == pytest.approx(ref.x_e_freeze_out, rel=0.25)


def test_self_regression(ours):
    # Values recorded on 2026-10-22 (results/key_results.json).
    assert ours.z_tau_one == pytest.approx(1089.44, abs=0.05)
    assert ours.z_peak == pytest.approx(1078.80, abs=0.05)
    assert ours.fwhm == pytest.approx(203.33, abs=0.1)
    assert ours.x_e_freeze_out == pytest.approx(3.849e-4, rel=1e-3)
