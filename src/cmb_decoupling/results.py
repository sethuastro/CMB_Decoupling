"""Key derived quantities of an ionisation history, computed in one place."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import NDArray

from .cosmology import Cosmology
from .visibility import XeFunction, decoupling_redshift, last_scattering

#: Redshift at which the "freeze-out" residual ionisation is quoted. Chosen
#: after recombination has frozen out but before the T_m = T_gamma
#: approximation of this package becomes poor (see docs/theory.md).
Z_FREEZE_OUT: float = 200.0

#: Redshift range over which x_e is compared with the reference (scripts/validate.py).
Z_VALIDATION_MIN: float = 200.0
Z_VALIDATION_MAX: float = 2000.0
#: Number of points of the uniform comparison grid (spacing 0.1 in z).
N_VALIDATION: int = 18001


@dataclass(frozen=True)
class KeyResults:
    """Summary numbers for one ionisation history.

    Attributes
    ----------
    z_decoupling : float
        Redshift at which Gamma_T = H.
    z_tau_one : float
        Redshift at which the optical depth equals 1 (Planck's z_*).
    z_peak : float
        Peak of g(z) per unit redshift.
    z_peak_conformal : float
        Peak of the visibility per unit conformal time.
    fwhm : float
        Full width at half maximum of g(z), in redshift.
    x_e_freeze_out : float
        x_e at ``Z_FREEZE_OUT``.
    x_e_today : float
        x_e at z = 0 with no reionisation.
    """

    z_decoupling: float
    z_tau_one: float
    z_peak: float
    z_peak_conformal: float
    fwhm: float
    x_e_freeze_out: float
    x_e_today: float

    def as_dict(self) -> dict[str, float]:
        """Return the results as a plain dictionary."""
        return asdict(self)


def key_results(x_e: XeFunction, cosmo: Cosmology) -> KeyResults:
    """Compute :class:`KeyResults` for an ionisation history.

    Parameters
    ----------
    x_e : callable
        Ionisation history x_e(z) = n_e / n_H.
    cosmo : Cosmology
        Background cosmology.
    """
    ls = last_scattering(x_e, cosmo)
    fwhm, _, _ = ls.fwhm()
    return KeyResults(
        z_decoupling=decoupling_redshift(x_e, cosmo),
        z_tau_one=ls.z_tau_one(),
        z_peak=ls.z_peak(),
        z_peak_conformal=ls.z_peak_conformal(),
        fwhm=fwhm,
        x_e_freeze_out=float(x_e(Z_FREEZE_OUT)),
        x_e_today=float(x_e(0.0)),
    )


def fractional_residual(x_e: XeFunction, x_e_ref: XeFunction, z: NDArray[np.float64]) -> NDArray[np.float64]:
    """Fractional residual Delta x_e / x_e = x_e(z) / x_e_ref(z) - 1 [dimensionless].

    Parameters
    ----------
    x_e, x_e_ref : callable
        Ionisation history under test and reference ionisation history, x_e(z) = n_e / n_H.
    z : ndarray
        Redshifts at which to compare.
    """
    return x_e(z) / x_e_ref(z) - 1.0


@dataclass(frozen=True)
class ResidualSummary:
    """Where and by how much an ionisation history departs from a reference.

    All residuals are fractional, Delta x_e / x_e = x_e / x_e_ref - 1 [dimensionless].

    Attributes
    ----------
    window_lo, window_hi : float
        Last-scattering window: the half-maximum redshifts of the *reference* g(z).
    max_abs_window : float
        Largest |Delta x_e / x_e| inside the last-scattering window.
    z_max_abs_window : float
        Redshift where ``max_abs_window`` occurs.
    min_residual, z_min_residual : float
        Most negative residual over the full comparison range, and its redshift.
    max_residual, z_max_residual : float
        Most positive residual over the full comparison range, and its redshift.
    z_min, z_max : float
        Comparison range.
    """

    window_lo: float
    window_hi: float
    max_abs_window: float
    z_max_abs_window: float
    min_residual: float
    z_min_residual: float
    max_residual: float
    z_max_residual: float
    z_min: float
    z_max: float

    def as_dict(self) -> dict[str, float]:
        """Return the summary as a plain dictionary."""
        return asdict(self)


def residual_summary(
    x_e: XeFunction,
    x_e_ref: XeFunction,
    cosmo: Cosmology,
    z_min: float = Z_VALIDATION_MIN,
    z_max: float = Z_VALIDATION_MAX,
    n: int = N_VALIDATION,
) -> ResidualSummary:
    """Summarise the fractional x_e residual against a reference.

    Extrema are taken on a uniform grid of ``n`` points over [``z_min``, ``z_max``]
    (no interpolation between grid points).

    Parameters
    ----------
    x_e, x_e_ref : callable
        Ionisation history under test and reference, x_e(z) = n_e / n_H.
    cosmo : Cosmology
        Background cosmology (used to compute the reference visibility function).
    z_min, z_max : float, optional
        Comparison range in redshift.
    n : int, optional
        Number of grid points.
    """
    z = np.linspace(z_min, z_max, n)
    r = fractional_residual(x_e, x_e_ref, z)
    _, w_lo, w_hi = last_scattering(x_e_ref, cosmo).fwhm()
    in_window = (z >= w_lo) & (z <= w_hi)
    i_w = int(np.argmax(np.abs(r[in_window])))
    i_min, i_max = int(np.argmin(r)), int(np.argmax(r))
    return ResidualSummary(
        window_lo=float(w_lo),
        window_hi=float(w_hi),
        max_abs_window=float(np.abs(r[in_window][i_w])),
        z_max_abs_window=float(z[in_window][i_w]),
        min_residual=float(r[i_min]),
        z_min_residual=float(z[i_min]),
        max_residual=float(r[i_max]),
        z_max_residual=float(z[i_max]),
        z_min=z_min,
        z_max=z_max,
    )
