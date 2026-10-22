"""Thomson scattering rate, optical depth and visibility function.

References
----------
Dodelson & Schmidt (2020), Modern Cosmology, 2nd ed.
Planck Collaboration (2020), A&A 641, A6 -- definition of z_* as tau(z_*) = 1.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from .constants import C_LIGHT, SIGMA_T
from .cosmology import Cosmology

XeFunction = Callable[[ArrayLike], NDArray[np.float64]]


def thomson_rate(z: ArrayLike, x_e: ArrayLike, cosmo: Cosmology) -> NDArray[np.float64]:
    """Thomson scattering rate per photon, Gamma_T = n_e sigma_T c [s^-1].

    Parameters
    ----------
    z : array_like
        Redshift.
    x_e : array_like
        Free-electron fraction n_e / n_H at ``z``.
    cosmo : Cosmology
        Background cosmology.
    """
    return np.asarray(x_e) * cosmo.n_H(z) * SIGMA_T * C_LIGHT


def decoupling_redshift(
    x_e: XeFunction, cosmo: Cosmology, z_lo: float = 200.0, z_hi: float = 1800.0
) -> float:
    """Redshift at which Gamma_T(z) = H(z).

    Parameters
    ----------
    x_e : callable
        Ionisation history x_e(z).
    cosmo : Cosmology
        Background cosmology.
    z_lo, z_hi : float, optional
        Bracketing redshifts.
    """

    def ratio_minus_one(z: float) -> float:
        return float(thomson_rate(z, x_e(z), cosmo) / cosmo.hubble(z)) - 1.0

    return brentq(ratio_minus_one, z_lo, z_hi, xtol=1e-6)


@dataclass
class LastScattering:
    """Optical depth and visibility function on a redshift grid.

    Attributes
    ----------
    z : ndarray
        Increasing redshift grid starting at 0.
    tau : ndarray
        Thomson optical depth from 0 to z (no reionisation).
    dtau_dz : ndarray
        dtau/dz.
    g : ndarray
        Visibility function g(z) = dtau/dz exp(-tau) [per unit z].
    hubble : ndarray
        H(z) on the grid [s^-1].
    """

    z: NDArray[np.float64]
    tau: NDArray[np.float64]
    dtau_dz: NDArray[np.float64]
    g: NDArray[np.float64]
    hubble: NDArray[np.float64]

    def z_tau_one(self) -> float:
        """Redshift where tau = 1 (Planck's definition of z_*)."""
        spline = CubicSpline(self.z, self.tau - 1.0)
        i = int(np.searchsorted(self.tau, 1.0))
        return brentq(spline, self.z[i - 1], self.z[i])

    def z_peak(self) -> float:
        """Redshift of the maximum of g(z) per unit z, refined with a cubic spline."""
        return _spline_peak(self.z, self.g)

    def z_peak_conformal(self) -> float:
        """Redshift of the maximum of the visibility per unit conformal time.

        Since d(eta)/dz = -c / H(z), g_eta = -dtau/deta exp(-tau) is
        proportional to g(z) H(z), which peaks at slightly higher z.
        """
        return _spline_peak(self.z, self.g * self.hubble)

    def fwhm(self) -> tuple[float, float, float]:
        """Full width at half maximum of g(z).

        Returns
        -------
        tuple of float
            (FWHM, z_low, z_high), the latter two being the half-maximum
            crossing redshifts.
        """
        spline = CubicSpline(self.z, self.g)
        z_pk = self.z_peak()
        half = 0.5 * float(spline(z_pk))
        i_pk = int(np.argmax(self.g))
        lo = np.nonzero(self.g[:i_pk] < half)[0][-1]
        hi = i_pk + np.nonzero(self.g[i_pk:] < half)[0][0]
        z_lo = brentq(lambda z: spline(z) - half, self.z[lo], z_pk)
        z_hi = brentq(lambda z: spline(z) - half, z_pk, self.z[hi])
        return z_hi - z_lo, z_lo, z_hi


def _spline_peak(z: NDArray[np.float64], f: NDArray[np.float64]) -> float:
    """Location of the maximum of sampled f(z), refined with a cubic spline."""
    spline = CubicSpline(z, f)
    i = int(np.argmax(f))
    return brentq(spline.derivative(), z[i - 1], z[i + 1])


def last_scattering(
    x_e: XeFunction, cosmo: Cosmology, z_max: float = 3000.0, n: int = 60001
) -> LastScattering:
    """Compute tau(z) and g(z) for an ionisation history.

    tau(z) = int_0^z n_e sigma_T c / [(1+z') H(z')] dz', integrated with the
    trapezoidal rule on a uniform grid (convergence is checked in the tests).

    Parameters
    ----------
    x_e : callable
        Ionisation history x_e(z).
    cosmo : Cosmology
        Background cosmology.
    z_max : float, optional
        Upper redshift of the grid.
    n : int, optional
        Number of grid points.
    """
    z = np.linspace(0.0, z_max, n)
    hubble = cosmo.hubble(z)
    dtau_dz = thomson_rate(z, x_e(z), cosmo) / ((1.0 + z) * hubble)
    tau = cumulative_trapezoid(dtau_dz, z, initial=0.0)
    g = dtau_dz * np.exp(-tau)
    return LastScattering(z=z, tau=tau, dtau_dz=dtau_dz, g=g, hubble=hubble)
