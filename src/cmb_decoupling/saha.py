"""Saha equilibrium ionisation of hydrogen and helium.

References
----------
Dodelson & Schmidt (2020), Modern Cosmology, 2nd ed. (hydrogen).
Seager, Sasselov & Scott (2000), ApJS 128, 407, eqs. 80-81 -- helium Saha
relations, including the statistical-weight factors 4 (He I <-> He II) and 1
(He II <-> He III).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq

from .constants import E_ION_H, E_ION_HEI, E_ION_HEII, HBAR, K_B, M_E
from .cosmology import Cosmology

#: Lower bracket for the charge-neutrality root find (x_e = n_e / n_H).
MIN_X_E: float = 1.0e-30


def thermal_density(T: ArrayLike) -> NDArray[np.float64]:
    """Electron thermal (quantum concentration) density (m_e k T / 2 pi hbar^2)^{3/2}.

    Parameters
    ----------
    T : array_like
        Temperature [K].

    Returns
    -------
    ndarray
        Quantum concentration [m^-3].
    """
    T = np.asarray(T, dtype=float)
    return (M_E * K_B * T / (2.0 * np.pi * HBAR**2)) ** 1.5


def saha_rhs_hydrogen(z: ArrayLike, cosmo: Cosmology) -> NDArray[np.float64]:
    """Right-hand side S(z) of the hydrogen Saha equation x^2/(1-x) = S.

    S = (1/n_H) (m_e k T / 2 pi hbar^2)^{3/2} exp(-E_ion / k T).

    Parameters
    ----------
    z : array_like
        Redshift.
    cosmo : Cosmology
        Background cosmology.

    Returns
    -------
    ndarray
        Dimensionless S(z).
    """
    T = cosmo.T_gamma(z)
    return thermal_density(T) / cosmo.n_H(z) * np.exp(-E_ION_H / (K_B * T))


def x_e_saha_hydrogen(z: ArrayLike, cosmo: Cosmology) -> NDArray[np.float64]:
    """Saha ionisation fraction of pure hydrogen, x_e = n_e / n_H.

    Solves x^2 / (1 - x) = S with the cancellation-free root
    x = 2 / (1 + sqrt(1 + 4/S)).

    Parameters
    ----------
    z : array_like
        Redshift.
    cosmo : Cosmology
        Background cosmology.

    Returns
    -------
    ndarray
        Free-electron fraction in [0, 1].
    """
    with np.errstate(divide="ignore", over="ignore"):
        S = saha_rhs_hydrogen(z, cosmo)
        return 2.0 / (1.0 + np.sqrt(1.0 + 4.0 / S))


def _saha_ratios(z: float, cosmo: Cosmology) -> tuple[float, float, float]:
    """Return Saha ratios (n_i+1 n_e / n_i) / n_H for H, He I and He II."""
    T = float(cosmo.T_gamma(z))
    pref = float(thermal_density(T)) / float(cosmo.n_H(z))
    kT = K_B * T
    r_H = pref * np.exp(-E_ION_H / kT)
    r_HeI = 4.0 * pref * np.exp(-E_ION_HEI / kT)
    r_HeII = pref * np.exp(-E_ION_HEII / kT)
    return r_H, r_HeI, r_HeII


def saha_species(z: float, cosmo: Cosmology) -> dict[str, float]:
    """Saha equilibrium for hydrogen and both helium ionisation stages.

    Solves charge neutrality
    x_e = x_HII + f_He (y_HeII + 2 y_HeIII)
    for x_e with Brent's method, where each species fraction follows from
    its Saha ratio at the given x_e.

    Parameters
    ----------
    z : float
        Redshift.
    cosmo : Cosmology
        Background cosmology.

    Returns
    -------
    dict
        ``x_e`` (n_e/n_H), ``x_HII`` (n_p/n_H), ``y_HeII`` and ``y_HeIII``
        (fractions of all helium nuclei).
    """
    r_H, r_HeI, r_HeII = _saha_ratios(z, cosmo)
    f_He = cosmo.f_He

    def fractions(x_e: float) -> tuple[float, float, float]:
        x_HII = r_H / (x_e + r_H)
        # log of n_HeI : n_HeII : n_HeIII = 1 : r_HeI/x_e : r_HeI r_HeII / x_e^2,
        # normalised in log space so that no ratio can overflow.
        log_a = np.log(r_HeI) - np.log(x_e)
        logs = np.array([0.0, log_a, log_a + np.log(r_HeII) - np.log(x_e)])
        w = np.exp(logs - logs.max())
        w /= w.sum()
        return float(x_HII), float(w[1]), float(w[2])

    def residual(log_x: float) -> float:
        x_e = np.exp(log_x)
        x_HII, y_HeII, y_HeIII = fractions(x_e)
        return x_HII + f_He * (y_HeII + 2.0 * y_HeIII) - x_e

    with np.errstate(divide="ignore", under="ignore"):
        if residual(np.log(MIN_X_E)) <= 0.0:
            # Ionisation is below MIN_X_E: helium is fully neutral, so the
            # hydrogen-only closed form is exact to machine precision.
            x_e = float(x_e_saha_hydrogen(z, cosmo))
            return {"x_e": x_e, "x_HII": x_e, "y_HeII": 0.0, "y_HeIII": 0.0}
        x_max = 1.0 + 2.0 * f_He
        log_x = brentq(residual, np.log(MIN_X_E), np.log(x_max), xtol=1e-14, rtol=1e-14)
        x_e = float(np.exp(log_x))
        x_HII, y_HeII, y_HeIII = fractions(x_e)
    return {"x_e": x_e, "x_HII": x_HII, "y_HeII": y_HeII, "y_HeIII": y_HeIII}


def x_e_saha_full(z: ArrayLike, cosmo: Cosmology) -> NDArray[np.float64]:
    """Saha x_e = n_e / n_H including hydrogen and helium (vectorised over z).

    Parameters
    ----------
    z : array_like
        Redshift.
    cosmo : Cosmology
        Background cosmology.

    Returns
    -------
    ndarray
        Free-electron fraction relative to hydrogen; can exceed 1
        (maximum 1 + 2 f_He when helium is doubly ionised).
    """
    z_arr = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.array([saha_species(zi, cosmo)["x_e"] for zi in z_arr])
    return out.reshape(np.shape(z))


def saha_redshift(cosmo: Cosmology, x_target: float = 0.5) -> float:
    """Redshift at which hydrogen-only Saha gives x_e = x_target.

    Parameters
    ----------
    cosmo : Cosmology
        Background cosmology.
    x_target : float, optional
        Target ionisation fraction in (0, 1).

    Returns
    -------
    float
        Redshift.
    """
    return brentq(lambda z: x_e_saha_hydrogen(z, cosmo) - x_target, 200.0, 5000.0)
