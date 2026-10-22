"""Peebles three-level-atom model of hydrogen recombination.

The free-electron fraction obeys (Peebles 1968; notation of Seager,
Sasselov & Scott 1999, eq. 1)

    dx_e/dz = C_r / [H(z)(1+z)] * [alpha_B n_H x_e^2
                                   - beta_B (1 - x_e) exp(-E_Lya / k T)],

    C_r = (Lambda_2s1s + Lambda_alpha) / (Lambda_2s1s + Lambda_alpha + beta_B),
    Lambda_alpha = 8 pi H / (lambda_Lya^3 n_1s),   n_1s = (1 - x_e) n_H,
    beta_B = alpha_B (m_e k T / 2 pi hbar^2)^{3/2} exp(-E_ion / 4 k T).

Approximations (documented in docs/theory.md):

* matter temperature equals radiation temperature, T_m = T_gamma;
* helium is fully neutral below the switch redshift, so n_e = n_p;
* above the switch redshift, H + He are in Saha equilibrium;
* no reionisation.

References
----------
Peebles (1968), ApJ 153, 1.
Seager, Sasselov & Scott (1999), ApJ 523, L1.
Seager, Sasselov & Scott (2000), ApJS 128, 407.
Pequignot, Petitjean & Boisson (1991), A&A 251, 680 -- case-B fit.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import OdeSolution, solve_ivp

from .constants import E_ION_H, K_B, LAMBDA_2S1S, LAMBDA_LYA
from .cosmology import Cosmology
from .saha import thermal_density, x_e_saha_full, x_e_saha_hydrogen

# Case-B recombination coefficient fit (Pequignot et al. 1991), as given in
# Seager et al. (1999), eq. 3:  alpha_B = 1e-19 a t^b / (1 + c t^d) m^3 s^-1,
# with t = T / 1e4 K.
ALPHA_B_NORM: float = 1.0e-19  # [m^3 s^-1]
ALPHA_B_A: float = 4.309
ALPHA_B_B: float = -0.6166
ALPHA_B_C: float = 0.6703
ALPHA_B_D: float = 0.5300
ALPHA_B_T_UNIT: float = 1.0e4  # [K]

#: RECFAST fudge factor that mimics the multi-level atom in the 3-level model
#: (Seager et al. 1999). Default in this package is 1 (pure Peebles model).
RECFAST_FUDGE: float = 1.14

#: Redshift at which the Saha (H + He) solution hands over to the Peebles ODE.
DEFAULT_Z_SWITCH: float = 1800.0

#: Integration tolerances for solve_ivp (Radau, stiff).
RTOL: float = 1.0e-8
ATOL: float = 1.0e-12

#: Floor on the neutral fraction to avoid division by zero in Lambda_alpha.
MIN_NEUTRAL_FRACTION: float = 1.0e-30


def alpha_B(T: ArrayLike, fudge: float = 1.0) -> NDArray[np.float64]:
    """Case-B recombination coefficient alpha_B(T).

    Parameters
    ----------
    T : array_like
        Temperature [K].
    fudge : float, optional
        Multiplicative factor F (1 for the pure Peebles model).

    Returns
    -------
    ndarray
        alpha_B [m^3 s^-1].
    """
    t = np.asarray(T, dtype=float) / ALPHA_B_T_UNIT
    return fudge * ALPHA_B_NORM * ALPHA_B_A * t**ALPHA_B_B / (1.0 + ALPHA_B_C * t**ALPHA_B_D)


def beta_B(T: ArrayLike, fudge: float = 1.0) -> NDArray[np.float64]:
    """Photo-ionisation rate from n=2 by detailed balance with alpha_B.

    Parameters
    ----------
    T : array_like
        Radiation temperature [K].
    fudge : float, optional
        Same factor applied to alpha_B.

    Returns
    -------
    ndarray
        beta_B [s^-1].
    """
    T = np.asarray(T, dtype=float)
    return alpha_B(T, fudge) * thermal_density(T) * np.exp(-E_ION_H / (4.0 * K_B * T))


def peebles_C(z: float, x_e: float, cosmo: Cosmology, fudge: float = 1.0) -> float:
    """Peebles suppression factor C_r (probability an n=2 atom reaches 1s).

    Parameters
    ----------
    z : float
        Redshift.
    x_e : float
        Ionised hydrogen fraction.
    cosmo : Cosmology
        Background cosmology.
    fudge : float, optional
        Factor applied to alpha_B and beta_B.

    Returns
    -------
    float
        C_r in (0, 1].
    """
    T = float(cosmo.T_gamma(z))
    n_1s = max(1.0 - x_e, MIN_NEUTRAL_FRACTION) * float(cosmo.n_H(z))
    lambda_alpha = 8.0 * np.pi * float(cosmo.hubble(z)) / (LAMBDA_LYA**3 * n_1s)
    escape = LAMBDA_2S1S + lambda_alpha
    return escape / (escape + float(beta_B(T, fudge)))


def peebles_rhs(z: float, y: NDArray[np.float64], cosmo: Cosmology, fudge: float = 1.0) -> list[float]:
    """Right-hand side dx_e/dz of the Peebles equation.

    Parameters
    ----------
    z : float
        Redshift.
    y : ndarray
        One-element array holding x_e.
    cosmo : Cosmology
        Background cosmology.
    fudge : float, optional
        Factor applied to alpha_B and beta_B.

    Returns
    -------
    list of float
        [dx_e/dz].
    """
    x_e = y[0]
    T = float(cosmo.T_gamma(z))
    n_H = float(cosmo.n_H(z))
    a_B = float(alpha_B(T, fudge))
    # beta_B exp(-E_Lya/kT) = alpha_B (m_e kT/2 pi hbar^2)^{3/2} exp(-E_ion/kT)
    ionisation = a_B * float(thermal_density(T)) * np.exp(-E_ION_H / (K_B * T))
    net = a_B * n_H * x_e**2 - ionisation * (1.0 - x_e)
    C_r = peebles_C(z, x_e, cosmo, fudge)
    return [C_r * net / (float(cosmo.hubble(z)) * (1.0 + z))]


@dataclass
class IonizationHistory:
    """Ionisation history: Saha (H + He) above ``z_switch``, Peebles below.

    Attributes
    ----------
    cosmo : Cosmology
        Background cosmology used.
    z_switch : float
        Hand-over redshift.
    fudge : float
        Factor applied to alpha_B and beta_B.
    solution : OdeSolution
        Dense output of the Peebles ODE on [z_end, z_switch].
    """

    cosmo: Cosmology
    z_switch: float
    fudge: float
    solution: OdeSolution

    def x_e(self, z: ArrayLike) -> NDArray[np.float64]:
        """Free-electron fraction n_e / n_H at redshift z.

        Parameters
        ----------
        z : array_like
            Redshift (>= the integration end point).

        Returns
        -------
        ndarray
            x_e(z).
        """
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        out = np.empty_like(z_arr)
        low = z_arr < self.z_switch
        if np.any(low):
            out[low] = self.solution(z_arr[low])[0]
        if np.any(~low):
            out[~low] = x_e_saha_full(z_arr[~low], self.cosmo)
        return out.reshape(np.shape(z))


def solve_peebles(
    cosmo: Cosmology,
    z_switch: float = DEFAULT_Z_SWITCH,
    z_end: float = 0.0,
    fudge: float = 1.0,
    rtol: float = RTOL,
    atol: float = ATOL,
) -> IonizationHistory:
    """Integrate the Peebles equation from ``z_switch`` down to ``z_end``.

    The initial condition is hydrogen Saha equilibrium at ``z_switch``.
    The ODE is stiff at high z (rates >> H), so the implicit Radau method is
    used.

    Parameters
    ----------
    cosmo : Cosmology
        Background cosmology.
    z_switch : float, optional
        Start redshift of the ODE integration.
    z_end : float, optional
        End redshift.
    fudge : float, optional
        Factor applied to alpha_B and beta_B (1 = pure Peebles).
    rtol, atol : float, optional
        Relative and absolute tolerances passed to ``solve_ivp``.

    Returns
    -------
    IonizationHistory
        Callable ionisation history.
    """
    x0 = float(x_e_saha_hydrogen(z_switch, cosmo))
    sol = solve_ivp(
        peebles_rhs,
        (z_switch, z_end),
        [x0],
        method="Radau",
        args=(cosmo, fudge),
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(f"Peebles integration failed: {sol.message}")
    return IonizationHistory(cosmo=cosmo, z_switch=z_switch, fudge=fudge, solution=sol.sol)
