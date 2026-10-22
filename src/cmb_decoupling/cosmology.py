"""Background cosmology: Planck 2018 parameters, H(z), temperatures, densities.

All cosmological parameters are defined in :data:`PLANCK2018` and nowhere else.

References
----------
Planck Collaboration (2020), "Planck 2018 results. VI. Cosmological
parameters", A&A 641, A6, Table 2, column "TT,TE,EE+lowE+lensing".
Fixsen (2009), ApJ 707, 916 -- T_CMB = 2.7255 K.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .constants import C_LIGHT, G_NEWTON, M_H_ATOM, M_HE4_ATOM, MPC, SIGMA_SB

#: Ratio of neutrino to photon temperature, (4/11)^(1/3), and the fermion
#: factor 7/8, giving rho_nu / rho_gamma = N_eff * 7/8 * (4/11)^(4/3)
#: (Dodelson & Schmidt 2020, Modern Cosmology, 2nd ed., Ch. 2).
NEUTRINO_DENSITY_FACTOR: float = 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0)


@dataclass(frozen=True)
class Cosmology:
    """Flat LCDM background with massless neutrinos.

    Parameters
    ----------
    H0_km_s_Mpc : float
        Hubble constant [km s^-1 Mpc^-1].
    omega_b : float
        Physical baryon density Omega_b h^2 [dimensionless].
    omega_c : float
        Physical cold dark matter density Omega_c h^2 [dimensionless].
    T_cmb : float
        CMB temperature today [K].
    N_eff : float
        Effective number of relativistic neutrino species.
    Y_p : float
        Primordial helium mass fraction.
    """

    H0_km_s_Mpc: float
    omega_b: float
    omega_c: float
    T_cmb: float
    N_eff: float
    Y_p: float

    @property
    def h(self) -> float:
        """Dimensionless Hubble parameter h = H0 / (100 km s^-1 Mpc^-1)."""
        return self.H0_km_s_Mpc / 100.0

    @property
    def H0(self) -> float:
        """Hubble constant [s^-1]."""
        return self.H0_km_s_Mpc * 1.0e3 / MPC

    @property
    def rho_crit0(self) -> float:
        """Critical density today, 3 H0^2 / (8 pi G) [kg m^-3]."""
        return 3.0 * self.H0**2 / (8.0 * np.pi * G_NEWTON)

    @property
    def Omega_b(self) -> float:
        """Baryon density parameter."""
        return self.omega_b / self.h**2

    @property
    def Omega_m(self) -> float:
        """Matter density parameter (baryons + CDM; neutrinos are massless here)."""
        return (self.omega_b + self.omega_c) / self.h**2

    @property
    def Omega_gamma(self) -> float:
        """Photon density parameter from rho_gamma c^2 = (4 sigma_SB / c) T^4."""
        rho_gamma = 4.0 * SIGMA_SB * self.T_cmb**4 / C_LIGHT**3
        return rho_gamma / self.rho_crit0

    @property
    def Omega_r(self) -> float:
        """Radiation density parameter: photons + N_eff massless neutrinos."""
        return self.Omega_gamma * (1.0 + NEUTRINO_DENSITY_FACTOR * self.N_eff)

    @property
    def Omega_lambda(self) -> float:
        """Dark-energy density parameter fixed by flatness."""
        return 1.0 - self.Omega_m - self.Omega_r

    @property
    def n_H0(self) -> float:
        """Total hydrogen number density today, (1 - Y_p) rho_b / m_H [m^-3]."""
        return (1.0 - self.Y_p) * self.Omega_b * self.rho_crit0 / M_H_ATOM

    @property
    def f_He(self) -> float:
        """Helium-to-hydrogen number ratio n_He / n_H."""
        return self.Y_p / (1.0 - self.Y_p) * M_H_ATOM / M_HE4_ATOM

    def hubble(self, z: ArrayLike) -> NDArray[np.float64]:
        """Hubble rate H(z) [s^-1].

        H(z) = H0 sqrt(Omega_m (1+z)^3 + Omega_r (1+z)^4 + Omega_Lambda).
        """
        zp1 = 1.0 + np.asarray(z, dtype=float)
        E2 = self.Omega_m * zp1**3 + self.Omega_r * zp1**4 + self.Omega_lambda
        return self.H0 * np.sqrt(E2)

    def T_gamma(self, z: ArrayLike) -> NDArray[np.float64]:
        """Photon temperature T = T_cmb (1+z) [K]."""
        return self.T_cmb * (1.0 + np.asarray(z, dtype=float))

    def n_H(self, z: ArrayLike) -> NDArray[np.float64]:
        """Hydrogen number density n_H(z) = n_H0 (1+z)^3 [m^-3]."""
        return self.n_H0 * (1.0 + np.asarray(z, dtype=float)) ** 3


#: Planck 2018 VI, Table 2, TT,TE,EE+lowE+lensing (mean values).
#: Y_p is the BBN-consistent helium fraction quoted for this fit.
#: TODO: verify the Y_p = 0.2454 value and its location in Planck 2018 VI.
#: N_eff = 3.046 (standard-model value used by Planck 2018; Mangano et al. 2005).
PLANCK2018 = Cosmology(
    H0_km_s_Mpc=67.36,
    omega_b=0.02237,
    omega_c=0.1200,
    T_cmb=2.7255,
    N_eff=3.046,
    Y_p=0.2454,
)
