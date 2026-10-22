"""Physical and atomic constants in SI units.

All fundamental constants are taken from :mod:`astropy.constants`
(CODATA 2018 in astropy >= 4.0) and converted once to plain floats so that
they can be used inside ODE right-hand sides without unit overhead.
Atomic data that astropy does not provide are given with their source.
"""

from __future__ import annotations

from astropy import constants as ac
from astropy import units as u

# --- Fundamental constants (CODATA 2018 via astropy) -----------------------
C_LIGHT: float = ac.c.si.value  # speed of light [m s^-1]
K_B: float = ac.k_B.si.value  # Boltzmann constant [J K^-1]
H_PLANCK: float = ac.h.si.value  # Planck constant [J s]
HBAR: float = ac.hbar.si.value  # reduced Planck constant [J s]
M_E: float = ac.m_e.si.value  # electron mass [kg]
M_P: float = ac.m_p.si.value  # proton mass [kg]
G_NEWTON: float = ac.G.si.value  # gravitational constant [m^3 kg^-1 s^-2]
SIGMA_T: float = ac.sigma_T.si.value  # Thomson cross-section [m^2]
SIGMA_SB: float = ac.sigma_sb.si.value  # Stefan-Boltzmann constant [W m^-2 K^-4]
EV: float = (1.0 * u.eV).to(u.J).value  # electron-volt [J]
MPC: float = (1.0 * u.Mpc).to(u.m).value  # megaparsec [m]
AMU: float = ac.u.si.value  # atomic mass unit [kg]

# --- Hydrogen atomic data ---------------------------------------------------
#: Rydberg energy for infinite nuclear mass, R_inf h c [J].
RYDBERG_ENERGY: float = (ac.Ryd * ac.h * ac.c).to(u.J).value
#: Hydrogen ground-state ionisation energy with the reduced-mass correction,
#: E_ion = R_inf h c * m_p / (m_e + m_p) ~ 13.598 eV [J].
E_ION_H: float = RYDBERG_ENERGY * M_P / (M_E + M_P)
#: Energy of the n=2 level above the ground state (Lyman-alpha), 3/4 E_ion [J].
E_LYA: float = 0.75 * E_ION_H
#: Binding energy of the n=2 level, E_ion / 4 [J].
E_ION_H_N2: float = 0.25 * E_ION_H
#: Lyman-alpha wavelength, lambda = h c / E_Lya [m].
LAMBDA_LYA: float = H_PLANCK * C_LIGHT / E_LYA
#: Two-photon 2s -> 1s decay rate [s^-1]. Value used in RECFAST
#: (Seager, Sasselov & Scott 1999, 2000), from Goldman (1989).
#: TODO: verify the Goldman (1989) citation details.
LAMBDA_2S1S: float = 8.22458
#: Mass of the 1H atom, 1.00782503223 u (Wang et al. 2017, AME2016).
#: TODO: verify AME citation.
M_H_ATOM: float = 1.00782503223 * AMU

# --- Helium atomic data -----------------------------------------------------
#: First ionisation energy of He I, 24.587387 eV (NIST Atomic Spectra Database).
E_ION_HEI: float = 24.587387 * EV
#: Ionisation energy of He II, 54.417760 eV (NIST Atomic Spectra Database).
E_ION_HEII: float = 54.417760 * EV
#: Mass of the 4He atom, 4.00260325413 u (Wang et al. 2017, AME2016).
M_HE4_ATOM: float = 4.00260325413 * AMU
