"""Reference ionisation history from CAMB (RECFAST recombination).

CAMB (Lewis, Challinor & Lasenby 2000, ApJ 538, 473) is run with *the same*
background parameters as :class:`~cmb_decoupling.cosmology.Cosmology`
(massless neutrinos, same Y_p and T_cmb) and with reionisation switched off,
so that residuals isolate differences in the recombination physics.
CAMB's default recombination module is RECFAST (Seager et al. 1999, 2000,
with later corrections; Wong, Moss & Scott 2008), which includes helium
recombination, matter-temperature evolution and fudge-factor corrections.
HyRec and CosmoRec are not compiled into the pip-distributed CAMB wheel.
"""

from __future__ import annotations

from dataclasses import dataclass

import camb
import numpy as np
from numpy.typing import ArrayLike, NDArray

from .cosmology import Cosmology


@dataclass
class CambReference:
    """CAMB background results for a given cosmology.

    Attributes
    ----------
    results : camb.CAMBdata
        CAMB background results.
    params : camb.CAMBparams
        CAMB parameters used.
    """

    results: camb.CAMBdata
    params: camb.CAMBparams

    def x_e(self, z: ArrayLike) -> NDArray[np.float64]:
        """CAMB free-electron fraction n_e / n_H (includes helium electrons)."""
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        xe = self.results.get_background_redshift_evolution(z_arr, ["x_e"], format="array")[:, 0]
        return xe.reshape(np.shape(z))

    def hubble(self, z: ArrayLike) -> NDArray[np.float64]:
        """CAMB Hubble rate [s^-1]."""
        km_per_mpc = 1.0e3 / camb.constants.Mpc
        return np.asarray(self.results.hubble_parameter(np.asarray(z, dtype=float))) * km_per_mpc

    @property
    def z_star(self) -> float:
        """CAMB z_* (tau = 1, excluding reionisation)."""
        return float(self.results.get_derived_params()["zstar"])

    @property
    def recombination_model(self) -> str:
        """Name of CAMB's recombination class."""
        return type(self.params.Recomb).__name__


def camb_reference(cosmo: Cosmology) -> CambReference:
    """Run CAMB (background only) matched to ``cosmo``.

    Parameters
    ----------
    cosmo : Cosmology
        Background cosmology to reproduce.

    Returns
    -------
    CambReference
        Wrapped CAMB results.
    """
    params = camb.set_params(
        H0=cosmo.H0_km_s_Mpc,
        ombh2=cosmo.omega_b,
        omch2=cosmo.omega_c,
        TCMB=cosmo.T_cmb,
        nnu=cosmo.N_eff,
        mnu=0.0,
        YHe=cosmo.Y_p,
    )
    params.Reion.Reionization = False
    results = camb.get_background(params)
    return CambReference(results=results, params=params)
