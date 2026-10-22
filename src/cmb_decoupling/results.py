"""Key derived quantities of an ionisation history, computed in one place."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .cosmology import Cosmology
from .visibility import XeFunction, decoupling_redshift, last_scattering

#: Redshift at which the "freeze-out" residual ionisation is quoted. Chosen
#: after recombination has frozen out but before the T_m = T_gamma
#: approximation of this package becomes poor (see docs/theory.md).
Z_FREEZE_OUT: float = 200.0


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
