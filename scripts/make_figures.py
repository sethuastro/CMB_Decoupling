"""Regenerate every figure and the key-results table used in the README.

Usage
-----
    python scripts/make_figures.py

Outputs go to ``figures/`` (PDF + PNG) and ``results/`` (JSON + Markdown).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cmb_decoupling import PLANCK2018, solve_peebles, x_e_saha_full, x_e_saha_hydrogen
from cmb_decoupling.peebles import RECFAST_FUDGE, peebles_C
from cmb_decoupling.plotting import GUIDE_COLOR, STYLES, save, use_style
from cmb_decoupling.reference import camb_reference
from cmb_decoupling.results import Z_FREEZE_OUT, key_results
from cmb_decoupling.visibility import last_scattering, thomson_rate

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
RES_DIR = ROOT / "results"

Z_PLOT_MIN, Z_PLOT_MAX = 100.0, 3000.0
N_PLOT = 2000
Z_TICKS = (100, 200, 300, 500, 700, 1000, 2000, 3000)


def no_label(key: str) -> dict[str, object]:
    """Style for ``key`` without its legend label."""
    return {k: v for k, v in STYLES[key].items() if k != "label"}


def set_log_z_axis(ax: plt.Axes, z_max: float = Z_PLOT_MAX, z_min: float = Z_PLOT_MIN) -> None:
    """Log redshift axis decreasing to the right (time increasing), with labelled ticks."""
    ax.set_xscale("log")
    ax.set_xlim(z_max, z_min)
    ticks = [t for t in Z_TICKS if z_min <= t <= z_max]
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel(r"Redshift $z$")


def build_models() -> dict[str, object]:
    """Compute every ionisation history shown in the figures."""
    cosmo = PLANCK2018
    peebles = solve_peebles(cosmo)
    fudge = solve_peebles(cosmo, fudge=RECFAST_FUDGE)
    ref = camb_reference(cosmo)
    return {
        "camb": ref.x_e,
        "peebles": peebles.x_e,
        "peebles_fudge": fudge.x_e,
        "saha_full": lambda z: x_e_saha_full(z, cosmo),
        "saha_h": lambda z: x_e_saha_hydrogen(z, cosmo),
        "_ref": ref,
        "_peebles": peebles,
    }


def fig_ionisation(models: dict[str, object]) -> None:
    """x_e(z) for all models with relative residuals against CAMB."""
    z = np.geomspace(Z_PLOT_MIN, Z_PLOT_MAX, N_PLOT)
    xe_ref = models["camb"](z)
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(6.4, 6.0), sharex=True, gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.06}
    )
    for key in ("camb", "saha_h", "saha_full", "peebles", "peebles_fudge"):
        ax.plot(z, models[key](z), **STYLES[key])
    for key in ("saha_full", "peebles", "peebles_fudge"):
        axr.plot(z, 100.0 * (models[key](z) / xe_ref - 1.0), **no_label(key))
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 3.0)
    ax.set_ylabel(r"Free-electron fraction $x_e = n_e/n_{\rm H}$")
    ax.legend(loc="upper right")
    axr.axhline(0.0, color=GUIDE_COLOR, lw=0.8)
    axr.set_ylim(-10, 20)
    axr.set_ylabel(r"$x_e/x_e^{\rm CAMB} - 1$ [%]")
    set_log_z_axis(axr)
    save(fig, FIG_DIR, "ionisation_history")


def fig_visibility(models: dict[str, object]) -> None:
    """Visibility function g(z) and optical depth tau(z)."""
    cosmo = PLANCK2018
    fig, (ax_g, ax_t) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    ls_p = last_scattering(models["peebles"], cosmo)
    _, z_lo, z_hi = ls_p.fwhm()
    ax_g.axvspan(z_lo, z_hi, color=STYLES["peebles"]["color"], alpha=0.10, lw=0, label="FWHM (Peebles)")
    for key in ("camb", "saha_full", "peebles"):
        ls = last_scattering(models[key], cosmo)
        ax_g.plot(ls.z, ls.g * 1e3, **STYLES[key])
        ax_t.plot(ls.z, ls.tau, **STYLES[key])
    ax_g.set_xlim(700, 1600)
    ax_g.set_ylim(0, 8)
    ax_g.set_xlabel(r"Redshift $z$")
    ax_g.set_ylabel(r"$g(z) = \frac{d\tau}{dz}\,e^{-\tau}$  [$10^{-3}$ per unit $z$]")
    ax_g.legend(loc="upper right")
    ax_t.axhline(1.0, color=GUIDE_COLOR, lw=0.8, ls="--")
    ax_t.set_yscale("log")
    ax_t.set_xlim(0, 2000)
    ax_t.set_ylim(1e-3, 1e3)
    ax_t.set_xlabel(r"Redshift $z$")
    ax_t.set_ylabel(r"Thomson optical depth $\tau(z)$")
    fig.tight_layout()
    save(fig, FIG_DIR, "visibility_function")


def fig_rates(models: dict[str, object]) -> None:
    """Gamma_T / H and the Peebles C_r factor."""
    cosmo = PLANCK2018
    z = np.geomspace(Z_PLOT_MIN, Z_PLOT_MAX, N_PLOT)
    fig, (ax, ax_c) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    for key in ("camb", "saha_full", "peebles"):
        ratio = thomson_rate(z, models[key](z), cosmo) / cosmo.hubble(z)
        ax.plot(z, ratio, **STYLES[key])
    ax.axhline(1.0, color=GUIDE_COLOR, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 1e5)
    set_log_z_axis(ax)
    ax.set_ylabel(r"$\Gamma_T / H = n_e \sigma_T c / H$")
    ax.legend(loc="upper right")
    hist = models["_peebles"]
    z_c = z[z < hist.z_switch]
    C = [peebles_C(zi, float(hist.x_e(zi)), cosmo) for zi in z_c]
    ax_c.plot(z_c, C, **no_label("peebles"))
    set_log_z_axis(ax_c, z_max=hist.z_switch)
    ax_c.set_ylim(0, 1.02)
    ax_c.set_ylabel(r"Peebles factor $C_r$")
    fig.tight_layout()
    save(fig, FIG_DIR, "rates")


def write_results(models: dict[str, object]) -> None:
    """Write the key-results table as JSON and Markdown."""
    cosmo = PLANCK2018
    rows = {key: key_results(models[key], cosmo).as_dict() for key in ("peebles", "peebles_fudge", "camb")}
    rows["camb"]["z_star_camb_internal"] = models["_ref"].z_star
    RES_DIR.mkdir(exist_ok=True)
    (RES_DIR / "key_results.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    labels = {
        "z_decoupling": r"$z_{\rm dec}$ where $\Gamma_T = H$",
        "z_tau_one": r"$z_*$ where $\tau = 1$",
        "z_peak": "peak of $g(z)$ (per unit $z$)",
        "z_peak_conformal": "peak of visibility (per conformal time)",
        "fwhm": r"FWHM of $g(z)$ [$\Delta z$]",
        "x_e_freeze_out": f"$x_e$ at $z = {Z_FREEZE_OUT:.0f}$",
        "x_e_today": "$x_e$ at $z = 0$ (no reionisation)",
    }
    lines = [
        "| Quantity | Peebles (this work) | Peebles, F=1.14 (this work) | CAMB (RECFAST) | Peebles − CAMB |",
        "|---|---|---|---|---|",
    ]
    for key, label in labels.items():
        p, f, c = rows["peebles"][key], rows["peebles_fudge"][key], rows["camb"][key]
        if key.startswith("x_e"):
            lines.append(f"| {label} | {p:.3e} | {f:.3e} | {c:.3e} | {100 * (p / c - 1):+.1f}% |")
        else:
            lines.append(f"| {label} | {p:.2f} | {f:.2f} | {c:.2f} | {p - c:+.2f} |")
    (RES_DIR / "key_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"CAMB internal z_* = {models['_ref'].z_star:.3f}")


def main() -> None:
    """Build all models, figures and tables."""
    # The summary table contains non-ASCII characters (e.g. the minus sign);
    # Windows consoles default to cp1252, which cannot print them.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    use_style()
    models = build_models()
    fig_ionisation(models)
    fig_visibility(models)
    fig_rates(models)
    write_results(models)


if __name__ == "__main__":
    main()
