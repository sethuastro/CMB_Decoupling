"""Validate the Peebles ionisation history against CAMB (RECFAST) on a matched background.

Usage
-----
    python scripts/validate.py

Both codes use the same Planck 2018 background (massless neutrinos, same Y_p
and T_cmb, no reionisation), so differences come from the recombination
physics alone. Nothing here is fitted or tuned to the reference.

Outputs
-------
figures/validation.{pdf,png}
    (a) x_e(z), (b) fractional residual Delta x_e / x_e against CAMB over
    z = 200-2000 for the pure Peebles model and the F = 1.14 literature variant,
    (c) the visibility function g(z).
results/validation.{json,md}
    Last-scattering quantities and the residual summary, also printed to the console.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cmb_decoupling import PLANCK2018, last_scattering, solve_peebles, x_e_saha_full
from cmb_decoupling.peebles import DEFAULT_Z_SWITCH, RECFAST_FUDGE
from cmb_decoupling.plotting import GUIDE_COLOR, STYLES, save, use_style
from cmb_decoupling.reference import camb_reference
from cmb_decoupling.results import (
    Z_FREEZE_OUT,
    Z_VALIDATION_MAX,
    Z_VALIDATION_MIN,
    fractional_residual,
    residual_summary,
)

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
RES_DIR = ROOT / "results"

N_PLOT = 2000
Z_TICKS = (200, 300, 500, 700, 1000, 1500, 2000)
#: Redshift range shown in the visibility-function panel.
G_PLOT_RANGE = (700.0, 1500.0)
#: Upper redshift for the helium-era residual (He I and He II recombination in Saha).
Z_HELIUM_MAX = 3000.0
#: Grid points on [z_switch, Z_HELIUM_MAX] (spacing 0.1 in z, as for the main range).
N_HELIUM = 12001

MODELS = ("peebles", "peebles_fudge")
LABELS = {"peebles": "Peebles (this work)", "peebles_fudge": f"Peebles, F = {RECFAST_FUDGE} (this work)"}


def compute() -> dict[str, object]:
    """Run both codes and collect every number reported by this script."""
    cosmo = PLANCK2018
    ref = camb_reference(cosmo)
    x_e = {
        "camb": ref.x_e,
        "peebles": solve_peebles(cosmo).x_e,
        "peebles_fudge": solve_peebles(cosmo, fudge=RECFAST_FUDGE).x_e,
    }
    ls = {k: last_scattering(f, cosmo) for k, f in x_e.items()}
    last_sc = {k: {"z_star": v.z_tau_one(), "z_peak": v.z_peak(), "fwhm": v.fwhm()[0]} for k, v in ls.items()}
    resid = {k: residual_summary(x_e[k], ref.x_e, cosmo).as_dict() for k in MODELS}
    for k in MODELS:
        resid[k]["residual_at_freeze_out"] = float(
            fractional_residual(x_e[k], ref.x_e, np.array(Z_FREEZE_OUT))
        )
    # Hydrogen-only era (Peebles ODE region) and the helium era beyond the main range.
    z_switch = DEFAULT_Z_SWITCH
    hydrogen_era = {
        k: residual_summary(x_e[k], ref.x_e, cosmo, z_min=Z_VALIDATION_MIN, z_max=z_switch).as_dict()
        for k in MODELS
    }
    z_he = np.linspace(z_switch, Z_HELIUM_MAX, N_HELIUM)
    r_he = fractional_residual(x_e["peebles"], ref.x_e, z_he)
    i_he = int(np.argmin(r_he))
    saha = last_scattering(lambda z: x_e_saha_full(z, cosmo), cosmo)
    rec = ref.params.Recomb
    recfast = {
        "fudge": rec.RECFAST_fudge,
        "Hswitch": bool(rec.RECFAST_Hswitch),
        "gauss_centres_z": [float(np.exp(rec.zGauss1) - 1.0), float(np.exp(rec.zGauss2) - 1.0)],
    }
    return {
        "x_e": x_e,
        "ls": ls,
        "numbers": {
            "camb_version": __import__("camb").__version__,
            "camb_recombination": ref.recombination_model,
            "last_scattering": last_sc,
            "camb_recfast_settings": recfast,
            "residuals_vs_camb": resid,
            "residuals_hydrogen_era": hydrogen_era,
            "residuals_helium_era": {
                "z_min": z_switch,
                "z_max": Z_HELIUM_MAX,
                "min_residual": float(r_he[i_he]),
                "z_min_residual": float(z_he[i_he]),
            },
            "tau_integrator_check": {
                "z_star_from_camb_xe": ls["camb"].z_tau_one(),
                "z_star_camb_internal": ref.z_star,
            },
            "saha": {"z_star": saha.z_tau_one(), "z_peak": saha.z_peak()},
        },
    }


def make_figure(data: dict[str, object]) -> None:
    """Three-panel validation figure: x_e, fractional residual, visibility function."""
    x_e, ls, num = data["x_e"], data["ls"], data["numbers"]
    window = (
        num["residuals_vs_camb"]["peebles"]["window_lo"],
        num["residuals_vs_camb"]["peebles"]["window_hi"],
    )
    z = np.geomspace(Z_VALIDATION_MIN, Z_VALIDATION_MAX, N_PLOT)

    fig = plt.figure(figsize=(6.4, 8.4))
    # Row 2 is an empty spacer: panel (c) has its own redshift axis.
    gs = fig.add_gridspec(4, 1, height_ratios=[2.0, 1.4, 0.55, 2.0], hspace=0.08)
    ax_x = fig.add_subplot(gs[0])
    ax_r = fig.add_subplot(gs[1], sharex=ax_x)
    ax_g = fig.add_subplot(gs[3])

    for k in ("camb", *MODELS):
        ax_x.plot(z, x_e[k](z), **STYLES[k])
    ax_x.set_yscale("log")
    ax_x.set_ylabel(r"Free-electron fraction $x_e = n_e/n_{\rm H}$")
    ax_x.legend(loc="lower left")
    ax_x.tick_params(labelbottom=False)

    for ax in (ax_x, ax_r):
        ax.axvspan(*window, color=GUIDE_COLOR, alpha=0.15, lw=0)
    ax_r.axhline(0.0, color=GUIDE_COLOR, lw=0.8)
    for k in MODELS:
        ax_r.plot(z, 100.0 * fractional_residual(x_e[k], x_e["camb"], z), **STYLES[k])
    ax_r.set_ylabel(r"$\Delta x_e / x_e^{\rm CAMB}$ [%]")
    ax_r.text(
        np.sqrt(window[0] * window[1]), 0.95, "last-scattering\nwindow", transform=ax_r.get_xaxis_transform(),
        ha="center", va="top", fontsize=8, color="#404040",
    )  # fmt: skip
    ax_r.set_xscale("log")
    ax_r.set_xlim(Z_VALIDATION_MAX, Z_VALIDATION_MIN)
    ax_r.set_xticks(Z_TICKS, [str(t) for t in Z_TICKS])
    ax_r.xaxis.set_minor_locator(plt.NullLocator())
    ax_r.set_xlabel(r"Redshift $z$")

    ax_g.axvspan(*window, color=GUIDE_COLOR, alpha=0.15, lw=0, label="FWHM of CAMB $g(z)$")
    for k in ("camb", *MODELS):
        ax_g.plot(ls[k].z, ls[k].g * 1e3, **STYLES[k])
    ax_g.set_xlim(*G_PLOT_RANGE[::-1])
    ax_g.set_ylim(bottom=0)
    ax_g.set_xlabel(r"Redshift $z$")
    ax_g.set_ylabel(r"$g(z) = \frac{d\tau}{dz}\,e^{-\tau}$  [$10^{-3}$ per unit $z$]")
    ax_g.legend(loc="upper left")

    # Panel tags placed where no curve passes.
    for ax, tag, y, va in ((ax_x, "a", 0.95, "top"), (ax_r, "b", 0.06, "bottom"), (ax_g, "c", 0.95, "top")):
        ax.text(0.98, y, f"({tag})", transform=ax.transAxes, ha="right", va=va, fontsize=11)
    save(fig, FIG_DIR, "validation")


def write_tables(num: dict[str, object]) -> str:
    """Write results/validation.{json,md} and return the Markdown."""
    ls, res = num["last_scattering"], num["residuals_vs_camb"]
    cols = ("peebles", "peebles_fudge", "camb")
    head = f"| | {LABELS['peebles']} | {LABELS['peebles_fudge']} | CAMB {num['camb_version']} (RECFAST) |"
    lines = [head, "|---|---|---|---|"]
    for key, label in (
        ("z_star", r"$z_\ast$ where $\tau = 1$"),
        ("z_peak", "Peak of $g(z)$ per unit $z$"),
        ("fwhm", "FWHM of $g(z)$"),
    ):
        lines.append(f"| {label} | " + " | ".join(f"{ls[c][key]:.2f}" for c in cols) + " |")

    def pct(x: float) -> str:
        return f"{100 * x:+.2f}%"

    p = res["peebles"]
    lines += [
        "",
        r"Residuals $\Delta x_e/x_e = x_e/x_e^{\rm CAMB} - 1$ "
        f"on $z = {p['z_min']:.0f}$–${p['z_max']:.0f}$ (grid spacing $\\Delta z = 0.1$). "
        "Last-scattering window = FWHM of CAMB's $g(z)$, "
        f"$z = {p['window_lo']:.1f}$–${p['window_hi']:.1f}$.",
        "",
        f"| | {LABELS['peebles']} | {LABELS['peebles_fudge']} |",
        "|---|---|---|",
    ]
    rows = (
        (
            r"(a) max $\lvert\Delta x_e/x_e\rvert$ in the last-scattering window",
            "max_abs_window",
            "z_max_abs_window",
            True,
        ),
        ("(b) most negative residual", "min_residual", "z_min_residual", False),
        ("(c) most positive residual", "max_residual", "z_max_residual", False),
    )
    for label, key, zkey, absolute in rows:
        cells = [
            (f"{100 * res[m][key]:.2f}%" if absolute else pct(res[m][key])) + f" at $z = {res[m][zkey]:.1f}$"
            for m in MODELS
        ]
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines.append(
        f"| residual at $z = {Z_FREEZE_OUT:.0f}$ | "
        + " | ".join(pct(res[m]["residual_at_freeze_out"]) for m in MODELS)
        + " |"
    )
    h_era = num["residuals_hydrogen_era"]
    he = num["residuals_helium_era"]
    lines.append(
        f"| range on $z = {Z_VALIDATION_MIN:.0f}$–${DEFAULT_Z_SWITCH:.0f}$ (hydrogen only) | "
        + " | ".join(f"{pct(h_era[m]['min_residual'])} to {pct(h_era[m]['max_residual'])}" for m in MODELS)
        + " |"
    )
    tc = num["tau_integrator_check"]
    rf = num["camb_recfast_settings"]
    lines += [
        "",
        f"Helium era, $z = {he['z_min']:.0f}$–${he['z_max']:.0f}$: most negative residual "
        f"{pct(he['min_residual'])} at $z = {he['z_min_residual']:.1f}$ (both models; helium is in Saha).",
        "",
        "$\\tau$-integrator check: CAMB's $x_e$ through this code gives "
        f"$z_\\ast = {tc['z_star_from_camb_xe']:.4f}$; "
        f"CAMB's internal value is {tc['z_star_camb_internal']:.4f} "
        f"(difference {tc['z_star_from_camb_xe'] - tc['z_star_camb_internal']:+.4f}).",
        "",
        f"CAMB RECFAST settings: $F = {rf['fudge']}$, Hswitch = {rf['Hswitch']} "
        f"(Gaussian corrections centred at $z \\approx {rf['gauss_centres_z'][0]:.0f}$ "
        f"and {rf['gauss_centres_z'][1]:.0f}).",
    ]
    saha = num["saha"]
    lines += [
        "",
        f"Saha (H + He) for comparison: $z_\\ast = {saha['z_star']:.2f}$, "
        f"peak of $g(z)$ at {saha['z_peak']:.2f}.",
    ]
    md = "\n".join(lines) + "\n"
    RES_DIR.mkdir(exist_ok=True)
    (RES_DIR / "validation.json").write_text(json.dumps(num, indent=2) + "\n", encoding="utf-8")
    (RES_DIR / "validation.md").write_text(md, encoding="utf-8")
    return md


def main() -> None:
    """Compute, plot and tabulate the validation against CAMB."""
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    use_style()
    data = compute()
    make_figure(data)
    print(write_tables(data["numbers"]))


if __name__ == "__main__":
    main()
