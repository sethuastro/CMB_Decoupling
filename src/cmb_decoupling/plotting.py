"""Plotting helpers: consistent, colour-blind-safe, publication-style figures.

Colours are from the Okabe & Ito (2008) palette, designed to remain
distinguishable under the common colour-vision deficiencies. Every series
also has its own line style, so identity never relies on colour alone.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

#: Okabe-Ito colours, assigned to models in a fixed order.
STYLES: dict[str, dict[str, object]] = {
    "camb": {"color": "#000000", "ls": "-", "lw": 2.2, "label": "CAMB 2 (RECFAST)"},
    "peebles": {"color": "#0072B2", "ls": "--", "lw": 1.8, "label": "Peebles (this work)"},
    "peebles_fudge": {"color": "#009E73", "ls": "-.", "lw": 1.6, "label": r"Peebles, $F=1.14$ (this work)"},
    "saha_full": {"color": "#D55E00", "ls": ":", "lw": 1.8, "label": "Saha H + He (this work)"},
    "saha_h": {"color": "#CC79A7", "ls": (0, (5, 1, 1, 1)), "lw": 1.4, "label": "Saha H only (this work)"},
}

GUIDE_COLOR = "#7f7f7f"


def use_style() -> None:
    """Apply the repository-wide matplotlib style."""
    mpl.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "legend.fontsize": 8.5,
            "legend.frameon": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.6,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.bbox": "tight",
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
        }
    )


def save(fig: Figure, outdir: Path, name: str) -> list[Path]:
    """Save ``fig`` as both PDF and PNG in ``outdir``.

    Parameters
    ----------
    fig : Figure
        Figure to save.
    outdir : Path
        Output directory (created if missing).
    name : str
        File stem.

    Returns
    -------
    list of Path
        Written file paths.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    paths = [outdir / f"{name}.{ext}" for ext in ("pdf", "png")]
    for p in paths:
        fig.savefig(p)
    plt.close(fig)
    return paths
