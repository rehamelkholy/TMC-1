"""
Shared figure style: vector PDF + PNG preview, 3.5in (one column) / 7.1in (two column)
widths, fonts >=8pt, no in-plot titles, colorblind-safe palette (Okabe & Ito 2008), same
color/label per method across every figure.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ONE_COL_W = 3.5
TWO_COL_W = 7.1

# Okabe-Ito colorblind-safe palette
OKABE_ITO = {
    "black": "#000000", "orange": "#E69F00", "sky_blue": "#56B4E9",
    "green": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
    "vermillion": "#D55E00", "purple": "#CC79A7",
}

METHOD_LABEL = {
    "combined": "Combined ranking",
    "detectability": "Detectability",
    "high_logcd": "Predicted abundance",
    "low_dist": "Similarity",
    "high_sigma": "Uncertainty sampling",
    "ucb": "UCB",
    "random": "Random",
}

METHOD_COLOR = {
    "combined": OKABE_ITO["black"],
    "detectability": OKABE_ITO["vermillion"],
    "high_logcd": OKABE_ITO["orange"],
    "low_dist": OKABE_ITO["sky_blue"],
    "high_sigma": OKABE_ITO["green"],
    "ucb": OKABE_ITO["purple"],
    "random": "#999999",
}

METHOD_LINESTYLE = {
    "combined": "-",
    "detectability": "-",
    "high_logcd": "-",
    "low_dist": "-",
    "high_sigma": "-",
    "ucb": "-",
    "random": "--",
}

REGIME_LABEL = {
    "released_2021": "Released 2021 (pre-registered)",
    "refit_2021": "Refit at 2021 (post hoc)",
    "lso_v5": "Updated leave-self-out (post hoc)",
}
REGIME_POSTHOC = {"released_2021": False, "refit_2021": True, "lso_v5": True}


def apply_rcparams():
    plt.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,   # embed as real fonts, not paths
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def save_fig(fig, out_stem, dpi=300):
    """Writes out_stem.pdf (vector) and out_stem.png (raster preview).
    No in-plot title is ever set here -- captions carry it."""
    fig.savefig(out_stem + ".pdf", bbox_inches="tight")
    fig.savefig(out_stem + ".png", dpi=dpi, bbox_inches="tight")
