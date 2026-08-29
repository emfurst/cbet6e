"""charts — the drawing craft shared by every property chart in the book.

Line weights, the two-tier major/minor treatment, chart-paper grids, and the
label-placement layer that writes a value *on* its own curve at a chosen height.
None of this is thermodynamics: it is the part of a property chart that is
draftsmanship, and it is identical whether the curves came from an equation of
state (`thermo.ph_chart`) or from a printed table (`thermo.steam_chart`).

WHY IT IS HERE. Figures 3.3-1(a), 3.3-1(b), 3.3-2 and 3.3-3 are chapter-3 figures,
but three more charts are *those charts with a process path drawn on them* and
belong to other chapters: Figure 5.1-3 (the LNG liquefaction path), `c05uf001`
(the Rankine cycle) and `c06uf002` (the cylinder-discharge window). A figure's
notebook lives in the figure's own chapter (author, 2026-08-08), so the machinery
they share cannot live in any one of them. Copied instead, it would be ~1,000
lines duplicated three ways, and the moment a parent chart is redrawn the copies
diverge silently -- which is the failure the recompute exists to prevent.

    from thermo.charts import chart_grid, label_at, label_end, use_book_style

`thermo` itself stays matplotlib-free: this module and the two chart modules are
loaded lazily (see `thermo/__init__.py`), so `from thermo import PengRobinson`
costs nothing extra.
"""
from __future__ import annotations

from shutil import which

import numpy as np


# --- line weights ----------------------------------------------------------
#
# The families must be distinguishable in black and white -- the book is printed
# that way -- so they differ in weight and dash pattern, not color. Each family is
# drawn in two tiers: a labeled MAJOR set at the weights below, and an unlabeled
# MINOR set between them at the `_minor` weights. THIS IS THE ONE PLACE TO CHANGE
# how heavy any family is on any chart in the book.
LW = {
    "sat": 1.7,                    # the saturation line: every chart's spine
    "P": 1.0, "P_minor": 0.25,     # isobars            (T-S, Mollier)
    "T": 1.0, "T_minor": 0.4,      # isotherms          (P-H)
    "S": 0.7, "S_minor": 0.3,      # isentropes         (P-H)
    "V": 0.7, "V_minor": 0.3,      # constant volume    (P-H)
    "x": 0.7, "x_minor": 0.3,      # constant quality
    "H": 0.7,                      # constant enthalpy  (T-S)
    # A MINOR TIER HAS TO BE TOLD APART FROM THE GRID, and the two tiers below cannot
    # do it on ink. They are drawn in GRAY_MINOR, 0.62, against a grid whose major
    # division is 0.70 -- lighter by almost nothing, and thicker at 0.45 pt -- so at
    # the 0.25-0.3 pt a minor tier would otherwise take, they read as grid. Both
    # therefore sit ABOVE the grid's 0.45 pt, and stay clear of their own labeled
    # tier. The minor isobars need none of this: at 0.25 gray they are near-black, and
    # separate from any grid tint on ink alone.
    "H_minor": 0.5,                # constant enthalpy  (T-S)
    "T_minor_mollier": 0.55,       # isotherms          (Mollier)
}

GRAY, GRAY_MINOR = "0.40", "0.62"

# The grid is read as much as the curves are, so it gets the same two-tier
# treatment: a labeled major division and five subdivisions inside it. Five is what
# makes the subdivisions land on round numbers -- 0.5/5 = 0.1 kJ/(kg K) and
# 100/5 = 20 kJ/kg -- which is the whole point of a grid you interpolate by eye.
GRID_MAJOR = dict(color="0.70", lw=0.45)
GRID_MINOR = dict(color="0.87", lw=0.25)


def use_book_style():
    """Typography for a figure that will be printed in the book.

    A real LaTeX installation when there is one, so the figure is set in Computer
    Modern like the rest of the book; matplotlib's own `mathtext` otherwise, so the
    notebook still runs in Colab, where there is no TeX. `fonttype 42` keeps text
    as text in the PDF rather than outlining it, which is what lets production
    re-set a label without regenerating the art.

    Call it once, before drawing. Every figure staged in
    `deliverables/figures/print-bw/` is set this way -- it is checked, and a figure
    that comes out in DejaVu Sans has not called this.
    """
    import matplotlib as mpl

    mpl.rcParams["text.usetex"] = which("latex") is not None
    mpl.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "cm",   # Computer Modern math for the no-TeX fallback
        "pdf.fonttype": 42,         # keep text as text in the PDF
        "ps.fonttype": 42,
    })
    return mpl.rcParams["text.usetex"]


def chart_grid(ax, x_major, y_major, minors=5, label_size=7):
    """Rule a linear-linear pair of axes like chart paper."""
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator

    ax.xaxis.set_major_locator(MultipleLocator(x_major))
    ax.yaxis.set_major_locator(MultipleLocator(y_major))
    if minors:
        ax.xaxis.set_minor_locator(AutoMinorLocator(minors))
        ax.yaxis.set_minor_locator(AutoMinorLocator(minors))
        ax.grid(True, which="minor", **GRID_MINOR)
    ax.grid(True, which="major", **GRID_MAJOR)
    ax.tick_params(labelsize=label_size)
    ax.tick_params(which="minor", length=2)
    ax.set_axisbelow(True)


def log_pressure_grid(ax, x_major, x_minors=5, label_size=7):
    """Rule a chart whose y axis is pressure on a log scale (the P-H charts).

    The decade subdivisions are the ones a reader actually interpolates between --
    1, 2, 3, 4, 5, 6, 8 get a labeled major line and the rest of the decade a minor
    one -- and the formatter prints `2` rather than `2 x 10^0`.
    """
    from matplotlib.ticker import (AutoMinorLocator, FuncFormatter, LogLocator,
                                   MultipleLocator, NullFormatter)

    ax.xaxis.set_major_locator(MultipleLocator(x_major))
    ax.xaxis.set_minor_locator(AutoMinorLocator(x_minors))
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 3, 4, 5, 6, 8)))
    ax.yaxis.set_minor_locator(
        LogLocator(base=10, subs=np.arange(1, 10) * 0.1, numticks=40))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.grid(True, which="minor", **GRID_MINOR)
    ax.grid(True, which="major", **GRID_MAJOR)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=label_size)


# --- labeling --------------------------------------------------------------

def angle_at(ax, x, y, i):
    """Local slope of a polyline at index `i`, in *display* degrees.

    Normalized to [-90, 90], so a curve stored right-to-left -- the saturation line
    is tabulated upward in temperature and so runs backwards in entropy -- does not
    come out with its label upside down. Display degrees, not data degrees: the
    label has to lie along the curve *as drawn*, which on a log axis or a
    non-square frame is not the same thing.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    j = [max(i - 1, 0), min(i + 1, len(x) - 1)]
    p = ax.transData.transform(np.column_stack([x[j], y[j]]))
    dx, dy = p[1] - p[0]
    a = float(np.degrees(np.arctan2(dy, dx)))
    return a - 180 if a > 90 else a + 180 if a < -90 else a


def label_at(ax, x, y, text, *, along, value, size=6.5, color="k", offset=(0, 0),
             pick="last", rotate=True, pad=0.6):
    """Write `text` on the polyline where `along` ('x'|'y') crosses `value`.

    Returns True if the label landed, False if the curve never crosses `value` or
    the crossing falls outside the axes -- so a caller can try a fallback height
    rather than silently dropping the label.

    A curve can cross `value` more than once -- a line of constant quality crosses a
    given enthalpy on both the rising and the falling side of the dome -- so `pick`
    ('first', 'last', 'max_x', 'min_x') chooses which crossing to label.

    `offset` is (x, y) in data units, EXCEPT that on a log y axis the y component is
    applied multiplicatively, as `y*(1 + offset[1])`. That is the only sane reading
    of "nudge it up a little" on a log scale, and it is what the P-H charts need;
    the linear charts get the additive form. Detected from the axis, not passed in,
    so a caller cannot get it wrong.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    c = x if along == "x" else y
    hits = np.nonzero((c[:-1] - value) * (c[1:] - value) <= 0)[0]
    if not len(hits):
        return False
    i = int({"first": hits[0], "last": hits[-1],
             "max_x": hits[np.argmax(x[hits])],
             "min_x": hits[np.argmin(x[hits])]}[pick])
    # Interpolate the crossing rather than snapping to the vertex before it. On a
    # 400-point curve the difference is sub-pixel, but the straight wet-region
    # segment of an isobar has only TWO points, and snapping puts its label on the
    # end of the line instead of where it was asked for.
    f = 0.0 if c[i + 1] == c[i] else (value - c[i]) / (c[i + 1] - c[i])
    xc = x[i] + f * (x[i + 1] - x[i])
    yc = y[i] + f * (y[i + 1] - y[i])
    (xlo, xhi), (ylo, yhi) = sorted(ax.get_xlim()), sorted(ax.get_ylim())
    if not (xlo <= xc <= xhi and ylo <= yc <= yhi):
        return False
    yc = yc * (1 + offset[1]) if ax.get_yscale() == "log" else yc + offset[1]
    ax.text(xc + offset[0], yc, text, fontsize=size, color=color,
            ha="center", va="center",
            rotation=angle_at(ax, x, y, i) if rotate else 0,
            rotation_mode="anchor", zorder=6,
            bbox=dict(fc="white", ec="none", pad=pad, alpha=0.85))
    return True


def label_end(ax, x, y, text, *, size=6.5, color="k", end="last", ha="left",
              pad=0.6):
    """Write `text` at the first or last point of the polyline inside the axes.

    The fallback when a curve leaves the frame before reaching any height worth
    labeling: put the value on the end of the line instead.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    (xlo, xhi), (ylo, yhi) = sorted(ax.get_xlim()), sorted(ax.get_ylim())
    inside = np.nonzero((x >= xlo) & (x <= xhi) & (y >= ylo) & (y <= yhi))[0]
    if not len(inside):
        return False
    i = int(inside[0] if end == "first" else inside[-1])
    ax.text(x[i], y[i], f" {text} ", fontsize=size, color=color, ha=ha,
            va="center", zorder=6,
            bbox=dict(fc="white", ec="none", pad=pad, alpha=0.85))
    return True


def fmt_P(P):
    """Format a pressure in MPa, dropping to kPa at the low end.

    The steam charts span five decades of pressure; `0.000611 MPa` on a label is
    unreadable where `0.611 kPa` is not.
    """
    return f"{P * 1e3:.3g} kPa" if P < 0.01 else f"{P:g} MPa"


__all__ = ["LW", "GRAY", "GRAY_MINOR", "GRID_MAJOR", "GRID_MINOR",
           "use_book_style", "chart_grid", "log_pressure_grid",
           "angle_at", "label_at", "label_end", "fmt_P"]
