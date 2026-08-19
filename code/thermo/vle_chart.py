"""vle_chart -- the drawing layer for Chapter 10's phase diagrams.

P-x-y, T-x-y and x-y diagrams, the tie lines and azeotrope marks that go on them, and
the McCabe-Thiele stage construction. The curves themselves come from `thermo.vle`
(`pxy`, `txy`) or from any solver with the same interface, so this module draws and
does no thermodynamics.

    from thermo.vle import Antoine, GammaPhi, pxy, txy
    from thermo.vle_chart import pxy_chart, txy_chart, xy_chart
    from thermo.charts import use_book_style
    from thermo.activity_models import VanLaar

    use_book_style()
    v = GammaPhi([Antoine(9.6830, 2842.2, -56.3209),
                  Antoine(9.3171, 2810.5, -51.2586)], VanLaar(1.15, 0.92))
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    pxy_chart(ax, *pxy(v, 348.15), species="ethyl acetate")

WHY IT IS A MODULE AND NOT CODE IN EACH NOTEBOOK. Chapter 10 draws this same picture
roughly thirteen times, and Section 10.3 draws it from a completely different physics
(the equation of state, not activity coefficients). The craft is identical either way --
which is the argument `charts.py` makes for its own existence -- and thirteen private
copies would be thirteen chances to disagree about tie-line weight, diagonal style and
where the region labels sit. The book is printed in black and white, so those choices
are the encoding, not decoration.

BLACK AND WHITE, AND THAT MEANS BLACK. Every line, marker and label this module draws is
pure black. Not dark gray -- black. Gray survives a laser printer and dies in offset
printing at small sizes, and the book is printed in black and white (author, 2026-08-13);
`charts.py` keeps a `GRAY` for the property charts' *grid rules*, which is a different
job, and nothing here uses it.

What separates one line from another is therefore weight, dash pattern and geometry:

    envelope (bubble + dew)   solid, LW["sat"]      the heaviest thing on the chart
    tie lines                 solid, 0.6            short horizontals inside the envelope
    x = y diagonal            dashed, 0.8           the only dashed line on an x-y chart
    operating lines           solid, 0.9            straight, against a curved equilibrium
    stage steps               solid, 0.7            a staircase, against everything else
    measured data             open markers          never filled, never colored

The bubble and dew curves are both solid, as the 5e art draws them: on a closed envelope
they meet at both pure ends, so position plus the `Liquid` / `Vapor` region labels
identify them without needing a dash or a tint.

GAPS ARE REAL. `pxy` and `txy` return NaN where the model has no solution, and matplotlib
renders NaN as a break in the line. That is deliberate: Figure 10.3-13 prints exactly such
a gap for acetone/water on the van der Waals one-fluid rule with k12 = 0. Do not filter
the NaNs out before plotting; the gap is the result.
"""
from __future__ import annotations

import numpy as np

from .charts import LW, chart_grid

__all__ = ["pxy_chart", "txy_chart", "xy_chart", "tie_lines", "mark_azeotrope",
           "mccabe_thiele", "total_reflux_steps", "BAR", "CELSIUS"]

BAR = ("bar", 1e-5, 0.0)          # (label, scale, offset) applied as v*scale + offset
CELSIUS = ("$^\\circ$C", 1.0, -273.15)
KELVIN = ("K", 1.0, 0.0)
PASCAL = ("Pa", 1.0, 0.0)


def _apply(unit, v):
    _, scale, offset = unit
    return np.asarray(v, dtype=float) * scale + offset


def _composition_axis(ax, species, phase_note=""):
    """Label the abscissa the way the chapter does -- by species, not by 'x'."""
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(f"mole fraction {species}{phase_note}")


# ---------------------------------------------------------------------------
# the two envelope diagrams
# ---------------------------------------------------------------------------
def _roomiest(comp, curve, bound, *, margin=0.12, n=81):
    """Where a single-phase region is widest, as (abscissa, ordinate) to label.

    A region here is bounded by one curve and one edge of the axes. Scanning the
    abscissa for the largest vertical gap between the two and taking its midpoint puts
    the label in the open, which a fixed offset from the axis limits cannot do: on a
    P-x-y at 50 C the bubble line runs from 0.19 to 1.56 bar, so "Liquid" placed a
    fixed 6 % below the top of the frame floats far above the curve at mid-composition
    while "Vapor" lands directly on the dew line. Both were wrong in the same figure,
    in opposite directions, which is what a geometry-blind rule buys.

    `margin` keeps the search away from the pure ends, where the region closes to
    nothing and a centered label would overhang the frame.
    """
    comp, curve = np.asarray(comp, float), np.asarray(curve, float)
    ok = np.isfinite(comp) & np.isfinite(curve)
    if ok.sum() < 2:
        return None
    comp, curve = comp[ok], curve[ok]
    order = np.argsort(comp)
    comp, curve = comp[order], curve[order]
    zs = np.linspace(margin, 1.0 - margin, n)
    edge = np.interp(zs, comp, curve)
    gap = np.abs(bound - edge)
    i = int(np.argmax(gap))
    if gap[i] <= 0:
        return None
    return float(zs[i]), float(0.5 * (bound + edge[i]))


def _widest_gap(x1, y1, curve, *, margin=0.10, n=97):
    """The abscissa at which the bubble and dew CURVES are furthest apart.

    Where to write "liquid" and "vapor" on the lines themselves, which is what an
    azeotropic envelope needs: there the single-phase regions are thin bands with no
    room for a region label (see `_envelope`), but the two curves are still widely
    separated somewhere between the pure end and the azeotrope, and that is the place
    to name them.

    Returns `(abscissa, v_bubble, v_dew)`, or None if the two never separate.
    """
    x1, y1, curve = (np.asarray(a, float) for a in (x1, y1, curve))
    ok = np.isfinite(x1) & np.isfinite(y1) & np.isfinite(curve)
    if ok.sum() < 3:
        return None
    xb, yb, cb = x1[ok], y1[ok], curve[ok]
    zs = np.linspace(margin, 1.0 - margin, n)

    def at(comp):                     # the curve's ordinate as a function of comp
        order = np.argsort(comp)
        return np.interp(zs, comp[order], cb[order])

    bub, dew = at(xb), at(yb)
    i = int(np.argmax(np.abs(bub - dew)))
    if not np.isfinite(bub[i] - dew[i]) or bub[i] == dew[i]:
        return None
    return float(zs[i]), float(bub[i]), float(dew[i])


def _outermost(ax, abscissa, start, *, toward, window=0.06):
    """How far out anything already drawn on `ax` reaches at this abscissa.

    Returns the ordinate of the artist nearest `toward` (the frame edge the label
    will be placed against), or `start` if nothing drawn goes past it. Every line
    and marker series already on the axes is considered, which is the point: a
    chart helper knows the curves it drew itself and nothing about what the
    notebook added afterwards.

    `window` is a tolerance in the abscissa, so a scatter of measured points near
    the label's composition counts even though none sits exactly on it.
    """
    best = float(start)
    up = toward > start
    for line in ax.get_lines():
        xs, ys = line.get_xdata(), line.get_ydata()
        if len(xs) == 0:
            continue
        xs, ys = np.asarray(xs, float), np.asarray(ys, float)
        near = np.abs(xs - abscissa) <= window
        vals = ys[near & np.isfinite(ys)]
        if vals.size == 0:
            continue
        edge = vals.max() if up else vals.min()
        if (edge > best) if up else (edge < best):
            best = float(edge)
    return best


def _envelope(ax, x1, y1, v, *, unit, species, ylabel, liquid_side,
              regions=True, region_labels=("Liquid", "Vapor"), region_pos=None,
              region_size=7, data=None, grid=None, lw=None, dew_style="-",
              curve_labels=None, curve_label_pos=None, curve_label_size=6.5):
    """Shared body of `pxy_chart` and `txy_chart`.

    The two diagrams are the same picture with the ordinate swapped and the phases on
    opposite sides -- `liquid_side` is 'upper' for P-x-y and 'lower' for T-x-y. That
    inversion is the point of Illustration 10.2-2's closing Comment, and putting it in
    one parameter rather than two functions is what keeps the two charts consistent.

    The region labels are placed where each region is actually widest (see
    `_roomiest`), which is right often enough to be the default and never a reason to
    accept a bad one: pass `region_pos=((xl, vl), (xv, vv))` in data units to put them
    exactly where you want, `region_labels` to rename them, or `regions=False` to drop
    them. A figure whose labels are worth arguing about is a figure worth tuning.

    On an AZEOTROPIC system, expect to use one of those. A maximum-pressure azeotrope
    puts the top of the frame just above the azeotropic point, so the liquid region is
    a thin band everywhere and no automatic placement will find real room in it. The
    5e's own azeotropic P-x-y (Illustration 10.2-2, `c10uf004`) carries no region
    labels at all, which is the honest answer for that shape: pass `regions=False`.

    `dew_style` dashes the dew line when the envelope is too narrow to label. Both
    curves solid is the default because on a normal envelope position identifies them;
    when the two lines nearly touch -- benzene / cyclohexane in Illustration 10.2-1
    spans 0.98 to 1.02 bar -- there is no room for a region label and no way to tell
    the curves apart, and the 5e itself dashes the vapor line there.

    `curve_labels=("liquid", "vapor")` names the two LINES instead of the two regions,
    which is the answer whenever `regions=False` was forced by the shape. It writes
    each name beside its own curve at the abscissa where the two are furthest apart,
    offset away from the other line, so nothing lands in the thin band or on top of
    the opposite curve. Author's request on Illustration 10.2-2, 2026-08-15: an
    unlabeled envelope leaves the reader to work out which line is which, and on an
    azeotropic diagram the usual cue -- liquid on top for P-x-y, on the bottom for
    T-x-y -- *reverses* on the far side of the azeotrope. `curve_label_pos` overrides
    placement with ((x, v), (x, v)) in display units, as `region_pos` does.
    """
    vv = _apply(unit, v)
    lw = LW["sat"] if lw is None else lw
    ax.plot(x1, vv, "-", color="k", lw=lw)              # bubble line: ordinate vs x
    ax.plot(y1, vv, dew_style, color="k", lw=lw)        # dew line:    ordinate vs y

    if data is not None:
        for comp, val, mk in _data_series(data, unit):
            ax.plot(comp, val, mk, mfc="none", color="k", ms=3.5, lw=0)

    if regions:
        ax.autoscale_view()
        lo, hi = ax.get_ylim()
        # the liquid is bounded by the BUBBLE curve, the vapor by the DEW curve, each
        # opening toward its own edge of the frame
        if liquid_side == "upper":
            liq, vap = _roomiest(x1, vv, hi), _roomiest(y1, vv, lo)
        else:
            liq, vap = _roomiest(x1, vv, lo), _roomiest(y1, vv, hi)
        if region_pos is not None:
            liq, vap = region_pos
        for pos, text in ((liq, region_labels[0]), (vap, region_labels[1])):
            if pos is not None and text:
                ax.text(pos[0], pos[1], text, fontsize=region_size,
                        ha="center", va="center")
        ax.set_ylim(lo, hi)

    if curve_labels:
        ax.autoscale_view()
        lo, hi = ax.get_ylim()
        spots = curve_label_pos
        if spots is None:
            found = _widest_gap(x1, y1, vv)
            if found is not None:
                a, v_bub, v_dew = found
                # Each label goes MIDWAY between its own curve and the frame edge that
                # curve opens toward -- the same idea as `_roomiest`, anchored at the
                # abscissa where the two curves are furthest apart. A fixed fractional
                # offset does not work: 3.5 % of the ordinate range put both labels
                # straight on top of their own lines on this figure.
                #
                # But the bubble and dew curves are not necessarily all that is on
                # the axes. Figure 10.2-8 adds a second model's curves and 26 data
                # points, and a label placed against the envelope alone landed on top
                # of them (author, 2026-08-15). So the outermost thing ALREADY DRAWN at
                # this abscissa is what the label has to clear, not just this envelope.
                l_edge, v_edge = (hi, lo) if liquid_side == "upper" else (lo, hi)
                v_bub = _outermost(ax, a, v_bub, toward=l_edge)
                v_dew = _outermost(ax, a, v_dew, toward=v_edge)
                spots = ((a, 0.5 * (v_bub + l_edge)), (a, 0.5 * (v_dew + v_edge)))
        if spots is not None:
            for pos, text in zip(spots, curve_labels):
                if pos is not None and text:
                    ax.text(pos[0], pos[1], text, fontsize=curve_label_size,
                            ha="center", va="center")
        ax.set_ylim(lo, hi)

    _composition_axis(ax, species)
    ax.set_ylabel(ylabel)
    if grid:
        chart_grid(ax, *grid)
    return ax


_MARKERS = ["o", "s", "^", "v", "D"]


def _data_series(data, unit):
    """Normalize `data` into a list of (composition, ordinate, marker) triples.

    Accepts one series as `(comp, value)` or `(comp, value, marker)`, or a list of
    those. An experimental VLE data set is normally two series against one ordinate --
    the measured liquid and the measured vapor at each pressure -- so overlaying both
    in one call is the common case, not the exotic one.
    """
    series = data if isinstance(data[0], (tuple, list)) else [data]
    out = []
    for i, s in enumerate(series):
        marker = s[2] if len(s) > 2 else _MARKERS[i % len(_MARKERS)]
        out.append((np.asarray(s[0], dtype=float), _apply(unit, s[1]), marker))
    return out


def pxy_chart(ax, x1, y1, P, *, species="species 1", unit=BAR, regions=True,
              region_labels=("Liquid", "Vapor"), region_pos=None, region_size=7,
              data=None, grid=None, lw=None, dew_style="-",
              curve_labels=None, curve_label_pos=None, curve_label_size=6.5):
    """A constant-temperature P-x-y diagram from a `pxy` sweep.

    `x1, y1, P` are exactly what `thermo.vle.pxy` returns. The liquid is the upper
    region: raising the pressure on a vapor condenses it.

    `region_pos=((x, P), (x, P))` overrides the automatic label placement, in the
    ordinate's *display* units (bar by default, not Pa).
    """
    return _envelope(ax, x1, y1, P, unit=unit, species=species,
                     ylabel=f"pressure ({unit[0]})", liquid_side="upper",
                     regions=regions, region_labels=region_labels,
                     region_pos=region_pos, region_size=region_size,
                     data=data, grid=grid, lw=lw, dew_style=dew_style,
                     curve_labels=curve_labels, curve_label_pos=curve_label_pos,
                     curve_label_size=curve_label_size)


def txy_chart(ax, x1, y1, T, *, species="species 1", unit=KELVIN, regions=True,
              region_labels=("Liquid", "Vapor"), region_pos=None, region_size=7,
              data=None, grid=None, lw=None, dew_style="-",
              curve_labels=None, curve_label_pos=None, curve_label_size=6.5):
    """A constant-pressure T-x-y diagram from a `txy` sweep.

    The liquid is the *lower* region here -- the inversion relative to `pxy_chart`,
    and the thing Illustration 10.2-2's Comment asks the reader to understand. Pass
    `unit=CELSIUS` for a Celsius ordinate.

    `region_pos=((x, T), (x, T))` overrides placement, in the ordinate's display units
    (so degrees Celsius when `unit=CELSIUS`).
    """
    return _envelope(ax, x1, y1, T, unit=unit, species=species,
                     ylabel=f"temperature ({unit[0]})", liquid_side="lower",
                     regions=regions, region_labels=region_labels,
                     region_pos=region_pos, region_size=region_size,
                     data=data, grid=grid, lw=lw, dew_style=dew_style,
                     curve_labels=curve_labels, curve_label_pos=curve_label_pos,
                     curve_label_size=curve_label_size)


# ---------------------------------------------------------------------------
# the x-y diagram
# ---------------------------------------------------------------------------
def xy_chart(ax, x1, y1, *, species="species 1", diagonal=True, data=None,
             grid=None, lw=None, label_diagonal=True, diagonal_label_at=0.62):
    """An x-y diagram, with the x = y diagonal the chapter always draws beside it.

    The diagonal is not decoration. Section 10.1 says the gap between the equilibrium
    line and x = y "is an indication of how easy or difficult it will be to separate
    the components by distillation," and an azeotrope is precisely where the two
    touch -- so a bare x-y curve omits the one comparison the diagram exists to make.

    `diagonal_label_at` is where along the diagonal the "$x = y$" label sits, in mole
    fraction. The default puts it in the upper half, clear of the curve on a system
    whose equilibrium line bows above the diagonal all the way across. It has to move
    on an AZEOTROPIC system, where the curve touches the diagonal and the label lands
    on the crossing: pass the composition of an open stretch instead.
    """
    lw = LW["sat"] if lw is None else lw
    if diagonal:
        ax.plot([0, 1], [0, 1], "--", color="k", lw=0.8)
        if label_diagonal:
            s = float(diagonal_label_at)
            ax.text(s, s - 0.04, "$x = y$", fontsize=6.5, color="k",
                    ha="center", va="center", rotation=45, rotation_mode="anchor")
    ax.plot(x1, y1, "-", color="k", lw=lw)
    if data is not None:
        d = np.asarray(data, dtype=float)
        ax.plot(d[0], d[1], "o", mfc="none", color="k", ms=3.5, lw=0)
    _composition_axis(ax, species, " in liquid")
    ax.set_ylabel(f"mole fraction {species} in vapor")
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal", adjustable="box")
    if grid:
        chart_grid(ax, *grid)
    return ax


# ---------------------------------------------------------------------------
# annotations
# ---------------------------------------------------------------------------
def _crossing(comp, ordinate, target):
    """Composition at which `ordinate` first crosses `target`, by interpolation.

    Deliberately NOT `np.interp`, which silently requires an increasing abscissa. On
    a T-x-y the temperature *falls* as the more volatile species is added, so the
    ordinate runs backwards and `np.interp` returns nonsense rather than an error --
    which is exactly how this went unnoticed until a T-x-y drew no tie lines at all.
    Crossing detection has no such requirement and also survives the non-monotonic
    ordinate of an azeotropic system.
    """
    comp, ordinate = np.asarray(comp, float), np.asarray(ordinate, float)
    ok = np.isfinite(comp) & np.isfinite(ordinate)
    comp, ordinate = comp[ok], ordinate[ok]
    d = ordinate - target
    hits = np.nonzero(d[:-1] * d[1:] <= 0)[0]
    if not len(hits):
        return None
    i = int(hits[0])
    span = ordinate[i + 1] - ordinate[i]
    f = 0.0 if span == 0 else (target - ordinate[i]) / span
    return float(comp[i] + f * (comp[i + 1] - comp[i]))


def tie_lines(ax, x1, y1, v, values, *, unit=BAR, color="k", lw=0.6):
    """Draw horizontal tie lines across an envelope at the given ordinate values.

    A tie line joins the liquid and vapor in equilibrium, so on a P-x-y it is a line
    of constant pressure and on a T-x-y one of constant temperature -- horizontal in
    both cases, which is the observation Figures 10.1-3 and 10.1-4 are drawn to make.

    `values` are in the ordinate's own SI units (Pa, K); `unit` converts them for the
    axes, exactly as it does for the envelope itself, so the caller never has to know
    which units the axes ended up in. A value outside the envelope is skipped.
    """
    vv = _apply(unit, v)
    for val in np.atleast_1d(values):
        target = float(_apply(unit, val))
        xa = _crossing(x1, vv, target)
        ya = _crossing(y1, vv, target)
        if xa is None or ya is None:
            continue
        ax.plot([xa, ya], [target, target], "-", color=color, lw=lw, zorder=1)
    return ax


def mark_azeotrope(ax, x_az, v_az, *, unit=BAR, text="azeotrope", size=6.5,
                   offset=(0.03, 0.0)):
    """Mark an azeotrope on an envelope or x-y diagram.

    Worth marking explicitly rather than leaving the reader to spot the maximum:
    Illustration 10.2-2 prints its azeotropic composition and pressure to three and
    four figures, and a printed number the figure does not point at is a number
    nobody checks. (In the 5e that number was wrong; see the chapter notes.)
    """
    vv = float(_apply(unit, v_az))
    ax.plot([x_az], [vv], "o", color="k", ms=3, zorder=5)
    ax.annotate(text, xy=(x_az, vv), xytext=(x_az + offset[0], vv + offset[1]),
                fontsize=size, ha="left", va="bottom")
    return ax


# ---------------------------------------------------------------------------
# McCabe-Thiele
# ---------------------------------------------------------------------------
def total_reflux_steps(ax, x1, y1, x_start, x_target=None, *, direction="down",
                       max_stages=80, draw=True, lw=0.7, clip=True, pinch_tol=1e-4):
    """Step off stages at TOTAL REFLUX, where the operating line is the diagonal.

    Returns `(stages, staircase, pinched)`: the stage count, the staircase as
    (x0, y0, x1, y1) segments, and whether the steps stalled before reaching
    `x_target`. Draws the staircase when `draw`.

    At total reflux the vapor rising from a tray equals the liquid falling from the
    one above, `y_{n+1} = x_n`, so each step is one application of the equilibrium
    curve or of its inverse -- and WHICH ONE IS THE DIRECTION OF TRAVEL:

    * `direction="down"` -- toward the reboiler. `x_{n+1} = f^-1(x_n)`, which moves
      AWAY from an azeotrope, toward whichever pure component lies on the starting
      composition's side of it. This is the branch that reaches a bottoms spec.
    * `direction="up"` -- toward the condenser. `x_{n+1} = f(x_n)`, which moves
      TOWARD an azeotrope. Past a certain point it cannot do anything else, which is
      what a pinch is: the steps shrink geometrically and `pinched` comes back True.

    **Both branches are needed to draw an azeotropic column, and asking for the
    wrong one silently walks the other way.** Before 2026-08-14 this function applied
    only the inverse, so a caller who wanted the rectifying section got the stripping
    section instead -- steps marching to the far pure end rather than piling up
    against the azeotrope. The two figures of §10.2's pressure-swing example were
    drawn that way, and the author caught it: the product and feed lines did not
    line up with the staircase, because the staircase was not the section the
    caption described. There is no way to infer the branch from `x_target` alone.

    **`x_start` must not be the azeotrope itself.** It is a fixed point of both
    maps, so a descent that begins there crawls out of the pinch over a dozen
    numerically meaningless steps before it separates anything, and the stage count
    then reports the arithmetic rather than the chemistry. Start from a real,
    specified composition -- the feed -- and step outward.

    `clip` truncates the last step at `x_target` so the staircase ends exactly on the
    specification line instead of overshooting it. The stage count is unaffected (the
    overshooting stage is still a stage, as it is in any hand construction); this
    only stops a figure from showing a staircase that misses its own product line.

    WHY THIS IS SEPARATE FROM `mccabe_thiele`. That function assumes the plotted
    species is the more volatile one, so the equilibrium curve lies above the diagonal
    and the steps march left-and-up from the distillate. On the far side of an
    azeotrope that assumption inverts: there y < x, the plotted species is the HEAVY
    one, and the steps march right-and-down toward the reboiler.

    The geometric invariant, which is what the drawing has to satisfy and what these
    doctests guard: **every corner of the staircase lies on the equilibrium curve**, in
    both directions, so no segment ever crosses the y = x diagonal.

        >>> import numpy as np
        >>> x = np.linspace(0, 1, 401)
        >>> y = 2.5 * x / (1 + 1.5 * x)            # constant relative volatility 2.5
        >>> f = lambda v: float(np.interp(v, x, y))
        >>> for d, target in (("up", 0.95), ("down", 0.05)):
        ...     n, segs, _ = total_reflux_steps(None, x, y, 0.5, target,
        ...                                     direction=d, draw=False, clip=False)
        ...     corners = max(abs(f(s[2]) - s[3]) for s in segs[::2])
        ...     crossed = [s for s in segs if (s[1] - s[0]) * (s[3] - s[2]) < -1e-12]
        ...     print(d, n, round(corners, 6), len(crossed))
        up 4 0.0 0
        down 4 0.0 0

    For a minimum-boiling azeotrope the azeotrope is the most volatile condition
    the mixture can reach, so the DISTILLATE is always driven toward the azeotropic
    composition and the BOTTOMS toward whichever pure component lies on the feed's
    side of it. Getting that backwards is erratum E17 in the 5e; see
    `revision_notes/c10.md`.
    """
    if direction not in ("up", "down"):
        raise ValueError("direction must be 'up' (toward the condenser) or "
                         "'down' (toward the reboiler)")
    xs, ys = np.asarray(x1, float), np.asarray(y1, float)
    ok = np.isfinite(xs) & np.isfinite(ys)
    xs, ys = xs[ok], ys[ok]
    order = np.argsort(xs)
    X, Y = xs[order], ys[order]
    # Descending inverts the curve: given the vapor leaving the tray above, find the
    # liquid on this one. np.interp needs an increasing abscissa, and y is monotone in
    # x across a single azeotrope, so sorting on Y gives a well-defined inverse.
    o2 = np.argsort(Y)
    Yi, Xi = Y[o2], X[o2]
    step = ((lambda v: float(np.interp(v, Yi, Xi))) if direction == "down"
            else (lambda v: float(np.interp(v, X, Y))))

    segs, x, n, pinched = [], float(x_start), 0, False
    for _ in range(max_stages):
        xn = step(x)
        if abs(xn - x) < pinch_tol:
            pinched = True
            break
        reached = x_target is not None and (
            (x_target > x_start and xn >= x_target) or
            (x_target < x_start and xn <= x_target))
        if reached and clip:
            xn = float(x_target)
        # THE CORNER OF EACH STEP MUST LIE ON THE EQUILIBRIUM CURVE, and which
        # corner that is depends on the direction. Descending, x_{n+1} = f^-1(x_n), so
        # the curve point is (x_{n+1}, x_n) -- go across, then down. Ascending,
        # x_{n+1} = f(x_n), so the curve point is (x_n, x_{n+1}) -- go UP first, then
        # across. Using the descending order for both draws the staircase reflected in
        # the diagonal: its corners land on the far side, the steps appear to cross
        # y = x, and they touch the equilibrium curve nowhere. The stage count is
        # unaffected (the recursion above is the same either way), which is exactly why
        # this survived a numerical check and had to be caught by eye -- author,
        # 2026-08-14, on the first printing of Figures 10.2-2 and 10.2-3.
        if direction == "down":
            segs.append((x, x, xn, x))      # across to the curve at (x_{n+1}, x_n)
            segs.append((xn, x, xn, xn))    # then down to the diagonal
        else:
            segs.append((x, x, x, xn))      # up to the curve at (x_n, x_{n+1})
            segs.append((x, xn, xn, xn))    # then across to the diagonal
        n += 1
        if reached:
            break
        x = xn
    else:
        pinched = x_target is not None
    if draw:
        for a, b, c, d in segs:
            ax.plot([a, c], [b, d], "-", color="k", lw=lw)
    return n, segs, pinched


def mccabe_thiele(ax, x1, y1, *, xD, xB, reflux, feed_ratio=None, xF=None,
                  species="species 1", max_stages=60, draw=True, lw=0.7):
    """The graphical stage-to-stage construction of Illustration 10.1-7.

    `reflux` is the book's `q = L/D`, the reflux ratio. The rectifying operating line
    runs from (xD, xD) with slope q/(q+1); when `feed_ratio` (F/D) and `xF` are given,
    the stripping line runs from (xB, xB) with slope (q + F/D)/(q + 1), as the
    illustration sets it up.

    Returns the number of stages stepped off, and draws the construction when `draw`.
    Returns `None` for the count if the steps stall -- which is not a failure but the
    answer: it is what a pinch point *is*, and Figure 10.1-11 exists to show it. The
    caller can sweep `reflux` downward and find the minimum as the value where this
    stops returning a number, rather than reading it off a plot.
    """
    x1 = np.asarray(x1, dtype=float)
    y1 = np.asarray(y1, dtype=float)
    ok = np.isfinite(x1) & np.isfinite(y1)
    xe, ye = x1[ok], y1[ok]
    order = np.argsort(xe)
    xe, ye = xe[order], ye[order]

    s_rect = reflux / (reflux + 1.0)
    b_rect = xD / (reflux + 1.0)
    rect = lambda x: s_rect * x + b_rect

    if feed_ratio is not None:
        s_strip = (reflux + feed_ratio) / (reflux + 1.0)
        b_strip = xB - s_strip * xB
        strip = lambda x: s_strip * x + b_strip
        # the operating lines cross at the feed composition; below it, strip governs
        x_cross = xF if xF is not None else (b_strip - b_rect) / (s_rect - s_strip)
        op = lambda x: rect(x) if x >= x_cross else strip(x)
    else:
        op, x_cross = rect, xB

    if draw:
        xr = np.linspace(x_cross, xD, 50)
        ax.plot(xr, rect(xr), "-", color="k", lw=lw + 0.2)
        if feed_ratio is not None:
            xs = np.linspace(xB, x_cross, 50)
            ax.plot(xs, strip(xs), "-", color="k", lw=lw + 0.2)

    # step: from (xD, xD) go horizontally to the equilibrium curve, then down to the
    # operating line, and repeat until the liquid is weaker than the bottoms spec.
    x, y = xD, xD
    stages = 0
    for _ in range(max_stages):
        x_eq = float(np.interp(y, ye, xe))          # horizontal to equilibrium
        if draw:
            ax.plot([x, x_eq], [y, y], "-", color="k", lw=lw)
        stages += 1
        if x_eq <= xB:
            return stages
        y_new = op(x_eq)                             # vertical to operating line
        if draw:
            ax.plot([x_eq, x_eq], [y, y_new], "-", color="k", lw=lw)
        if abs(x_eq - x) < 1e-6 and abs(y_new - y) < 1e-6:
            return None                              # pinched: steps stopped advancing
        x, y = x_eq, y_new
    return None
