"""ternary — triangular composition diagrams, the drawing half of Sec. 11.2.

Five of chapter 11's figures are triangular diagrams (11.2-7 through 11.2-11) and
nothing else in the book draws one: `charts.py` rules rectangular chart paper and
`vle_chart.py` draws binary phase envelopes. This module is the missing piece. It
is draftsmanship only -- the thermodynamics is in `thermo.lle` (`lle_flash`,
`tie_line_split`, `mix_streams`), and nothing here solves anything.

    from thermo.ternary import ternary_axes, plot, tie_line, point, to_xy

**WHY THIS IS BESPOKE AND NOT `mpltern`.** The author's practice files
(`code/test code/ternary.ipynb`) use `mpltern`, which is a good library, and it was
still not the right choice here:

  - It is a **new student-facing dependency**, and `code/pyproject.toml` is
    deliberately a short list of floors with no lock file, so that a student in
    2032 still gets an environment that resolves. A triangle is about a hundred lines
    of geometry, which is not worth a dependency.
  - The house rules for printed art are specific and this module has to obey them:
    **pure black ink**, no tints ([[black-not-gray]]); **Computer Modern**, set by
    `charts.use_book_style` ([[figure-typography-cm]]); and **nothing below 7 pt**
    ([[no-figure-text-below-7pt]]). Bending a third-party projection to all three is
    more work than owning the geometry.
  - Figure 11.2-7 is not a data plot at all. It is the schematic that teaches a
    reader **how to read** a triangular diagram, so it needs annotation that no
    plotting library provides.

THE COORDINATE CONVENTION, stated once because everything here depends on it. The
triangle is equilateral with the LEFT corner at (0, 0), the RIGHT corner at (1, 0)
and the APEX at (1/2, sqrt(3)/2). A composition is given as `(top, left, right)` --
the fraction at each named corner, summing to one:

              top
              /\\
             /  \\        x = right + top/2
            /    \\       y = top * sqrt(3)/2
           /______\\
        left        right

so `(1, 0, 0)` is the apex, `(0, 1, 0)` the left corner and `(0, 0, 1)` the right.
Each corner is the pure component named there.

**Two different edges belong to each species, and confusing them is the easy
mistake.** A species' fraction is proportional to the perpendicular distance from
the edge **opposite** its corner -- that is the geometry, and it is why the three
fractions sum to one. But the printed **ruler** for that species is drawn on an
**adjacent** edge, the one ending at its own corner, so that the scale runs toward
the thing it measures. Fig. 11.2-7 is the schematic that says this, `_edge_points`
implements the ruler placement, and `read_construction` draws the reading.

**Fraction-agnostic, like `thermo.lle`.** Section 11.2's ternary illustrations are
worked in WEIGHT fractions and weights; the mathematics is identical either way, so
nothing here assumes moles. Pass `percent=True` to label the axes 0-100.
"""
from __future__ import annotations

import numpy as np

__all__ = ["to_xy", "from_xy", "ternary_axes", "plot", "scatter", "tie_line",
           "point", "text", "region_label", "lever_arm", "read_construction",
           "check_labels"]

_H = np.sqrt(3.0) / 2.0          # height of a unit equilateral triangle

# Corner positions, in the (top, left, right) order used throughout.
APEX, LEFT, RIGHT = np.array([0.5, _H]), np.array([0.0, 0.0]), np.array([1.0, 0.0])


# --- geometry --------------------------------------------------------------

def _norm(comp):
    """(top, left, right) -> a normalized float array, with the sum checked.

    Compositions read off a printed diagram rarely sum to exactly one -- Sec. 11.2's
    own illustrations are the reason `tie_line_split` solves by least squares -- so a
    small shortfall is renormalized silently and a large one is an error, because at
    that size it is a transcription mistake rather than rounding.
    """
    c = np.asarray(comp, dtype=float)
    if c.shape[-1] != 3:
        raise ValueError("a ternary composition needs three components "
                         "(top, left, right)")
    s = c.sum(axis=-1, keepdims=True)
    if np.any(np.abs(s - 1.0) > 0.02):
        bad = float(np.max(np.abs(s - 1.0)))
        raise ValueError(f"composition sums differ from 1 by up to {bad:.3f}; "
                         "that is too much to be rounding -- check the transcription")
    return c / s


def to_xy(comp):
    """(top, left, right) -> (x, y) on the unit triangle. Accepts an (N, 3) array."""
    c = _norm(comp)
    top, right = c[..., 0], c[..., 2]
    return np.stack([right + top / 2.0, top * _H], axis=-1)


def from_xy(xy):
    """(x, y) -> (top, left, right). The inverse of `to_xy`, for reading a diagram.

    Used when a composition has been picked off a figure by eye -- which is how every
    ternary answer in Sec. 11.2 is obtained, and why Illustration 11.2-9's tie line
    sums to 86 % instead of 100.
    """
    p = np.asarray(xy, dtype=float)
    top = p[..., 1] / _H
    right = p[..., 0] - top / 2.0
    return np.stack([top, 1.0 - top - right, right], axis=-1)


def _unpack(*args):
    """Accept either one (N, 3) array or three separate sequences."""
    if len(args) == 1:
        return np.atleast_2d(np.asarray(args[0], dtype=float))
    if len(args) == 3:
        return np.stack([np.atleast_1d(np.asarray(a, dtype=float)) for a in args],
                        axis=-1)
    raise TypeError("give one (N, 3) array of compositions, or three sequences")


# --- the axes --------------------------------------------------------------

def _edge_points(f):
    """The three tick anchors at fraction `f`, one per species, on the BOOK's edges.

    **Which edge carries which species is not arbitrary, and a plausible-looking
    choice is wrong.** The first draft of this module put the apex species' scale on
    the left edge; the 5e's own Figs. 11.2-7 and 11.2-8 put it on the right. The
    rule, read off those two figures: **each species is scaled on the edge that ENDS
    AT ITS OWN CORNER, going clockwise** --

        apex species  -> RIGHT edge,  increasing upward   (x_B in Fig. 11.2-7)
        left species  -> LEFT edge,   increasing downward (x_C)
        right species -> BOTTOM edge, increasing rightward (x_A)

    so every scale runs *toward* the corner of the thing it measures, which is what
    makes the diagram readable without a legend.
    """
    return (to_xy([f, 0.0, 1 - f]),        # top species, on the right edge
            to_xy([1 - f, f, 0.0]),        # left species, on the left edge
            to_xy([0.0, 1 - f, f]))        # right species, on the bottom edge


def ternary_axes(ax, *, top, left, right, ticks=0.2, minors=2, grid=True,
                 percent=False, symbol=None, corner_values=True,
                 label_size=9, tick_size=7, lw=0.9, pad=0.040, gray_grid=False):
    """Draw the triangle, its scales, its labels and (optionally) its grid.

    `top`, `left`, `right` are the species names at the three corners. Pass
    `symbol="w"` or `"x"` to set each edge's own axis label the way the book does --
    `w_MIK` along the left edge of Fig. 11.2-8, and so on -- built from the corner
    names. `symbol=None` leaves the edges unlabeled, which is right only when the
    corner names alone are unambiguous.

    The grid gets the two-tier major/minor treatment `charts.chart_grid` uses,
    because a triangular diagram is read by interpolating between rulings exactly as
    chart paper is -- but in **pure black hairlines, not gray**.

    **An earlier version of this docstring claimed a triangular grid was "the one
    place tint is allowed." It is not, and `tools/check_print_art.py` said so**: it
    rejected the first Figure 11.2-8 on `0.7 G` and `0.87 G` in the content stream.
    That gate carries an allowlist of author-granted gray waivers and this was not on
    it. [[black-not-gray]] is the rule -- separate lines by **weight, dash and
    geometry, never by tint** -- and it is satisfied here the way the rule says: the
    grid is 0.20 pt black where the binodal is 1.1 pt black, a factor of five in
    weight. Pass `gray_grid=True` for a screen figure that will never be staged.

    Nothing else in the book stages gray art, and `charts.GRID_MAJOR`/`GRID_MINOR`
    are unused by any staged figure -- this module was about to be the first, which is
    what makes the gate worth having.

    **Label sizes are floored at 7 pt.** A long species name on a 3.3 in figure is
    the failure recorded in [[no-figure-text-below-7pt]]: the fix is a SHORTER NAME,
    not smaller type. "Methyl isobutyl ketone" does not fit; "MIK" does, and the
    caption carries the expansion -- which is exactly what the 5e's Fig. 11.2-8 does.
    """
    if min(label_size, tick_size) < 7:
        raise ValueError("nothing in a printed figure goes below 7 pt; shorten the "
                         "label instead [[no-figure-text-below-7pt]]")
    if gray_grid:
        from .charts import GRID_MAJOR, GRID_MINOR
    else:
        GRID_MAJOR = dict(color="k", lw=0.30)
        GRID_MINOR = dict(color="k", lw=0.18)

    scale = 100.0 if percent else 1.0
    fmt = (lambda v: f"{v * scale:.0f}") if percent else (lambda v: f"{v:g}")

    # --- grid, drawn first so the triangle and the data sit on top of it -------
    if grid and ticks:
        step = ticks / max(int(minors), 1)
        n = int(round(1.0 / step))
        for i in range(1, n):
            f = i * step
            major = abs(f / ticks - round(f / ticks)) < 1e-9
            style = GRID_MAJOR if major else GRID_MINOR
            for p, q in ((to_xy([f, 1 - f, 0.0]), to_xy([f, 0.0, 1 - f])),
                         (to_xy([0.0, f, 1 - f]), to_xy([1 - f, f, 0.0])),
                         (to_xy([0.0, 1 - f, f]), to_xy([1 - f, 0.0, f]))):
                ax.plot([p[0], q[0]], [p[1], q[1]], zorder=0, **style)

    corners = np.array([APEX, LEFT, RIGHT, APEX])
    ax.plot(corners[:, 0], corners[:, 1], "-", color="k", lw=lw, zorder=3,
            solid_joinstyle="miter")

    # --- the three scales -----------------------------------------------------
    if ticks:
        for f in np.arange(ticks, 1.0 - 1e-9, ticks):
            r, l, b = _edge_points(f)
            r, l, b = r.ravel(), l.ravel(), b.ravel()
            ax.text(r[0] + pad * 0.55, r[1], fmt(f), size=tick_size, color="k",
                    ha="left", va="center")
            ax.text(l[0] - pad * 0.55, l[1], fmt(f), size=tick_size, color="k",
                    ha="right", va="center")
            ax.text(b[0], b[1] - pad * 0.50, fmt(f), size=tick_size, color="k",
                    ha="center", va="top")

    # --- corner names outside, scale ends just inside -------------------------
    # The `1.0` sits INSIDE the triangle, on the inward bisector, rather than
    # beside the species name. Outside it collides with the name at every corner --
    # "100W" instead of "100  W" -- and the collision is worse the longer the name.
    one = fmt(1.0)
    ax.text(APEX[0], APEX[1] + pad, top, size=label_size, color="k",
            ha="center", va="bottom")
    ax.text(LEFT[0] - pad * 0.55, LEFT[1] - pad * 0.55, left, size=label_size,
            color="k", ha="right", va="top")
    ax.text(RIGHT[0] + pad * 0.55, RIGHT[1] - pad * 0.55, right, size=label_size,
            color="k", ha="left", va="top")
    if corner_values:
        # The apex needs a longer inward offset than the base corners: it is the
        # acute corner, its two walls converge on the label, and a wide value like
        # "100" crowds them where "1" does not.
        for corner, inward, reach in ((APEX, np.array([0.0, -1.0]), 2.9),
                                      (LEFT, np.array([_H, 0.5]), 1.9),
                                      (RIGHT, np.array([-_H, 0.5]), 1.9)):
            p = corner + inward * pad * reach
            ax.text(p[0], p[1], one, size=tick_size, color="k",
                    ha="center", va="center")

    # --- one axis label per edge, set outside it and turned to match ----------
    # Set well clear of the tick numbers (which sit at 0.55 pad) and of anything
    # `read_construction` places (1.7 pad), because all three live on the same
    # outward normals.
    if symbol:
        def _sub(name):
            return rf"${symbol}_{{\mathrm{{{name}}}}}$"
        for mid, normal, rot, s in (
                ((APEX + RIGHT) / 2.0, np.array([_H, 0.5]), -60, top),
                ((APEX + LEFT) / 2.0, np.array([-_H, 0.5]), 60, left),
                ((LEFT + RIGHT) / 2.0, np.array([0.0, -1.0]), 0, right)):
            p = mid + normal * pad * 3.4
            ax.text(p[0], p[1], _sub(s), size=label_size, color="k",
                    ha="center", va="center", rotation=rot,
                    rotation_mode="anchor")

    ax.set_xlim(-0.22, 1.22)
    ax.set_ylim(-0.20, _H + 0.15)
    ax.set_aspect("equal")
    ax.axis("off")
    return ax


# --- drawing on it ---------------------------------------------------------

def plot(ax, *comp, **kw):
    """Plot a curve. `plot(ax, arr)` or `plot(ax, top, left, right)`.

    The binodal of Fig. 11.2-8 is one call. Defaults to the heaviest black line in
    the module, because on these diagrams the phase boundary is the spine the way
    the saturation line is on a property chart.
    """
    xy = to_xy(_unpack(*comp))
    kw.setdefault("color", "k")
    kw.setdefault("lw", 1.1)
    kw.setdefault("zorder", 4)
    return ax.plot(xy[:, 0], xy[:, 1], **kw)


def scatter(ax, *comp, **kw):
    """Measured points. Open circles by default, so they read against a black curve."""
    xy = to_xy(_unpack(*comp))
    kw.setdefault("marker", "o")
    kw.setdefault("s", 12)
    kw.setdefault("facecolors", "none")
    kw.setdefault("edgecolors", "k")
    kw.setdefault("linewidths", 0.7)
    kw.setdefault("zorder", 5)
    return ax.scatter(xy[:, 0], xy[:, 1], **kw)


def tie_line(ax, cI, cII, *, ends=True, **kw):
    """One tie line, between two coexisting phases.

    Lighter and dashed by default so a fan of them does not compete with the
    binodal -- weight and dash, never color. `ends=True` marks both ends, which is
    what makes a tie line readable as *two phases* rather than as a chord.
    """
    p, q = to_xy([cI, cII])
    kw.setdefault("color", "k")
    kw.setdefault("lw", 0.6)
    kw.setdefault("zorder", 3)
    out = ax.plot([p[0], q[0]], [p[1], q[1]], **kw)
    if ends:
        ax.scatter([p[0], q[0]], [p[1], q[1]], s=7, c="k", zorder=5,
                   linewidths=0)
    return out


def point(ax, comp, label=None, *, size=7, offset=(0.018, 0.018), marker="o",
          ms=3.2, **kw):
    """A single labeled composition -- a feed point, a plait point, a mixture."""
    xy = to_xy(comp).ravel()
    kw.setdefault("color", "k")
    kw.setdefault("zorder", 6)
    ax.plot([xy[0]], [xy[1]], marker=marker, ms=ms, **kw)
    if label:
        if size < 7:
            raise ValueError("nothing below 7 pt [[no-figure-text-below-7pt]]")
        ax.text(xy[0] + offset[0], xy[1] + offset[1], label, size=size, color="k",
                ha="left", va="bottom")
    return xy


def text(ax, comp, s, *, size=7, **kw):
    """Free text placed at a composition rather than at an (x, y)."""
    if size < 7:
        raise ValueError("nothing below 7 pt [[no-figure-text-below-7pt]]")
    xy = to_xy(comp).ravel()
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    kw.setdefault("color", "k")
    return ax.text(xy[0], xy[1], s, size=size, **kw)


def region_label(ax, comps, s, *, size=7, **kw):
    """Label a region, placed at the centroid of the compositions bounding it.

    The `2L` / `3L` labels of Fig. 11.2-11(e) -- and the reason it takes a list
    rather than a point is that a region's visual center is rarely a composition
    anybody has a number for.
    """
    c = _norm(np.atleast_2d(np.asarray(comps, dtype=float))).mean(axis=0)
    return text(ax, c, s, size=size, **kw)


def read_construction(ax, comp, *, symbol="x", names=None, size=7, lw=0.9,
                      arrows=True, which="tlr"):
    """Show how a point is READ: one construction line per species, to its own edge.

    This is Figure 11.2-7, which is not a data plot but the schematic that teaches
    the diagram. For each species it draws the constant-composition line through
    `comp` -- horizontal for the apex species, parallel to an edge for the other two
    -- out to the edge carrying that species' scale, and labels it with the value.

    `which` selects which of the three to draw ("tlr" for all; Fig. 11.2-7 draws all
    three, and a figure that only needs to make the point once can ask for one).

    The three lines meet at the point *by construction*, and that is the lesson:
    two of the three fractions fix the composition and the third is not free.
    """
    c = _norm(np.asarray(comp, dtype=float)).ravel()
    p = to_xy(c).ravel()
    if names is None:
        names = ("A", "B", "C")
    ends = _edge_points_for(c)
    order = {"t": 0, "l": 1, "r": 2}
    for ch in which:
        i = order[ch]
        e = ends[i].ravel()
        ax.plot([p[0], e[0]], [p[1], e[1]], "-", color="k", lw=lw, zorder=4)
        if arrows:
            ax.annotate("", xy=e, xytext=p, zorder=4,
                        arrowprops=dict(arrowstyle="-|>", color="k", lw=lw,
                                        shrinkA=0, shrinkB=0, mutation_scale=7))
        lbl = rf"${symbol}_{{\mathrm{{{names[i]}}}}} = {c[i]:g}$"
        if size < 7:
            raise ValueError("nothing below 7 pt [[no-figure-text-below-7pt]]")
        # placed on the same outward normal as that edge's tick numbers, but
        # further out, so the value clears the scale it is being read against
        normal = (np.array([_H, 0.5]), np.array([-_H, 0.5]),
                  np.array([0.0, -1.0]))[i]
        q = e + normal * 0.095
        ha = ("left", "right", "center")[i]
        va = ("center", "center", "top")[i]
        ax.text(q[0], q[1], lbl, size=size, color="k", ha=ha, va=va)
    ax.plot([p[0]], [p[1]], marker="o", ms=3.0, color="k", zorder=6)
    return p


def _edge_points_for(c):
    """Where each species' constant-composition line through `c` meets its own edge.

    Same edge assignment as `_edge_points`: apex species to the right edge, left
    species to the left edge, right species to the bottom edge.
    """
    t, l, r = float(c[0]), float(c[1]), float(c[2])
    return (to_xy([t, 0.0, 1 - t]),        # constant top, out to the right edge
            to_xy([1 - l, l, 0.0]),        # constant left, out to the left edge
            to_xy([0.0, 1 - r, r]))        # constant right, down to the bottom edge


def lever_arm(ax, z, cI, cII, *, label_I=None, label_II=None, size=7):
    """Draw the lever rule: the feed `z` on the tie line joining `cI` and `cII`.

    Illustrations 11.2-7, 11.2-8 and 11.2-9 are this picture. The feed is marked
    with a filled square so it cannot be mistaken for a phase, and the two arms are
    drawn at different weights only if labeled -- the *lengths* carry the meaning,
    and they are already drawn to scale by construction.

    It does NOT check that `z` lies on the segment. `thermo.lle.tie_line_split`
    does that properly, by least squares over all species, and returns the residual
    that says how far off the data are. Drawing is not validation.
    """
    tie_line(ax, cI, cII, ends=True, lw=0.8, ls="-")
    p = to_xy(z).ravel()
    ax.plot([p[0]], [p[1]], marker="s", ms=3.6, color="k", zorder=6)
    a, b = to_xy([cI, cII])
    if label_I:
        m = (a + p) / 2.0
        ax.text(m[0], m[1] - 0.030, label_I, size=size, color="k",
                ha="center", va="top")
    if label_II:
        m = (p + b) / 2.0
        ax.text(m[0], m[1] - 0.030, label_II, size=size, color="k",
                ha="center", va="top")
    return p


def _renderer(fig):
    """A renderer that works under `jupyter nbconvert --execute` too.

    `fig.canvas.get_renderer()` exists on the Agg canvas and NOT on the plain
    `FigureCanvasBase` a headless kernel can hand you, so a check that passes in
    JupyterLab can die with AttributeError in a batch re-execution -- which is
    exactly where a figure gate needs to work. Fall back to attaching an Agg canvas.
    """
    try:
        return fig.canvas.get_renderer()
    except AttributeError:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        return canvas.get_renderer()


def check_labels(ax, *, tol=0.0, verbose=True):
    """Fail if any two text labels on `ax` overlap. Call it before saving.

    **This is the ch6/ch11 lesson made reusable.** `tools/check_print_art.py`
    gates staged art on color and font family; it cannot see that two labels are
    sitting on top of each other, and a triangular diagram has three scales, three
    corner names, three edge symbols and whatever the figure itself adds -- all
    crowded around the same small perimeter. The first render of Fig. 11.2-7 put
    `x_C = 0.25` straight through the `0.2` tick and `100W` through the corner name,
    and both looked deliberate at thumbnail size.

    Returns the list of overlapping pairs. Raises when there are any, unless `tol`
    is set to allow a small overlap in points.
    """
    fig = ax.figure
    fig.canvas.draw()                      # extents are only real after a draw
    rend = _renderer(fig)
    # TEXT box only. `Annotation.get_window_extent` includes the leader, so a
    # label that correctly points at a curve would be reported as colliding with
    # every other label its leader happens to pass. `read_construction` draws
    # leaders, so this matters here and not only in the notebooks.
    import matplotlib.text as mtext
    texts = [t for t in ax.texts if (t.get_text() or "").strip()]
    boxes = [(t, mtext.Text.get_window_extent(t, renderer=rend)) for t in texts]
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            (ti, bi), (tj, bj) = boxes[i], boxes[j]
            dx = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
            dy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            if dx > tol and dy > tol:
                bad.append((ti.get_text(), tj.get_text(), round(float(dx), 1),
                            round(float(dy), 1)))
    if bad:
        if verbose:
            for a, b, dx, dy in bad:
                print(f"  label overlap: {a!r} and {b!r}  ({dx} x {dy} px)")
        raise AssertionError(f"{len(bad)} overlapping label pair(s) on this axes; "
                             "shorten a label, do not shrink the type "
                             "[[no-figure-text-below-7pt]]")
    return bad
