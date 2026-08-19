"""ph_chart — pressure-enthalpy charts from a cubic equation of state.

The machinery behind **Figures 3.3-2 (methane) and 3.3-3 (nitrogen)**, and behind
the two derived figures that redraw them with a process path: **Figure 5.1-3** (the
LNG liquefaction path of Illustration 5.1-1) and **`c06uf002`** (the
cylinder-discharge window of Illustration 6.5-1).

    from thermo import PengRobinson, TABLE_6_6_1, APPENDIX_A2_CP_CRYO
    from thermo.ph_chart import ChartFluid, ph_chart
    from thermo.charts import use_book_style

    use_book_style()
    ch4 = ChartFluid("methane", M=16.043, T_triple=90.69)
    fig, ax = plt.subplots(figsize=(6.5, 8.0))
    ph_chart(ax, ch4, H_lim=(0, 1300), isotherms=range(120, 401, 40), ...)

TWO THINGS THIS MODULE IS PARTICULAR ABOUT, both of which cost a visible defect
when they are got wrong, and both of which are commented at the point of use:

1. **The dome closes on the ANALYTIC critical point.** The vapor-pressure loop dies
   a hair below `Tc`, and the cubic degenerates to a triple root exactly at it, so
   neither the loop nor the root finder can supply that point. `ChartFluid.crit`
   evaluates it in closed form instead.
2. **Constant-volume curves are pinned to the dome.** An isochore that crosses the
   two-phase region leaves it at one exact temperature; without solving for that
   temperature the polyline steps across the corner and the curve visibly bends away
   before it reaches the envelope.

**Peng-Robinson puts saturated-liquid density about 12 percent high**, which is
the textbook failing of a cubic. It lands on exactly one family: the constant-volume
curves in the compressed-liquid region, left of the saturation curve. Those should
not be read quantitatively there. Nothing else on the chart is affected -- the dome,
the isotherms and the isentropes all come through within chart-reading accuracy.
"""
from __future__ import annotations

import numpy as np
from scipy.constants import R
from scipy.optimize import brentq

from .charts import GRAY, GRAY_MINOR, LW, label_at, log_pressure_grid
from .data import APPENDIX_A2_CP_CRYO, TABLE_6_6_1
from .peng_robinson import PengRobinson

TREF = 298.15          # K, where the ideal-gas integrals start (arbitrary, cancels)
PREF = 1e5             # Pa


class ChartFluid:
    """One fluid's Peng-Robinson model, per kilogram, on the chart's own datum.

    The charts are read in kJ/kg and m^3/kg, and their datum is the one the
    cryogenic charts they replace used: **saturated liquid at the triple point has
    H = S = 0**. Both of those are chart conventions rather than thermodynamics, so
    they live here rather than in `PengRobinson`.

    `cp` defaults to the compound's `APPENDIX_A2_CP_CRYO` row and `Tc/Pc/omega` to
    its `TABLE_6_6_1` entry -- the book's own constants, so a reader can reproduce
    the chart from the printed appendix.
    """

    T_TOP = 1 - 1e-4            # the spinodal loop vanishes numerically before Tc

    # The critical compressibility of the Peng-Robinson equation. It is ANALYTIC,
    # and it has to be, because the vapor-pressure loop cannot reach it: at
    # (Tc, Pc) the cubic collapses to the triple root Zc = 0.30740. Do not ask the
    # root finder for it -- a triple root spreads under rounding, and
    # `pr.compressibility(Tc, Pc)` comes back 0.32138, wrong in the second digit.
    ZC = 0.30740130869693604

    def __init__(self, key, M, T_triple, cp=None, constants=None):
        self.key, self.M, self.T_triple = key, M, T_triple
        self.cp = cp if cp is not None else APPENDIX_A2_CP_CRYO[key]
        self.pr = PengRobinson(**(constants or TABLE_6_6_1[key]), cp=self.cp)
        self.Tc, self.Pc = self.pr.Tc, self.pr.Pc
        self.H0 = self.S0 = 0.0
        P = self.Psat(T_triple)
        self.H0 = self.H(T_triple, P, "liquid")
        self.S0 = self.S(T_triple, P, "liquid")

    # --- ideal gas ---------------------------------------------------------
    def _H_ig(self, T):
        a, b, c, d = self.cp
        f = lambda t: a * t + b * t**2 / 2 + c * t**3 / 3 + d * t**4 / 4
        return f(T) - f(TREF)

    def _S_ig(self, T, P):
        a, b, c, d = self.cp
        f = lambda t: a * np.log(t) + b * t + c * t**2 / 2 + d * t**3 / 3
        return f(T) - f(TREF) - R * np.log(P / PREF)

    # --- per kilogram ------------------------------------------------------
    def H(self, T, P, phase):
        """kJ/kg  (J/mol divided by g/mol is J/g, which is kJ/kg)."""
        return (self._H_ig(T) + self.pr.departure_H(T, P, phase)) / self.M - self.H0

    def S(self, T, P, phase):
        """kJ/(kg K)."""
        return (self._S_ig(T, P) + self.pr.departure_S(T, P, phase)) / self.M - self.S0

    def V(self, T, P, phase):
        """m^3/kg."""
        return self.pr.molar_volume(T, P, phase) * 1e3 / self.M

    # --- saturation --------------------------------------------------------
    def Psat(self, T):
        try:
            return self.pr.vapor_pressure(T)
        except (ValueError, TypeError):
            return None

    def sat(self, T):
        """(P, HL, HV, SL, SV, VL, VV) at temperature T."""
        P = self.Psat(T)
        VL, VV = self.pr.saturation_volumes(T, P)
        return (P, self.H(T, P, "liquid"), self.H(T, P, "vapor"),
                self.S(T, P, "liquid"), self.S(T, P, "vapor"),
                VL * 1e3 / self.M, VV * 1e3 / self.M)

    def crit(self):
        """(P, H, S, V) at the critical point, from the analytic triple root."""
        T, P, Z = self.Tc, self.Pc, self.ZC
        B = self.pr.b * P / (R * T)
        s2 = np.sqrt(2.0)
        lt = np.log((Z + (1 + s2) * B) / (Z + (1 - s2) * B))
        k = 1 / (2 * s2 * self.pr.b)
        depH = R * T * (Z - 1) + (T * self.pr.dadT(T) - self.pr.a(T)) * k * lt
        depS = R * np.log(Z - B) + self.pr.dadT(T) * k * lt
        return (P, (self._H_ig(T) + depH) / self.M - self.H0,
                (self._S_ig(T, P) + depS) / self.M - self.S0,
                Z * R * T / P * 1e3 / self.M)

    def phase_at(self, T, P):
        """Which root to take at (T, P): 'liquid', 'vapor', or the one that exists."""
        if T >= self.Tc:
            return "vapor" if P < self.Pc else "liquid"   # one root; label is cosmetic
        Ps = self.Psat(T)
        if Ps is None:
            return "liquid" if P > 1e7 else "vapor"
        return "liquid" if P > Ps else "vapor"


# --- the curve families ----------------------------------------------------

def sat_T(f, n):
    """Saturation temperatures, bunched toward Tc.

    The dome's whole shoulder -- where the two branches turn over and meet -- lies in
    the last percent of the temperature range, so a grid uniform in T spends its
    points where the curve is straight and starves the part that has the corner.
    """
    u = np.linspace(0.0, 1.0, n)
    return f.T_triple + (f.Tc * f.T_TOP - f.T_triple) * (1 - (1 - u)**2)


def dome(f, n=400):
    """Saturation dome, CLOSED on the critical point.

    Returns (T, P_MPa, HL, HV, SL, SV).

    The vapor-pressure loop dies a hair below Tc, and stopping there leaves the two
    branches blunt-ended with open air between them -- 15 kJ/kg of it for nitrogen.
    The critical point does not need the loop (`ChartFluid.crit`), so append it:
    both branches then end on the one point they are required to share.
    """
    T = sat_T(f, n)
    rows = [f.sat(t) for t in T]
    Pc, Hc, Sc, _ = f.crit()
    col = lambda j, last: np.append([r[j] for r in rows], last)
    return (np.append(T, f.Tc), col(0, Pc) / 1e6, col(1, Hc), col(2, Hc),
            col(3, Sc), col(4, Sc))


def isotherm(f, T, P_lim, n=260):
    """[(H, P_MPa), ...] segments along one isotherm; a sub-critical one is 3 pieces."""
    Plo, Phi = P_lim[0] * 1e6, P_lim[1] * 1e6
    segs = []
    if T < f.Tc:
        Ps = f.Psat(T)
        if Ps is not None and Plo < Ps < Phi:
            Pv = np.geomspace(Plo, Ps, n)
            segs.append((np.array([f.H(T, p, "vapor") for p in Pv]), Pv / 1e6))
            _, HL, HV, *_ = f.sat(T)
            segs.append((np.array([HL, HV]), np.array([Ps, Ps]) / 1e6))
            Pl = np.geomspace(Ps, Phi, n)
            segs.append((np.array([f.H(T, p, "liquid") for p in Pl]), Pl / 1e6))
            return segs
    P = np.geomspace(Plo, Phi, 2 * n)
    H = np.array([f.H(T, p, f.phase_at(T, p)) for p in P])
    return [(H, P / 1e6)]


def isentrope(f, S0, P_lim, n=220):
    """Constant-entropy line: at each pressure, find the state with entropy S0."""
    P = np.geomspace(P_lim[0] * 1e6, P_lim[1] * 1e6, n)
    H, keep = [], []
    for p in P:
        Ts = None
        # A pressure just under Pc can still be above the highest vapor pressure the
        # loop survives to, which is a hair below Pc -- so test the bracket, do not
        # assume p < Pc means a saturation temperature exists.
        if p < f.Pc and f.Psat(f.Tc * f.T_TOP) > p:
            Ts = brentq(lambda t: f.Psat(t) - p, f.T_triple, f.Tc * f.T_TOP, xtol=1e-8)
            _, HL, HV, SL, SV, _, _ = f.sat(Ts)
            if SL <= S0 <= SV:                       # inside the dome: lever rule
                x = (S0 - SL) / (SV - SL)
                H.append(HL + x * (HV - HL)); keep.append(p); continue
        lo, hi = (f.T_triple, Ts - 1e-6) if (Ts and S0 < f.sat(Ts)[3]) else \
                 ((Ts + 1e-6, 1200.0) if Ts else (f.T_triple, 1200.0))
        ph = "liquid" if (Ts and S0 < f.sat(Ts)[3]) else "vapor"
        try:
            g = lambda t: f.S(t, p, f.phase_at(t, p) if not Ts else ph) - S0
            if g(lo) * g(hi) > 0:
                continue
            T = brentq(g, lo, hi, xtol=1e-6)
            H.append(f.H(T, p, f.phase_at(T, p) if not Ts else ph)); keep.append(p)
        except (ValueError, RuntimeError):
            continue
    return np.array(H), np.array(keep) / 1e6


def isochore_T(f, Vhat, n):
    """Temperatures for one constant-volume curve, PINNED TO THE DOME.

    An isochore that crosses the two-phase region leaves it at the one temperature
    where Vhat is exactly a saturated volume. Solve for that temperature and put
    points on both sides of it. Without this the polyline steps straight across the
    corner and the curve appears to bend away before it reaches the envelope -- a
    visible discontinuity, worst on the curves whose junction is a long way from a
    grid point.

    Dense below Tc, where the curve both crosses the dome and turns hardest; sparse
    above it, where it is close to straight.
    """
    T = np.concatenate([np.linspace(f.T_triple, f.Tc * f.T_TOP, n),
                        np.linspace(f.Tc, 1000.0, n)])
    pins = []
    for idx in (5, 6):                                # sat(T) -> (..., VL, VV)
        try:
            Tj = brentq(lambda t: f.sat(t)[idx] - Vhat,
                        f.T_triple, f.Tc * f.T_TOP, xtol=1e-12)
        except ValueError:
            continue                                  # this branch never reaches Vhat
        pins += list(Tj + np.array([-.4, -.1, -.02, -1e-6, 1e-6, .02, .1, .4]))
    if pins:
        T = np.concatenate([T, np.clip(pins, f.T_triple, 1000.0)])
    return np.unique(T)


def isochore(f, Vhat, P_lim, n=300):
    """Constant specific volume, m^3/kg. Parametrised by T: PR is explicit in P."""
    Vm = Vhat * f.M / 1e3                             # m^3/mol
    H, P = [], []
    for t in isochore_T(f, Vhat, n):
        if t < f.Tc:
            _, HL, HV, _, _, VL, VV = f.sat(t)
            if VL <= Vhat <= VV:                      # inside the dome
                x = (Vhat - VL) / (VV - VL)
                H.append(HL + x * (HV - HL)); P.append(f.Psat(t)); continue
        p = f.pr.pressure(Vm, t)
        if p <= 0:
            continue
        H.append(f.H(t, p, f.phase_at(t, p))); P.append(p)
    P = np.array(P) / 1e6; H = np.array(H)
    m = (P >= P_lim[0]) & (P <= P_lim[1])
    return H[m], P[m]


def quality(f, x, n=250):
    """The lever rule up the dome.

    Every quality line ends at the critical point too -- at Tc the lever has zero
    length, so all of them converge there.
    """
    T = sat_T(f, n)
    rows = [f.sat(t) for t in T]
    Pc, Hc, _, _ = f.crit()
    P = np.append([r[0] for r in rows], Pc) / 1e6
    H = np.append([r[1] + x * (r[2] - r[1]) for r in rows], Hc)
    return H, P


# --- the chart -------------------------------------------------------------

def ph_chart(ax, f, *, H_lim, P_lim=(0.1, 20), isotherms=(), isotherms_minor=(),
             isentropes=(), volumes=(), qualities=(), H_grid=50, H_minors=5,
             T_label_P=None, S_label_P=None, V_label_P=None, x_label_P=None,
             family_labels=None, label_fmt=None, mark_critical=False, title=None):
    """Draw one P-H chart on `ax`.

    Each family is labeled at its OWN pressure -- isotherms high, isentropes in the
    middle, constant-volume lines low -- because all three run steeply through the
    same corner and a shared label height puts them on top of one another. A
    `*_label_P` may also be a `{value: pressure}` mapping when one height will not
    do for the whole family, since the liquid-side members run through a different
    part of the frame from the vapor-side ones.

    `label_fmt` gives a family's labels their units, e.g.
    `{"T": "{:g} K", "S": "{:g} kJ/(kg K)", "V": "{:g} m$^3$/kg"}`. The full-size
    charts leave it off and name each family once instead (`family_labels`), because
    30 isotherms cannot each carry a unit. A small window with four of them can, and
    should.

    `mark_critical=True` puts a dot on the critical point. Off by default, which is
    how the printed figures are drawn: the dome closing on it, and the quality lines
    converging there, already say where it is. Turn it on when the critical point is
    what is being discussed rather than something to be read past.
    """
    ax.set_xlim(*H_lim)
    ax.set_yscale("log")
    ax.set_ylim(*P_lim)
    H_right = H_lim[1] - 0.04 * (H_lim[1] - H_lim[0])

    for T in isotherms_minor:
        for H, P in isotherm(f, T, P_lim):
            ax.plot(H, P, "-", color=GRAY_MINOR, lw=LW["T_minor"], zorder=1)

    for x in qualities:
        H, P = quality(f, x)
        ax.plot(H, P, "--", color=GRAY, lw=LW["x"], zorder=2, dashes=(2.6, 2.2))
        label_at(ax, H, P, f"{x:g}", along="y", value=x_label_P, color=GRAY,
                 size=5.5, pad=0.5)

    def fam_text(kind, v):
        spec = (label_fmt or {}).get(kind)
        return spec.format(v) if spec else f"{v:g}"

    def label_family(H, P, txt, spec, key, color, size=5.5):
        """Label one curve, falling back to the right-hand edge then to the bottom."""
        value = spec.get(key, spec.get(None)) if isinstance(spec, dict) else spec
        if value is not None and label_at(ax, H, P, txt, along="y", value=value,
                                          color=color, size=size, pad=0.5):
            return True
        if label_at(ax, H, P, txt, along="x", value=H_right, color=color, size=size,
                    pad=0.5):
            return True
        return label_at(ax, H, P, txt, along="y", value=P_lim[0] * 1.12, color=color,
                        size=size, pad=0.5)

    for V in volumes:
        H, P = isochore(f, V, P_lim)
        if len(H) > 1:
            ax.plot(H, P, "--", color="0.25", lw=LW["V"], zorder=3, dashes=(5, 2))
            label_family(H, P, fam_text("V", V), V_label_P, V, "0.25")

    for S in isentropes:
        H, P = isentrope(f, S, P_lim)
        if len(H) > 1:
            ax.plot(H, P, "-", color=GRAY, lw=LW["S"], zorder=2)
            label_family(H, P, fam_text("S", S), S_label_P, S, GRAY)

    for T in isotherms:
        segs = isotherm(f, T, P_lim)
        for H, P in segs:
            ax.plot(H, P, "-", color="k", lw=LW["T"], zorder=4)
        H, P = max(segs, key=lambda s: float(np.nanmax(s[0])))
        label_family(H, P, fam_text("T", T), T_label_P, T, "k")

    # one member of each family carries the family's name and its units
    for spec in (family_labels or []):
        kind, value, at_P, text = spec
        if kind == "T":
            segs = isotherm(f, value, P_lim)
            H, P = max(segs, key=lambda s: float(np.nanmax(s[0])))
        elif kind == "S":
            H, P = isentrope(f, value, P_lim)
        else:
            H, P = isochore(f, value, P_lim)
        label_at(ax, H, P, text, along="y", value=at_P, size=6, pad=0.5,
                 color="k" if kind == "T" else (GRAY if kind == "S" else "0.25"))

    # The dome is ONE continuous path: up the liquid branch, around the critical
    # point, back down the vapor branch. Drawn as two separate curves it cannot
    # close the corner, and the top of the envelope is left standing open.
    _, Pd, HL, HV, _, _ = dome(f)
    ax.plot(np.concatenate([HL, HV[::-1]]), np.concatenate([Pd, Pd[::-1]]),
            "-", color="k", lw=LW["sat"], zorder=5, solid_joinstyle="round")
    if mark_critical:
        Pc, Hc, _, _ = f.crit()
        ax.plot([Hc], [Pc / 1e6], "o", color="k", ms=2.6, zorder=6)

    log_pressure_grid(ax, H_grid, H_minors)
    ax.set_xlabel(r"enthalpy  $\hat H$  (kJ/kg)")
    ax.set_ylabel(r"pressure  $P$  (MPa)")
    if title:
        ax.set_title(title, fontsize=9)
    return ax


__all__ = ["ChartFluid", "ph_chart", "dome", "isotherm", "isentrope", "isochore",
           "quality", "sat_T", "isochore_T", "TREF", "PREF"]
