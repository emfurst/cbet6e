"""steam_chart — the Mollier and temperature-entropy charts, from Appendix A.III.

The machinery behind **Figures 3.3-1(a) and 3.3-1(b)**, and behind **`c05uf001`**,
the Rankine cycle of Illustration 5.2-1 drawn on the $T$-$\\hat S$ chart.

Unlike `thermo.ph_chart`, none of this comes from an equation of state. The curves
are built from the book's own **Appendix A.III**, digitized to CSV in `code/data/`:
the saturation table is the dome, the superheat table gives the isobars, and the
compressed-liquid table gives the liquid branches. That is why steam gets its own
module rather than a `ChartFluid` -- for water the book *tabulates* the properties,
and a cubic would be a step backwards from data the reader already has.

The CSVs are the appendix as the 6e prints it, five corrected cells included, so
nothing here has to patch the data before it draws with it. How those five were
found and what each was corrected to is
`code/ch3/validation/Appendix_A3_corrections.ipynb`.

    from thermo.steam_chart import SteamTables, temperature_entropy
    from thermo.charts import use_book_style

    use_book_style()
    st = SteamTables()                       # finds code/data/ automatically
    fig, ax = plt.subplots(figsize=(7.5, 9.5))
    temperature_entropy(ax, st)

THE HARD PART IS NOT THE THERMODYNAMICS, IT IS WHERE A CURVE STOPS. Three families
are built from two different tables and have to *meet* on the saturation curve;
below 0.01 MPa the superheat table simply ends and the curves must be continued on
the ideal-gas limit; and a line of constant enthalpy can exist over two disjoint
pressure ranges, so joining its points into one polyline draws a chord clear across
the chart. Each is handled at the point it arises and commented there.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from .charts import (GRAY, GRAY_MINOR, LW, chart_grid, fmt_P, label_at,
                     label_end)

# The critical point: the last row of the saturation table.
TC, PC = 374.14, 22.09
HC, SC = 2099.3, 4.4298

# R/M for water, kJ/(kg K) -- the ideal-gas slope the low-pressure continuations use.
R_M = 8.314462618 / 18.015268e-3 / 1e3

# MPa. The triple point: the first row of the saturation pressure table, and the
# lowest pressure at which there is a saturated LIQUID at all. Below it the dome does
# not exist, so it is the floor for an isobar on either chart.
P_TRIPLE = 0.0006113



def _default_data_dir():
    """`code/data/`, whether the caller is in code/, code/chN/, or a Colab copy."""
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, os.pardir, "data"),
                 os.path.join(os.getcwd(), "data"),
                 os.path.join(os.getcwd(), os.pardir, "data")):
        if os.path.isdir(cand):
            return os.path.abspath(cand)
    raise FileNotFoundError(
        "could not find code/data/ — pass SteamTables(data_dir=...) explicitly")


class SteamTables:
    """Appendix A.III, digitized, with the curve families the charts are drawn from.

    Holds the four tables and the interpolators built from them. Every curve builder
    is a method, so one instance carries the whole chart's data and its caches.
    """

    def __init__(self, data_dir=None):
        d = data_dir or _default_data_dir()
        self.data_dir = d
        self.sat = pd.read_csv(f"{d}/steam_saturation_T.csv")      # the dome
        self.sh = pd.read_csv(f"{d}/steam_superheat.csv")          # superheated vapor
        self.cl = pd.read_csv(f"{d}/steam_compressed_liquid.csv")  # compressed liquid
        self.solid = pd.read_csv(f"{d}/steam_solid_vapor.csv")     # solid-vapor

        self.ALL_P = tuple(sorted(self.sh.P_MPa.unique()))
        # the dome, as smooth functions of temperature
        self.dome = {c: PchipInterpolator(self.sat.T_C, self.sat[c]) for c in
                     ("P_MPa", "Hl", "Hv", "Sl", "Sv", "dH", "dS")}
        self.Tsat_of_P = PchipInterpolator(self.sat.P_MPa, self.sat.T_C)

        self.CL = self.cl[~self.cl["sat"]].copy()

        self._isobar_cache = {}
        self._merged_cache = {}

    # --- the dome ----------------------------------------------------------
    def quality_line(self, x, T=None):
        """(S, H) along a line of constant quality x, swept up the dome (lever rule)."""
        T = np.linspace(self.sat.T_C.min(), TC, 400) if T is None else T
        return (self.dome["Sl"](T) + x * self.dome["dS"](T),
                self.dome["Hl"](T) + x * self.dome["dH"](T))

    # --- isobars -----------------------------------------------------------
    def isobar(self, P, table=None, raw=False):
        """(S, H, T) up one isobar: its saturated-vapor anchor, then the superheat.

        `raw=True` returns the tabulated knots instead of the interpolated curve, for
        a caller that needs to join them to knots from another table before smoothing.
        """
        if table is None and (P, raw) in self._isobar_cache:
            return self._isobar_cache[(P, raw)]
        if table is None and P < self.ALL_P[0]:        # below the superheat table
            out = self.low_pressure_isobar(P)
            self._isobar_cache[(P, raw)] = out
            return out
        t = self.sh if table is None else table
        g = t[t.P_MPa == P]
        if g.empty:
            raise KeyError(f"no isobar at P = {P} MPa; have {sorted(t.P_MPa.unique())}")
        anchor, g = g[g["sat"]], g[~g["sat"]].sort_values("T_C")
        T, S, H = g.T_C.to_numpy(), g.S.to_numpy(), g.H.to_numpy()
        if len(anchor):                        # start the curve on the dome
            T = np.r_[float(anchor.Tsat_C.iloc[0]), T]
            S = np.r_[float(anchor.S.iloc[0]), S]
            H = np.r_[float(anchor.H.iloc[0]), H]
        if raw:
            out = (S, H, T)
        else:
            fine = np.linspace(T.min(), T.max(), 400)
            out = (PchipInterpolator(T, S)(fine), PchipInterpolator(T, H)(fine), fine)
        if table is None:
            self._isobar_cache[(P, raw)] = out
        return out

    def low_pressure_isobar(self, P):
        """(S, H, T) for an isobar below 0.01 MPa, the lowest the superheat table has.

        Not an extrapolation of the chart's own curves. At these pressures steam is
        very nearly an ideal gas, for which, at fixed temperature,

            S(P) = S(P0) - (R/M) ln(P/P0)        and        H(P) = H(P0)

        exactly. Every point is anchored on the nearest state the appendix does
        tabulate: below 45.81 C on the saturated vapor at that temperature, whose
        pressure is already at or below 0.01 MPa, and above it on the 0.01 MPa isobar.
        The two agree where they meet, because Psat(45.81 C) IS 0.01 MPa -- so the
        curve leaves the saturation line exactly at Tsat(P) and joins the tabulated
        data without a step.

        What this neglects is the residual pressure dependence of H, under
        1.5 kJ/kg here: well below the width of the printed line.
        """
        if P < P_TRIPLE:
            raise ValueError(f"{P} MPa is below the triple point, {P_TRIPLE} MPa: "
                             "there is no saturated liquid, so no isobar to anchor "
                             "on the dome")
        S0, H0, T0 = self.isobar(0.01)
        T_join = float(T0.min())                   # 45.81 C: where 0.01 MPa saturates

        Ta = np.linspace(float(self.Tsat_of_P(P)), T_join, 150)   # anchored on the dome
        Sa = self.dome["Sv"](Ta) + R_M * np.log(self.dome["P_MPa"](Ta) / P)
        keep = T0 >= T_join                                       # anchored on 0.01 MPa
        return (np.r_[Sa, S0[keep] + R_M * np.log(0.01 / P)],
                np.r_[self.dome["Hv"](Ta), H0[keep]],
                np.r_[Ta, T0[keep]])

    def liquid_branch(self, P):
        """(S, H, T) up the liquid side of an isobar, from the compressed-liquid table.

        An isobar at a pressure the compressed-liquid table does not carry --
        including 60 MPa, the highest on the chart -- would otherwise begin abruptly
        at 375 C, where the superheat table starts. Compressed-liquid properties are
        very nearly linear in pressure at fixed temperature, so the missing branches
        are filled in across pressure: an interpolation for 25, 35 and 40 MPa, and a
        one-column extrapolation for 60.
        """
        cols, T = {"S": [], "H": []}, []
        for t, g in self.CL.groupby("T_C"):
            g = g.sort_values("P_MPa")
            p = g.P_MPa.to_numpy()
            if len(p) < 2 or P < p[0]:
                continue                       # no clamping: skip, rather than guess
            for name, out in cols.items():
                v = g[name].to_numpy()
                if P <= p[-1]:
                    out.append(float(np.interp(P, p, v)))
                else:                          # linear in P beyond the last column
                    out.append(float(v[-1] + (v[-1] - v[-2]) / (p[-1] - p[-2])
                                     * (P - p[-1])))
            T.append(t)
        if len(T) < 2:
            return None
        S, H, T = np.array(cols["S"]), np.array(cols["H"]), np.array(T)
        o = np.argsort(T)
        S, H, T = S[o], H[o], T[o]
        if P < PC:
            # The table's last row is 2-11 C below saturation, so the branch is closed
            # on the saturated liquid itself. That is not only tidier to draw: without
            # it the branch stops short in ENTHALPY too, and a line of constant
            # enthalpy that ends between the two never finds this isobar at all.
            Ts = float(self.Tsat_of_P(P))
            keep = T <= Ts
            S = np.r_[S[keep], float(self.dome["Sl"](Ts))]
            H = np.r_[H[keep], float(self.dome["Hl"](Ts))]
            T = np.r_[T[keep], Ts]
        return (S, H, T) if len(T) > 1 else None

    def merged_isobar(self, P):
        """(S, H, T) for a SUPERCRITICAL isobar: compressed liquid and superheat joined.

        Above the critical pressure there is no phase change, so the two tables
        describe one curve and it must be drawn -- and searched -- as one. They
        overlap, liquid running to 380 C and superheat starting at 375, and the liquid
        knots within 1 C of the superheat table are dropped: that is where the pressure
        interpolation of the liquid side is weakest, since the isobars fan out near
        the critical point.

        Joining them also closes the 360-375 C hole at 25 MPa, where the
        compressed-liquid table has no row and the superheat table has not started.
        """
        if P in self._merged_cache:
            return self._merged_cache[P]
        Sk, Hk, Tk = self.isobar(P, raw=True)
        lb = self.liquid_branch(P)
        if lb is not None:
            keep = lb[2] < Tk.min() - 1
            Tk = np.r_[lb[2][keep], Tk]
            Sk = np.r_[lb[0][keep], Sk]
            Hk = np.r_[lb[1][keep], Hk]
        fine = np.linspace(Tk.min(), Tk.max(), 600)
        out = (PchipInterpolator(Tk, Sk)(fine), PchipInterpolator(Tk, Hk)(fine), fine)
        self._merged_cache[P] = out
        return out

    def wet_isobar(self, P, Hmin):
        """The wet-region part of an isobar on an H-S chart: a straight line.

        Inside the dome an isobar is also an isotherm, so dH = T dS with T constant.

        Returns None above the critical pressure, where there is no wet region at all.
        Without that guard `Tsat_of_P` extrapolates past the end of the saturation
        table -- it reports a saturation temperature of 2025 C at 60 MPa -- and the
        dome interpolators are then evaluated far outside their range, which produced
        a segment reaching H = 5e9 kJ/kg.
        """
        if P >= PC:
            return None
        T = float(self.Tsat_of_P(P))
        Sl, Sv = float(self.dome["Sl"](T)), float(self.dome["Sv"](T))
        Hl, Hv = float(self.dome["Hl"](T)), float(self.dome["Hv"](T))
        if Hv <= Hmin:
            return None
        x = np.array([max(0.0, (Hmin - Hl) / (Hv - Hl)), 1.0])
        return Sl + x * (Sv - Sl), Hl + x * (Hv - Hl)

    # --- isotherms ---------------------------------------------------------
    def extend_to_low_P(self, P, S, H, S_end):
        """Continue an isotherm below the lowest pressure the table carries, 0.01 MPa.

        Needed because the table stops there, which on the Mollier frame is short of
        the right-hand edge for every isotherm below about 350 C. The continuation is
        not a fit: at these pressures steam is nearly an ideal gas, and for an ideal
        gas at fixed temperature S(P) = S(P0) - (R/M) ln(P/P0) and H(P) = H(P0),
        exactly. The small residual in H is kept by extrapolating it linearly in P
        from the two lowest tabulated isobars; the residual in S is smaller than the
        line width.
        """
        if S_end <= S[0]:
            return P, S, H
        dHdP = (H[1] - H[0]) / (P[1] - P[0])
        P_end = P[0] * np.exp(-(S_end - S[0]) / R_M)
        P_new = np.geomspace(P_end, P[0], 60)[:-1]
        return (np.r_[P_new, P],
                np.r_[S[0] - R_M * np.log(P_new / P[0]), S],
                np.r_[H[0] + dHdP * (P_new - P[0]), H])

    def isotherm(self, T, extend_to=None):
        """(S, H, P) along one isotherm, from the dome out to low pressure.

        Built by evaluating every isobar at T, not by reading a row of the table. The
        table lists the intermediate temperature levels only at higher pressures --
        350 C first appears at 0.5 MPa, 450 C at 2.5, 550 C at 6.0, 650 C at 9.0 --
        so a row-based isotherm at those temperatures stops in the middle of the
        chart. Interpolating along each isobar instead carries every isotherm down
        to 0.01 MPa.
        """
        pts = [(P, float(np.interp(T, Tl, S)), float(np.interp(T, Tl, H)))
               for P in self.ALL_P
               for S, H, Tl in [self.isobar(P)]
               if Tl.min() <= T <= Tl.max()]
        P, S, H = (np.array(a) for a in zip(*sorted(pts)))

        if T < TC:                             # anchor it on the dome
            P = np.r_[P, float(self.dome["P_MPa"](T))]
            S = np.r_[S, float(self.dome["Sv"](T))]
            H = np.r_[H, float(self.dome["Hv"](T))]
            o = np.argsort(P)
            P, S, H = P[o], S[o], H[o]

        lp = np.log10(P)                       # smooth in log P: an isotherm is
        fine = np.linspace(lp.min(), lp.max(), 400)   # nearly straight there
        S, H, P = (PchipInterpolator(lp, S)(fine), PchipInterpolator(lp, H)(fine),
                   10 ** fine)
        if extend_to is not None:
            P, S, H = self.extend_to_low_P(P, S, H, extend_to)
        return S, H, P

    # --- lines of constant enthalpy ----------------------------------------
    def _dome_crossings(self, H0):
        """(P, S, T) where the line of constant enthalpy H0 meets the saturation curve.

        It can meet it more than once: the saturated-vapor enthalpy rises to a maximum
        of 2804 kJ/kg at 235 C and falls again, so H0 = 2600 leaves the dome at 54 C
        and re-enters it at 344 C.
        """
        T = np.linspace(self.sat.T_C.min(), TC, 4000)
        out = []
        for Hkey, Skey in (("Hv", "Sv"), ("Hl", "Sl")):
            f = self.dome[Hkey](T) - H0
            for i in np.nonzero(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]:
                t = T[i] - f[i] * (T[i + 1] - T[i]) / (f[i + 1] - f[i])
                out.append((float(self.dome["P_MPa"](t)), float(self.dome[Skey](t)),
                            float(t)))
        return out

    def isenthalp(self, H0, isobars=None):
        """Segments of the line of constant enthalpy H0, as [(S, T, where), ...].

        Two things make this the fiddliest family on the chart, and both are about
        where a segment STOPS.

        A line of constant enthalpy runs through the two-phase region and out of it,
        and the two parts are computed from different tables -- the lever rule inside,
        the isobars outside -- so each part is closed on the exact point where the line
        meets the saturation curve. Without that they stop a few degrees short of it
        and the line visibly fails to cross the phase envelope.

        And the part outside the dome is sampled one point per isobar, so a pressure
        range where no tabulated isobar carries H0 breaks the line in two. At
        2600 kJ/kg the line exists at 0.01 MPa and again from 17.5 MPa upward, and
        joining those into one polyline draws a chord clear across the chart. The
        points are therefore grouped into runs of CONSECUTIVE isobars, and each run is
        a separate segment.
        """
        isobars = tuple(sorted(self.ALL_P if isobars is None else isobars))
        crossings = self._dome_crossings(H0)
        segments = []

        # --- inside the dome: the quality that gives H0 at each temperature ---
        T = np.linspace(self.sat.T_C.min(), TC, 300)
        inside = (self.dome["Hl"](T) <= H0) & (H0 <= self.dome["Hv"](T))
        if inside.sum() > 1:
            t = T[inside]
            x = (H0 - self.dome["Hl"](t)) / self.dome["dH"](t)
            S = self.dome["Sl"](t) + x * self.dome["dS"](t)
            for _, sc, tc in crossings:          # close it onto the saturation curve
                if t.min() - 3 <= tc <= t.max() + 3:
                    t, S = np.r_[t, tc], np.r_[S, sc]
            o = np.argsort(t)
            segments.append((S[o], t[o], "wet"))

        # --- outside the dome: one point per isobar, grouped into consecutive runs ---
        def point_on(P):
            """(S, T) where the isobar at P has enthalpy H0, on either side of the dome.

            BOTH sides have to be searched. A line of constant enthalpy that leaves
            the dome on the LIQUID side continues as a compressed liquid, and at
            1800 kJ/kg that is exactly where it spends the 19-30 MPa stretch:
            invisible to the superheat table, so the line came out in two pieces with
            a gap between them. Above the critical pressure the tabulated branch is
            preferred over the pressure-interpolated liquid one, which is the weaker
            of the two there.
            """
            pieces = ((self.merged_isobar(P),) if P >= PC
                      else (self.liquid_branch(P), self.isobar(P)))
            for piece in pieces:
                if piece is None:
                    continue
                S, H, T = piece
                if H.min() <= H0 <= H.max():
                    o = np.argsort(H)
                    return (float(np.interp(H0, H[o], S[o])),
                            float(np.interp(H0, H[o], T[o])))
            return None

        runs, run = [], []
        for P in isobars:
            pt = point_on(P)
            if pt is not None:
                run.append((P, pt[0], pt[1]))
            elif run:
                runs.append(run)
                run = []
        if run:
            runs.append(run)

        # a run that ends next to the dome -- with no tabulated isobar between it and
        # the crossing -- gets the crossing itself, so the line arrives ON the
        # saturation curve and meets the segment inside it
        for pc, sc, tc in crossings:
            for r in runs:
                if pc < r[0][0] and not any(pc < P < r[0][0] for P in isobars):
                    r.insert(0, (pc, sc, tc))
                elif pc > r[-1][0] and not any(r[-1][0] < P < pc for P in isobars):
                    r.append((pc, sc, tc))

        for r in runs:                           # already ordered in pressure
            if len(r) > 1:
                _, S, T = (np.array(a) for a in zip(*r))
                segments.append((S, T, "superheat"))
        return segments


# --- the charts ------------------------------------------------------------
#
# Two tiers per family: a labeled major set, and an unlabeled minor set between them.
# Every set is a parameter -- together with `charts.LW`, this is where the chart's
# density is decided.

# Below 0.01 MPa the superheat table stops, but the chart does not: without these the
# whole lower-right corner is empty and the constant-quality lines stop in mid air.
# They end on the triple-point isobar, the lowest pressure at which steam has a
# saturated liquid at all.
MOLLIER_P = (P_TRIPLE, 0.001, 0.002, 0.005,
             0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0)
MOLLIER_T = (100, 200, 300, 400, 500, 600, 700, 800)
MOLLIER_T_MINOR = (150, 250, 350, 450, 550, 650, 750)
MOLLIER_X = (0.75, 0.80, 0.85, 0.90, 0.95)
MOLLIER_X_MINOR = (0.725, 0.775, 0.825, 0.875, 0.925, 0.975)

TS_P = (0.01, 0.1, 1.0, 5.0, 10.0, 20.0, 40.0, 60.0)
TS_H = (400, 800, 1200, 1600, 2000, 2400, 2800, 3200, 3600, 4000)
TS_H_MINOR = (200, 600, 1000, 1400, 1800, 2200, 2600, 3000, 3400, 3800)
TS_X = (0.2, 0.4, 0.6, 0.8)
TS_X_MINOR = (0.1, 0.3, 0.5, 0.7, 0.9)


def minor_isobars(st, major):
    """Every tabulated isobar the caller did not list as major."""
    return tuple(P for P in st.ALL_P if P not in major)


def mollier(ax, st, *, S_lim=(4.5, 9.5), H_lim=(2000, 4200),
            isobars=MOLLIER_P, isobars_minor=None,
            isotherms=MOLLIER_T, isotherms_minor=MOLLIER_T_MINOR,
            qualities=MOLLIER_X, qualities_minor=MOLLIER_X_MINOR,
            label_H=4080, T_key=True, T_key_title=True, grid=True,
            S_grid=0.5, H_grid=100, grid_minors=5, x_label_between=(0.01, 0.05)):
    """Draw the enthalpy-entropy (Mollier) chart, Figure 3.3-1(a), on `ax`."""
    if isobars_minor is None:
        isobars_minor = minor_isobars(st, isobars)
    dome, sat, solid = st.dome, st.sat, st.solid
    ax.set_xlim(*S_lim)          # limits first: labels are rotated in display space
    ax.set_ylim(*H_lim)
    S_edge = S_lim[1]

    # --- minor lines first, so the labeled majors sit on top of them ---
    for Tv in isotherms_minor:
        S, H, _ = st.isotherm(Tv, extend_to=S_edge)
        ax.plot(S, H, "-", color=GRAY_MINOR, lw=LW["T_minor"], zorder=1)
    for x in qualities_minor:
        S, H = st.quality_line(x)
        ax.plot(S, H, "--", color=GRAY_MINOR, lw=LW["x_minor"], zorder=1,
                dashes=(4, 2))
    for P in isobars_minor:
        S, H, _ = st.isobar(P)
        ax.plot(S, H, "-", color="0.25", lw=LW["P_minor"], zorder=3)
        seg = st.wet_isobar(P, H_lim[0])
        if seg:
            ax.plot(*seg, "-", color="0.25", lw=LW["P_minor"], zorder=3)

    # --- major isotherms, carried out to the right-hand edge ---
    S_key = S_lim[1] - 0.42            # x of the temperature key column
    for Tv in isotherms:
        S, H, _ = st.isotherm(Tv, extend_to=S_edge)
        ax.plot(S, H, "-", color=GRAY, lw=LW["T"], zorder=2)
        if T_key and S.max() >= S_key:
            # the isotherm is stored from high pressure to low, so entropy runs
            # DOWN the array: flip it for np.interp, which needs increasing x
            ax.text(S_key, float(np.interp(S_key, S[::-1], H[::-1])),
                    f"{Tv}$^\\circ$", fontsize=6.5, color=GRAY, ha="left",
                    va="center", zorder=6,
                    bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))
        else:
            label_end(ax, S, H, f"{Tv}$^\\circ$C", color=GRAY)
    if T_key and T_key_title:      # off for small panels: it lands in the margin
        ax.text(1.012, 0.5, "temperature, $^\\circ$C", transform=ax.transAxes,
                rotation=90, ha="left", va="center", fontsize=7.5, color=GRAY)

    # --- major constant-quality lines ---
    # A point of a constant-quality line at temperature T lies exactly ON the isobar
    # whose saturation temperature is T -- inside the dome the isobar IS the T tie
    # line. So labeling the whole family at one temperature puts every label on one
    # isobar, and labeling at the mean of two pressures runs them up the band between
    # those two isobars, which is where the ASME chart puts its moisture labels: out
    # of the crowded left-hand fan and clear of the isobar labels.
    T_xlab = float(st.Tsat_of_P(np.sqrt(x_label_between[0] * x_label_between[1])))
    for x in qualities:
        S, H = st.quality_line(x)
        ax.plot(S, H, "--", color=GRAY, lw=LW["x"], zorder=2, dashes=(4, 2))
        H_lab = float(dome["Hl"](T_xlab) + x * dome["dH"](T_xlab))
        label_at(ax, S, H, f"$x$ = {x:g}", along="y", value=H_lab, color=GRAY,
                 pick="max_x")

    # --- major isobars, labeled in a fan across the top ---
    low = [P for P in isobars if P < st.ALL_P[0]]      # the sub-0.01 MPa group
    for P in isobars:
        S, H, _ = st.isobar(P)
        ax.plot(S, H, "-", color="k", lw=LW["P"], zorder=4)
        seg = st.wet_isobar(P, H_lim[0])
        if seg:
            ax.plot(*seg, "-", color="k", lw=LW["P"], zorder=4)
        # three places to put the label, in order of preference: the fan across the
        # top; on the line itself for an isobar that leaves through the RIGHT edge,
        # clear of the temperature key column; and, for the lowest pressures of all,
        # whose superheat branch starts too far right to reach either, on the straight
        # wet-region segment -- staggered in enthalpy so they do not collide
        if label_at(ax, S, H, fmt_P(P), along="y", value=label_H):
            continue
        if label_at(ax, S, H, fmt_P(P), along="x", value=S_key - 0.6):
            continue
        if seg and P in low:
            label_at(ax, seg[0], seg[1], fmt_P(P), along="y",
                     value=2080 + 180 * low.index(P))

    # --- the saturation line, heaviest ---
    T = np.linspace(sat.T_C.min(), TC, 500)
    Ssat, Hsat = dome["Sv"](T), dome["Hv"](T)
    ax.plot(Ssat, Hsat, "-", color="k", lw=LW["sat"], zorder=5)
    # Past the triple point the vapor is saturated with respect to ICE, not liquid,
    # and A.III tabulates that too. Drawing it carries the line off the edge of the
    # frame instead of letting it stop in mid-air at the triple point.
    sub = solid.sort_values("T_C")
    ax.plot(sub.Sv, sub.Hv, "-", color="k", lw=LW["sat"], zorder=5)
    label_at(ax, Ssat, Hsat, "saturated vapor", along="x", value=7.6, size=7,
             offset=(0, -22))

    if grid:
        chart_grid(ax, S_grid, H_grid, grid_minors)
    ax.set_xlabel(r"entropy  $\hat S$  (kJ/kg K)")
    ax.set_ylabel(r"enthalpy  $\hat H$  (kJ/kg)")
    return ax


def temperature_entropy(ax, st, *, S_lim=(0, 10), T_lim=(0, 800),
                        isobars=TS_P, isobars_minor=None,
                        isenthalps=TS_H, isenthalps_minor=TS_H_MINOR,
                        qualities=TS_X, qualities_minor=TS_X_MINOR,
                        grid=True, liquid_isobars=True, mark_critical=False,
                        S_grid=0.5, T_grid=50, grid_minors=5,
                        H_label_T=60, H_family_label=(3275, 6.0)):
    """Draw the temperature-entropy chart (Figure 3.3-1b) on `ax`.

    This is also the frame `c05uf001` draws the Rankine cycle of Illustration 5.2-1
    on: call it, then plot the cycle's four states over the top.
    """
    if isobars_minor is None:
        isobars_minor = minor_isobars(st, isobars)
    dome, sat = st.dome, st.sat
    ax.set_xlim(*S_lim)
    ax.set_ylim(*T_lim)

    def draw_isobar(P, lw, color, wet=True):
        """Draw one isobar, and return the curve its label should sit on.

        An isobar is drawn differently on the two sides of the critical pressure,
        because it IS a different object there. Above P_c there is no phase change:
        the compressed-liquid table and the superheat table describe one continuous
        line, and it has to be drawn as one. Below P_c the line is genuinely three
        pieces -- liquid, a horizontal tie line across the dome, vapor -- and the
        three have to MEET, exactly, on the saturation curve.
        """
        lb = st.liquid_branch(P) if liquid_isobars else None

        if P >= PC:
            S, _, T = st.merged_isobar(P)
            ax.plot(S, T, "-", color=color, lw=lw, zorder=4)
            return S, T

        Ts = float(st.Tsat_of_P(P))
        Sl, Sv = float(dome["Sl"](Ts)), float(dome["Sv"](Ts))
        S, T = st.isobar(P)[0], st.isobar(P)[2]      # the vapor branch
        ax.plot(S, T, "-", color=color, lw=lw, zorder=4)
        # The tie line is horizontal and spans the whole dome, so drawing it for all
        # 36 isobars turns the dome into a ladder and buries the quality lines. Major
        # isobars only -- one is enough to make the point that boiling at constant
        # pressure is boiling at constant temperature.
        if wet:
            ax.plot([Sl, Sv], [Ts, Ts], "-", color=color, lw=lw, zorder=4)
        if lb is not None:                           # already closed on the dome
            fine = np.linspace(lb[2].min(), lb[2].max(), 200)
            ax.plot(PchipInterpolator(lb[2], lb[0])(fine), fine, "-", color=color,
                    lw=lw, zorder=4)
        return S, T

    # --- minor lines first ---
    for H0 in isenthalps_minor:
        for S, T, _ in st.isenthalp(H0):
            ax.plot(S, T, "-", color=GRAY_MINOR, lw=LW["H_minor"], zorder=1)
    # A 0.3 pt line in GRAY_MINOR antialiases its own dashes away: the minor quality
    # lines came out looking SOLID, which is the one thing they must not look like
    # here, because solid inside the dome means constant enthalpy. Two changes fix it
    # -- the full GRAY, and a tighter dash than the major tier's (4, 2), since a fine
    # line needs a short dash to read as broken. The tiers are still told apart by
    # weight, 0.7 against 0.3.
    for x in qualities_minor:
        Tq = np.linspace(sat.T_C.min(), TC, 400)
        S, _ = st.quality_line(x, Tq)
        ax.plot(S, Tq, "--", color=GRAY, lw=LW["x_minor"], zorder=1,
                dashes=(2.6, 2.2))
    for P in isobars_minor:
        draw_isobar(P, LW["P_minor"], "0.25", wet=False)

    # --- major constant-enthalpy lines, labeled once each ---
    for H0 in isenthalps:
        segs = st.isenthalp(H0)
        for S, T, where in segs:
            ax.plot(S, T, "-", color=GRAY, lw=LW["H"], zorder=2)
        # Label once, on whichever segment reaches furthest right.
        S, T, where = max(segs, key=lambda seg: float(np.nanmax(seg[0])))
        if S[0] > S[-1]:
            S, T = S[::-1], T[::-1]
        if where == "wet":
            # A line that never leaves the dome runs DOWN to the bottom of the frame,
            # so its right-hand end is ON the temperature axis and label_end prints
            # the number over the axis labels. Put it on the line instead, at one
            # temperature for the whole family.
            if not label_at(ax, S, T, f"{H0}", along="y", value=H_label_T,
                            color=GRAY):
                label_at(ax, S, T, f"{H0}", along="y",
                         value=0.5 * (T.min() + T.max()), color=GRAY)
        else:
            label_end(ax, S, T, f"{H0}", color=GRAY, end="last")

    # One member of the family is named as well as numbered, the way the saturation
    # line is named on the Mollier chart. `H_family_label` is (which line, at what
    # entropy); the height is not free, because the text sits ON the curve and
    # rotates to it. Watch two things when moving it: the text is CENTERED and about
    # 1.3 kJ/(kg K) wide, so keep that much clear of the segment's ends and of the
    # line's own number at its right-hand end; and `label_at` draws NOTHING,
    # silently, if the curve never reaches the entropy asked for.
    if H_family_label and isenthalps:
        S, T, _ = max(st.isenthalp(H_family_label[0]),
                      key=lambda seg: float(np.nanmax(seg[0])))
        if S[0] > S[-1]:
            S, T = S[::-1], T[::-1]
        label_at(ax, S, T, "constant enthalpy, kJ/kg", along="x",
                 value=H_family_label[1], color=GRAY, size=7)

    # --- major constant-quality lines ---
    Tq = np.linspace(sat.T_C.min(), TC, 400)
    for x in qualities:
        S, _ = st.quality_line(x, Tq)
        ax.plot(S, Tq, "--", color=GRAY, lw=LW["x"], zorder=2, dashes=(4, 2))
        label_at(ax, S, Tq, f"{x:g}", along="y", value=150, color=GRAY)

    # --- major isobars ---
    for P in isobars:
        S, T = draw_isobar(P, LW["P"], "k")
        # as on the Mollier chart: the fan across the top first, then on the line
        # itself for an isobar that leaves through the RIGHT edge instead
        if label_at(ax, S, T, fmt_P(P), along="y", value=T_lim[1] - 60):
            continue
        if label_at(ax, S, T, fmt_P(P), along="x", value=S_lim[1] - 0.7):
            continue
        label_end(ax, S, T, fmt_P(P), ha="right")

    # --- the dome, heaviest ---
    T = np.linspace(sat.T_C.min(), TC, 500)
    ax.plot(dome["Sl"](T), T, "-", color="k", lw=LW["sat"], zorder=5)
    ax.plot(dome["Sv"](T), T, "-", color="k", lw=LW["sat"], zorder=5)
    if mark_critical:
        ax.plot([SC], [TC], "o", color="k", ms=3.5, zorder=6)
        ax.annotate("critical point", xy=(SC, TC), xytext=(SC - 1.4, TC + 55),
                    fontsize=7, ha="center", zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.5))

    if grid:
        chart_grid(ax, S_grid, T_grid, grid_minors)
    ax.set_xlabel(r"entropy  $\hat S$  (kJ/kg K)")
    ax.set_ylabel(r"temperature  $T$  ($^\circ$C)")
    return ax


__all__ = ["SteamTables", "mollier", "temperature_entropy", "minor_isobars",
           "TC", "PC", "HC", "SC", "R_M", "P_TRIPLE",
           "MOLLIER_P", "MOLLIER_T", "MOLLIER_T_MINOR", "MOLLIER_X",
           "MOLLIER_X_MINOR", "TS_P", "TS_H", "TS_H_MINOR", "TS_X", "TS_X_MINOR"]
