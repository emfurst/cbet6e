"""ph_chart_blend -- pressure-enthalpy charts for a refrigerant BLEND.

`thermo.ph_chart` draws a pure fluid, which has one saturation temperature at each
pressure. A zeotropic blend such as R-448A does not: at fixed pressure it starts to
boil at its BUBBLE point and finishes at its DEW point, several kelvin higher, and the
difference (the temperature glide) is the feature its chart exists to show. Every
isotherm therefore SLOPES across the two-phase region instead of running flat.

This module supplies a `BlendChartFluid` that `thermo.ph_chart.ph_chart` draws exactly
as it draws a pure fluid -- same line weights, labels and grid -- because the fluid
object provides its own curve families:

    from thermo import PengRobinson
    from thermo.pr_mixture import PRMixture
    from thermo.ph_chart import ph_chart
    from thermo.ph_chart_blend import BlendChartFluid

    mix = PRMixture([PengRobinson(...), ...], kij=K)
    r448a = BlendChartFluid(mix, z, M)            # z mole fractions, M in g/mol
    ph_chart(ax, r448a, H_lim=(50, 650), P_lim=(0.01, 20),
             isotherms=range(-80, 221, 20), T_celsius=True, ...)

THE THREE PIECES OF NUMERICS, each commented at the point of use:

1. **The phase envelope is TRACED, not solved point by point.** The bubble- and
   dew-point drivers in `thermo.phi_phi` fall to the trivial root (K = 1 for every
   species) above about half the critical pressure. `envelope()` follows the curve by
   continuation instead (Michelsen, Fluid Phase Equilib. 4, 1 (1980)): Newton on
   (ln K, ln T, ln P), stepping in whichever variable changes fastest. It runs up the
   bubble branch, through the mixture critical point, and down the dew branch in one
   pass.
2. **The dome interior is one grid.** Each row is a pressure; along a row the vapor
   fraction runs from 0 (bubble point) to 1 (dew point). Isotherms, isentropes,
   quality lines and isochores inside the dome are all interpolated along rows of
   that one grid, so they agree with one another and with the dome itself.
3. **The critical point closes the dome.** The rows stop a fraction of a percent
   below it, where the two phases are too alike to converge; every family that
   reaches the top is closed onto the critical point, as `ph_chart` closes a pure
   fluid's dome.

DATUM. The chart uses the IIR convention of refrigeration practice, so it can be read
against manufacturers' tables: **saturated (bubble-point) liquid at 0 deg C has
H = 200 kJ/kg and S = 1.00 kJ/(kg K).** The pure-fluid charts in Chapter 3 put their
zero at the triple point instead; a blend has no single triple point.

Eric M. Furst
November 2026
"""
from __future__ import annotations

import numpy as np
from scipy.constants import R
from scipy.optimize import brentq

TREF = 298.15          # K, where the ideal-gas integrals start (arbitrary, cancels)
PREF = 1e5             # Pa
T_LOW, T_HIGH = 140.0, 900.0     # K, the search bracket for single-phase states


class BlendChartFluid:
    """A fixed-composition blend on a Peng-Robinson mixture model, per kilogram.

    Parameters
    ----------
    mix : PRMixture
        The mixture equation of state. Each component must carry its ideal-gas
        heat capacity `cp = (a, b, c, d)`, Cp* = a + bT + cT^2 + dT^3 in J/(mol K).
    z : sequence of float
        Mole fractions of the blend.
    M : sequence of float
        Molar masses of the components, g/mol.
    datum : (T, H, S)
        Bubble-point liquid at temperature T (K) is assigned H (kJ/kg) and
        S (kJ/(kg K)). The default is the IIR convention.
    P_floor : float
        Lowest pressure (Pa) the two-phase grid has to cover. Set it a little below
        the bottom of the chart.
    """

    def __init__(self, mix, z, M, datum=(273.15, 200.0, 1.00), P_floor=8e3,
                 n_rows=160, n_beta=61):
        self.mix = mix
        self.z = np.asarray(z, float)
        self.Mi = np.asarray(M, float)
        self.n = len(self.z)
        self.M = float(self.z @ self.Mi)                 # g/mol of the blend
        self.cp = np.array([c.cp for c in mix.components], float)
        if self.cp.shape != (self.n, 4):
            raise ValueError("every component needs cp = (a, b, c, d)")
        self.H0 = self.S0 = 0.0
        self.envelope(P_start=0.5 * P_floor)
        self.two_phase_grid(P_floor, n_rows, n_beta)
        # Set the datum last: everything above was computed with H0 = S0 = 0, and the
        # grid is shifted afterwards so that the two agree.
        Td, Hd, Sd = datum
        P = self.P_bubble(Td)
        h, s, _ = self.feed(Td, P, "liquid")
        self.H0, self.S0 = h - Hd, s - Sd
        for k, off in (("H", self.H0), ("S", self.S0)):
            self.grid[k] = self.grid[k] - off
        self.Hc -= self.H0
        self.Sc -= self.S0

    # --- ideal gas, per mole of each component --------------------------------
    def _H_ig(self, T):
        a, b, c, d = self.cp.T
        f = lambda t: a * t + b * t**2 / 2 + c * t**3 / 3 + d * t**4 / 4
        return f(T) - f(TREF)

    def _S_ig(self, T, P):
        a, b, c, d = self.cp.T
        f = lambda t: a * np.log(t) + b * t + c * t**2 / 2 + d * t**3 / 3
        return f(T) - f(TREF) - R * np.log(P / PREF)

    # --- one phase ------------------------------------------------------------
    def phase_molar(self, x, T, P, phase):
        """(H J/mol, S J/(mol K), V m^3/mol) of one phase of composition x.

        The departure functions in `PRMixture` are measured from the ideal-gas
        MIXTURE, which already holds the ideal entropy of mixing, so that term is
        added here, once.
        """
        x = np.asarray(x, float)
        xs = x[x > 0]
        H = x @ self._H_ig(T) + self.mix.departure_H(x, T, P, phase)
        S = (x @ self._S_ig(T, P) - R * np.sum(xs * np.log(xs))
             + self.mix.departure_S(x, T, P, phase))
        return H, S, self.mix.molar_volume(x, T, P, phase)

    def feed(self, T, P, phase):
        """(H kJ/kg, S kJ/(kg K), V m^3/kg) of the blend as ONE phase of its own
        composition. J/mol divided by g/mol is J/g, which is kJ/kg."""
        H, S, V = self.phase_molar(self.z, T, P, phase)
        return H / self.M - self.H0, S / self.M - self.S0, V * 1e3 / self.M

    def single(self, T, P):
        """Single-phase (H, S, V) at (T, P), taking the root of lower Gibbs energy.

        Only for states known to lie outside the dome: above the critical
        temperature, or on the correct side of the bubble or dew curve.
        """
        # A root with Z < B has a negative free volume and no physical meaning; left
        # in, its log(Z - B) is NaN and the Gibbs-energy comparison silently fails.
        B = self.mix.b_mix(self.z) * P / (R * T)
        zs = [Z for Z in self.mix.compressibility(self.z, T, P) if Z > B]
        phase = "vapor"
        if len(zs) > 1:
            gL, gV = (np.sum(self.z * self.mix.ln_phi(self.z, T, P, ph))
                      for ph in ("liquid", "vapor"))
            phase = "liquid" if gL < gV else "vapor"
        return self.feed(T, P, phase)

    # --- 1. the phase envelope ------------------------------------------------
    def _wilson_K(self, T, P):
        return np.array([c.Pc / P * np.exp(5.373 * (1 + c.omega) * (1 - c.Tc / T))
                         for c in self.mix.components])

    def _F_env(self, X, s, S, z_phase):
        """Incipient-phase equations at vapor fraction 0, plus the specification."""
        n = self.n
        T, P = np.exp(X[n]), np.exp(X[n + 1])
        y = np.exp(X[:n]) * self.z
        y_phase = "vapor" if z_phase == "liquid" else "liquid"
        F = np.empty(n + 2)
        F[:n] = (X[:n] + self.mix.ln_phi(y, T, P, y_phase)
                 - self.mix.ln_phi(self.z, T, P, z_phase))
        F[n] = y.sum() - 1.0
        F[n + 1] = X[s] - S
        return F

    @staticmethod
    def _jacobian(F, X, *args):
        F0 = F(X, *args)
        J = np.empty((len(F0), len(X)))
        for j in range(len(X)):
            h = 1e-7 * max(1.0, abs(X[j]))
            dX = np.zeros_like(X)
            dX[j] = h
            J[:, j] = (F(X + dX, *args) - F0) / h
        return F0, J

    def _newton_env(self, X, s, S, z_phase, it=30):
        for _ in range(it):
            F, J = self._jacobian(self._F_env, X, s, S, z_phase)
            dx = np.linalg.solve(J, -F)
            X = X + dx
            if np.max(np.abs(dx)) < 1e-10:
                return X, J
        raise RuntimeError("envelope Newton did not converge")

    def envelope(self, P_start=4e3, max_pts=3000):
        """Trace the bubble curve up, through the critical point, and the dew curve down.

        The equations are those of an incipient phase y = K z in equilibrium with the
        feed z. Below the critical point the feed is the liquid and y the first bubble
        of vapor. Past it, ln K of every species changes sign and the roles swap: the
        feed is the vapor and y the first drop of dew. The trace watches the most
        volatile species' ln K and swaps the two cubic roots when it crosses zero.
        That crossing is also how the critical point is located.
        """
        with np.errstate(over="ignore", invalid="ignore"):   # rejected trial steps
            return self._trace(P_start, max_pts)

    def _trace(self, P_start, max_pts):
        n = self.n
        T0 = brentq(lambda T: (self._wilson_K(T, P_start) * self.z).sum() - 1, 80, 600)
        X = np.concatenate([np.log(self._wilson_K(T0, P_start)),
                            [np.log(T0), np.log(P_start)]])
        s = n + 1                                        # first specify ln P
        X, J = self._newton_env(X, s, X[s], "liquid")
        light = int(np.argmax(X[:n]))                    # the most volatile species
        pts, z_phase, ds = [X.copy()], "liquid", 0.05
        while len(pts) < max_pts:
            # dX/dS from the converged Jacobian; specify the fastest-moving variable
            e = np.zeros(n + 2)
            e[-1] = 1.0
            dXdS = np.linalg.solve(J, e)
            k = int(np.argmax(np.abs(dXdS)))
            if k != s:
                dXdS, s = dXdS / dXdS[k], k
            sgn = np.sign((pts[-1] - pts[-2])[s]) if len(pts) > 1 else 1.0
            Xg = X + dXdS * (sgn or 1.0) * ds
            phase = z_phase
            if np.sign(Xg[light]) != np.sign(X[light]):
                phase = "vapor" if z_phase == "liquid" else "liquid"
            try:
                Xn, Jn = self._newton_env(Xg, s, Xg[s], phase)
            except (RuntimeError, np.linalg.LinAlgError, ValueError):
                ds *= 0.5
                if ds < 1e-5:
                    raise
                continue
            if np.sign(Xn[light]) != np.sign(X[light]):
                z_phase = "vapor" if z_phase == "liquid" else "liquid"
            X, J = Xn, Jn
            pts.append(X.copy())
            # Do NOT shorten the step near the critical point. At ln K = 0 the trivial
            # solution K = 1 satisfies the equations at every (T, P), the Jacobian is
            # singular, and a short step lands in it. Step OVER it, and locate it
            # afterwards by interpolation.
            ds = min(1.3 * ds, 0.1)
            if z_phase == "vapor" and np.exp(X[n + 1]) < P_start:
                break
        pts = np.array(pts)
        self.env_lnK = pts[:, :n]
        self.env_T, self.env_P = np.exp(pts[:, n]), np.exp(pts[:, n + 1])
        self.env_bubble = self.env_lnK[:, light] > 0
        cross = np.nonzero(np.diff(np.sign(self.env_lnK[:, light])))[0]
        if len(cross) != 1:
            raise RuntimeError("the envelope trace did not pass the critical point once")
        # A cubic through the two traced points on each side, evaluated at ln K = 0.
        # A straight line through the nearest pair alone puts Tc 0.2 K and Pc 0.4 %
        # off, depending on where the steps happened to fall.
        i = int(cross[0])
        k = slice(max(i - 1, 0), i + 3)
        lnK = self.env_lnK[k, light]
        self.Tc = float(np.polyval(np.polyfit(lnK, self.env_T[k], len(lnK) - 1), 0.0))
        self.Pc = float(np.exp(np.polyval(np.polyfit(lnK, np.log(self.env_P[k]),
                                                     len(lnK) - 1), 0.0)))
        # The maxima of the traced points can fall short of the true ones, which lie
        # between steps. Neither can be below the critical point.
        self.cricondentherm = max(float(self.env_T.max()), self.Tc)
        self.cricondenbar = max(float(self.env_P.max()), self.Pc)
        return self

    # --- 2. the dome interior -------------------------------------------------
    def _F_row(self, X, beta, lnP):
        """Flash at fixed P and vapor fraction beta: unknowns (ln K, ln T)."""
        n = self.n
        K, T, P = np.exp(X[:n]), np.exp(X[n]), np.exp(lnP)
        x = self.z / (1 + beta * (K - 1))
        y = K * x
        F = np.empty(n + 1)
        F[:n] = (X[:n] + self.mix.ln_phi(y, T, P, "vapor")
                 - self.mix.ln_phi(x, T, P, "liquid"))
        F[n] = y.sum() - x.sum()
        return F

    def _newton_row(self, X, beta, lnP, it=60):
        for _ in range(it):
            F, J = self._jacobian(self._F_row, X, beta, lnP)
            if np.max(np.abs(F)) < 1e-12:
                return X
            dx = np.linalg.solve(J, -F)
            # Damped: ln T may move 2 % and ln K 0.5 per step. Undamped, a step near
            # the critical point walks to the trivial root, which satisfies these
            # equations too and is the nearer one there.
            dx *= min(1.0, 0.5 / max(np.max(np.abs(dx[:-1])), 1e-30),
                      0.02 / max(abs(dx[-1]), 1e-30))
            X = X + dx
            if np.max(np.abs(dx)) < 1e-9:
                return X
        raise RuntimeError("two-phase Newton did not converge")

    def _bubble_start(self, lnP):
        """A converged bubble point at ln P, started from the traced envelope."""
        m = self.env_bubble
        o = np.argsort(self.env_P[m])
        lnPe = np.log(self.env_P[m][o])
        lnK = np.array([np.interp(lnP, lnPe, self.env_lnK[m][o][:, i])
                        for i in range(self.n)])
        T = np.interp(lnP, lnPe, self.env_T[m][o])
        X = self._newton_row(np.append(lnK, np.log(T)), 0.0, lnP)
        if abs(np.exp(X[-1]) - T) > 2.0 or np.max(np.abs(X[:-1])) < 1e-6:
            raise RuntimeError("bubble point fell to the trivial root")
        return X

    def two_phase_grid(self, P_floor, n_rows, n_beta, P_top=0.9995):
        """Rows at fixed P from P_floor to just below Pc, vapor fraction 0 to 1."""
        n = self.n
        u = np.linspace(0, 1, n_rows)
        # bunched toward Pc, where the dome's shoulder turns over
        lnP = np.log(P_floor) + np.log(self.Pc * P_top / P_floor) * (1 - (1 - u)**1.6)
        beta = np.linspace(0, 1, n_beta)
        G = {k: [] for k in ("T", "H", "S", "V", "q", "P")}
        for lp in lnP:
            P = np.exp(lp)
            try:
                X = self._bubble_start(lp)
                row = []
                for c, b in enumerate(beta):
                    if c:
                        X = self._newton_row(X, b, lp)
                    row.append(X)
            except RuntimeError:
                # The last fraction of a percent below Pc may fail to converge; it
                # is closed onto the critical point instead. Anything lower is a
                # real failure and is raised.
                if P < 0.99 * self.Pc:
                    raise
                break
            cols = {k: [] for k in ("T", "H", "S", "V", "q")}
            for b, X in zip(beta, row):
                K, T = np.exp(X[:n]), np.exp(X[n])
                x = self.z / (1 + b * (K - 1))
                y = K * x
                x, y = x / x.sum(), y / y.sum()
                HL, SL, VL = self.phase_molar(x, T, P, "liquid")
                HV, SV, VV = self.phase_molar(y, T, P, "vapor")
                cols["T"].append(T)
                cols["H"].append(((1 - b) * HL + b * HV) / self.M)
                cols["S"].append(((1 - b) * SL + b * SV) / self.M)
                cols["V"].append(((1 - b) * VL + b * VV) * 1e3 / self.M)
                cols["q"].append(b * float(y @ self.Mi) / self.M)   # MASS quality
            for k in cols:
                G[k].append(cols[k])
            G["P"].append(P)
        self.grid = {k: np.array(v) for k, v in G.items()}
        self.grid["beta"] = beta
        # At the mixture critical point the cubic has one root, so liquid and vapor
        # agree; average them anyway in case rounding separates the two.
        HL, SL, VL = self.phase_molar(self.z, self.Tc, self.Pc, "liquid")
        HV, SV, VV = self.phase_molar(self.z, self.Tc, self.Pc, "vapor")
        self.Hc, self.Sc = (HL + HV) / 2 / self.M, (SL + SV) / 2 / self.M
        self.Vc = (VL + VV) / 2 * 1e3 / self.M
        return self.grid

    # --- boundaries, from the grid's first and last columns -------------------
    def _edge(self, col):
        g = self.grid
        return g["T"][:, col], np.log(g["P"])

    def P_bubble(self, T):
        Tb, lnP = self._edge(0)
        return float(np.exp(np.interp(T, Tb, lnP)))

    def P_dew(self, T):
        Td, lnP = self._edge(-1)
        return float(np.exp(np.interp(T, Td, lnP)))

    def T_bubble(self, P):
        Tb, lnP = self._edge(0)
        return float(np.interp(np.log(P), lnP, Tb))

    def T_dew(self, P):
        Td, lnP = self._edge(-1)
        return float(np.interp(np.log(P), lnP, Td))

    def glide(self, P):
        """Dew minus bubble temperature at pressure P, K."""
        return self.T_dew(P) - self.T_bubble(P)

    def crit(self):
        """(P, H, S, V) at the mixture critical point, as `ChartFluid.crit` returns."""
        return self.Pc, self.Hc, self.Sc, self.Vc

    def _in_rows(self, P):
        return self.grid["P"][0] <= P <= self.grid["P"][-1]

    # --- 3. the curve families, in the form `ph_chart` draws ------------------
    def ph_dome(self):
        """(T, P_MPa, HL, HV, SL, SV), closed on the critical point."""
        g = self.grid
        ap = lambda a, v: np.append(a, v)
        return (ap(g["T"][:, 0], self.Tc), ap(g["P"], self.Pc) / 1e6,
                ap(g["H"][:, 0], self.Hc), ap(g["H"][:, -1], self.Hc),
                ap(g["S"][:, 0], self.Sc), ap(g["S"][:, -1], self.Sc))

    def _row_interp(self, key, value, out="H"):
        """Along every row where `key` brackets `value`: (out, P)."""
        g = self.grid
        res, P = [], []
        for r in range(len(g["P"])):
            a = g[key][r]
            if a[0] <= value <= a[-1]:
                res.append(np.interp(value, a, g[out][r]))
                P.append(g["P"][r])
        return np.array(res), np.array(P)

    def ph_isotherm(self, T, P_lim, n=200):
        """[(H, P_MPa), ...]: vapor, then across the dome, then liquid."""
        Plo, Phi = P_lim[0] * 1e6, P_lim[1] * 1e6
        g = self.grid
        if T < g["T"][-1, 0]:                            # crosses the gridded dome
            segs = []
            Pd, Pb = self.P_dew(T), self.P_bubble(T)
            # Bunched toward the dome on both sides: next to the critical point an
            # isotherm turns hardest just outside the dew and bubble curves, and
            # evenly log-spaced points (3 % apart) show there as a visible corner.
            near = np.geomspace(1e-5, 0.3, 40)
            if Pd > Plo:
                Pv = np.unique(np.concatenate([np.geomspace(Plo, Pd, n), Pd * (1 - near)]))
                Pv = Pv[Pv >= Plo]
                segs.append((np.array([self.feed(T, p, "vapor")[0] for p in Pv]), Pv))
            H2, P2 = self._row_interp("T", T)
            H2 = np.concatenate([[self.feed(T, Pd, "vapor")[0]], H2[::-1],
                                 [self.feed(T, Pb, "liquid")[0]]])
            P2 = np.concatenate([[Pd], P2[::-1], [Pb]])
            o = np.argsort(P2)
            segs.append((H2[o], P2[o]))
            if Pb < Phi:
                Pl = np.unique(np.concatenate([np.geomspace(max(Pb, Plo), Phi, n),
                                               Pb * (1 + near)]))
                Pl = Pl[(Pl >= Plo) & (Pl <= Phi)]
                segs.append((np.array([self.feed(T, p, "liquid")[0] for p in Pl]), Pl))
            return [(H, P / 1e6) for H, P in segs if len(H) > 1 and P[-1] > Plo]
        if T < self.cricondentherm + 0.05:
            raise ValueError(f"T = {T:.2f} K is inside the last {self.cricondentherm - g['T'][-1, 0]:.2f} K "
                             "below the cricondentherm, which the grid does not resolve")
        P = np.geomspace(Plo, Phi, 2 * n)
        return [(np.array([self.single(T, p)[0] for p in P]), P / 1e6)]

    def _solve_T(self, prop, target, P, lo, hi, phase):
        """Temperature in [lo, hi] at which property `prop` (1 = S, 2 = V) = target."""
        f = (lambda t: self.single(t, P)[prop] - target) if phase is None else \
            (lambda t: self.feed(t, P, phase)[prop] - target)
        try:
            if f(lo) * f(hi) > 0:
                return None
            return brentq(f, lo, hi, xtol=1e-6)
        except (ValueError, RuntimeError):
            return None

    def _family(self, prop, key, value, P_lim, n):
        """A constant-S or constant-V line: inside the dome from the grid rows,
        outside it by solving for T at each pressure."""
        H, P = [], []
        for p in np.geomspace(P_lim[0] * 1e6, P_lim[1] * 1e6, n):
            if self._in_rows(p):
                r = np.searchsorted(self.grid["P"], p)
                r = min(max(r, 1), len(self.grid["P"]) - 1)
                w = (np.log(p) - np.log(self.grid["P"][r - 1])) / \
                    np.log(self.grid["P"][r] / self.grid["P"][r - 1])
                lo_v = (1 - w) * self.grid[key][r - 1, 0] + w * self.grid[key][r, 0]
                hi_v = (1 - w) * self.grid[key][r - 1, -1] + w * self.grid[key][r, -1]
                Tb, Td = self.T_bubble(p), self.T_dew(p)
                if lo_v <= value <= hi_v:
                    continue                             # inside: from the rows below
                if value > hi_v:
                    T = self._solve_T(prop, value, p, Td + 1e-6, T_HIGH, "vapor")
                    ph = "vapor"
                else:
                    T = self._solve_T(prop, value, p, T_LOW, Tb - 1e-6, "liquid")
                    ph = "liquid"
                if T is not None:
                    H.append(self.feed(T, p, ph)[0]); P.append(p)
            else:
                lo = T_LOW if p > self.Pc else self.T_dew(p) + 1e-6
                T = self._solve_T(prop, value, p, lo, T_HIGH, None)
                if T is not None:
                    H.append(self.single(T, p)[0]); P.append(p)
        Hd, Pd = self._row_interp(key, value)
        # the exact crossings of the bubble and dew curves pin the corners
        Hx, Px = [], []
        for col in (0, -1):
            a = self.grid[key][:, col]
            s = np.nonzero(np.diff(np.sign(a - value)))[0]
            for i in s:
                f = (value - a[i]) / (a[i + 1] - a[i])
                Px.append(np.exp(np.log(self.grid["P"][i]) + f *
                                 np.log(self.grid["P"][i + 1] / self.grid["P"][i])))
                Hx.append(self.grid["H"][i, col] + f *
                          (self.grid["H"][i + 1, col] - self.grid["H"][i, col]))
        H = np.concatenate([H, Hd, Hx])
        P = np.concatenate([P, Pd, Px])
        o = np.argsort(P)
        H, P = H[o], P[o] / 1e6
        m = (P >= P_lim[0]) & (P <= P_lim[1])
        return H[m], P[m]

    def ph_isentrope(self, S, P_lim, n=220):
        return self._family(1, "S", S, P_lim, n)

    def ph_isochore(self, V, P_lim, n=260):
        return self._family(2, "V", V, P_lim, n)

    def ph_quality(self, x):
        """Constant MASS quality across the dome, closed on the critical point."""
        H, P = self._row_interp("q", x)
        return np.append(H, self.Hc), np.append(P, self.Pc) / 1e6


__all__ = ["BlendChartFluid"]
