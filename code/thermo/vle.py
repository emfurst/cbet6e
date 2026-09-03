"""vle -- low-pressure vapor-liquid equilibrium by the gamma-phi method, SIS Chapter 10.

The liquid is described by an activity coefficient model and the vapor is treated as
an ideal gas, so the equilibrium relation is SIS Eq. 10.2-1,

    y_i P = x_i gamma_i(x, T) Pvap_i(T)

and Sections 10.1 and 10.2 are **one procedure, not two**: Raoult's law is this with
gamma = 1. That is why `Ideal` is a model here rather than a separate code path -- the
book presents the ideal and nonideal calculations as different procedures, and they are
not. Section 10.3 is the same procedure again with both fugacities from an equation of
state; that is in `pr_mixture.py`, and see WHAT MAKES THE DIAGRAMS SHARED below.

    from thermo import PengRobinson
    from thermo.vle import Antoine, ClausiusClapeyron, GammaPhi, Ideal, pxy, txy
    from thermo.activity_models import VanLaar

    EA = Antoine(9.6830, 2842.2, -56.3209)          # Illustration 10.2-2, as printed
    BZ = Antoine(9.3171, 2810.5, -51.2586)
    vle = GammaPhi([EA, BZ], VanLaar(1.15, 0.92))
    P, y = vle.bubble_pressure([0.4, 0.6], 348.15)  # Pa, and the incipient vapor

**WHAT MAKES THE DIAGRAMS SHARED.** Every solver below has the same name and the same
signature as its counterpart on `PRMixture` -- `bubble_pressure(x, T) -> (P, y)`,
`flash(z, T, P) -> (beta, x, y)`, and so on. That is deliberate and it is the whole
reason `pxy` and `txy` at the bottom of this module take a *model* rather than being
methods: they never learn whether they were handed a `GammaPhi` or a `PRMixture`, so
Chapter 10's low-pressure diagrams and Section 10.3's high-pressure diagrams come out
of one generator. It is the same trick as `ActivityModel.gamma(x, T)` in Chapter 9.

 **The generators live here, in the physics tier, and not in a chart module** -- which
is the opposite of where `ph_chart.py` keeps `dome` and `isotherm`, and the difference is
deliberate. A saturation dome exists in order to be drawn. A VLE curve does not:
Chapter 10 prints these numbers as *tables* (Illustration 10.1-1's x, y and P; the
L = 1.0 -> 0.0 sweep of Illustration 10.1-5; Table 10.2-3, which Illustration 10.2-6
generates). Notebooks that want the numbers and draw nothing must not have to import
matplotlib to get them. The drawing is in `thermo/vle_chart.py`, which is lazy.

## Units

SI, as everywhere else in the package: **T in K and P in Pa**, and every solver returns
Pa. The vapor-pressure correlations are the one place the book's own units intrude --
Chapter 10 prints its constants for P in bar -- so `Antoine` and `ClausiusClapeyron`
**take the constants exactly as printed and convert internally**, the same bargain
`RegularSolution.from_table_9_6_1` strikes with cal and cc. Type in what the book shows;
get Pa back.

## The three correlation forms Chapter 10 actually uses

The chapter does not use one vapor-pressure equation, it uses three, and a module that
offers only "Antoine" cannot reproduce its illustrations:

1. `Antoine(A, B, C)`    -- ln(Pvap/bar) = A - B/(T + C)          Illustration 10.2-2
2. `ClausiusClapeyron`   -- ln(Pvap/bar) = A - dHvap/(RT)          Illustration 10.1-1
3. `Antoine.log10(A, B)` -- log10(Pvap/bar) = -A/T + B             Problem 10.1-1

 **Note the sign convention on C.** The book writes Illustration 10.2-2's denominator
as `T(K) - 56.3209`, so the constant passed here is **-56.3209**: `Antoine` is defined
with `T + C` throughout, which is the usual form, and negating at the call site is
visible where a silent internal flip would not be.

## When the book gives no constants: `psat_from_database`

Some figures need a vapor pressure the chapter never prints -- Figures 10.1-4 and 10.1-5
are at constant *pressure*, so they need hexane and triethylamine as functions of
temperature, and Section 10.1 gives only two numbers at 60 C.
`psat_from_database(name)` reads `code/data/pure_property.csv` and returns the right
correlation object, dispatching on that file's own `Eq` selector: 1 Wagner, 2 Riedel
(implicit, and solved), 3 Antoine. Checked against every species in the file that has
both a fit and a normal boiling point -- `Psat(Tb)` lands within a median 0.24 % of
1 atm for Wagner (n = 299), 0.22 % for Antoine (n = 170) and 0.79 % for Riedel (n = 49).

**Prefer the book's own numbers where the book prints them.** Where a printed value
and the database disagree, the printed one is what the surrounding text was computed
from, and switching sources silently moves numbers the prose quotes. See the
hexane/triethylamine notebook, which uses both on purpose and says which where.
"""
from __future__ import annotations

import numpy as np
from scipy import constants, optimize

from .activity_models import ActivityModel
from .pr_mixture import PRMixture

R = constants.R
BAR = constants.bar          # 1e5 Pa

__all__ = ["VaporPressure", "Antoine", "ClausiusClapeyron", "Wagner", "Riedel",
           "psat_from_database", "TabulatedPsat",
           "Ideal", "GammaPhi", "pxy", "txy", "azeotrope"]


# ---------------------------------------------------------------------------
# vapor pressure correlations
# ---------------------------------------------------------------------------
class VaporPressure:
    """Base class: a correlation supplies `P(T)` in Pa and inherits `T_sat`."""

    def P(self, T):
        """Vapor pressure at T (K), in Pa. Scalar in, scalar out."""
        raise NotImplementedError

    @staticmethod
    def _shape(T, value):
        """Return a plain float for scalar T, an array for array T -- so that a
        one-temperature call in a notebook prints `0.945` and not `np.float64(0.945)`."""
        return value if np.ndim(T) else float(value)

    #: hard upper limit on T for this correlation, or None. `Wagner` sets it to Tc,
    #: because its (1 - T/Tc)^1.5 term is NaN above the critical point rather than
    #: merely inaccurate -- and a NaN silently defeats a sign test, so a bracketing
    #: search that widened past Tc would hunt forever instead of failing.
    T_upper = None

    def T_sat(self, P, bracket=(100.0, 1000.0)):
        """Inverse: the temperature at which `P(T)` equals P (Pa).

        Used to seed and bracket the bubble- and dew-temperature solves, where a
        pure-component boiling point is the only guess available that is guaranteed
        to be in range.
        """
        lo, hi = bracket
        if self.T_upper is not None:
            hi = min(hi, self.T_upper * (1.0 - 1e-9))
            lo = min(lo, hi * 0.5)

        def f(T):
            v = self.P(T)
            if not np.isfinite(v) or v <= 0:
                raise ValueError("correlation undefined here")
            return float(np.log(v) - np.log(P))

        for _ in range(80):
            try:
                if f(lo) * f(hi) <= 0:
                    return float(optimize.brentq(f, lo, hi, xtol=1e-10))
            except (ValueError, FloatingPointError, ZeroDivisionError):
                pass
            lo *= 0.9
            if self.T_upper is None:
                hi *= 1.1
            if lo < 1.0:
                break
        raise ValueError(f"no saturation temperature for P = {P:.6g} Pa in "
                         f"({lo:.1f}, {hi:.1f}) K -- check the correlation's units "
                         f"and range")

    def __repr__(self):
        ps = ", ".join(f"{k}={v!r}" for k, v in self._params().items())
        return f"{type(self).__name__}({ps})"

    def _params(self):
        return {}


class Antoine(VaporPressure):
    """ln(Pvap/unit) = A - B/(T + C), with T in K. Default unit is the bar the book
    prints its constants in; the returned pressure is Pa.

    Illustration 10.2-2's ethyl acetate, typed in exactly as printed except for the
    sign of C (see the module docstring):

        >>> EA = Antoine(9.6830, 2842.2, -56.3209)
        >>> round(EA.P(348.15) / 1e5, 4)      # the book rounds this to 0.946 bar
        0.9453
    """

    def __init__(self, A, B, C=0.0, *, unit=BAR, base10=False):
        self.A, self.B, self.C = float(A), float(B), float(C)
        self.unit, self.base10 = float(unit), bool(base10)

    @classmethod
    def log10(cls, A, B, *, unit=BAR):
        """Problem 10.1-1's form, log10(Pvap/bar) = -A/T + B.

         Note the argument order and the sign: the book tabulates A as the
        *numerator over T* with a minus sign already in the equation, so
        `Antoine.log10(817.08, 4.402229)` is ethane as printed.
        """
        return cls(B, A, 0.0, unit=unit, base10=True)

    def _params(self):
        d = {"A": self.A, "B": self.B, "C": self.C}
        if self.base10:
            d["base10"] = True
        return d

    def P(self, T):
        T = np.asarray(T, dtype=float)
        e = self.A - self.B / (T + self.C)
        return self._shape(T, (10.0 ** e if self.base10 else np.exp(e)) * self.unit)


class ClausiusClapeyron(VaporPressure):
    """ln(Pvap/unit) = A - dHvap/(R T) -- the two-constant integrated form.

    This is what Illustration 10.1-1 prints, and the reason `B` is named `dHvap`
    is that in this form it *is* the enthalpy of vaporization in J/mol (SIS
    Eq. 7.7-5a), which is worth having visible rather than buried as a fitted
    constant:

        >>> C5 = ClausiusClapeyron(10.422, 26799.0)     # n-pentane
        >>> round(C5.P(323.15) / 1e5, 4)
        1.5648
    """

    def __init__(self, A, dHvap, *, unit=BAR):
        self.A, self.dHvap, self.unit = float(A), float(dHvap), float(unit)

    def _params(self):
        return {"A": self.A, "dHvap": self.dHvap}

    def P(self, T):
        T = np.asarray(T, dtype=float)
        return self._shape(T, np.exp(self.A - self.dHvap / (R * T)) * self.unit)


class Wagner(VaporPressure):
    """ln(Pvap/Pc) = (A t + B t^1.5 + C t^3 + D t^6)/(1 - t), with t = 1 - T/Tc.

    `Eq=1` in `code/data/pure_property.csv`, and the form most of that table uses.
    Reduced, so it needs Tc and Pc as well as the four constants; Pc carries the unit,
    which for that file is bar.
    """

    def __init__(self, A, B, C, D, Tc, Pc, *, unit=BAR):
        self.A, self.B, self.C, self.D = map(float, (A, B, C, D))
        self.Tc, self.Pc, self.unit = float(Tc), float(Pc), float(unit)
        self.T_upper = self.Tc

    def _params(self):
        return {"A": self.A, "B": self.B, "C": self.C, "D": self.D,
                "Tc": self.Tc, "Pc": self.Pc}

    def P(self, T):
        T = np.asarray(T, dtype=float)
        t = 1.0 - T / self.Tc
        # t^1.5 is NaN for t < 0, i.e. above Tc, where there is no vapor pressure to
        # report anyway. Return NaN deliberately rather than letting the power raise.
        with np.errstate(invalid="ignore"):
            lnPr = (self.A * t + self.B * np.where(t >= 0, t, np.nan) ** 1.5
                    + self.C * t ** 3 + self.D * t ** 6) / (1.0 - t)
        return self._shape(T, self.Pc * np.exp(lnPr) * self.unit)


class Riedel(VaporPressure):
    """ln Pvap = A - B/T + C ln T + D Pvap/T^2 -- `Eq=2`, and implicit in Pvap.

    The last term contains the answer, so this one is solved rather than evaluated.
    Fixed-point iteration from the D = 0 estimate converges in a few passes over the
    range the table covers; it is the only one of the three forms that cannot be
    written down explicitly, which is worth knowing before trusting a vectorized call.
    """

    def __init__(self, A, B, C, D, *, unit=BAR):
        self.A, self.B, self.C, self.D = map(float, (A, B, C, D))
        self.unit = float(unit)

    def _params(self):
        return {"A": self.A, "B": self.B, "C": self.C, "D": self.D}

    def P(self, T):
        T = np.asarray(T, dtype=float)
        # The iteration is exponential in its own output, so a starting temperature
        # well outside the fit range can run away before it converges. Overflow is
        # allowed to produce inf quietly and the caller sees a non-finite pressure --
        # which `T_sat` already treats as "undefined here" -- rather than a screenful
        # of RuntimeWarnings in the middle of a notebook.
        with np.errstate(over="ignore", invalid="ignore"):
            P = np.exp(self.A - self.B / T + self.C * np.log(T))   # D = 0 seed
            for _ in range(200):
                P_new = np.exp(self.A - self.B / T + self.C * np.log(T)
                               + self.D * P / T ** 2)
                if not np.all(np.isfinite(P_new)):
                    P = P_new
                    break
                if np.all(np.abs(P_new - P) <= 1e-12 * np.abs(P_new)):
                    P = P_new
                    break
                P = P_new
        return self._shape(T, P * self.unit)


def psat_from_database(name, *, strict=True):
    """The vapor-pressure correlation for `name` from `code/data/pure_property.csv`.

    Dispatches on that file's `Eq` selector -- 1 Wagner, 2 Riedel, 3 Antoine -- so a
    caller asks for a species and gets back something with `.P(T)` in Pa, without
    having to know which of the three forms the table happened to use for it. With
    `strict`, a species whose fit does not cover the temperature you later ask for is
    still returned; the range is in `.Tmin`/`.Tmax` for the caller to check.

        >>> hx = psat_from_database("n-hexane")
        >>> round(hx.P(333.15) / 1e5, 4)          # the book says 0.7583 bar at 60 C
        0.7632
    """
    from .data import load_pure_properties

    df = load_pure_properties()
    hit = df[df["Name"].str.lower() == str(name).lower()]
    if hit.empty:
        raise KeyError(f"{name!r} is not in pure_property.csv")
    r = hit.iloc[0]
    eq = int(r["Eq"])
    if eq == 1:
        vp = Wagner(r["VpA"], r["VpB"], r["VpC"], r["VpD"], r["Tc"], r["Pc"])
    elif eq == 2:
        vp = Riedel(r["VpA"], r["VpB"], r["VpC"], r["VpD"])
    elif eq == 3:
        vp = Antoine(r["VpA"], r["VpB"], r["VpC"])
    else:
        raise ValueError(f"{name!r} has no vapor-pressure fit (Eq = {eq})")
    vp.name, vp.Tmin, vp.Tmax = str(r["Name"]), float(r["Tmin"]), float(r["Tmax"])
    return vp


class TabulatedPsat(VaporPressure):
    """A single measured vapor pressure held constant -- for the illustrations that
    give one number at one temperature and nothing else.

    Illustration 10.2-1 is the case: benzene and cyclohexane at 77.6 C are given as
    0.993 and 0.980 bar with no correlation, and the whole illustration is isothermal,
    so no temperature dependence is needed or available.  `T_sat` therefore cannot
    work on this class and raises rather than returning a wrong answer quietly.
    """

    def __init__(self, P_value, T=None):
        self.P_value, self.T = float(P_value), (None if T is None else float(T))

    def _params(self):
        return {"P_value": self.P_value, "T": self.T}

    def P(self, T):
        return (np.full_like(np.asarray(T, dtype=float), self.P_value, dtype=float)
                if np.ndim(T) else self.P_value)

    def T_sat(self, P, bracket=(100.0, 1000.0)):
        raise TypeError("TabulatedPsat has no temperature dependence -- a constant "
                        "vapor pressure cannot support a bubble- or dew-temperature "
                        "calculation. Supply an Antoine or ClausiusClapeyron model.")


# ---------------------------------------------------------------------------
# the ideal liquid
# ---------------------------------------------------------------------------
class Ideal(ActivityModel):
    """gamma_i = 1 for every species -- Raoult's law, SIS Eq. 10.1-3.

    This exists so that Section 10.1 is not a separate code path from Section 10.2.
    The book presents the ideal-mixture calculation and the nonideal one as two
    procedures; they are one procedure and one substitution, and a reader who can
    swap `Ideal()` for `VanLaar(...)` in a single line sees that immediately.
    """

    n = None

    def lngamma(self, x, T):
        return np.zeros(np.asarray(x, dtype=float).size)

    def gex_over_RT(self, x, T):
        return 0.0

    def _params(self):
        return {}


# ---------------------------------------------------------------------------
# the gamma-phi solver
# ---------------------------------------------------------------------------
class GammaPhi:
    """Low-pressure VLE: activity coefficients in the liquid, ideal gas in the vapor.

    `psats` is one `VaporPressure` per species, in the same order as the activity
    model's species. `model` defaults to `Ideal()`.

    **The bubble-point pressure is explicit here, and that is the whole
    simplification of low-pressure VLE.** Because gamma does not depend on pressure,
    P = sum_i x_i gamma_i Pvap_i needs no iteration at all -- contrast
    `PRMixture.bubble_pressure`, which must iterate because the fugacity coefficients
    depend on the pressure being solved for. Everything that *does* iterate below
    (dew point, both temperature solves, the flash) iterates for a different reason:
    the unknown composition or temperature enters gamma or Pvap, not the pressure.
    """

    def __init__(self, psats, model=None, names=None):
        self.psats = list(psats)
        self.model = Ideal() if model is None else model
        self.names = list(names) if names is not None else None
        n = getattr(self.model, "n", None)
        if n is not None and n != len(self.psats):
            raise ValueError(f"{type(self.model).__name__} is written for {n} species "
                             f"but {len(self.psats)} vapor pressures were given")
        if self.names is not None and len(self.names) != len(self.psats):
            raise ValueError("names and psats must be the same length")

    @property
    def n(self):
        return len(self.psats)

    def __repr__(self):
        return (f"GammaPhi({self.n} species, model={type(self.model).__name__}"
                + (f", names={self.names}" if self.names else "") + ")")

    # --- the pieces -------------------------------------------------------
    def Psat(self, T):
        """Pure-component vapor pressures at T (K), Pa."""
        return np.array([float(p.P(T)) for p in self.psats])

    def gamma(self, x, T):
        """Activity coefficients -- just the model's, exposed for the illustrations
        that tabulate them (Illustration 10.2-1 prints a gamma column)."""
        return np.asarray(self.model.gamma(self._x(x), T), dtype=float)

    def K(self, x, T, P):
        """K_i = y_i/x_i = gamma_i Pvap_i / P, SIS Eq. 10.1-5."""
        return self.gamma(x, T) * self.Psat(T) / P

    @staticmethod
    def _x(x):
        x = np.asarray(x, dtype=float)
        s = x.sum()
        return x / s if s > 0 else x

    # --- bubble / dew points ---------------------------------------------
    def bubble_pressure(self, x, T, P_guess=None, max_iter=100, tol=1e-9):
        """Bubble-point pressure and incipient vapor composition at (x, T).

        Returns (P, y), P in Pa. No iteration: SIS Eq. 10.2-2b directly.
        `P_guess`, `max_iter` and `tol` are accepted and ignored, so that this
        signature matches `PRMixture.bubble_pressure` exactly.
        """
        x = self._x(x)
        p = x * self.gamma(x, T) * self.Psat(T)
        P = float(p.sum())
        return P, p / P

    def dew_pressure(self, y, T, P_guess=None, max_iter=200, tol=1e-12):
        """Dew-point pressure and incipient liquid composition at (y, T).

        Returns (P, x). Iterates because gamma depends on the unknown liquid x.
        """
        y = self._x(y)
        Ps = self.Psat(T)
        x = y.copy()                       # ideal-solution seed
        P = float(1.0 / np.sum(y / Ps))
        for _ in range(max_iter):
            g = self.gamma(x, T)
            P = float(1.0 / np.sum(y / (g * Ps)))
            x_new = y * P / (g * Ps)
            x_new = x_new / x_new.sum()
            if np.max(np.abs(x_new - x)) < tol:
                return P, x_new
            x = x_new
        return P, x

    def bubble_temperature(self, x, P, T_guess=None, max_iter=100, tol=1e-9):
        """Bubble-point temperature and vapor composition at (x, P).

        Returns (T, y). Solves sum_i x_i gamma_i(x, T) Pvap_i(T) = P for T.
        """
        x = self._x(x)
        f = lambda T: float(np.sum(x * self.gamma(x, T) * self.Psat(T))) - P
        T = self._solve_T(f, x, P, T_guess)
        return T, self.bubble_pressure(x, T)[1]

    def dew_temperature(self, y, P, T_guess=None, max_iter=100, tol=1e-9):
        """Dew-point temperature and liquid composition at (y, P).

        Returns (T, x). The inner loop is `dew_pressure`'s; the outer solves for T.
        """
        y = self._x(y)
        f = lambda T: self.dew_pressure(y, T)[0] - P
        T = self._solve_T(f, y, P, T_guess)
        return T, self.dew_pressure(y, T)[1]

    def _solve_T(self, f, z, P, T_guess):
        """Bracket and solve a scalar equilibrium-temperature equation.

         **The bracket cannot be the range of pure boiling points.** For a
        minimum-boiling azeotrope -- ethyl acetate/benzene in Illustration 10.2-2 is
        exactly one -- the bubble temperature falls *below* both pure boiling points,
        so the obvious bracket excludes the root. The pure boiling points are used
        only as a starting scale, then the interval is widened until it brackets.
        """
        Tb = []
        for p in self.psats:
            try:
                Tb.append(p.T_sat(P))
            except (TypeError, ValueError):
                pass
        if not Tb:
            raise ValueError("no vapor-pressure correlation can supply a boiling "
                             "point to seed the temperature solve")
        lo, hi = min(Tb), max(Tb)
        if T_guess is not None:
            lo, hi = min(lo, T_guess), max(hi, T_guess)
        span = max(hi - lo, 10.0)
        for _ in range(80):
            try:
                if f(lo) * f(hi) <= 0:
                    return float(optimize.brentq(f, lo, hi, xtol=1e-10))
            except (ValueError, FloatingPointError, ZeroDivisionError):
                pass
            lo, hi = lo - 0.25 * span, hi + 0.25 * span
            if lo <= 1.0:
                lo = 1.0
        raise ValueError(f"could not bracket the equilibrium temperature at "
                         f"P = {P:.6g} Pa")

    # --- isothermal flash -------------------------------------------------
    def flash(self, z, T, P, max_iter=500, tol=1e-12):
        """Isothermal (T, P) flash of a feed of composition z.

        Returns (beta, x, y): beta is the vapor molar fraction, and beta == 0 or 1
        means the feed is single-phase at (T, P). Same contract as
        `PRMixture.flash`, and it reuses that class's Rachford-Rice solver rather
        than carrying a second copy -- the equation is the same one (SIS Eq. 10.1-7);
        only the source of the K values differs.
        """
        z = self._x(z)
        Ps = self.Psat(T)
        x = z.copy()
        K = self.gamma(x, T) * Ps / P
        beta = 0.5
        for _ in range(max_iter):
            beta = PRMixture._rachford_rice(z, K, beta)
            x = z / (1 + beta * (K - 1))
            x = x / x.sum()
            K_new = self.gamma(x, T) * Ps / P
            if np.max(np.abs(K_new / K - 1)) < tol:
                K = K_new
                break
            K = K_new
        beta = float(min(max(PRMixture._rachford_rice(z, K, beta), 0.0), 1.0))
        x = z / (1 + beta * (K - 1))
        y = K * x
        return beta, x / x.sum(), y / y.sum()


# ---------------------------------------------------------------------------
# diagram generators -- method-agnostic by construction
# ---------------------------------------------------------------------------
def pxy(model, T, n=201, x1=None):
    """The P-x-y curve of a **binary** at fixed T. Returns (x1, y1, P), P in Pa.

    Works with any object exposing `bubble_pressure(x, T) -> (P, y)`, which means
    `GammaPhi` and `PRMixture` alike -- Illustration 10.2-2's ethyl acetate/benzene
    and Section 10.3's acetone/water come out of this one function.

    **One sweep gives both curves.** The bubble calculation returns the incipient
    vapor with the pressure, so P-vs-x1 is the bubble line and P-vs-y1 is the dew
    line of the same diagram. There is no second, separate dew sweep to run, and
    running one would be a different (and wrong) picture.

    **Points with no solution come back as NaN rather than raising**, so one bad
    composition does not destroy a sweep -- and matplotlib draws NaN as a gap, which
    is the correct picture. This is not defensive padding: **Figure 10.3-13 prints
    exactly such a gap**, labeled "region of nonconvergence with k12 = 0", for
    acetone/water on the van der Waals one-fluid mixing rule. The equation of state
    genuinely has no bubble point there, the band closes as k12 is fitted upward, and
    a generator that raised instead of returning NaN could not reproduce the figure.
    """
    x1 = np.linspace(0.0, 1.0, n) if x1 is None else np.asarray(x1, dtype=float)
    P, y1 = np.full(x1.size, np.nan), np.full(x1.size, np.nan)
    for i, a in enumerate(_clip01(x1)):
        try:
            Pi, yi = model.bubble_pressure([a, 1.0 - a], T)
        except Exception:
            continue                      # no solution here; leave the NaN
        if np.isfinite(Pi) and Pi > 0:
            P[i], y1[i] = Pi, yi[0]
    return x1, y1, P


def txy(model, P, n=201, x1=None):
    """The T-x-y curve of a **binary** at fixed P. Returns (x1, y1, T), T in K.

    The constant-pressure counterpart of `pxy`, and the expensive one: every point
    is a root find, because temperature enters both gamma and every vapor pressure.
    """
    x1 = np.linspace(0.0, 1.0, n) if x1 is None else np.asarray(x1, dtype=float)
    T, y1 = np.full(x1.size, np.nan), np.full(x1.size, np.nan)
    T_guess = None
    for i, a in enumerate(_clip01(x1)):
        try:
            Ti, yi = model.bubble_temperature([a, 1.0 - a], P, T_guess=T_guess)
        except Exception:
            continue                      # no solution here; see the note in `pxy`
        if np.isfinite(Ti) and Ti > 0:
            T[i], y1[i] = Ti, yi[0]
            T_guess = Ti                  # march along; each point seeds the next
    return x1, y1, T


def _clip01(x1, eps=1e-9):
    """Pure-component ends make gamma and Rachford-Rice degenerate; step just inside."""
    return np.clip(x1, eps, 1.0 - eps)


def azeotrope(model, T=None, P=None, bracket=(1e-6, 1.0 - 1e-6)):
    """Find a binary azeotrope -- the composition where y1 = x1.

    Give `T` for the isothermal azeotrope (returns `(x1, P)`, P in Pa) or `P` for the
    isobaric one (returns `(x1, T)`, T in K). Returns None when the residual y1 - x1
    does not change sign, i.e. there is no azeotrope in range.

    Worth having as a function rather than read off a plot: Chapter 10 prints
    azeotropic compositions and pressures to three and four figures, and those printed
    values are checkable only if the crossing can be located numerically.
    It finds **one** crossing. The hexafluorobenzene-benzene system of Table 10.2-4
    has *two* (a minimum and a maximum boiling azeotrope); scan sub-intervals for that.
    """
    if (T is None) == (P is None):
        raise ValueError("give exactly one of T or P")

    def resid(a):
        x = [a, 1.0 - a]
        if T is not None:
            return model.bubble_pressure(x, T)[1][0] - a
        return model.bubble_temperature(x, P)[1][0] - a

    lo, hi = bracket
    try:
        if resid(lo) * resid(hi) > 0:
            return None
        a = float(optimize.brentq(resid, lo, hi, xtol=1e-12))
    except (ValueError, FloatingPointError):
        return None
    if T is not None:
        return a, float(model.bubble_pressure([a, 1.0 - a], T)[0])
    return a, float(model.bubble_temperature([a, 1.0 - a], P)[0])
