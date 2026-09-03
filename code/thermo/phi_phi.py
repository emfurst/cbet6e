"""The phi-phi VLE drivers: bubble point, dew point and isothermal flash.

Section 10.3's method. The equilibrium condition is the same one the whole book uses,

    f_i^L = f_i^V   ->   x_i phi_i^L(x, T, P) = y_i phi_i^V(y, T, P)

and everything here follows from it. What makes the drivers reusable is that they
touch the mixture model through exactly two things: `self.components` (for the Wilson
K-value seed, which needs Tc, Pc and omega) and `self.phi(x, T, P, phase)`. **Any
mixture model that can return a fugacity coefficient can be driven by them.**

WHY THIS IS ITS OWN MODULE. These solvers were written inside `PRMixture`, where the
mixing rule is van der Waals one-fluid. Chapter 10 needs the identical solvers driven by
the **Wong-Sandler** rule (Figs. 10.3-9 to 10.3-13), and `wong_sandler.py` importing them
from `pr_mixture.py` would say the Wong-Sandler rule depends on the van der Waals one --
which is false, and is the kind of import that quietly turns a package into a tangle.
Both mixing rules now inherit from here instead, which is also the chapter's own thesis:
**one procedure, three ways of getting the fugacities.** Sections 10.1 and 10.2 reach it
through `thermo.vle.GammaPhi`, whose solvers carry these same names and signatures, so
`pxy`/`txy` draw a diagram without knowing which of the three produced it.

Deliberately Tier 1: no matplotlib, no drawing. Import cost is numpy + scipy only.

Eric M. Furst
August 2026
"""
import numpy as np
from scipy import optimize

_TRIVIAL_TOL = 1e-6      # |y - x| below this is one phase, not two

# How much the residual must move across a +-0.1% probe in the outer variable for its
# zero to mean anything. MEASURED, not chosen: successive substitution on K = phi_L/phi_V
# has a noise floor near a mixture critical point of about 1e-7 in sum(K x) - 1 -- that is
# what the residual reads once the iteration has drifted onto the trivial branch. A
# crossing that moves the residual by less than ten times that has its location set by
# round-off rather than by thermodynamics. Along the 377.65 K carbon dioxide/isopentane
# envelope this quantity falls smoothly from 9.6e-4 at x = 0.40 to 6.6e-8 at x = 0.70, so
# the cut is not near anything the figure needs.
_RESOLVABLE = 1e-6


class PhiPhiVLE:
    """Bubble/dew/flash for any mixture model providing `components` and `phi`.

    A mixin, not a model: it defines no thermodynamics of its own. The host class
    supplies

        self.components           sequence with .Tc, .Pc, .omega  (for the K seed)
        self.phi(x, T, P, phase)  fugacity coefficients, phase in {"liquid","vapor"}

    and inherits `bubble_pressure`, `dew_pressure`, `bubble_temperature`,
    `dew_temperature` and `flash`.

    **THE TRIVIAL ROOT, and why these solvers are not a plain Newton.** The
    equilibrium condition x_i phi_i^L = y_i phi_i^V is satisfied *identically* by
    y = x, K = 1 -- "equilibrium" between a phase and itself. It is a real root of
    the equations and a physically empty one, and successive substitution runs
    straight into it wherever the trial pressure lies outside the range in which
    two phases exist: there the cubic has one real root, so `Z(..., "liquid")` and
    `Z(..., "vapor")` return the same number, K = 1 exactly, and the outer residual
    sum(K x) - 1 is **zero at every pressure**. A Newton solver then declares
    convergence wherever it happens to be standing.

    This is not a corner case in Section 10.3. Above the critical temperature of the
    more volatile species -- 377.65 K for carbon dioxide in Figure 10.3-7 -- the
    upper half of every isotherm is near-critical, and the trivial root is what a
    naive solver returns for all of it. So the outer solve here **rejects any
    iterate whose incipient phase has collapsed onto the fixed one** (`_is_trivial`),
    brackets the residual on a ladder of pressures where two phases do exist, and
    only then converges. See `_solve_outer`.

    Convergence is successive substitution on the inner loop and a bracketed
    Newton on the outer. Where a model has no two-phase solution -- inside a region
    the equation of state makes single-phase -- these return **NaN** rather than
    raising, and rather than the last iterate. That is deliberate: `thermo.vle.pxy`
    passes NaN through to the plot, and matplotlib draws it as a gap. **Figure
    10.3-13 prints exactly such a gap** ("region of nonconvergence with k_ij = 0"),
    so a driver that raised could not draw the figure the book has -- and one that
    returned its last iterate, as this module did until 2026-08-15, silently drew
    the trivial root as if it were data.
    """

    # --- K-value initial guess ------------------------------------------
    def _wilson_K(self, T, P):
        """Wilson-correlation K_i estimate, used only to start the iterations."""
        Tc = np.array([c.Tc for c in self.components])
        Pc = np.array([c.Pc for c in self.components])
        w = np.array([c.omega for c in self.components])
        return Pc / P * np.exp(5.373 * (1 + w) * (1 - Tc / T))

    # --- bubble / dew points --------------------------------------------
    def _equilibrate_y(self, x, T, P, y0, max_iter=200, tol=1e-10):
        """At fixed (x, T, P), the incipient vapor y. Returns (sum K x, y)."""
        y = np.array(y0, dtype=float)
        y /= y.sum()
        s = 1.0
        for _ in range(max_iter):
            phiL = self.phi(x, T, P, "liquid")
            phiV = self.phi(y, T, P, "vapor")
            K = phiL / phiV
            y_new = K * np.asarray(x, dtype=float)
            s = y_new.sum()
            y_new = y_new / s
            if np.max(np.abs(y_new - y)) < tol:
                return s, y_new
            y = y_new
        return s, y

    def _equilibrate_x(self, y, T, P, x0, max_iter=200, tol=1e-10):
        """At fixed (y, T, P), the incipient liquid x. Returns (sum y/K, x)."""
        x = np.array(x0, dtype=float)
        x /= x.sum()
        s = 1.0
        for _ in range(max_iter):
            phiV = self.phi(y, T, P, "vapor")
            phiL = self.phi(x, T, P, "liquid")
            K = phiL / phiV
            x_new = np.asarray(y, dtype=float) / K
            s = x_new.sum()
            x_new = x_new / s
            if np.max(np.abs(x_new - x)) < tol:
                return s, x_new
            x = x_new
        return s, x

    # --- the outer solve, and the trivial root -----------------------------
    def _is_trivial(self, incipient, fixed, T, tol=_TRIVIAL_TOL):
        """Has the incipient phase collapsed onto the phase held fixed?

        The test is on composition, not on the number of roots of the cubic. A
        single real root does **not** by itself mean the trivial root: the liquid
        and vapor fugacities are evaluated at *different* compositions, so K is
        still meaningful. What identifies the trivial root is that the two
        compositions have become the same one, and when they have, they agree to
        machine precision rather than merely closely -- the collapse is a fixed
        point of the iteration, not an approach to one.

        **A SUBCRITICAL pure component is exempt, and must be.** At x = (1, 0)
        the vapor really is the liquid's composition, so y = x is the answer rather
        than the trivial root -- the condition phi^L = phi^V it satisfies is the
        ordinary vapor-pressure criterion. Without this exemption every P-x-y sweep
        loses both of its endpoints.

        **Above that species' own critical temperature the exemption must not
        apply.** A supercritical pure fluid has no vapor pressure at all, so there
        the collapse is the trivial root again and returning a pressure for it is
        worse than returning nothing: carbon dioxide at Figure 10.3-7's 377.65 K is
        73 K above its critical temperature, and a bubble point at x_CO2 = 1 would
        be a number with no physical referent anchoring the end of the curve. The
        figure's own closing question asks the student *why* the envelope cannot
        reach x = 1 there.

        A genuine mixture critical point also has y -> x, so points within
        `tol` of it are rejected too. That is the honest outcome: there the two
        are numerically indistinguishable. Figure 10.3-7's critical points are
        therefore *estimated* by extrapolating y - x to zero, which is what its
        caption says they are.
        """
        fixed = np.asarray(fixed, dtype=float)
        pure = int(np.argmax(fixed))
        if float(fixed[pure]) > 1.0 - tol:
            return T >= self.components[pure].Tc      # pure: only if supercritical
        return float(np.max(np.abs(np.asarray(incipient, dtype=float)
                                   - fixed))) < tol

    def _P_bounds(self):
        """Ladder limits for a bracketing search in pressure (Pa)."""
        Pc = max(c.Pc for c in self.components)
        return 1e2, 5.0 * Pc          # a mixture critical P can exceed both pure ones

    def _T_bounds(self):
        """Ladder limits for a bracketing search in temperature (K)."""
        Tc = [c.Tc for c in self.components]
        return 0.25 * min(Tc), 2.0 * max(Tc)

    @staticmethod
    def _bracket(f, guess, lo_b, hi_b, ratio=1.1, n_walk=80, n_creep=12):
        """Bracket a sign change of `f` by walking outward from `guess`.

        `f(v) -> float`, NaN where the trivial root is all there is. Returns
        `(a, b)` with f(a) and f(b) of opposite sign, or None.

        A blind ladder over the whole pressure range does not work here, and
        that is worth stating because it is the obvious thing to try. The window
        in which two phases exist can be **narrower than one rung** near a mixture
        critical point -- for carbon dioxide/isopentane at 377.65 K and
        x = 0.59 it is a few bar wide, on a range that runs to hundreds -- so a
        ladder coarse enough to be affordable steps straight over the answer.
        Walking outward from the ideal-solution seed finds it instead, because
        that seed is the right order of magnitude even where it is a factor of
        two high.
        """
        # the nearest point to the seed at which two phases exist
        anchor = next((v for i in range(n_walk)
                       for v in ((guess,) if i == 0
                                 else (guess * ratio ** -i, guess * ratio ** i))
                       if lo_b < v < hi_b and np.isfinite(f(v))), None)
        if anchor is None:
            return None

        f_a = f(anchor)
        step = ratio if f_a > 0.0 else 1.0 / ratio   # march toward the crossing
        v, f_v = anchor, f_a
        for _ in range(n_walk):
            nxt = v * step
            if not (lo_b < nxt < hi_b):
                return None
            f_n = f(nxt)
            for _ in range(n_creep):                 # the window may end first;
                if np.isfinite(f_n):                 # creep up on its edge
                    break
                nxt = np.sqrt(v * nxt)
                f_n = f(nxt)
            else:
                return None
            if f_v * f_n <= 0.0:
                return (v, nxt) if v < nxt else (nxt, v)
            v, f_v = nxt, f_n
        return None

    def _solve_outer(self, residual, fixed, T_at, guess, state, bounds,
                     max_iter, tol):
        """Solve `residual(v, state) -> (sum, state)` for sum = 1, rejecting the
        trivial root.

        `fixed` is the composition of the phase held fixed, against which each
        iterate's incipient composition is tested; `T_at(v)` gives the temperature
        of an iterate, which is constant for a pressure solve and the variable
        itself for a temperature one. Returns `(v, incipient)`, or `(nan, fixed)`
        where the model has no two-phase solution.

        Newton from the caller's guess is the fast path and handles the ordinary
        subcritical points. When it walks into a single-phase region -- which is
        where the trivial root appears -- the residual there is NaN rather than
        zero, Newton gives up, and `_bracket` walks outward from the seed to find
        two points at which two phases *do* exist and the residual changes sign.
        Bisection inside that bracket cannot escape it, so the near-critical
        branch converges to the physical root instead of the empty one.

        Every evaluation is seeded from the same `state`, never warm-started
        from the previous iterate. Warm starting is faster and makes the answer
        depend on the path taken to it; a figure has to be reproducible.
        """
        lo_b, hi_b = bounds
        cache = {}

        def f(v):
            """(sum - 1, incipient composition) at v; NaN at the trivial root."""
            key = float(v)
            if key not in cache:
                s, st = residual(key, state)
                trivial = self._is_trivial(st, fixed, T_at(key))
                cache[key] = (np.nan if trivial else s - 1.0, st)
            return cache[key]

        def crosses(v, rel=1e-3):
            """Does the residual change SIGN through v, or is v on a flat plateau?

            `|f| < tol` is not enough to identify a bubble point near a mixture
            critical point, and this is the subtlest failure in the module. The
            collapse onto the trivial root is not sudden: over a *range* of
            pressures above the true bubble point the iteration drifts onto the
            trivial branch, where the residual is ~0 at every pressure. The whole
            plateau therefore passes `|f| < tol`, and the composition test does not
            catch it either, because y - x decays smoothly through the tolerance
            rather than jumping to zero -- for carbon dioxide/isopentane at 377.65 K
            and x = 0.655 the plateau runs from about 93.5 bar upward with
            |y - x| ~ 1e-4, a hundred times `_TRIVIAL_TOL`.

            A real bubble point is a **transversal** zero: the residual is positive
            below it and negative above. A plateau point has the same sign on both
            sides. That distinction needs no new tolerance, and it is what stops
            Newton from reporting 94.0 bar where the answer is 92.6.

            Sign change alone is still not enough *very* close to the critical
            point, where the residual flattens toward the noise floor and its zero
            wanders. The crossing must also be RESOLVABLE -- see `_RESOLVABLE`.
            Points that fail this come back as NaN, which is the honest answer:
            the bubble and dew branches are no longer distinguishable there, and
            the critical point is properly placed by extrapolating y - x to zero
            from the part of the envelope that is resolved.
            """
            below, _ = f(v * (1.0 - rel))
            above, _ = f(v * (1.0 + rel))
            if not (np.isfinite(below) and np.isfinite(above)):
                return False
            return below * above < 0.0 and abs(below - above) > _RESOLVABLE

        def newton(v0, lo, hi):
            v = float(v0)
            for _ in range(max_iter):
                fv, st = f(v)
                if not np.isfinite(fv):
                    return None                      # walked into a single phase
                if abs(fv) < tol:
                    return (float(v), st) if crosses(v) else None
                dv = v * 1e-6
                f2, _ = f(v + dv)
                if not np.isfinite(f2) or f2 == fv:
                    return None
                v = v - fv * dv / (f2 - fv)
                if not (lo < v < hi):
                    return None
            return None

        hit = newton(guess, lo_b, hi_b)
        if hit is not None:
            return hit

        # --- the seed sat in a single-phase region: bracket by walking out
        bracket = self._bracket(lambda v: f(v)[0], guess, lo_b, hi_b)
        if bracket is None:
            return float("nan"), np.asarray(fixed, dtype=float)

        lo, hi = bracket
        hit = newton(np.sqrt(lo * hi), lo, hi)
        if hit is not None:
            return hit

        f_lo = f(lo)[0]                              # safeguarded fallback
        for _ in range(200):
            mid = np.sqrt(lo * hi)
            f_mid, st = f(mid)
            if not np.isfinite(f_mid):
                break
            if hi / lo - 1.0 < 1e-13:
                return (float(mid), st) if crosses(mid) else (
                    float("nan"), np.asarray(fixed, dtype=float))
            # A near-zero residual is accepted only where it CROSSES. On the
            # trivial plateau f is ~0 without changing sign, so an early return on
            # |f| < tol alone would stop here -- and the plateau always lies ABOVE
            # the bubble point, so the search continues downward instead.
            if abs(f_mid) < tol:
                if crosses(mid):
                    return float(mid), st
                hi = mid
                continue
            if f_mid * f_lo > 0.0:
                lo, f_lo = mid, f_mid
            else:
                hi = mid
        return float("nan"), np.asarray(fixed, dtype=float)

    def bubble_pressure(self, x, T, P_guess=None, max_iter=100, tol=1e-9):
        """Bubble-point pressure and incipient vapor composition at (x, T) -> (P, y).

        The seed is the ideal-solution bubble pressure sum(x_i K_i) P_ref, which
        is the Raoult's-law estimate the Wilson correlation was built to give. It is
        *not* 1 / sum(x_i / K_i) -- that is the dew-point form, and using it here
        seeded the solve six orders of magnitude low.

        Figure 10.3-7's near-critical branch is the regression test, because it is
        where both faults showed. Carbon dioxide/isopentane at 377.65 K, k12 = 0.121:

        >>> import numpy as np
        >>> from thermo.peng_robinson import PengRobinson
        >>> from thermo.pr_mixture import PRMixture
        >>> m = PRMixture([PengRobinson.from_database("carbon dioxide"),
        ...                PengRobinson.from_database("2-methyl butane")],
        ...               kij=[[0.0, 0.121], [0.121, 0.0]])

        A composition well up the near-critical branch, where the trivial root used
        to be returned as P = 0 (Besserer and Robinson measure 59.16 bar here):

        >>> P, y = m.bubble_pressure([0.3481, 0.6519], 377.65)
        >>> round(P / 1e5, 1)
        57.7

        The vapor must be richer in carbon dioxide than the liquid -- the trivial
        root's signature is that it is not:

        >>> bool(y[0] > 0.3481 + 0.1)
        True

        And pure carbon dioxide has no bubble point at 377.65 K at all, being 73 K
        above its critical temperature:

        >>> bool(np.isnan(m.bubble_pressure([1.0, 0.0], 377.65)[0]))
        True

        while at 277.59 K it has one, and it is the measured 38.96 bar:

        >>> round(m.bubble_pressure([1.0, 0.0], 277.59)[0] / 1e5, 1)
        38.9

        The trivial PLATEAU, which is a different failure from the trivial root and was
        the second one found. At this composition the residual is ~0 at every pressure
        above about 93.5 bar, so `|f| < tol` alone accepted 94.0 bar; the answer is the
        transversal crossing at 92.5, and |y - x| there is 100x `_TRIVIAL_TOL`, so the
        composition test cannot separate the two:

        >>> P, y = m.bubble_pressure([0.6547, 0.3453], 377.65)
        >>> round(P / 1e5, 1)
        92.5

        And past the point where the branches are resolvable the driver returns NaN
        rather than a crossing whose location is round-off:

        >>> bool(np.isnan(m.bubble_pressure([0.75, 0.25], 377.65)[0]))
        True
        """
        x = np.asarray(x, dtype=float)
        K = self._wilson_K(T, 1e5)
        if P_guess is None:
            P_guess = float(np.sum(x * K)) * 1e5 if np.all(K > 0) else 1e5
        y = K * x
        y /= y.sum()
        return self._solve_outer(lambda P, y0: self._equilibrate_y(x, T, P, y0),
                                 x, lambda P: T, P_guess, y,
                                 self._P_bounds(), max_iter, tol)

    def dew_pressure(self, y, T, P_guess=None, max_iter=100, tol=1e-9):
        """Dew-point pressure and incipient liquid composition at (y, T) -> (P, x)."""
        y = np.asarray(y, dtype=float)
        K = self._wilson_K(T, 1e5)
        if P_guess is None:
            P_guess = float(1.0 / np.sum(y / K)) * 1e5 if np.all(K > 0) else 1e5
        x = y / K
        x /= x.sum()
        return self._solve_outer(lambda P, x0: self._equilibrate_x(y, T, P, x0),
                                 y, lambda P: T, P_guess, x,
                                 self._P_bounds(), max_iter, tol)

    def bubble_temperature(self, x, P, T_guess=None, max_iter=100, tol=1e-9):
        """Bubble-point temperature and vapor composition at (x, P) -> (T, y)."""
        x = np.asarray(x, dtype=float)
        if T_guess is None:
            T_guess = float(x @ np.array([c.Tc * 0.7 for c in self.components]))
        y = self._wilson_K(T_guess, P) * x
        y /= y.sum()
        return self._solve_outer(lambda T, y0: self._equilibrate_y(x, T, P, y0),
                                 x, lambda T: T, T_guess, y,
                                 self._T_bounds(), max_iter, tol)

    def dew_temperature(self, y, P, T_guess=None, max_iter=100, tol=1e-9):
        """Dew-point temperature and liquid composition at (y, P) -> (T, x)."""
        y = np.asarray(y, dtype=float)
        if T_guess is None:
            T_guess = float(y @ np.array([c.Tc * 0.7 for c in self.components]))
        x = y / self._wilson_K(T_guess, P)
        x /= x.sum()
        return self._solve_outer(lambda T, x0: self._equilibrate_x(y, T, P, x0),
                                 y, lambda T: T, T_guess, x,
                                 self._T_bounds(), max_iter, tol)

    # --- isothermal flash ------------------------------------------------
    def flash(self, z, T, P, max_iter=200, tol=1e-10):
        """Isothermal (T, P) flash of a feed z -> (beta, x, y).

        beta is the vapor molar fraction; 0 is a subcooled liquid and 1 a superheated
        vapor, i.e. the feed does not split at (T, P).

        **The trivial root reaches here too, and is not guarded against.** If the
        K values collapse to 1, Rachford-Rice returns beta = 0 and the answer reads
        "a subcooled liquid" -- which is the right answer whenever the feed really
        is single-phase, and the wrong one near a mixture critical point, where it
        is indistinguishable from a genuine no-split. Chapter 10's flash work is all
        `thermo.vle.GammaPhi` at low pressure, so nothing in the book currently
        depends on this path; a near-critical flash would need the same bracketing
        treatment `_solve_outer` gives the bubble and dew points.
        """
        z = np.asarray(z, dtype=float)
        K = self._wilson_K(T, P)
        beta = 0.5
        for _ in range(max_iter):
            beta = self._rachford_rice(z, K, beta)
            x = z / (1 + beta * (K - 1))
            y = K * x
            x = x / x.sum()
            y = y / y.sum()
            K_new = self.phi(x, T, P, "liquid") / self.phi(y, T, P, "vapor")
            if np.max(np.abs(K_new / K - 1)) < tol:
                K = K_new
                break
            K = K_new
        beta = min(max(self._rachford_rice(z, K, beta), 0.0), 1.0)
        x = z / (1 + beta * (K - 1))
        y = K * x
        return beta, x / x.sum(), y / y.sum()

    @staticmethod
    def _rachford_rice(z, K, beta0=0.5):
        """Solve sum_i z_i (K_i - 1) / (1 + beta (K_i - 1)) = 0 for beta in [0, 1],
        falling back to the single-phase edges when the feed does not split."""
        z, K = np.asarray(z, dtype=float), np.asarray(K, dtype=float)

        def g(beta):
            return np.sum(z * (K - 1) / (1 + beta * (K - 1)))

        if g(0.0) <= 0:          # bubble-point test
            return 0.0
        if g(1.0) >= 0:          # dew-point test
            return 1.0
        lo = 1.0 / (1.0 - K.max()) + 1e-10
        hi = 1.0 / (1.0 - K.min()) - 1e-10
        return float(optimize.brentq(g, max(lo, 0.0), min(hi, 1.0), xtol=1e-12))
