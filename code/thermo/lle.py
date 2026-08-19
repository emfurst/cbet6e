"""lle -- liquid-liquid and vapor-liquid-liquid equilibrium, SIS Sections 11.2 and 11.3.

Chapter 11 opens by saying that what separates its five cases from Chapter 10 "is how
the fugacity of each species is computed," and Section 11.2 makes that literal. The
equilibrium condition is the same equality of species fugacities used since Chapter 8;
because both phases are liquids, the pure-component fugacity cancels from both sides
and what is left is SIS Eq. 11.2-2,

    x_i^I gamma_i(x^I, T) = x_i^II gamma_i(x^II, T)        i = 1, 2, ..., C

with no vapor pressure in it anywhere. That is the whole module. Everything below is
either that equation solved, that equation differentiated (stability, the spinodal, the
consolute temperature), or that equation with SIS Chapter 10's bubble point bolted on
top of it (Section 11.3's vapor-liquid-liquid equilibrium).

    from thermo import VanLaar
    from thermo.lle import binary_lle

    xI, xII = binary_lle(VanLaar(2.62, 3.02), 310.95)   # Illustration 11.2-2
    xI[0], xII[0]                                       # 0.1128, 0.9284

## Why this is not a solver you can point at the equations and walk away from

**Equation 11.2-2 has a trivial root.** x^I = x^II satisfies it at every temperature and
every composition, for every model, and it is the root a naive solver finds. Worse, it
is the *correct* root wherever the mixture does not split, so a solver that returns it
silently cannot be distinguished from one that works. Everything in the binary path
below exists to avoid that: the seeds do not come from a guess but from the **lower
convex hull** of the Gibbs energy of mixing, which is the common-tangent construction
of SIS Fig. 11.2-5 done numerically, and which cannot produce a trivial seed because a
hull that never leaves the curve reports *no split* instead. After the root solve,
`binary_lle` checks that the common tangent it found actually lies below the curve
everywhere -- a local tangent that does not is a wrong answer that satisfies
Eq. 11.2-2 exactly.

That is deliberate belt-and-braces. Two liquid compositions that satisfy the equilibrium
equation to 1e-12 and are still wrong is not hypothetical here; it is the ordinary
failure of this calculation.

## What the book's illustrations are, and where each is solved

| SIS | what | entry point |
|-----|------|-------------|
| Ill. 11.2-1 | amounts of two phases from a read phase diagram | `tie_line_split` |
| Ill. 11.2-2 | binary LLE from van Laar | `binary_lle` |
| Ill. 11.2-3, 11.2-4 | activity coefficients from solubility | `gamma_from_solubility` |
| Ill. 11.2-6 | polymer compatibility, Flory-Huggins, T-dependent chi | `binary_lle_envelope` |
| Ill. 11.2-7 | mass balance on a triangular diagram | `mix_streams` |
| Ill. 11.2-8, 11.2-9 | one- and two-stage extraction | `tie_line_split` |
| Ill. 11.2-5 | LLE from an equation of state | `eos_binary_lle` |
| Ill. 11.3-1 | VLLE from an activity coefficient model | `vlle_binary` |
| Ill. 11.3-2 | three-phase equilibrium from an equation of state | `eos_vlle_binary` |
| Ill. 11.3-3 | steam distillation of turpentine | `steam_distillation` |
| Fig. 11.2-5 | G of mixing, ideal and not, with the tangent | `gmix_over_RT`, `common_tangent` |
| Fig. 11.3-1, 11.3-2 | P-x with and without LLE | `pxy_lle` |
| Fig. 11.3-4 | T-x at fixed P, heterogeneous azeotrope | `txy_lle` |

Illustrations 11.2-5 and 11.3-2 are the equation-of-state route, and the seam is at the
fugacity rather than at the phase: `lle_flash` accepts anything with a `gamma(x, T)`,
and `eos_binary_lle`/`eos_vlle_binary` accept anything with an
`ln_phi(x, T, P, phase=...)` -- so `PRMixture` and `wong_sandler` both work, with the
liquid root taken in both phases, which is the change SIS p. 630 describes.

⚠️ **The two routes do not have the same signature, and the difference is physical.**
The activity-coefficient functions take T and no P; the equation-of-state ones take
both, because a cubic has no incompressible-liquid shortcut. `eos_vlle_binary` therefore
solves all four unknowns at once rather than solving LLE first and taking a bubble point
afterwards, as `vlle_binary` legitimately can.

## Models whose parameters depend on temperature

Illustration 11.2-6 sweeps 25 C to 600 C with chi = 1473/T, so the model is different at
every temperature. Every function here that varies T therefore takes **either** an
`ActivityModel` (parameters fixed, the Table 9.5-1 case) **or a callable `T -> model`**:

    envelope = binary_lle_envelope(lambda T: FloryHuggins(1473.0 / T, 0.943), T_grid)

## Units

SI, as everywhere else in the package: T in K, P in Pa, energies in J/mol. Compositions
are mole fractions unless a function says otherwise -- `mix_streams` and
`tie_line_split` are fraction-agnostic and are used with the *weight* fractions the
ternary illustrations are printed in.

Eric M. Furst
August 2026
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy import constants, optimize

from .activity_models import ActivityModel

R = constants.R

__all__ = ["as_model_of_T", "gmix_over_RT", "d2gmix_over_RT", "is_stable",
           "spinodal", "common_tangent", "binary_lle", "binary_lle_envelope",
           "consolute_temperature", "lle_flash", "tie_line_split", "mix_streams",
           "gamma_from_solubility", "solubility_at_T",
           "vlle_binary", "vlle_temperature", "immiscible_pressure",
           "steam_distillation", "eos_binary_lle", "eos_vlle_binary",
           "pxy_lle", "txy_lle"]


# ---------------------------------------------------------------------------
# the model-of-temperature convention
# ---------------------------------------------------------------------------
def as_model_of_T(model):
    """Normalize `model` to a callable `T -> ActivityModel`.

    An `ActivityModel` whose parameters are already fixed becomes a callable that
    ignores T; a callable is returned unchanged. This is the one concession the
    module makes to Illustration 11.2-6, where chi = 1473/T means there is a
    different Flory-Huggins model at every point of the phase boundary.
    """
    if isinstance(model, ActivityModel):
        return lambda T, _m=model: _m
    if callable(model):
        return model
    raise TypeError("model must be an ActivityModel or a callable T -> ActivityModel")


# ---------------------------------------------------------------------------
# Sec. 11.2 -- the Gibbs energy of mixing and its curvature
# ---------------------------------------------------------------------------
def gmix_over_RT(model, x1, T):
    """Delta_mix G / RT for a binary at mole fraction `x1`, SIS Eq. 11.2-16.

        Delta_mix G / RT = x1 ln(x1 gamma_1) + x2 ln(x2 gamma_2)

    Accepts a scalar or an array of x1 and returns the same shape. The endpoints
    are 0 by the limit x ln x -> 0, taken directly rather than divided by zero.

    This is the quantity plotted in SIS Fig. 11.2-5, whose whole message is that
    the ideal curve (gamma = 1) is convex and cannot support a tangent touching it
    twice, while a curve with A/RT = 3 can.
    """
    m = as_model_of_T(model)(T)
    scalar = np.isscalar(x1)
    out = []
    for a in np.atleast_1d(np.asarray(x1, dtype=float)):
        if a <= 0.0 or a >= 1.0:
            out.append(0.0)
            continue
        x = np.array([a, 1.0 - a])
        lg = np.asarray(m.lngamma(x, T), dtype=float)
        out.append(float(np.sum(x * (np.log(x) + lg))))
    out = np.array(out)
    return float(out[0]) if scalar else out


def _gex_over_RT_from_gamma(m, x1, T):
    """SUM_i x_i ln gamma_i at x1 -- G^ex/RT *as the activity coefficients imply it*.

    Not `model.gex_over_RT`, and the difference is deliberate. For most models in
    `activity_models` those two are the same function written twice, which is the
    whole point of `check_gibbs_duhem`. For a model that is **thermodynamically
    inconsistent** -- the printed variants of Flory-Huggins Eq. 9.5-18, which exist in
    that module precisely so the book's own numbers can be reproduced -- they are not,
    and then only this one is the right choice: Eq. 11.2-2 is written in activity
    coefficients, so every curvature, spinodal and consolute temperature here has to
    follow the same gamma the equilibrium condition does. Reading `gex_over_RT`
    instead silently reported the *consistent* model's critical temperature for an
    inconsistent one.
    """
    x = np.array([x1, 1.0 - x1])
    return float(x @ np.asarray(m.lngamma(x, T), dtype=float))


def d2gmix_over_RT(model, x1, T, h=1e-4):
    """Second composition derivative of `gmix_over_RT`, SIS Eq. 11.2-11's left side.

    Split as

        d2(Delta_mix G/RT)/dx1^2 = 1/(x1 x2) + d2(G^ex/RT)/dx1^2

    so the ideal part is exact and only the excess part is differenced. That matters:
    1/(x1 x2) diverges at both ends, and differencing it numerically there is what
    makes a naive stability test report nonsense in the dilute corners where the
    answer is in fact unambiguous.

    For the one-constant Margules model this returns RT/(x1 x2) - 2A, divided by RT,
    which is SIS Eq. 11.2-11 exactly.
    """
    m = as_model_of_T(model)(T)
    scalar = np.isscalar(x1)
    out = []
    for a in np.atleast_1d(np.asarray(x1, dtype=float)):
        step = min(h, 0.25 * a, 0.25 * (1.0 - a))
        if step <= 0:
            out.append(np.inf)
            continue
        g = [_gex_over_RT_from_gamma(m, b, T) for b in (a - step, a, a + step)]
        out.append(1.0 / (a * (1.0 - a)) + (g[0] - 2.0 * g[1] + g[2]) / step**2)
    out = np.array(out)
    return float(out[0]) if scalar else out


def is_stable(model, x1, T):
    """True where a single liquid phase is intrinsically stable, SIS Eq. 11.2-9.

    The test is (d2 G / dx1^2)_{T,P} > 0. **A True here does not mean one phase is
    the equilibrium state** -- between the binodal and the spinodal the single phase
    is metastable, stable against small fluctuations and unstable against the finite
    one that nucleates the second phase. SIS says this in the footnote to
    Eq. 11.2-10 by analogy with the van der Waals loop of Sec. 7.3: the limit of
    stability is not the equilibrium composition. Use `binary_lle` for that.
    """
    d2 = d2gmix_over_RT(model, x1, T)
    return d2 > 0.0 if np.isscalar(x1) else np.asarray(d2) > 0.0


def spinodal(model, T, n=401):
    """The limits of stability at T: the two roots of d2(Delta_mix G)/dx1^2 = 0.

    Returns `(x1_lo, x1_hi)`, or **None** when the curvature is positive everywhere
    and the single phase is stable at every composition. Because 1/(x1 x2) diverges
    at both ends the curvature is always positive near the pure components, so there
    are either zero roots or two.

    These are the inflection points of SIS Fig. 11.2-5, and the dashed spinodal curve
    of Problem 11.2-2. They are *not* the coexisting compositions.
    """
    xs = np.linspace(1e-6, 1.0 - 1e-6, n)
    d2 = d2gmix_over_RT(model, xs, T)
    if np.min(d2) > 0.0:
        return None
    i = int(np.argmin(d2))
    f = lambda a: d2gmix_over_RT(model, float(a), T)
    lo_hi = []
    for side in (slice(None, i + 1), slice(i, None)):
        seg_x, seg_d = xs[side], d2[side]
        sign = np.sign(seg_d)
        cross = np.nonzero(np.diff(sign) != 0)[0]
        if cross.size == 0:
            return None
        k = cross[0] if side.start is None else cross[-1]
        lo_hi.append(float(optimize.brentq(f, seg_x[k], seg_x[k + 1], xtol=1e-12)))
    return tuple(lo_hi)


def _hull_gaps(x, g):
    """Grid points spanned by the lower convex hull of the graph (x, g).

    Returns a list of `(i, j)` index pairs such that the hull goes straight from
    `x[i]` to `x[j]` while the curve dips above it in between -- one pair per
    two-phase region. An empty list means the curve is already convex, that is, no
    phase split at all.

    This is the common-tangent construction of SIS Fig. 11.2-5 done by machine, and
    it is the seed for every binary solve in this module. Monotone chain, O(n).
    """
    hull = []
    for k in range(len(x)):
        while (len(hull) >= 2
               and ((g[k] - g[hull[-2]]) * (x[hull[-1]] - x[hull[-2]])
                    <= (g[hull[-1]] - g[hull[-2]]) * (x[k] - x[hull[-2]]))):
            hull.pop()
        hull.append(k)
    return [(a, b) for a, b in zip(hull, hull[1:]) if b - a > 1]


def common_tangent(model, T, n=1001):
    """Bracketing compositions of every two-phase region at T, from the hull.

    Returns a list of `(x1_lo, x1_hi)` grid brackets -- **approximate**, at the
    resolution of the grid; `binary_lle` refines one of them to the equilibrium
    compositions. The list is normally empty (one phase) or one pair (two phases);
    more than one pair is the signature of a curve with several distinct two-phase
    regions, which is what SIS Fig. 11.2-11d shows for a ternary and what
    Eq. 11.2-25 generalizes the equations to.
    """
    x = np.linspace(1e-9, 1.0 - 1e-9, n)
    g = gmix_over_RT(model, x, T)
    return [(float(x[a]), float(x[b])) for a, b in _hull_gaps(x, g)]


# ---------------------------------------------------------------------------
# Sec. 11.2 -- binary liquid-liquid equilibrium
# ---------------------------------------------------------------------------
def binary_lle(model, T, guess=None, n=1001, tol=1e-10, min_gap=1e-6,
               check_tangent=True):
    """Coexisting liquid compositions of a binary at T, SIS Eq. 11.2-2.

    Returns `(xI, xII)`, each a 2-vector of mole fractions, **ordered so that
    xI[0] < xII[0]** -- phase I is the one lean in species 1, which is the order
    Illustration 11.2-2 prints (isobutane x1^I = 0.1128, x1^II = 0.9284). Returns
    **None** when the mixture does not split at T.

    `model` is an `ActivityModel` or a callable `T -> ActivityModel`; `guess` is an
    optional `(x1^I, x1^II)` pair that replaces the hull seed.

    The solve is on `ln(x_i gamma_i)` rather than `x_i gamma_i`, in a logit variable
    so the iterates cannot leave (0, 1). Three things are then checked, and any of
    them failing returns None rather than a plausible pair of numbers:

    1. the two compositions differ by more than `min_gap` (not the trivial root);
    2. the residuals of Eq. 11.2-2 are below `tol`;
    3. the common tangent through the two points lies below the Gibbs energy of
       mixing everywhere (`check_tangent`) -- the guard against converging onto a
       *local* tangent, which satisfies (1) and (2) and is still not the equilibrium
       state.

    ⚠️ **Check 3 has no meaning for a thermodynamically inconsistent model**, because
    such a model has no Gibbs energy for the tangent to lie under. Reproducing SIS
    Illustration 11.2-6's printed table needs `FloryHuggins(..., printed_chi=True)`,
    which is inconsistent by construction, and this function rejects it -- with a
    warning saying so, rather than silently. Pass `check_tangent=False` there, and
    read the answer as "what the book's equations give," not as an equilibrium state.
    """
    if guess is None:
        brackets = common_tangent(model, T, n=n)
        if not brackets:
            return None
        a0, b0 = max(brackets, key=lambda p: p[1] - p[0])
    else:
        a0, b0 = float(guess[0]), float(guess[1])
        if not (0.0 < a0 < 1.0 and 0.0 < b0 < 1.0):
            raise ValueError("guess components must lie strictly inside (0, 1)")

    m = as_model_of_T(model)(T)

    def ln_a(a):
        """ln(x_i gamma_i) for both species at x1 = a."""
        x = np.array([a, 1.0 - a])
        return np.log(x) + np.asarray(m.lngamma(x, T), dtype=float)

    def residual(u):
        a, b = 1.0 / (1.0 + np.exp(-u))
        return ln_a(a) - ln_a(b)

    u0 = np.log(np.array([a0, b0]) / (1.0 - np.array([a0, b0])))
    sol = optimize.root(residual, u0, method="hybr", tol=1e-14)
    a, b = 1.0 / (1.0 + np.exp(-sol.x))

    if abs(a - b) < min_gap:                      # the trivial root
        return None
    if np.max(np.abs(residual(sol.x))) > tol:
        return None
    if a > b:
        a, b = b, a
    if check_tangent and not _tangent_is_global(model, T, a, b):
        warnings.warn(
            f"binary_lle found x1 = {a:.6g} and {b:.6g} at T = {T:.6g} K satisfying "
            f"Eq. 11.2-2, but the common tangent through them does not lie below the "
            f"Gibbs energy of mixing, so they are not an equilibrium state; returning "
            f"None. A thermodynamically inconsistent activity coefficient model will "
            f"always fail this test -- pass check_tangent=False if that is what you "
            f"have.", RuntimeWarning, stacklevel=2)
        return None
    return np.array([a, 1.0 - a]), np.array([b, 1.0 - b])


def _tangent_is_global(model, T, a, b, n=401, slack=1e-9):
    """Does the chord from a to b lie on or below Delta_mix G everywhere?"""
    x = np.linspace(1e-9, 1.0 - 1e-9, n)
    g = gmix_over_RT(model, x, T)
    ga, gb = gmix_over_RT(model, a, T), gmix_over_RT(model, b, T)
    line = ga + (gb - ga) * (x - a) / (b - a)
    return bool(np.min(g - line) > -slack - 1e-12 * max(1.0, abs(ga), abs(gb)))


def binary_lle_envelope(model, T, **kw):
    """The liquid-liquid coexistence (binodal) curve over an array of temperatures.

    Returns `(T_split, x1_I, x1_II)` -- only the temperatures at which two phases
    exist, so the arrays are shorter than `T` whenever part of the range is single
    phase. Pass a callable for `model` when its parameters depend on temperature.

    This is Illustration 11.2-6's table and the phase boundary of its figure, and
    the LLE curve of Figs. 11.2-1, 11.2-2 and 11.3-3.
    """
    Ts, lo, hi = [], [], []
    for Ti in np.atleast_1d(np.asarray(T, dtype=float)):
        got = binary_lle(model, float(Ti), **kw)
        if got is None:
            continue
        Ts.append(float(Ti))
        lo.append(float(got[0][0]))
        hi.append(float(got[1][0]))
    return np.array(Ts), np.array(lo), np.array(hi)


def consolute_temperature(model, T_bounds, upper=True, xtol=1e-6):
    """The consolute (critical solution) temperature, SIS Eqs. 11.2-10a and 11.2-10b.

    Bisects on the *existence of a spinodal* rather than solving
    d2G/dx1^2 = d3G/dx1^3 = 0, because the third derivative of a numerically
    differenced excess Gibbs energy is noise at the level the root needs. Existence
    is a clean yes/no at every temperature and the bisection inherits its accuracy.

    `T_bounds` must bracket the transition: for an upper consolute temperature
    (`upper=True`) the low end must split and the high end must not. Returns the
    temperature; `spinodal` at that temperature gives the critical composition.

    For a one-constant Margules mixture the answer is A/2R, SIS Eq. 11.2-14, and for
    a regular solution it is Problem 11.2-1's formula -- both checked in validation.
    """
    lo, hi = (float(t) for t in T_bounds)
    splits = lambda T: spinodal(model, T) is not None
    if upper:
        if not splits(lo) or splits(hi):
            raise ValueError("T_bounds must bracket an upper consolute temperature: "
                             "the low end must phase-split and the high end must not")
    else:
        if splits(lo) or not splits(hi):
            raise ValueError("T_bounds must bracket a lower consolute temperature: "
                             "the low end must not phase-split and the high end must")
    while hi - lo > xtol:
        mid = 0.5 * (lo + hi)
        if splits(mid) == upper:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Sec. 11.2 -- multicomponent LLE, and the mass balances
# ---------------------------------------------------------------------------
def lle_flash(model, z, T, K0=None, max_iter=500, tol=1e-12, min_gap=1e-4):
    """Two-liquid-phase flash of a feed `z` at T: SIS Eqs. 11.2-2 and 11.2-24 together.

    Returns `(psi, xI, xII)` with `psi` the mole fraction of the feed in phase II, or
    **None** if the feed does not split. Successive substitution on

        K_i = x_i^II / x_i^I = gamma_i^I / gamma_i^II

    with the Rachford-Rice equation for `psi`, which is the same machinery as the
    vapor-liquid flash of SIS Eq. 10.1-7 -- only the source of K differs, which is
    the point Sec. 11.2 makes when it says the two-phase-liquid problem is "six
    coupled equations" for a ternary and is best done on a computer.

    **The initial K cannot come from vapor pressures.** SIS p. 630 is explicit about
    this: in vapor-liquid equilibrium the seed is Pvap_i/P, and here there is no such
    thing, so K is "chosen arbitrarily (for example, K1 = 10 and K2 = 0.1)." That is
    the default. Pass `K0` when a data-based seed is available -- the book's other
    suggestion, and much the better one for a correlation.

    ⚠️ **This solver can converge to the trivial root** where the two phases become
    identical, and it reports None when it does. Successive substitution has no
    defense against it beyond that; `binary_lle` does, and for two species it is the
    function to use. The check here is `min_gap` on the largest |K_i - 1|.
    """
    m = as_model_of_T(model)(T)
    z = np.asarray(z, dtype=float)
    z = z / z.sum()
    nc = z.size
    if K0 is None:
        K = np.full(nc, 0.1)
        K[0] = 10.0
    else:
        K = np.asarray(K0, dtype=float).copy()

    psi = 0.5
    xI = z.copy()
    for _ in range(max_iter):
        psi = _rachford_rice(z, K, psi)
        if psi <= 0.0 or psi >= 1.0:
            return None
        xI = z / (1.0 + psi * (K - 1.0))
        xI = xI / xI.sum()
        xII = K * xI
        xII = xII / xII.sum()
        K_new = (np.asarray(m.gamma(xI, T), dtype=float)
                 / np.asarray(m.gamma(xII, T), dtype=float))
        if np.max(np.abs(K_new / K - 1.0)) < tol:
            K = K_new
            break
        K = K_new

    if np.max(np.abs(K - 1.0)) < min_gap:         # collapsed onto one phase
        return None
    psi = _rachford_rice(z, K, psi)
    if not (0.0 < psi < 1.0):
        return None
    xI = z / (1.0 + psi * (K - 1.0))
    xI = xI / xI.sum()
    xII = K * xI
    return psi, xI, xII / xII.sum()


def _rachford_rice(z, K, beta0=0.5):
    """sum_i z_i (K_i - 1)/(1 + beta (K_i - 1)) = 0, clamped to [0, 1].

    The same equation `phi_phi.PhiPhi._rachford_rice` solves; kept local because
    importing it would pull the equation-of-state module in for a calculation that
    has no equation of state in it.
    """
    z, K = np.asarray(z, dtype=float), np.asarray(K, dtype=float)
    g = lambda b: float(np.sum(z * (K - 1.0) / (1.0 + b * (K - 1.0))))
    if g(0.0) <= 0.0:
        return 0.0
    if g(1.0) >= 0.0:
        return 1.0
    lo = 1.0 / (1.0 - K.max()) + 1e-10
    hi = 1.0 / (1.0 - K.min()) - 1e-10
    return float(optimize.brentq(g, max(lo, 0.0), min(hi, 1.0), xtol=1e-12))


def tie_line_split(z, xI, xII, total=1.0):
    """Amounts of two phases of known composition that make up a feed `z`.

    Returns `(amount_I, amount_II)` summing to `total`. Works in mass or mole
    fractions and for any number of components -- it is SIS Eq. 11.2-1b, the lever
    rule, solved by least squares over **all** the species rather than by picking
    one and hoping the data are consistent.

    That choice is the point. Illustration 11.2-8 balances on water alone and
    Illustration 11.2-1 on both species; when the compositions have been read off a
    triangular diagram by eye, as they are throughout Sec. 11.2, the balances on
    different species do not agree, and which one you pick moves the answer. The
    residual returned by `numpy.linalg.lstsq` is the size of that disagreement, so
    it is worth looking at rather than discarding.
    """
    z = np.asarray(z, dtype=float)
    xI, xII = np.asarray(xI, dtype=float), np.asarray(xII, dtype=float)
    if not (z.shape == xI.shape == xII.shape):
        raise ValueError("z, xI and xII must have the same number of components")
    z, xI, xII = z / z.sum(), xI / xI.sum(), xII / xII.sum()
    A = (xI - xII).reshape(-1, 1)
    b = z - xII
    f = float(np.linalg.lstsq(A, b, rcond=None)[0][0])
    return f * total, (1.0 - f) * total


def mix_streams(amounts, compositions):
    """Combine streams: returns `(total, composition)` of the mixture.

    Illustration 11.2-7 -- and the first step of Illustrations 11.2-8 and 11.2-9,
    where the feed point that gets located on the triangular diagram is the mixture
    of the solvent and the feed. Fraction-agnostic; the illustrations are in weight
    fractions and weight.
    """
    amounts = np.asarray(amounts, dtype=float)
    comps = np.atleast_2d(np.asarray(compositions, dtype=float))
    if comps.shape[0] != amounts.size:
        raise ValueError("one composition vector per stream is required")
    comps = comps / comps.sum(axis=1, keepdims=True)
    total = float(amounts.sum())
    return total, (amounts @ comps) / total


# ---------------------------------------------------------------------------
# Sec. 11.2 -- the sparingly-soluble limit
# ---------------------------------------------------------------------------
def gamma_from_solubility(x_saturated):
    """Activity coefficient of a sparingly soluble species from its solubility.

    SIS Eq. 11.2-18, gamma = 1/x, which follows from Eq. 11.2-2 when the *other*
    phase is essentially the pure species so that x gamma there is 1. Illustrations
    11.2-3 and 11.2-4 are both this one line, and because the solubilities involved
    are so small the numbers it returns are effectively infinite-dilution values.

    ⚠️ Illustration 11.2-4 makes the approximation explicit and then says it is an
    over-estimate: the product x gamma in the *saturated* phase is a little below 1,
    not equal to it, so 1/x is a little high. Slightly, for a solubility of 0.7e-4;
    not at all slightly for a phase that is only 74 mol % pure, which is the same
    illustration's octanol phase.
    """
    return 1.0 / np.asarray(x_saturated, dtype=float)


def solubility_at_T(x1, T1, T2, hex_):
    """Shift a solubility from T1 to T2 at constant excess enthalpy, Eq. 11.2-22.

        x(T2) = x(T1) exp[-(H^ex/R)(1/T2 - 1/T1)]

    `hex_` is the partial molar excess enthalpy of the solute in J/mol, taken
    independent of temperature and composition. Positive H^ex gives solubility
    rising with temperature and negative H^ex gives it falling, which is the whole
    content of the equation.
    """
    return np.asarray(x1, dtype=float) * np.exp(-hex_ / R * (1.0 / T2 - 1.0 / T1))


# ---------------------------------------------------------------------------
# Sec. 11.3 -- vapor-liquid-liquid equilibrium
# ---------------------------------------------------------------------------
def vlle_binary(model, T, gp=None, rtol=1e-6, **kw):
    """Three-phase equilibrium of a binary at T: SIS Eqs. 11.3-2 and 11.3-3.

    Returns `(P, y, xI, xII)` with P in Pa, or **None** if the liquids do not split.
    `model` is a `GammaPhi` (from which the activity model is taken), or an activity
    model with the `GammaPhi` passed separately as `gp`.

    Sec. 11.3's own method, and one worth not automating away: solve the
    liquid-liquid problem, then take a **bubble point on either liquid**, because a
    vapor in equilibrium with one coexisting liquid is in equilibrium with the other.
    This function does both and requires them to agree to `rtol` -- which is not
    decoration, it is the identity x_i^I gamma_i^I = x_i^II gamma_i^II re-emerging
    as a pressure, and Illustration 11.3-1 prints both to make exactly that point.

    ⓘ The Gibbs phase rule leaves one degree of freedom here (C = 2, P = 3), so at
    a given T there is a single three-phase pressure and all three compositions are
    fixed. Changing the overall feed moves mass between the phases and nothing else.
    """
    from .vle import GammaPhi

    if isinstance(model, GammaPhi):
        gp, act = model, model.model
    elif gp is None:
        raise TypeError("pass a GammaPhi, or an activity model together with gp=")
    else:
        act = model

    split = binary_lle(act, T, **kw)
    if split is None:
        return None
    xI, xII = split
    P_I, y_I = gp.bubble_pressure(xI, T)
    P_II, y_II = gp.bubble_pressure(xII, T)
    if abs(P_I - P_II) > rtol * max(P_I, P_II):
        raise RuntimeError(
            f"the two liquids give different bubble pressures, {P_I:.6g} and "
            f"{P_II:.6g} Pa; they must agree, so the liquid-liquid solution is wrong")
    return 0.5 * (P_I + P_II), 0.5 * (y_I + y_II), xI, xII


def vlle_temperature(model, P, T_bounds, gp=None, xtol=1e-8, **kw):
    """The three-phase *temperature* at fixed pressure, for T-x diagrams.

    Returns `(T, y, xI, xII)`. `T_bounds` must bracket the root, and both ends must
    be temperatures at which the liquids still split -- above the consolute
    temperature there is no three-phase state to find, so `consolute_temperature`
    is the natural upper bound.

    This is the VLLE line of SIS Fig. 11.3-4, and the composition `y` at that
    temperature is the heterogeneous azeotrope: the single vapor that both dew-point
    curves run into.
    """
    def f(T):
        got = vlle_binary(model, float(T), gp=gp, **kw)
        if got is None:
            raise ValueError(f"the liquids do not split at T = {T:.6g} K, so there "
                             f"is no three-phase state to bracket")
        return got[0] - P

    T = float(optimize.brentq(f, float(T_bounds[0]), float(T_bounds[1]), xtol=xtol))
    return (T,) + vlle_binary(model, T, gp=gp, **kw)[1:]


def immiscible_pressure(psats, T):
    """Equilibrium pressure and vapor composition over two *immiscible* liquids.

    SIS Eq. 11.3-5: P = sum_i Pvap_i(T), each species behaving as though the other
    were not there. Returns `(P, y)` with P in Pa. `psats` is a list of
    `VaporPressure` objects, one per species -- the same objects `GammaPhi` takes.

    The result is the limit of Eq. 11.3-3 as each phase becomes pure, and it is why
    steam distillation works: the mixture boils below the boiling point of *either*
    pure component, so a compound that would decompose at its own boiling point can
    still be distilled.
    """
    p = np.array([float(ps.P(T)) for ps in psats])
    P = float(p.sum())
    return P, p / P


def steam_distillation(psat_organic, psat_water, T, mw_organic, mw_water=18.015):
    """Kilograms of organic recovered per kilogram of steam condensed.

    Illustration 11.3-3, turpentine at 100 C: `steam_distillation(0.177 bar as a
    correlation, water, 373.15, 140)` returns 1.36. Molecular weights in g/mol; the
    ratio is dimensionless so the units only have to match each other.

    The vapor leaving the kettle is at the composition `immiscible_pressure` gives,
    and the mass ratio is that mole ratio weighted by molecular weight.
    """
    P, y = immiscible_pressure([psat_organic, psat_water], T)
    return float(y[0] * mw_organic / (y[1] * mw_water))


# ---------------------------------------------------------------------------
# Secs. 11.2 and 11.3 -- the equation-of-state route
#
# Everything above solves Eq. 11.2-2, in which the pure-component fugacity has
# cancelled and only activity coefficients remain. An equation of state does not
# work that way: it computes the fugacity of a species in a phase directly, so the
# equilibrium condition stays in its original form,
#
#     x_i^I  phi_i(x^I,  T, P)  =  x_i^II phi_i(x^II, T, P)
#
# and the *only* change from the vapor-liquid calculation of Chapter 9 is that the
# liquid root of the cubic is taken in both phases -- the one sentence SIS p. 630
# spends on it. These functions take any object with `ln_phi(x, T, P, phase=...)`,
# so `PRMixture` and the Wong-Sandler mixing rule of Sec. 9.9 both work.
# ---------------------------------------------------------------------------
def _eos_g_over_RT(eos, x1, T, P, phase="liquid"):
    """SUM_i x_i ln(x_i phi_i) for a binary at T and P.

    This is the molar Gibbs energy of the phase over RT **to within a linear
    function of x1** -- the pure-component reference terms, which are linear and
    therefore invisible to a common-tangent construction. Returns `nan` where the
    cubic has no usable root, so those points can be dropped rather than poison the
    hull.
    """
    out = []
    for a in np.atleast_1d(np.asarray(x1, dtype=float)):
        if a <= 0.0 or a >= 1.0:
            out.append(np.nan)
            continue
        x = np.array([a, 1.0 - a])
        try:
            lnphi = np.asarray(eos.ln_phi(x, T, P, phase=phase), dtype=float)
        except (ValueError, FloatingPointError):
            out.append(np.nan)
            continue
        g = float(np.sum(x * (np.log(x) + lnphi)))
        out.append(g if np.isfinite(g) else np.nan)
    return np.array(out)


def eos_binary_lle(eos, T, P, guess=None, n=201, tol=1e-11, min_gap=1e-4,
                   check_tangent=True):
    """Coexisting liquid compositions of a binary from an equation of state.

    Returns `(xI, xII)` ordered so that `xI[0] < xII[0]`, or **None** where the
    liquid does not split at this T and P. This is Illustration 11.2-5's
    calculation -- CO2 + n-decane from Peng-Robinson with k12 = 0.114.

    ⚠️ **Pressure is an argument here and is not in `binary_lle`.** With activity
    coefficients the liquid is taken incompressible and P drops out; an equation of
    state has no such shortcut, and near a liquid-liquid critical point -- which for
    CO2 + n-decane is only a few kelvin above the range of interest -- the
    compositions move with pressure. Ask for the three-phase state instead of
    choosing P yourself with `eos_vlle_binary`.

    The seed comes from the lower convex hull of the phase's own Gibbs energy, the
    same construction `binary_lle` uses and for the same reason: x^I = x^II solves
    the equilibrium condition at every T and P, so a solver started from a guess
    finds it and reports success.
    """
    if guess is None:
        xs = np.linspace(1e-6, 1.0 - 1e-6, n)
        g = _eos_g_over_RT(eos, xs, T, P)
        ok = np.isfinite(g)
        if ok.sum() < 3:
            return None
        xs, g = xs[ok], g[ok]
        gaps = _hull_gaps(xs, g)
        if not gaps:
            return None
        i, j = max(gaps, key=lambda p: xs[p[1]] - xs[p[0]])
        a0, b0 = float(xs[i]), float(xs[j])
    else:
        a0, b0 = float(guess[0]), float(guess[1])

    def ln_f(a):
        x = np.array([a, 1.0 - a])
        return np.log(x) + np.asarray(eos.ln_phi(x, T, P, phase="liquid"), dtype=float)

    def residual(u):
        a, b = 1.0 / (1.0 + np.exp(-u))
        return ln_f(a) - ln_f(b)

    u0 = np.log(np.array([a0, b0]) / (1.0 - np.array([a0, b0])))
    sol = optimize.root(residual, u0, method="hybr", tol=1e-14)
    a, b = 1.0 / (1.0 + np.exp(-sol.x))

    if abs(a - b) < min_gap:
        return None
    if np.max(np.abs(residual(sol.x))) > tol:
        return None
    if a > b:
        a, b = b, a
    if check_tangent:
        xs = np.linspace(1e-6, 1.0 - 1e-6, 401)
        g = _eos_g_over_RT(eos, xs, T, P)
        ok = np.isfinite(g)
        ga, gb = float(_eos_g_over_RT(eos, a, T, P)[0]), float(_eos_g_over_RT(eos, b, T, P)[0])
        line = ga + (gb - ga) * (xs[ok] - a) / (b - a)
        if np.min(g[ok] - line) < -1e-9:
            warnings.warn(
                f"eos_binary_lle found x1 = {a:.6g} and {b:.6g} at T = {T:.6g} K, "
                f"P = {P:.6g} Pa satisfying equal fugacities, but the common tangent "
                f"through them does not lie below the Gibbs energy of the phase, so "
                f"they are not an equilibrium state; returning None.",
                RuntimeWarning, stacklevel=2)
            return None
    return np.array([a, 1.0 - a]), np.array([b, 1.0 - b])


def eos_vlle_binary(eos, T, P_guess=None, guess=None, tol=1e-10, min_gap=1e-4):
    """Three-phase equilibrium of a binary from an equation of state, at fixed T.

    Returns `(P, y, xI, xII)` with P in Pa, or **None** if the liquids do not split.
    This is Illustration 11.3-2 -- CO2 + n-decane, Peng-Robinson, k12 = 0.114 --
    and it is the calculation the 5e sends the reader to a Visual Basic program, a
    DOS program or a Mathcad worksheet to do.

    ⚠️ **Not two two-phase calculations.** `vlle_binary` can solve the liquid-liquid
    problem first and take a bubble point afterwards, because with activity
    coefficients the liquid compositions do not depend on pressure. Here they do, so
    all four unknowns -- both liquid compositions, the vapor composition and the
    pressure -- are solved together:

        x_i^I phi_i(x^I) = x_i^II phi_i(x^II)      both liquid roots
        x_i^I phi_i(x^I) = y_i    phi_i(y)         vapor root

    Four equations, four unknowns, and the Gibbs phase rule agrees: C = 2 and
    P = 3 leave one degree of freedom, which T has just used.

    ⓘ Doing it as two separate calculations instead -- LLE at a guessed pressure,
    then a bubble point on one liquid -- is what the seeding below does, and for
    CO2 + n-decane it lands within 0.001 bar of the simultaneous answer because the
    liquid compositions barely move over the quarter-bar between the seed and the
    result. That is a fact about this system, not about the method.

    ⛔ **There is a trivial root here too, and it is not the one `binary_lle`
    guards against.** `y = x^II` satisfies the vapor equation *exactly* wherever
    the cubic has a single real root, because then the "vapor" root and the liquid
    root are the same number and phi^V(y) = phi^L(y) identically. Above the
    three-phase pressure that is the situation, so a solve seeded at, say, 30 bar
    converges to P = 30 bar with a "vapor" that is the CO2-rich liquid -- a
    perfectly convergent answer with three liquids in it. The check below rejects
    it. ⚠️ **`PhiPhiVLE.bubble_pressure` has the same hole and is not guarded**:
    given an LLE composition and a seed above the three-phase pressure it returns
    the liquid-liquid state as a bubble point. Its default seed is far below, so
    the ordinary path is unaffected, and that is why this function does not lean on
    it. See `revision_notes/c11.md` Section 13.
    """
    if P_guess is None:
        P_guess = max(float(c.vapor_pressure(T)) for c in eos.components)

    if guess is None:
        split = eos_binary_lle(eos, T, P_guess)
        if split is None:
            return None
        xI0, xII0 = split
        P0, y0 = eos.bubble_pressure(xI0, T, P_guess=P_guess)
    else:
        xI0, xII0, y0, P0 = guess

    def unpack(u):
        a, b, c = 1.0 / (1.0 + np.exp(-u[:3]))
        return (np.array([a, 1.0 - a]), np.array([b, 1.0 - b]),
                np.array([c, 1.0 - c]), float(np.exp(u[3])))

    def residual(u):
        xI, xII, y, P = unpack(u)
        lnf = lambda x, ph: np.log(x) + np.asarray(
            eos.ln_phi(x, T, P, phase=ph), dtype=float)
        fI = lnf(xI, "liquid")
        return np.concatenate([fI - lnf(xII, "liquid"), fI - lnf(y, "vapor")])

    logit = lambda v: np.log(v / (1.0 - v))
    u0 = np.array([logit(xI0[0]), logit(xII0[0]), logit(y0[0]), np.log(P0)])
    sol = optimize.root(residual, u0, method="hybr", tol=1e-14)
    if np.max(np.abs(residual(sol.x))) > tol:
        return None
    xI, xII, y, P = unpack(sol.x)
    if xI[0] > xII[0]:
        xI, xII = xII, xI
    if abs(xI[0] - xII[0]) < min_gap:
        return None                                   # the two liquids collapsed
    for phase, x in (("I", xI), ("II", xII)):         # ... and the vapor onto one
        if np.max(np.abs(y - x)) < min_gap:
            warnings.warn(
                f"eos_vlle_binary converged at T = {T:.6g} K to P = {P:.6g} Pa with "
                f"a vapor indistinguishable from liquid {phase}; that is the trivial "
                f"root of the vapor equation, which every pressure above the "
                f"three-phase pressure admits. Returning None -- try a lower "
                f"P_guess, the default being the largest pure-component vapor "
                f"pressure.", RuntimeWarning, stacklevel=2)
            return None
    return P, y, xI, xII


# ---------------------------------------------------------------------------
# Sec. 11.3 -- the phase diagrams
# ---------------------------------------------------------------------------
def pxy_lle(gp, T, n=201, **kw):
    """P-x-y for a binary at fixed T, **with** the liquid-liquid region resolved.

    Returns a dict with

        x1, y1, P       the bubble-point curve ignoring LLE  (SIS Fig. 11.3-2, solid)
        xI, xII, P_vlle the three-phase tie line, or None     (Fig. 11.3-2, dashed)
        y_vlle          the vapor composition on that line

    Both curves, from one call, because the figure's argument needs both: computed
    as vapor-liquid equilibrium alone, the isobutane-furfural system shows a
    maximum *and* a minimum in pressure and looks like a double azeotrope, which
    would be extraordinary. It is not one. It is two liquid phases, and the true
    pressure across that composition range is the flat dashed line. SIS Sec. 11.3
    turns that into a working rule: **a computed P-x curve with both a maximum and
    a minimum is a signal to go and do a liquid-liquid calculation.**
    """
    from .vle import pxy

    x1, y1, P = pxy(gp, T, n=n)
    out = {"x1": x1, "y1": y1, "P": P,
           "xI": None, "xII": None, "P_vlle": None, "y_vlle": None}
    got = vlle_binary(gp, T, **kw)
    if got is not None:
        P3, y3, xI, xII = got
        out.update(xI=xI, xII=xII, P_vlle=P3, y_vlle=y3)
    return out


def txy_lle(gp, P, T_bounds, n=201, n_lle=101, **kw):
    """T-x-y for a binary at fixed P, with the LLE envelope and the three-phase line.

    Returns a dict with

        x1, y1, T           the bubble- and dew-point curves
        T_lle, xI, xII      the liquid-liquid envelope over `T_bounds`
        T_vlle, y_vlle      the three-phase temperature and vapor, or None

    SIS Fig. 11.3-4. `T_bounds` sets the temperature window the LLE envelope is
    swept over and is also the bracket for the three-phase temperature; the LLE
    curve is pressure-independent (liquids are not compressible, SIS Sec. 11.3) so
    that window is a drawing choice, while the VLE curves are not.
    """
    from .vle import txy

    x1, y1, T = txy(gp, P, n=n)
    T_grid = np.linspace(float(T_bounds[0]), float(T_bounds[1]), n_lle)
    T_lle, xI, xII = binary_lle_envelope(gp.model, T_grid, **kw)
    out = {"x1": x1, "y1": y1, "T": T,
           "T_lle": T_lle, "xI": xI, "xII": xII,
           "T_vlle": None, "y_vlle": None}
    try:
        T3, y3, _, _ = vlle_temperature(gp, P, T_bounds, **kw)
        out.update(T_vlle=T3, y_vlle=y3)
    except (ValueError, RuntimeError):
        pass
    return out
