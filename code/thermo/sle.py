"""sle -- solid-liquid and solid-fluid equilibrium, SIS Secs. 12.1 and 12.3.

**The whole of Sec. 12.1 and Sec. 12.3 is one equation solved for two different
unknowns.** Sec. 12.1 asks *how much solid dissolves at this temperature* and solves
Eq. 12.1-7 for x; Sec. 12.3 asks *at what temperature does solid first appear from this
liquid* and solves Eq. 12.3-2 for T. Those two printed equations are the same function
of (T, T_m, dH_fus, dCp), with the triple-point temperature of one playing the role of
the normal melting point of the other -- which Sec. 12.1 says outright in the sentence
that introduces Eq. 12.1-7. So there is exactly one of them in here, `ln_x_gamma`, and
everything else is that function inverted one way or the other:

    ln (x_1 gamma_1) = - dH_fus(T_m)/R * (T_m - T)/(T_m T)
                       - dCp/R * [1 - T_m/T + ln(T_m/T)]

    from thermo.sle import ideal_solubility, solubility, freezing_point

    ideal_solubility(293.15, 353.35, 18804.0)        # Ill. 12.1-1, the ideal answer
    solubility(293.15, 353.35, 18804.0, gamma)       # Ill. 12.1-1 / 12.1-2, with gamma(x)
    freezing_point(0.9, 178.16, 6610.7, dCp=48.6)    # Ill. 12.3-2, toluene as the solid

## What Secs. 12.1 and 12.3 ask for, and where each piece is

| SIS | what | entry point |
|-----|------|-------------|
| Eqs. 12.1-6, 12.1-7, 12.3-2 | the one equation, as ln(x_1 gamma_1) | `ln_x_gamma` |
| Eq. 12.1-8 | ideal solubility of a solid in a liquid (dCp = 0) | `ideal_solubility` |
| Eq. 12.1-7 with gamma(x) | the real solubility -- a fixed point, not a formula | `solubility` |
| Eq. 12.1-6 rearranged | gamma_1 from a *measured* solubility | `activity_coefficient` |
| Eqs. 12.1-11 to 12.1-13 | delta_1 of a subcooled liquid, from sublimation data | `solubility_parameter` |
| Eq. 7.7-5a | dH_sub from the slope of log10 P^sub vs 1/T | `heat_of_sublimation` |
| Eq. 12.1-19b | dH_fus from two solubilities | `heat_of_fusion` |
| Eq. 12.1-16 | ideal solubility of a solid in a *gas* | `ideal_solubility_in_gas` |
| Eq. 12.1-17 | solubility in a gas or SCF, given phi-bar | `solubility_in_gas` |
| Eq. 12.1-18a | the enhancement factor E, and its Poynting half | `enhancement_factor`, `poynting` |
| Eq. 12.3-2 solved for T | the freezing point of a liquid mixture | `freezing_point` |
| Eq. 12.3-5 | the dilute-solute freezing-point depression | `freezing_point_depression` |
| Sec. 12.3 | the eutectic, as the crossing of two freezing curves | `eutectic` |

## The trap: gamma depends on x, so the solubility is a fixed point

This is what the highlighted sentence in Illustration 12.1-2 is actually about. UNIFAC
(or regular solution theory, or anything else) returns gamma *given* a composition, but
the composition is what is being solved for. So Eq. 12.1-7 is not evaluated, it is
iterated:

    x <- exp(ln_x_gamma(...)) / gamma(x)

`solubility` does that and **returns the whole iterate history**, because the book
prints it -- Illustration 12.1-1's "the results of the next two iterations are
x_1 = 0.0768 and x_1 = 0.0772" and Illustration 12.1-2's 0.07 -> 0.0821 -> 0.0856 are
the sequence, not just its limit, and a reader checking the book needs to see the same
steps. **On failure to converge it raises** rather than returning the last iterate:
a fixed point that has not converged is not an answer, and returning it silently is how
a wrong number reaches print.

## Units

`ln_x_gamma` and everything built on it are in J/mol and K, so R is 8.314 and dCp is
J/(mol K) -- the units Secs. 12.1 and 12.3 print. The gas-side functions
(`ideal_solubility_in_gas`, `poynting`, `enhancement_factor`, `solubility_in_gas`) are
unit-agnostic in pressure as long as P and P^sat agree, except that `poynting` takes a
molar volume and so needs pressure in Pa and volume in m^3/mol; pass `R=83.14` and
cc/mol and bar if you would rather work the way Illustration 12.1-5 prints.

Eric Furst
August 2026
"""
import numpy as np
from scipy import constants, optimize

R = constants.R

__all__ = ["ln_x_gamma", "ideal_solubility", "solubility", "activity_coefficient",
           "heat_of_sublimation", "solubility_parameter", "heat_of_fusion",
           "ideal_solubility_in_gas", "poynting", "enhancement_factor",
           "solubility_in_gas", "freezing_point", "freezing_point_depression",
           "eutectic", "Iteration"]


class Iteration:
    """The result of a fixed-point solve: the answer, and how it got there.

    `x` is the converged value, `history` every iterate including the initial guess,
    and `gammas` the activity coefficient evaluated at each of them. `direction` is
    "rising", "falling" or "non-monotone" -- recorded because a sweep's *direction* is
    part of what makes it checkable, and a solve that walked the wrong way while
    landing on a plausible number is a bug that agreement hides.
    """

    def __init__(self, x, history, gammas):
        self.x = float(x)
        self.history = list(map(float, history))
        self.gammas = list(map(float, gammas))
        d = np.diff(self.history)
        self.direction = ("rising" if np.all(d > 0) else
                          "falling" if np.all(d < 0) else "non-monotone")

    @property
    def iterations(self):
        return len(self.history) - 1

    def __float__(self):
        return self.x

    def __repr__(self):
        return (f"Iteration(x={self.x:.6g}, gamma={self.gammas[-1]:.6g}, "
                f"iterations={self.iterations}, direction={self.direction!r})")


# --- the one equation -------------------------------------------------------

def ln_x_gamma(T, T_m, dH_fus, dCp=0.0):
    """ln(x_1 gamma_1) at equilibrium with the pure solid -- SIS Eq. 12.1-6 / 12.3-2.

        ln(x_1 gamma_1) = - dH_fus/R * (T_m - T)/(T_m T) - dCp/R [1 - T_m/T + ln(T_m/T)]

    Sec. 12.1 writes the first term as -dH_fus/(R T) (1 - T/T_m), which is the same
    thing, and calls T_m the triple-point temperature T_t; Eq. 12.1-7's own
    introduction is the argument that the two are interchangeable here. Setting
    dCp = 0 gives Eq. 12.1-8.

    Parameters
    ----------
    T : temperature of the mixture, K.
    T_m : normal melting (or triple-point) temperature of the solute, K.
    dH_fus : heat of fusion at T_m, J/mol.
    dCp : Cp(liquid) - Cp(solid), J/(mol K), taken independent of temperature.
        *That subtraction order matters and is easy to get backwards.*
        Illustration 12.3-2 prints 48.6 for toluene, which is its Cp^L of 135.6 less
        its Cp^S of 87.0 -- liquid minus solid.
    """
    T = np.asarray(T, float)
    fusion = -dH_fus / R * (T_m - T) / (T_m * T)
    if dCp == 0.0:
        return fusion
    return fusion - dCp / R * (1.0 - T_m / T + np.log(T_m / T))


# --- Sec. 12.1: solved for the composition ----------------------------------

def ideal_solubility(T, T_m, dH_fus, dCp=0.0):
    """Saturation mole fraction of a solid in a liquid for gamma_1 = 1 -- Eq. 12.1-8.

    Illustration 12.1-1's Comment is the reason this is worth a function of its own:
    for naphthalene in n-hexane at 20 C it returns 0.269 against a measured 0.09, "a
    factor of 3 too large," and that failure is the whole motivation for Sec. 12.1's
    activity-coefficient machinery.
    """
    return np.exp(ln_x_gamma(T, T_m, dH_fus, dCp))


def solubility(T, T_m, dH_fus, gamma, dCp=0.0, x0=None, tol=1e-10, max_iter=200):
    """Saturation mole fraction with a composition-dependent activity coefficient.

    Eq. 12.1-7 solved as the fixed point x <- x_ideal / gamma(x), which is the
    iteration Illustrations 12.1-1 and 12.1-2 both print.

    Parameters
    ----------
    gamma : callable x -> gamma_1(x). Anything that maps the solute mole fraction to
        its activity coefficient: a regular-solution model, a UNIFAC call, a lambda
        closing over a fitted correlation.
    x0 : first guess. Defaults to the ideal solubility, which is the natural start and
        is what "assume that x_1 will be small" amounts to; pass 0.07 to reproduce
        Illustration 12.1-2's printed sequence exactly.

    Returns
    -------
    Iteration -- `.x` is the answer, `.history` the printed sequence.

    Raises
    ------
    RuntimeError if the iteration has not converged in `max_iter` steps. It does not
        return the last iterate; see the module docstring.
    """
    x_ideal = float(ideal_solubility(T, T_m, dH_fus, dCp))
    x = x_ideal if x0 is None else float(x0)
    history, gammas = [x], []
    for _ in range(max_iter):
        g = float(gamma(x))
        if not np.isfinite(g) or g <= 0.0:
            raise RuntimeError(f"gamma({x:.6g}) returned {g!r}; the activity "
                               f"coefficient must be finite and positive")
        gammas.append(g)
        x_new = x_ideal / g
        history.append(x_new)
        if abs(x_new - x) < tol * max(1.0, abs(x_new)):
            return Iteration(x_new, history, gammas + [float(gamma(x_new))])
        x = x_new
    raise RuntimeError(
        f"the solubility fixed point did not converge in {max_iter} iterations; "
        f"last iterates {history[-4:]}. Successive substitution diverges when "
        f"|d ln gamma / d ln x| > 1 near the root -- damp the step or bracket the "
        f"root instead. The last iterate is NOT the answer.")


def activity_coefficient(x_sat, T, T_m, dH_fus, dCp=0.0):
    """gamma_1 of a dissolved solid from its *measured* saturation solubility.

    Eq. 12.1-6 rearranged, and Illustration 12.1-3 is the point of it: for
    benzo[a]pyrene in water at 25 C, from x = 3.37e-10, it returns 3.74e8. Because the
    solubility is so small this is also the infinite-dilution value, which is what
    Sec. 12.5 then uses to get an air-water partition coefficient.
    """
    return np.exp(ln_x_gamma(T, T_m, dH_fus, dCp)) / np.asarray(x_sat, float)


def heat_of_sublimation(slope_log10, T=None):
    """dH_sub from the slope of log10 P^sub against 1/T -- Clausius-Clapeyron, Eq. 7.7-5a.

    A sublimation pressure reported as log10 P^sub = A - B/T has
    dH_sub = 2.303 R B, since d ln P/d(1/T) = -dH_sub/R. Pass B (a positive number,
    K) as `slope_log10`.

    Illustration 12.1-1's data line is log10 P^sub(bar) = 8.722 - 3783/T, so
    `heat_of_sublimation(3783.0)` is its 72.4 kJ/mol.

    `T` is accepted and ignored; it is there so a caller can record the temperature the
    slope was taken at without keeping it in a separate variable.
    """
    return np.log(10.0) * R * np.asarray(slope_log10, float)


def solubility_parameter(dH_sub, dH_fus, V_liquid, T):
    """delta_1 of the *subcooled liquid* solute -- Eqs. 12.1-11 to 12.1-13.

        dH_vap(subcooled liquid) = dH_sub - dH_fus         (Eq. 12.1-12)
        dU_vap = dH_vap - RT
        delta_1 = sqrt(dU_vap / V_liquid)                  (Eq. 12.1-11)

    **This is the awkward step in regular solution theory and Sec. 12.1 says so**:
    the solubility parameter wanted is that of a liquid which, at the temperature of
    interest, does not exist. It is reached from the solid's sublimation pressure and
    its heat of fusion, which is why `heat_of_sublimation` sits next to it.

    Units: dH in J/mol and V_liquid in m^3/mol gives delta in Pa^(1/2). To get the
    (cal/cc)^(1/2) the book prints, divide by `thermo.activity_models.CAL_CC_HALF`.
    """
    dU_vap = (np.asarray(dH_sub, float) - np.asarray(dH_fus, float)
              - R * np.asarray(T, float))
    return np.sqrt(dU_vap / np.asarray(V_liquid, float))


def heat_of_fusion(x1, T1, x2, T2):
    """Apparent dH_fus from two solubilities at two temperatures -- Eq. 12.1-19b.

        dH_fus = R T1 T2 / (T1 - T2) * ln(x1/x2)

    Illustration 12.1-7 is the reason the word *apparent* is in Sec. 12.1's own
    sentence: the activity coefficients have been assumed to cancel in the ratio, so
    what comes out carries their temperature and composition dependence with it.

    **Solubilities in any units that are proportional to mole fraction will do.**
    Illustration 12.1-7 hands it mg/mL directly, and its own justification is that at
    a molecular weight of 34 800 the mole fraction is 9.4e-8, where x is linear in S.
    """
    T1, T2 = float(T1), float(T2)
    return R * T1 * T2 / (T1 - T2) * np.log(np.asarray(x1, float)
                                            / np.asarray(x2, float))


# --- Sec. 12.1: the solid in a gas or a supercritical fluid -----------------

def ideal_solubility_in_gas(P_sub, P):
    """y_i^ID = P_i^sub(T) / P -- Eq. 12.1-16.

    The low-pressure limit, in which the Poynting factor and every fugacity
    coefficient is unity. Illustration 12.1-4 is this and nothing else.
    """
    return np.asarray(P_sub, float) / np.asarray(P, float)


def poynting(V_solid, P, P_sat, T, R=R):
    """exp[V^S (P - P^sat) / RT], the Poynting factor of Eqs. 12.1-15 and 12.1-17.

    Pass SI (m^3/mol, Pa) or the book's units with `R=83.14` (cc bar / mol K), cc/mol
    and bar. Illustration 12.1-6's fourth Comment is worth reading before dropping it:
    the factor reaches 3.5 there, so it is not a correction that can be neglected at
    supercritical pressures.
    """
    V_solid, P, P_sat, T = (np.asarray(v, float) for v in (V_solid, P, P_sat, T))
    return np.exp(V_solid * (P - P_sat) / (R * T))


def enhancement_factor(phi_bar, poynting_factor, f_over_P=1.0):
    """E = (f/P)_sat,T * Poynting / phi-bar_i^V -- Eq. 12.1-18a.

    E -> 1 as P -> P^sat, and Sec. 12.1's own note is that it collects two distinct
    effects: the Poynting factor and the vapor-phase fugacity coefficient. Illustration
    12.1-6's second Comment reports E > 17 700, "among the larger nonideal corrections
    encountered in chemical engineering thermodynamics."

    `f_over_P` is the fugacity coefficient of the *pure saturated solid*, which
    Sec. 12.1 says is usually unity because the sublimation pressure is small.
    """
    return (np.asarray(f_over_P, float) * np.asarray(poynting_factor, float)
            / np.asarray(phi_bar, float))


def solubility_in_gas(P_sat, P, phi_bar, V_solid, T, f_over_P=1.0, y0=None,
                      R_gas=R, tol=1e-12, max_iter=200, method="bracket",
                      y_max=0.8, n_scan=300, all_roots=False):
    """y_i of a solid in a gas or supercritical fluid -- Eq. 12.1-17.

        y_i = P_i^sat (f/P)_sat exp[V^S(P - P_i^sat)/RT] / (P phi-bar_i^V(T, P, y))

    **phi-bar depends on y, so this has to be solved, not evaluated** -- Sec. 12.1
    says so in the sentence after Eq. 12.1-17, and Illustration 12.1-5's procedure
    paragraph is the hand version. `phi_bar` is a callable y -> phi-bar_i^V.

    **AND THE EQUATION HAS MORE THAN ONE ROOT.** For CO2 + naphthalene at 60.4 C
    and 133.8 bar it has *three*: 1.05e-2, 0.124 and 0.542. Only the smallest is the
    solubility of a trace solid in a supercritical fluid; the others are the same
    equality of fugacities satisfied by a dense, naphthalene-rich phase, which is not
    what Illustration 12.1-6 is asking about. A bracketed solve over (0, 1) finds the
    LARGEST and reports it as success -- so:

      - ``method="bracket"`` (default) scans **upward** from y ~ 0 and returns the
        FIRST crossing. Robust, and right.
      - ``method="substitution"`` is the book's hand loop, y <- y^ID E(y). It
        reproduces Illustration 12.1-5's printed procedure and stays near the physical
        root by construction, but it diverges once the two lowest roots approach each
        other -- which they do above about 190 bar at 60.4 C.

    Pass ``all_roots=True`` to get every crossing instead of the first, which is how to
    see the branches merge as pressure rises: when the list drops from three entries to
    one, the physical branch has ceased to exist for these parameters.

    Returns an `Iteration` (or a list of roots if `all_roots`); raises rather than
    returning an unconverged iterate.

    Units: pass `R_gas=83.14` with cc/mol, bar and K to work the way Sec. 12.1 prints.
    """
    E_pressure = float(poynting(V_solid, P, P_sat, T, R=R_gas))
    y_ideal = float(np.asarray(P_sat, float) / np.asarray(P, float))
    f_solid = y_ideal * float(f_over_P) * E_pressure         # y^ID (f/P) Poynting

    def residual(y):
        phi = float(phi_bar(y))
        if not np.isfinite(phi) or phi <= 0.0:
            raise RuntimeError(f"phi_bar({y:.6g}) returned {phi!r}; a fugacity "
                               f"coefficient must be finite and positive")
        return y * phi - f_solid

    if method == "bracket":
        grid = np.logspace(-14, np.log10(y_max), n_scan)
        f = np.array([residual(y) for y in grid])
        cross = np.nonzero(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
        if cross.size == 0:
            raise RuntimeError(
                f"Eq. 12.1-17 has no root in (0, {y_max}] at P = {float(P):.4g}: the "
                f"residual does not change sign. For a near-critical mixture this is "
                f"real, not a numerical failure -- the low-solubility branch can cease "
                f"to exist as the branches merge with rising pressure. Call again with "
                f"all_roots=True at a lower pressure to watch them approach.")
        roots = [float(optimize.brentq(residual, grid[i], grid[i + 1],
                                       xtol=1e-16, rtol=1e-13)) for i in cross]
        if all_roots:
            return roots
        y = roots[0]
        return Iteration(y, [grid[cross[0]], y], [float(phi_bar(y))])

    if method != "substitution":
        raise ValueError("method must be 'bracket' or 'substitution'")

    y = f_solid if y0 is None else float(y0)
    history, phis = [y], []
    for _ in range(max_iter):
        phi = float(phi_bar(y))
        if not np.isfinite(phi) or phi <= 0.0:
            raise RuntimeError(f"phi_bar({y:.6g}) returned {phi!r}; a fugacity "
                               f"coefficient must be finite and positive")
        phis.append(phi)
        y_new = float(y_ideal * enhancement_factor(phi, E_pressure, f_over_P))
        history.append(y_new)
        if abs(y_new - y) < tol * max(1.0, abs(y_new)):
            return Iteration(y_new, history, phis + [float(phi_bar(y_new))])
        y = y_new
    raise RuntimeError(
        f"the vapor-phase solubility did not converge in {max_iter} iterations; "
        f"last iterates {history[-4:]}. Successive substitution fails where the two "
        f"lowest roots of Eq. 12.1-17 approach each other; use method='bracket'. "
        f"The last iterate is NOT the answer.")


# --- Sec. 12.3: the same equation solved for the temperature -----------------

def freezing_point(x1, T_m, dH_fus, dCp=0.0, gamma=1.0, T_min=None, n_scan=400):
    """The temperature at which pure solid 1 first appears from a liquid -- Eq. 12.3-2.

    `ln_x_gamma(T_f, ...) = ln(x_1 gamma_1)` solved for T_f. Same equation as
    `solubility`, other unknown.

    Parameters
    ----------
    x1 : liquid mole fraction of the species that freezes out (the solvent).
    gamma : its activity coefficient -- a number, or a callable T -> gamma_1(T) if the
        liquid model is temperature-dependent. Illustration 12.3-2 uses 1.0, and says
        why: toluene and ethyl benzene have solubility parameters of 8.9 and 8.8, so
        the liquid mixture is very nearly ideal.
    T_min : lower end of the search. Defaults to 0.2 T_m, which spans every freezing
        point Illustration 12.3-2 reports (its lowest is 0.56 T_m).

    **WITH dCp != 0 THIS EQUATION HAS TWO ROOTS, AND ONLY THE UPPER ONE IS
    PHYSICAL.** As T -> 0 the dCp term goes as +dCp T_m/(R T) while the fusion term
    goes as -dH_fus/(R T), so whenever dCp T_m > dH_fus -- which is the case for *both*
    species in Illustration 12.3-2 (48.6 x 178.16 = 8659 against 6611) -- the residual
    turns around and diverges to +infinity, giving a second, spurious root at low
    temperature. A plain bracket of (0, T_m] therefore either finds the wrong root or
    reports no sign change at all, depending on which end it starts from. This function
    scans **down from T_m** and takes the first crossing, which is the physical one:
    the freezing point is the highest temperature at which solid can appear.

    That is also why Sec. 12.3 warns "since dCp is so large, none of the terms in this
    equation can be neglected" -- the dCp term is not a small correction here, it
    changes the shape of the equation.
    """
    x1 = float(x1)
    if not 0.0 < x1 <= 1.0:
        raise ValueError(f"x1 must be in (0, 1]; got {x1}")
    T_m = float(T_m)
    if x1 == 1.0:
        return T_m
    lo = 0.2 * T_m if T_min is None else float(T_min)

    def residual(T):
        g = gamma(T) if callable(gamma) else gamma
        return float(ln_x_gamma(T, T_m, dH_fus, dCp) - np.log(x1 * g))

    grid = np.linspace(T_m, lo, n_scan)              # downward, so the first crossing wins
    f = np.array([residual(T) for T in grid])
    crossings = np.nonzero(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
    if crossings.size == 0:
        raise ValueError(
            f"Eq. 12.3-2 has no root on [{lo:.1f}, {T_m:.1f}] K for x1 = {x1}: the "
            f"residual runs from {f[0]:.4g} to {f[-1]:.4g} without changing sign. "
            f"Lower T_min, or check that x1 is the mole fraction of the species assumed "
            f"to freeze out -- passing the *other* component's mole fraction is the "
            f"usual cause.")
    i = int(crossings[0])
    return float(optimize.brentq(residual, grid[i + 1], grid[i], xtol=1e-10))


def freezing_point_depression(x2, T_m, dH_fus):
    """dT = T_m - T_f = R T_m^2 x_2 / dH_fus -- Eq. 12.3-5.

    The dilute limit, in which gamma_1 -> 1 and ln x_1 = ln(1 - x_2) -> -x_2.
    **It contains nothing about the solute but its mole fraction**, which is
    Sec. 12.3's own observation and the reason it can be compared term by term with
    the osmotic pressure of Eq. 11.5-4. `x2` may be a sum over several solutes, which
    is Eq. 12.3-6.
    """
    return R * float(T_m) ** 2 * np.asarray(x2, float) / float(dH_fus)


def eutectic(x_grid, T_a, T_b):
    """The eutectic, as the crossing of two freezing curves.

    Sec. 12.3's own procedure, stated as an algorithm: compute the freezing point on
    the assumption that species A freezes out, then again for B, and at each
    composition the real freezing point is **the higher of the two**, because that is
    the solid that appears first. Where the two curves cross, both solids appear
    together -- the eutectic.

    Parameters
    ----------
    x_grid : compositions, as the mole fraction of the species whose curve is `T_a`.
    T_a, T_b : the two freezing curves on that grid. NaN is allowed and is skipped;
        `freezing_point` legitimately has no solution at the ends.

    Returns
    -------
    (x_eutectic, T_eutectic) -- found by linear interpolation on the difference
    T_a - T_b across its sign change, so the answer is not restricted to a grid point.

    **The eutectic is the MINIMUM of the freezing curve, not of either branch.**
    Reporting min(T_a) or min(T_b) gives a temperature far below the real one; on
    Illustration 12.3-2's grid min(T_a) is 99.8 K against a eutectic of 155.9 K.
    """
    x = np.asarray(x_grid, float)
    d = np.asarray(T_a, float) - np.asarray(T_b, float)
    ok = np.isfinite(d)
    x, d = x[ok], d[ok]
    Ta = np.asarray(T_a, float)[ok]
    sign_change = np.nonzero(np.diff(np.sign(d)))[0]
    if sign_change.size == 0:
        raise ValueError("the two freezing curves do not cross on this grid, so there "
                         "is no eutectic in range; widen x_grid")
    i = int(sign_change[0])
    t = d[i] / (d[i] - d[i + 1])                      # linear in the difference
    return (float(x[i] + t * (x[i + 1] - x[i])),
            float(Ta[i] + t * (Ta[i + 1] - Ta[i])))
