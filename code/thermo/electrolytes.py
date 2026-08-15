"""Electrolyte solutions -- SIS Section 9.10, Debye-Huckel and its extensions.

The mean ionic activity coefficient of a dissolved salt, from the three equations
the section gives:

    ln gamma_pm = -alpha |z+ z-| sqrt(I)                                (Eq. 9.10-15)
    ln gamma_pm = -alpha |z+ z-| sqrt(I) / (1 + beta a sqrt(I))         (Eq. 9.10-17)
    ln gamma_pm = -alpha |z+ z-| sqrt(I) / (1 + beta a sqrt(I)) + delta I
                                                                        (Eq. 9.10-18)

behind one class, because the book presents the last two as *modifications of the
first* rather than as separate models:

    DebyeHuckel(NaCl)                             -> the limiting law, 9.10-15
    DebyeHuckel(NaCl, beta_a=1)                   -> 9.10-17
    DebyeHuckel(NaCl, beta_a=1, delta=0.1)        -> 9.10-18

`model.equation` reports which of the three you built, so a notebook that sweeps the
parameters can label its own curves.

## Why this module is separate from `activity_models`

Everything in Sec. 9.5 is a function of *mole fraction* and returns one activity
coefficient per species. An electrolyte model is a function of *ionic strength* and
returns a single mean ionic coefficient for the salt -- a different independent
variable and a different dependent one, so `ActivityModel.gamma(x, T)` is the wrong
interface and inheriting it would be a lie. The reason for the difference is the
physics of Sec. 9.10 itself: the individual ion activity coefficients are **not
separately measurable**, because electroneutrality (Eq. 9.10-3) forbids varying the
cations while holding the anions fixed. Only gamma_pm is measurable, so only
gamma_pm is what a model can be asked for.

## Units: molality, and the one place the book is loose

Concentrations here are **molalities**, mol per kg of solvent -- not mole fractions
and not molarities. `alpha` is tabulated in (mol/kg)^(-1/2) and `beta` in
[(mol/kg)^(1/2) angstrom]^(-1), so `beta_a` (with a in angstroms) is
(mol/kg)^(-1/2) and beta_a*sqrt(I) is dimensionless.

⚠️ **The book says it is ignoring the molality/molarity distinction** ("in the
application of these formulas the distinction between molality ... and ... molarity
will sometimes be ignored"), and Figs. 9.10-1 and 9.10-2 do exactly that -- 9.10-1's
abscissa is labeled molarity M and 9.10-2's caption calls gamma_pm the mean *molar*
activity coefficient while its abscissa is the square root of the ionic strength
computed from a *molality* table. In dilute aqueous solution the two agree to within
a percent or so and the point of the section survives; at the 16 molal end of
Illustration 9.10-1 they do not. This module is molal throughout and says so, which
is the only way the arithmetic can be checked.

## Reused by Chapter 15

Ionic strength, mean ionic molality and gamma_pm are the machinery of Sec. 15.2
(ionic strength effects), and the ionization equilibria of Secs. 15.1 and 15.6 need
activity coefficients from the same source. This module is written to be imported
there, not copied.

Eric M. Furst
August 2026
"""
import numpy as np
from scipy import interpolate, optimize

__all__ = ["Electrolyte", "DebyeHuckel", "ionic_strength", "water_parameters",
           "TABLE_9_10_1", "ELECTROLYTES"]


# ---------------------------------------------------------------------------
# Table 9.10-1 -- the solvent parameters for water
# ---------------------------------------------------------------------------
# SIS Table 9.10-1, "Values of the Parameters in the Equations for gamma_pm for
# Aqueous Solutions". T in degrees C, alpha in (mol/kg)^(-1/2), beta in
# [(mol/kg)^(1/2) angstrom]^(-1).
#
# ⚠️ THE TABLE CARRIES NO SOURCE. Its footnote (5e note 21) defines the angstrom
# and the liter and says nothing about where the numbers came from. They are the
# standard Debye-Huckel constants for water -- alpha is 2.303 A_gamma, with A_gamma
# the base-10 constant of the physical-chemistry literature -- and they agree with
# the accepted values to the digits printed (see `check_alpha_against_theory` in
# code/ch9/validation/). ⬜ The 6e should give the table a credit line; Robinson and
# Stokes, *Electrolyte Solutions*, 2nd ed. (1959) is already cited in Fig. 9.10-1's
# caption and tabulates them.
TABLE_9_10_1 = np.array([
    # T (C),  alpha,   beta
    (0.0,   1.129,   0.3245),
    (5.0,   1.137,   0.3253),
    (10.0,  1.146,   0.3261),
    (15.0,  1.155,   0.3269),
    (20.0,  1.164,   0.3276),
    (25.0,  1.175,   0.3284),
    (30.0,  1.184,   0.3292),
    (40.0,  1.206,   0.3309),
    (50.0,  1.230,   0.3326),
    (60.0,  1.255,   0.3343),
    (70.0,  1.283,   0.3361),
    (80.0,  1.313,   0.3380),
    (90.0,  1.345,   0.3400),
    (100.0, 1.379,   0.3420),
])

_T_C = TABLE_9_10_1[:, 0]
# Cubic splines through the table. alpha rises by 22% over the range with visible
# curvature (it goes as (eps T)^(-3/2)), so linear interpolation between the 10-degree
# gaps above 30 C would be a poor reading of the book's own table; beta is nearly
# linear and the spline costs nothing there.
_ALPHA = interpolate.CubicSpline(_T_C, TABLE_9_10_1[:, 1])
_BETA = interpolate.CubicSpline(_T_C, TABLE_9_10_1[:, 2])


def water_parameters(T):
    """(alpha, beta) for water at temperature `T` in **kelvin**, from Table 9.10-1.

    T is in K because everything else in this package is; the table is printed in
    degrees C and converted here rather than in the caller.

    Raises outside 273.15-373.15 K. Extrapolating alpha past the table is not a
    small sin -- it is the coefficient the whole limiting law hangs on, and above
    100 C water is not liquid at 1 bar anyway.
    """
    T_C = np.asarray(T, dtype=float) - 273.15
    if np.any(T_C < _T_C[0] - 1e-9) or np.any(T_C > _T_C[-1] + 1e-9):
        raise ValueError(
            f"Table 9.10-1 covers 0-100 C (273.15-373.15 K); got {np.min(T_C):.1f}"
            f"-{np.max(T_C):.1f} C. Pass alpha= and beta= explicitly for another "
            f"solvent or temperature.")
    return float(_ALPHA(T_C)), float(_BETA(T_C))


# ---------------------------------------------------------------------------
# the salt
# ---------------------------------------------------------------------------
class Electrolyte:
    """An electrically neutral electrolyte that dissociates as SIS Eq. 9.10-1.

        A_(nu+) B_(nu-)  ->  nu+ A^(z+)  +  nu- B^(z-)

    >>> Electrolyte("CaCl2", nu_plus=1, z_plus=2, nu_minus=2, z_minus=-1).nu
    3

    Electroneutrality (Eq. 9.10-2, nu+ z+ + nu- z- = 0) is checked in the
    constructor, not left to the caller: a salt that fails it is not a salt, and
    every quantity below would silently return a number anyway.
    """

    def __init__(self, name, nu_plus, z_plus, nu_minus, z_minus):
        self.name = str(name)
        self.nu_plus = int(nu_plus)
        self.nu_minus = int(nu_minus)
        self.z_plus = int(z_plus)
        self.z_minus = int(z_minus)
        if self.nu_plus < 1 or self.nu_minus < 1:
            raise ValueError("nu+ and nu- are numbers of ions and must be >= 1")
        if self.z_plus <= 0 or self.z_minus >= 0:
            raise ValueError(
                f"{name}: z+ must be positive and z- negative (SIS writes the anion "
                f"charge signed, e.g. z- = -1 for chloride); got z+ = {z_plus}, "
                f"z- = {z_minus}")
        net = self.nu_plus * self.z_plus + self.nu_minus * self.z_minus
        if net != 0:
            raise ValueError(
                f"{name} is not electrically neutral: nu+ z+ + nu- z- = {net}, "
                f"which violates Eq. 9.10-2")

    @property
    def nu(self):
        """nu = nu+ + nu-, the total ions per formula unit (used in Eq. 9.10-14)."""
        return self.nu_plus + self.nu_minus

    @property
    def z_product(self):
        """|z+ z-|, the bracketed term of Eqs. 9.10-15, 9.10-17 and 9.10-18."""
        return abs(self.z_plus * self.z_minus)

    def ionic_strength(self, M):
        """Ionic strength I of a solution of this salt alone, at molality `M`.

        Eq. 9.10-16, I = (1/2) SUM_i z_i^2 M_i, with the sum over the ions this
        salt provides: M_+ = nu+ M and M_- = nu- M.

        The result is a multiple of M fixed by the salt's stoichiometry -- I = M
        for a 1:1 salt like NaCl, 3M for CaCl2, 4M for CuSO4, the three cases
        Fig. 9.10-1's caption states.
        """
        M = np.asarray(M, dtype=float)
        return 0.5 * (self.nu_plus * self.z_plus**2
                      + self.nu_minus * self.z_minus**2) * M

    def mean_ionic_molality(self, M):
        """Mean ionic molality M_pm, SIS Eqs. 9.10-12 and 9.10-13.

            M_pm^nu = M_+^(nu+) M_-^(nu-)     with M_+ = nu+ M, M_- = nu- M

        so M_pm = (nu+^nu+ nu-^nu-)^(1/nu) M -- the geometric mean of the two ion
        molalities, weighted by how many of each the salt provides.
        """
        M = np.asarray(M, dtype=float)
        return (self.nu_plus**self.nu_plus
                * self.nu_minus**self.nu_minus) ** (1.0 / self.nu) * M

    def __repr__(self):
        return (f"Electrolyte({self.name!r}, nu_plus={self.nu_plus}, "
                f"z_plus={self.z_plus}, nu_minus={self.nu_minus}, "
                f"z_minus={self.z_minus})")


#: The salts Sec. 9.10 uses by name -- Illustrations 9.10-1 and 9.10-2 and
#: Fig. 9.10-1. Build any other with `Electrolyte(...)` directly; there is
#: deliberately no formula parser here, because guessing the dissociation of an
#: arbitrary formula string is exactly the kind of silent wrong answer this
#: package tries not to produce.
ELECTROLYTES = {
    "HCl":   Electrolyte("HCl", nu_plus=1, z_plus=1, nu_minus=1, z_minus=-1),
    "NaCl":  Electrolyte("NaCl", nu_plus=1, z_plus=1, nu_minus=1, z_minus=-1),
    "CaCl2": Electrolyte("CaCl2", nu_plus=1, z_plus=2, nu_minus=2, z_minus=-1),
    "CuSO4": Electrolyte("CuSO4", nu_plus=1, z_plus=2, nu_minus=1, z_minus=-2),
}


def ionic_strength(molalities, charges):
    """Ionic strength of a **mixed** electrolyte solution, Eq. 9.10-16.

        I = (1/2) SUM_i z_i^2 M_i,     summed over every ion present

    >>> ionic_strength([0.1, 0.1], [1, -1])        # 0.1 molal NaCl
    0.1

    The section's closing paragraph is the reason this exists as a free function:
    "the equations in this section are valid for ... mixed electrolyte solutions",
    where I is computed once by summing over *all* ions and is then the same for
    every salt in the solution. `Electrolyte.ionic_strength` is the single-salt
    shortcut; pass the ion list here when there is more than one salt.
    """
    M = np.asarray(molalities, dtype=float)
    z = np.asarray(charges, dtype=float)
    if M.shape != z.shape:
        raise ValueError(f"got {M.size} molalities and {z.size} charges")
    if np.any(M < 0):
        raise ValueError("molalities must be non-negative")
    net = float(M @ z)
    if not np.isclose(net, 0.0, atol=1e-8 * max(1.0, float(M.sum()))):
        raise ValueError(
            f"the ion list is not electrically neutral: SUM_i z_i M_i = {net:.6g}, "
            f"which Eq. 9.10-3 forbids")
    return 0.5 * float(z**2 @ M)


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------
class DebyeHuckel:
    """Mean ionic activity coefficient from Eq. 9.10-15, 9.10-17 or 9.10-18.

    Parameters
    ----------
    salt : Electrolyte or str
        The dissolved salt; a string is looked up in `ELECTROLYTES`.
    T : float
        Temperature in K, default 298.15. Sets alpha and beta from Table 9.10-1
        unless they are given explicitly.
    beta_a : float, optional
        The product beta*a of Eq. 9.10-17, in (mol/kg)^(-1/2). `None` (default)
        selects the limiting law. The book notes that a is "a constant related to
        the average hydrated radius of ions, usually about 4 angstrom", but that
        "in practice the product beta a is sometimes set equal to unity or treated
        as an adjustable parameter" -- and both illustrations set it to 1.
    delta : float, optional
        The linear term of Eq. 9.10-18, in (mol/kg)^(-1). `None` (default) omits
        it. Requires `beta_a`: the book prints no form with a delta term and no
        denominator, and the two together are what make 9.10-18 work.
    alpha, beta : float, optional
        Override Table 9.10-1 -- for a solvent other than water, or to reproduce a
        printed number that used a different constant (Illustration 9.10-2 does).

    Notes
    -----
    `a` is not stored separately from `beta_a` because the book never uses one
    without the other, and carrying an angstrom radius around invites a unit error
    at the one multiplication where it matters.
    """

    def __init__(self, salt, T=298.15, beta_a=None, delta=None,
                 alpha=None, beta=None):
        if isinstance(salt, str):
            try:
                salt = ELECTROLYTES[salt]
            except KeyError:
                raise KeyError(
                    f"{salt!r} is not one of the salts Sec. 9.10 names "
                    f"({sorted(ELECTROLYTES)}); build it with Electrolyte(...)"
                ) from None
        self.salt = salt
        self.T = float(T)
        if delta is not None and beta_a is None:
            raise ValueError(
                "delta is the correction term of Eq. 9.10-18, which also has the "
                "1 + beta a sqrt(I) denominator; pass beta_a as well (both "
                "illustrations use beta_a = 1)")
        self.beta_a = None if beta_a is None else float(beta_a)
        self.delta = None if delta is None else float(delta)
        if alpha is None or beta is None:
            a_tab, b_tab = water_parameters(self.T)
            alpha = a_tab if alpha is None else alpha
            beta = b_tab if beta is None else beta
        self.alpha = float(alpha)
        self.beta = float(beta)

    # -- what equation is this -------------------------------------------
    @property
    def equation(self):
        """The book's number for the equation this instance evaluates."""
        if self.beta_a is None:
            return "9.10-15"
        return "9.10-17" if self.delta is None else "9.10-18"

    @property
    def slope(self):
        """The limiting-law slope alpha|z+ z-|, i.e. -d ln gamma_pm / d sqrt(I)
        as I -> 0. Every one of the three equations has it, which is the sense in
        which 9.10-17 and 9.10-18 are corrections rather than rival models."""
        return self.alpha * self.salt.z_product

    # -- the two public answers -------------------------------------------
    def ln_gamma_pm(self, M):
        """ln gamma_pm at molality `M` of the salt (mol/kg)."""
        I = self.salt.ionic_strength(M)
        return self.ln_gamma_pm_from_I(I)

    def gamma_pm(self, M):
        """Mean ionic activity coefficient at molality `M` (mol/kg)."""
        return np.exp(self.ln_gamma_pm(M))

    def ln_gamma_pm_from_I(self, I):
        """ln gamma_pm as a function of ionic strength directly.

        The form the equations are actually written in, and the one a mixed
        electrolyte solution needs -- there I comes from summing over every ion
        present, not from this salt's own molality.
        """
        I = np.asarray(I, dtype=float)
        if np.any(I < 0):
            raise ValueError("ionic strength cannot be negative")
        rootI = np.sqrt(I)
        ln_g = -self.slope * rootI                                # Eq. 9.10-15
        if self.beta_a is not None:
            ln_g = ln_g / (1.0 + self.beta_a * rootI)             # Eq. 9.10-17
        if self.delta is not None:
            ln_g = ln_g + self.delta * I                          # Eq. 9.10-18
        return ln_g

    # -- fitting -----------------------------------------------------------
    def fit_delta(self, M, gamma_pm, weight_log=True):
        """Least-squares delta for Eq. 9.10-18 against measured gamma_pm.

        Returns a **new** model carrying the fitted delta, leaving this one alone,
        so a notebook can plot the fitted and unfitted curves side by side.

        `weight_log=True` (the default) fits ln gamma_pm rather than gamma_pm.
        That is the right choice here and not a detail: over Illustration 9.10-1's
        range gamma_pm runs from 0.76 to 42, so a fit in gamma_pm is a fit to the
        four most concentrated points and ignores the dilute end where the theory
        is supposed to be exact.

        Only delta is fitted. alpha comes from Table 9.10-1 and beta_a is the
        book's chosen 1; making all three adjustable would fit anything and
        demonstrate nothing.
        """
        if self.beta_a is None:
            raise ValueError("fit_delta needs beta_a; delta belongs to Eq. 9.10-18")
        M = np.asarray(M, dtype=float)
        g = np.asarray(gamma_pm, dtype=float)
        if M.shape != g.shape:
            raise ValueError(f"got {M.size} molalities and {g.size} gamma values")
        if np.any(g <= 0):
            raise ValueError("activity coefficients must be positive")
        I = self.salt.ionic_strength(M)
        base = DebyeHuckel(self.salt, self.T, beta_a=self.beta_a,
                           alpha=self.alpha, beta=self.beta)

        def residual(p):
            model = base.ln_gamma_pm_from_I(I) + p[0] * I
            return model - np.log(g) if weight_log else np.exp(model) - g

        delta = float(optimize.least_squares(residual, [0.1]).x[0])
        return DebyeHuckel(self.salt, self.T, beta_a=self.beta_a, delta=delta,
                           alpha=self.alpha, beta=self.beta)

    def rms(self, M, gamma_pm, in_log=True):
        """RMS deviation from measured gamma_pm -- in ln gamma_pm by default.

        Reported in the log because that is the quantity every one of the three
        equations is linear-ish in, and because it is the only way the dilute
        points carry any weight (see `fit_delta`).
        """
        M = np.asarray(M, dtype=float)
        g = np.asarray(gamma_pm, dtype=float)
        if in_log:
            return float(np.sqrt(np.mean((self.ln_gamma_pm(M) - np.log(g))**2)))
        return float(np.sqrt(np.mean((self.gamma_pm(M) - g)**2)))

    def __repr__(self):
        bits = [f"{self.salt.name!r}", f"T={self.T:g}"]
        if self.beta_a is not None:
            bits.append(f"beta_a={self.beta_a:g}")
        if self.delta is not None:
            bits.append(f"delta={self.delta:g}")
        return f"DebyeHuckel({', '.join(bits)})  # Eq. {self.equation}"


# ---------------------------------------------------------------------------
# the book's own data
# ---------------------------------------------------------------------------
# Illustration 9.10-1: mean ionic activity coefficients of HCl in water at 25 C,
# as printed. The abscissa of Fig. 9.10-2 is sqrt(M_HCl), and for a 1:1 salt
# I = M, so these molalities are also the ionic strengths.
ILLUSTRATION_9_10_1_HCL = np.array([
    # M (mol/kg), gamma_pm
    (0.0005, 0.975), (0.001, 0.965), (0.005, 0.928), (0.01, 0.904),
    (0.05, 0.830), (0.1, 0.796), (0.5, 0.757), (1.0, 0.809),
    (3.0, 1.316), (5.0, 2.38), (8.0, 5.90), (10.0, 10.44),
    (12.0, 17.25), (14.0, 27.3), (16.0, 42.4),
])

# Illustration 9.10-2: mean ionic activity coefficients of NaCl in water at 25 C.
# ⚠️ The printed table's first row is M = 0, gamma_pm = 1.00 -- the infinite-dilution
# limit, not a measurement. It is dropped here: every model returns exactly 1 there
# by construction, so including it in a fit or an RMS flatters all three equally.
ILLUSTRATION_9_10_2_NACL = np.array([
    (0.1, 0.778), (0.25, 0.720), (0.5, 0.681), (0.75, 0.665), (1.0, 0.657),
    (2.0, 0.669), (3.0, 0.714), (4.0, 0.782), (5.0, 0.873), (6.0, 0.987),
])
