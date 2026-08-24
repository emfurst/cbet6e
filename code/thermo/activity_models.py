"""Liquid-phase activity coefficient models -- SIS Sections 9.5 and 9.6.

The correlative models of Sec. 9.5 (Margules, Redlich-Kister, van Laar, Wilson,
NRTL, Flory-Huggins, UNIQUAC) and the predictive regular solution model of Sec. 9.6,
behind **one interface**:

    model.gamma(x, T)   ->  activity coefficients, dimensionless
    model.gex(x, T)     ->  molar excess Gibbs energy, J/mol

That uniformity is the point of the module, not an accident of style. Section 9.11 --
the payoff of the chapter -- is a pair of tables recommending one model over another
for a given kind of mixture. With a common interface those recommendations can be
*demonstrated*: the same data, the same call, seven curves. Without it they can only
be asserted. It is also what lets the Chapter 10 bubble- and dew-point drivers accept
any of these models without knowing which one they were handed.

The Aspen-optional substitute for activity-coefficient calculations, and the Python
replacement for the 4th-edition ACTCOEFF worksheet.

**Every model here is temperature-aware at the interface even when its parameters are
not.** `gamma(x, T)` requires T throughout. For van Laar, Wilson, NRTL and UNIQUAC the
book's parameters already absorb T, so T only sets the scale of `gex`; for Margules,
Redlich-Kister and regular solution T enters the activity coefficient itself. Requiring
it everywhere costs one argument and removes the class of error where a Margules A in
J/mol is silently used as A/RT.

## Reusing a model at a second temperature -- read this before extrapolating

"Absorb T" means **at the temperature the parameters were fitted at.** Three of these
models are defined by an energy divided by RT, and the book prints the law:

    Wilson    Lambda_ij = (V_j/V_i) exp[-(lambda_ij - lambda_ii)/RT]   Eq. 9.5-11
    NRTL      tau_ij    = (g_ij - g_jj)/RT                             Eq. 9.5-13
    UNIQUAC   ln tau_ij = -(u_ij - u_jj)/RT                            para. 1303

So the **temperature-independent parameter is the energy**, and Lambda or tau carries a
1/T. A model built with `NRTL(tau12=..., tau21=...)` stores the ratio, not the energy;
it has no way to move it, and its `lngamma(x, T)` therefore uses the same numbers at
every T. That is correct for one isotherm and wrong for an extrapolation, and the
failure is silent: the fit at the reference temperature is identical either way, because
the map between the energy and the ratio is one-to-one there.

**Reuse across temperature is therefore an explicit choice, and this module makes you
make it.** Evaluating one instance at a second temperature emits a `UserWarning`, once
per instance. Three ways to say what you meant:

    NRTL(tau12=dg12 / (R * T), tau21=dg21 / (R * T), alpha=0.3)   # rebuild at each T
    NRTL(t12, t21, alpha=0.3).at_T(298.15)                        # declare the origin
    m.rescaled_to(523.15)                                         # move it, by the law

`at_T(T)` records the temperature the parameters were fitted at and returns the model,
so the warning can name it. `rescaled_to(T)` returns a **new** model with the parameters
carried to T by the law above -- NRTL and UNIQUAC only, because Wilson's Lambda bundles
the volume ratio V_j/V_i in with the exponential and the energy cannot be recovered from
Lambda alone (build Wilson from lambda_ij and the volumes instead).

Rebuilding by hand is the clearest of the three in a notebook, because it puts the 1/T
on the page next to the equation. The other two are for code that is handed a model it
did not build.

**Not guarded, on purpose:** van Laar (Table 9.5-1 tabulates alpha and beta per
temperature *range* -- the book gives no law to apply, so there is nothing to correct)
and Flory-Huggins (chi's temperature dependence is substance-specific and the book fits
it per system, so the notebooks already build a fresh model at each T). Margules,
Redlich-Kister and regular solution hold energies in J/mol and divide by RT themselves.

This guard exists because the class of error cost a session on Fig. 11.2-6, where the
Wong-Sandler curve missed the measured liquid-liquid dome badly with tau frozen and
matched it to 0.004 in mole fraction with the energies frozen.

## Units

SI, as everywhere else in the package: T in K, `gex` in J/mol, and the Margules /
Redlich-Kister constants in J/mol. Two exceptions are the book's, not ours:

- **van Laar alpha, beta, Wilson Lambda_ij, NRTL tau_ij and UNIQUAC tau_ij are
  dimensionless** -- they are defined as contributions to ln gamma, so there is no
  choice to make.
- **Regular solution theory is tabulated in cal and cc.** Table 9.6-1 gives molar
  volumes in cc/mol and solubility parameters in (cal/cc)^(1/2), and says in its own
  note that it keeps those units by tradition. `RegularSolution` stores SI and
  `RegularSolution.from_table_9_6_1` does the conversion, so the traditional numbers
  can be typed in as printed and the arithmetic still happens in J/mol.

## The consistency check is not decoration

For most of these models the book prints *both* an expression for G^ex and expressions
for the activity coefficients, derived from it. Those two must satisfy

    G^ex / RT = SUM_i x_i ln gamma_i                                    (Eq. 9.3-12)

exactly, at every composition -- and it is an *independent* check on the algebra,
because the two printed expressions look nothing alike. `check_gibbs_duhem` runs it.
It is how the UNIQUAC residual term in this module was verified against SIS
Eqs. 9.5-21 and 9.5-23b, and it is worth running in a notebook rather than trusting:
the same identity, applied to the 5e's own hand-drawn tangent construction in Table
8.6-4, misses by up to 16 J/mol (see `fitting.py`).

Eric M. Furst
August 2026
"""
import functools
import warnings

import numpy as np
from scipy import constants, optimize

R = constants.R

# Table 9.6-1 is printed in cc/mol and (cal/cc)^(1/2). One (cal/cc)^(1/2) is
# sqrt(4.184 J / 1e-6 m^3) = 2045.48 Pa^(1/2); one cc/mol is 1e-6 m^3/mol. Then
# V (delta - deltabar)^2 comes out in m^3/mol * Pa = J/mol.
CC_PER_MOL = 1e-6                                  # cc/mol   -> m^3/mol
CAL_CC_HALF = np.sqrt(constants.calorie / 1e-6)     # (cal/cc)^(1/2) -> Pa^(1/2)

__all__ = ["ActivityModel", "FrozenParameterWarning",
           "OneConstantMargules", "TwoConstantMargules",
           "RedlichKisterGex", "VanLaar", "Wilson", "NRTL", "FloryHuggins",
           "UNIQUAC", "RegularSolution", "fit_binary",
           "TABLE_9_5_1_VAN_LAAR", "TABLE_9_6_1", "TABLE_9_6_1_BLOCK",
           "CC_PER_MOL", "CAL_CC_HALF"]


# ---------------------------------------------------------------------------
# the temperature guard -- see "Reusing a model at a second temperature" above
# ---------------------------------------------------------------------------
class FrozenParameterWarning(UserWarning):
    """One model instance was evaluated at two temperatures.

    Its stored parameter is a ratio that hides a 1/T (Eq. 9.5-11, 9.5-13, para.
    1303), so the second evaluation is a different model from the one that was
    fitted -- silently, because both give the same answer at the first temperature.
    """


def _guard_T(func):
    """Wrap `lngamma` / `gex_over_RT` so the temperature they are called at is seen.

    Both are wrapped, because they are independent entry points: `gamma` routes
    through `lngamma`, but `gex` routes through `gex_over_RT`, and the Sec. 9.5
    models override the latter with the book's own printed expression for G^ex.
    `WongSandler` calls `gex(x, T)`, so guarding only `lngamma` would miss the case
    this guard was written for.
    """
    @functools.wraps(func)
    def wrapper(self, x, T, *args, **kwargs):
        self._note_T(T)
        return func(self, x, T, *args, **kwargs)

    wrapper._guards_T = True
    return wrapper


# ---------------------------------------------------------------------------
# the interface
# ---------------------------------------------------------------------------
class ActivityModel:
    """Base class: a model supplies `lngamma(x, T)` and inherits the rest.

    Subclasses may also override `gex_over_RT` when the book prints an independent
    expression for G^ex -- that is what makes `check_gibbs_duhem` a real test rather
    than a tautology.
    """

    #: number of species the model is written for; None means any
    n = None

    #: The book's law for the stored parameter, for the models whose parameter is an
    #: energy divided by RT. Not None means one instance describes **one**
    #: temperature; see the module docstring. None means the parameters are energies
    #: in J/mol and `lngamma` divides by RT itself, so the model travels freely.
    _T_LAW = None

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        for name in ("lngamma", "gex_over_RT"):
            f = cls.__dict__.get(name)
            if f is not None and not getattr(f, "_guards_T", False):
                setattr(cls, name, _guard_T(f))

    # -- temperature bookkeeping -------------------------------------------
    def at_T(self, T):
        """Declare the temperature these parameters were fitted at; returns `self`.

        Only bookkeeping -- it changes no number. It lets the warning name the
        reference temperature, and it is the honest way to record a fit whose
        parameters arrived as bare ratios from a table or a paper.
        """
        self._T_ref = float(T)
        return self

    def rescaled_to(self, T):
        """A **new** model with the parameters carried to `T` by the book's law.

        Defined where the law can be inverted from what is stored. `self` is left
        alone, because the fitted model is a record of a measurement and should not
        change under the reader's feet.
        """
        if self._T_LAW is None:
            return self                     # parameters are energies; nothing to move
        raise NotImplementedError(
            f"{type(self).__name__} has no rescaling: its stored parameter cannot be "
            f"inverted for the energy on its own ({self._T_LAW}). Build the model "
            f"from the energies at the temperature you want.")

    def _note_T(self, T):
        """Warn the first time this instance is used at a second temperature."""
        if self._T_LAW is None or T is None:
            return
        T = float(T)
        seen = self.__dict__.get("_T_ref")
        if seen is None:
            self.__dict__["_T_ref"] = T
        elif abs(T - seen) > 1e-6 and not self.__dict__.get("_T_warned"):
            self.__dict__["_T_warned"] = True
            warnings.warn(
                f"{type(self).__name__} built for T = {seen:g} K is being evaluated "
                f"at T = {T:g} K, and its parameters were not moved. The book defines "
                f"{self._T_LAW}, so the temperature-independent quantity is the "
                f"energy, not what is stored here -- the two agree at {seen:g} K and "
                f"nowhere else. Rebuild the model at each temperature, or call "
                f"`.rescaled_to(T)`. If holding it fixed is deliberate, say so with "
                f"`warnings.simplefilter('ignore', FrozenParameterWarning)`.",
                FrozenParameterWarning, stacklevel=3)

    # -- the two public answers -------------------------------------------
    def gamma(self, x, T):
        """Activity coefficients at composition `x` and temperature `T` (K)."""
        return np.exp(self.lngamma(x, T))

    def gex(self, x, T):
        """Molar excess Gibbs energy, J/mol."""
        return R * T * self.gex_over_RT(x, T)

    # -- what subclasses provide ------------------------------------------
    def lngamma(self, x, T):
        raise NotImplementedError

    def gex_over_RT(self, x, T):
        """G^ex/RT. Default is the summation, Eq. 9.3-12; override with the
        book's own printed expression where there is one."""
        x = self._x(x)
        return float(x @ self.lngamma(x, T))

    # -- checks ------------------------------------------------------------
    def check_gibbs_duhem(self, T, x1=None, rtol=1e-10):
        """Largest |SUM_i x_i ln gamma_i - G^ex/RT| over a composition scan.

        Zero to machine precision means the activity coefficient expressions and
        the G^ex expression are the same thermodynamics. A model that does not
        override `gex_over_RT` returns 0 trivially, and says so.
        """
        if type(self).gex_over_RT is ActivityModel.gex_over_RT:
            return 0.0  # nothing independent to compare against
        if x1 is None:
            x1 = np.linspace(0.001, 0.999, 199)
        worst = 0.0
        for a in np.atleast_1d(x1):
            x = np.array([a, 1.0 - a]) if (self.n in (2, None)) else None
            if x is None:
                raise ValueError("pass x1 explicitly for multicomponent models")
            lhs = float(x @ self.lngamma(x, T))
            rhs = self.gex_over_RT(x, T)
            worst = max(worst, abs(lhs - rhs))
        return worst

    # -- helpers -----------------------------------------------------------
    def _x(self, x):
        x = np.asarray(x, dtype=float)
        if self.n is not None and x.size != self.n:
            raise ValueError(f"{type(self).__name__} is written for {self.n} "
                             f"species, got {x.size}")
        if np.any(x < 0):
            raise ValueError("mole fractions must be non-negative")
        s = x.sum()
        if not np.isclose(s, 1.0, atol=1e-8):
            raise ValueError(f"mole fractions must sum to 1, got {s:.6g}")
        return x

    def __repr__(self):
        ps = ", ".join(f"{k}={v!r}" for k, v in self._params().items())
        return f"{type(self).__name__}({ps})"

    def _params(self):
        return {}


# ---------------------------------------------------------------------------
# Sec. 9.5 -- the polynomial G^ex models
# ---------------------------------------------------------------------------
class OneConstantMargules(ActivityModel):
    """One-constant Margules, SIS Eqs. 9.5-1 and 9.5-5.

        G^ex = A x1 x2          RT ln gamma_1 = A x2^2    RT ln gamma_2 = A x1^2

    Parameters
    ----------
    A : J/mol. Positive gives gamma > 1, negative gamma < 1 (SIS p. 449).

    The two activity coefficients are mirror images -- a consequence of choosing a
    symmetric G^ex, not a general result (SIS Fig. 9.5-3).
    """

    n = 2

    def __init__(self, A):
        self.A = float(A)

    def _params(self):
        return {"A": self.A}

    def lngamma(self, x, T):
        x1, x2 = self._x(x)
        return np.array([self.A * x2**2, self.A * x1**2]) / (R * T)

    def gex_over_RT(self, x, T):
        x1, x2 = self._x(x)
        return self.A * x1 * x2 / (R * T)

    @classmethod
    def fit(cls, x1, gamma1, gamma2, T):
        """Least squares on both species' RT ln gamma at once.

        One parameter against two curves, so the fit cannot follow both unless the
        data really are symmetric -- which is the whole content of Illustration
        9.5-1's answer for this model.
        """
        x1 = np.asarray(x1, float)
        x2 = 1.0 - x1
        design = np.concatenate([x2**2, x1**2])
        obs = R * T * np.concatenate([np.log(gamma1), np.log(gamma2)])
        A = float(design @ obs / (design @ design))
        return cls(A)


class TwoConstantMargules(ActivityModel):
    """Two-constant Margules, SIS Eqs. 9.5-6 (truncated) and 9.5-7.

        G^ex = x1 x2 [A + B(x1 - x2)]
        RT ln gamma_1 = alpha_1 x2^2 + beta_1 x2^3
        RT ln gamma_2 = alpha_2 x1^2 + beta_2 x1^3

    with alpha_i = A + 3(-1)^(i+1) B and beta_i = 4(-1)^i B, so

        alpha_1 = A + 3B, beta_1 = -4B      alpha_2 = A - 3B, beta_2 = +4B

    A and B in J/mol. Unlike the one-constant form, G^ex is not symmetric and the two
    activity coefficients are not mirror images.
    """

    n = 2

    def __init__(self, A, B):
        self.A, self.B = float(A), float(B)

    def _params(self):
        return {"A": self.A, "B": self.B}

    def lngamma(self, x, T):
        x1, x2 = self._x(x)
        A, B = self.A, self.B
        return np.array([(A + 3 * B) * x2**2 - 4 * B * x2**3,
                         (A - 3 * B) * x1**2 + 4 * B * x1**3]) / (R * T)

    def gex_over_RT(self, x, T):
        x1, x2 = self._x(x)
        return x1 * x2 * (self.A + self.B * (x1 - x2)) / (R * T)

    @classmethod
    def fit(cls, x1, gamma1, gamma2, T):
        """Least squares on both species' RT ln gamma at once."""
        x1 = np.asarray(x1, float)
        x2 = 1.0 - x1
        # columns are the coefficients of A and of B in each species' RT ln gamma
        col_A = np.concatenate([x2**2, x1**2])
        col_B = np.concatenate([3 * x2**2 - 4 * x2**3, -3 * x1**2 + 4 * x1**3])
        obs = R * T * np.concatenate([np.log(gamma1), np.log(gamma2)])
        (A, B), *_ = np.linalg.lstsq(np.column_stack([col_A, col_B]), obs, rcond=None)
        return cls(A, B)


class RedlichKisterGex(ActivityModel):
    """The Redlich-Kister expansion of G^ex, SIS Eq. 9.5-6, to any order.

        G^ex = x1 x2 SUM_i a_i (x1 - x2)^i

    Coefficients in J/mol. This is the *same expansion* that `fitting.RedlichKister`
    applies to a property change on mixing in Sec. 8.6; it is a separate class
    because the quantity of interest here is the activity coefficient, and the
    derivative is taken with respect to mole number rather than reported as a partial
    molar property.

    Differentiating (SIS Illustration 10.2-4, Eqs. 2 and 3) gives, with u = x1 - x2,

        RT ln gamma_1 = x2^2 S(u) + 2 x1 x2^2 S'(u)
        RT ln gamma_2 = x1^2 S(u) - 2 x1^2 x2 S'(u)

    where S(u) = SUM_i a_i u^i. Illustration 10.2-4's three-term fit to the
    benzene / 2,2,4-trimethyl pentane data at 55 C is the book's own worked case:
    a = (1389.0, 419.45, 109.83) J/mol. Those coefficients reproduce all fourteen
    smoothed activity coefficients and all seven smoothed G^ex values printed in
    Illustration 10.2-4, to the last digit:

        >>> rk = RedlichKisterGex([1389.0, 419.45, 109.83])
        >>> g1, g2 = rk.gamma([0.5256, 0.4744], 328.15)
        >>> round(float(g1), 3), round(float(g2), 3)
        (1.166, 1.107)

    **They are not recoverable from the seven data points the book prints.**
    `fit_gex` on the printed G^ex column gives (1392.5, 430.4, 96.6) -- a 12% shift in
    the third coefficient. That is expected and must not be tuned away.

    **They are also NOT Weissman and Wood's own constants**, though this docstring
    said so until the paper was read (J. Chem. Phys. 32, 1153 (1960), 2026-08-15).
    Their Eqs. (8)-(10) give the coefficients as functions of temperature, and at
    328.15 K those come to (1354.7, 424.5, 101.0) J/mol -- not these. Nor does the
    difference come from the corrections the book drops: applying their Eq. (2)
    virial correction to the seven printed points moves a *further* away, to 1347.
    The printed constants are within 0.25% in a of a plain least-squares fit to the
    seven points as the book reduces them, so they are the book's own fit to its own
    reduction, and the c that is 12% adrift is simply the coefficient the seven
    points constrain least. Use the printed coefficients to reproduce the book; use
    `fit_gex` to fit your own data, and say which you did.
    """

    n = 2

    def __init__(self, a):
        self.a = np.asarray(a, dtype=float).ravel()
        if self.a.size == 0:
            raise ValueError("need at least one coefficient")

    def _params(self):
        return {"a": list(self.a)}

    def _S(self, u):
        i = np.arange(self.a.size)
        return float(np.sum(self.a * u**i))

    def _dS(self, u):
        i = np.arange(1, self.a.size)
        if i.size == 0:
            return 0.0
        return float(np.sum(self.a[1:] * i * u ** (i - 1)))

    def lngamma(self, x, T):
        x1, x2 = self._x(x)
        u = x1 - x2
        S, dS = self._S(u), self._dS(u)
        return np.array([x2**2 * S + 2 * x1 * x2**2 * dS,
                         x1**2 * S - 2 * x1**2 * x2 * dS]) / (R * T)

    def gex_over_RT(self, x, T):
        x1, x2 = self._x(x)
        return x1 * x2 * self._S(x1 - x2) / (R * T)

    @classmethod
    def fit_gex(cls, x1, gex, order=2):
        """Least squares of G^ex (J/mol) directly -- the route of Illustration 10.2-4.

        Linear in the coefficients, so one `lstsq` and no starting guess. `order` is
        the highest power of (x1 - x2), so order=2 gives three coefficients a, b, c.
        """
        x1 = np.asarray(x1, float)
        gex = np.asarray(gex, float)
        x2 = 1.0 - x1
        u = x1 - x2
        design = np.column_stack([x1 * x2 * u**i for i in range(order + 1)])
        a, *_ = np.linalg.lstsq(design, gex, rcond=None)
        return cls(a)


class VanLaar(ActivityModel):
    """van Laar, SIS Eq. 9.5-9.

        ln gamma_1 = alpha / [1 + (alpha/beta)(x1/x2)]^2
        ln gamma_2 = beta  / [1 + (beta/alpha)(x2/x1)]^2

    alpha and beta are dimensionless (SIS Table 9.5-1 tabulates them for 30 binaries
    with the temperature range over which each pair applies).

    When alpha = beta the model collapses to the one-constant Margules form
    ln gamma_1 = alpha x2^2 -- the note printed under Table 9.5-1, and worth checking
    in a notebook rather than believing:

        >>> vl, m = VanLaar(0.5, 0.5), OneConstantMargules(0.5 * R * 300.0)
        >>> bool(np.allclose(vl.gamma([0.3, 0.7], 300.0), m.gamma([0.3, 0.7], 300.0)))
        True
    """

    n = 2

    def __init__(self, alpha, beta):
        self.alpha, self.beta = float(alpha), float(beta)
        if self.alpha == 0.0 or self.beta == 0.0:
            raise ValueError("van Laar alpha and beta must both be nonzero")

    def _params(self):
        return {"alpha": self.alpha, "beta": self.beta}

    def lngamma(self, x, T):
        x1, x2 = self._x(x)
        a, b = self.alpha, self.beta
        # the limits are finite; take them directly rather than dividing by zero
        ln1 = a if x1 == 0.0 else a / (1.0 + (a / b) * (x1 / x2)) ** 2 if x2 else 0.0
        ln2 = b if x2 == 0.0 else b / (1.0 + (b / a) * (x2 / x1)) ** 2 if x1 else 0.0
        return np.array([ln1, ln2])

    def gex_over_RT(self, x, T):
        """SIS p. 450: G^ex/RT = alpha beta x1 x2 / (alpha x1 + beta x2).

        Printed there as 2 a12 x1 q1 x2 q2 / (x1 q1 + x2 q2) with alpha = 2 q1 a12
        and beta = 2 q2 a12; the form above is that with the q's eliminated.
        """
        x1, x2 = self._x(x)
        a, b = self.alpha, self.beta
        return a * b * x1 * x2 / (a * x1 + b * x2)

    @classmethod
    def from_single_point(cls, x1, gamma1, gamma2):
        """SIS Eq. 9.5-10 -- both constants from activity coefficients at *one*
        composition.

        This is the route Chapter 10 uses, and the first half of Illustration
        9.5-1(b): at x1 = 0.6 the benzene / 2,2,4-trimethyl pentane data give
        alpha = 0.415 and beta = 0.706.

        **It is exact algebra and exquisitely sensitive to the point you read.**
        The round trip is exact -- feeding back the activity coefficients that
        alpha = 0.415, beta = 0.706 imply at x1 = 0.6 (gamma_1 = 1.1243,
        gamma_2 = 1.1677) returns those constants to six figures. But reading the
        smoothed curve 0.2% differently, at gamma_1 = 1.1266 and gamma_2 = 1.1572,
        gives alpha = 0.393 and beta = 0.723 -- a 5% shift in alpha from a change in
        the third decimal place of an input. That is the argument for `fit`, which
        uses all the data, and it is why Illustration 9.5-1 does both.
        """
        l1, l2 = np.log(float(gamma1)), np.log(float(gamma2))
        x2 = 1.0 - float(x1)
        alpha = (1.0 + (x2 * l2) / (float(x1) * l1)) ** 2 * l1
        beta = (1.0 + (float(x1) * l1) / (x2 * l2)) ** 2 * l2
        return cls(alpha, beta)

    @classmethod
    def fit(cls, x1, gamma1, gamma2, p0=(1.0, 1.0)):
        """Least squares on ln gamma over all the data -- the other half of
        Illustration 9.5-1(b). Compare with `from_single_point`; the book does both
        and the difference is the lesson."""
        return fit_binary(cls, x1, gamma1, gamma2, T=None, p0=p0)


class FloryHuggins(ActivityModel):
    """Flory-Huggins, SIS Eqs. 9.5-17 and 9.5-18 -- for molecules of very different
    size, including polymer solutions.

        G^ex/RT    = [x1 ln(phi_1/x1) + x2 ln(phi_2/x2)] + chi (x1 + m x2) phi_1 phi_2
        ln gamma_1 = ln(phi_1/x1) + (1 - 1/m) phi_2 + chi phi_2^2
        ln gamma_2 = ln(phi_2/x2) + (1 - m) phi_1   + m chi phi_1^2

    with phi_1 = x1/(x1 + m x2), phi_2 = m x2/(x1 + m x2), m = v2/v1 the ratio of
    molecular volumes, and chi the Flory interaction parameter (dimensionless). The
    first bracket in G^ex is the combinatorial (entropic) contribution, the second the
    residual (enthalpic) one.

    **ERRATUM, SIS Eq. 9.5-18.** The book prints the second equation with
    **+(m - 1) phi_1**. That sign is inconsistent with the book's own Eq. 9.5-17: the
    two disagree by 2(m - 1) x1 x2 / (x1 + m x2), which for m = 3, chi = 0.5 reaches
    **0.535 in ln gamma**. With (1 - m) phi_1 the identity G^ex/RT = SUM_i x_i ln
    gamma_i holds to machine precision, and that is also the standard Flory-Huggins
    result. This module implements the consistent form; `check_gibbs_duhem` is what
    found the discrepancy, and returns ~1e-16 as written and ~5e-1 with the printed
    sign. Problem 9.29 asks the reader to derive Eq. 9.5-18, so a student who does the
    derivation correctly will disagree with the printed answer.

    **A SECOND ERRATUM, SIS Illustration 11.2-6.** That illustration restates these
    equations "from Eqs. 9.5-18" and drops the **m** in front of chi phi_1^2. Its
    printed table of PS/PMMA coexistence compositions was computed with the m missing
    *and* with the sign the way this class writes it -- neither Eq. 9.5-18 as printed
    nor the consistent form reproduces the table, and the combination
    `printed_sign=False, printed_chi=True` reproduces every row of it to 0.0004 in
    mole fraction. So the sign in Illustration 11.2-6 is a typesetting slip that did
    not reach the arithmetic, and the missing m did. Set both flags to reproduce the
    printed table; leave both alone to get the physics. The two answers differ by up
    to 0.10 in mole fraction near the critical point and they disagree about where
    the two polymers become miscible.
    """

    n = 2

    def __init__(self, chi, m):
        self.chi, self.m = float(chi), float(m)
        if self.m <= 0:
            raise ValueError("m = v2/v1 must be positive")

    def _params(self):
        return {"chi": self.chi, "m": self.m}

    def _phi(self, x1, x2):
        d = x1 + self.m * x2
        return x1 / d, self.m * x2 / d

    def lngamma(self, x, T, printed_sign=False, printed_chi=False):
        """Both flags select a printed variant of Eq. 9.5-18 over the consistent one.

        `printed_sign=True` gives the +(m - 1) phi_1 of Eq. 9.5-18 as typeset;
        `printed_chi=True` drops the m from m chi phi_1^2, as SIS Illustration 11.2-6
        restates it. Either is thermodynamically inconsistent -- see the class
        docstring for what each is for.
        """
        x1, x2 = self._x(x)
        p1, p2 = self._phi(x1, x2)
        m, chi = self.m, self.chi
        s = +1.0 if printed_sign else -1.0
        mc = 1.0 if printed_chi else m
        with np.errstate(divide="ignore", invalid="ignore"):
            ln1 = (np.log(p1 / x1) if x1 else np.log(1.0 / m)) + (1 - 1 / m) * p2 + chi * p2**2
            ln2 = (np.log(p2 / x2) if x2 else np.log(m)) + s * (m - 1) * p1 + mc * chi * p1**2
        return np.array([ln1, ln2])

    def gex_over_RT(self, x, T):
        x1, x2 = self._x(x)
        p1, p2 = self._phi(x1, x2)
        comb = (x1 * np.log(p1 / x1) if x1 else 0.0) + (x2 * np.log(p2 / x2) if x2 else 0.0)
        return comb + self.chi * (x1 + self.m * x2) * p1 * p2


# ---------------------------------------------------------------------------
# Sec. 9.5 -- the local-composition models
# ---------------------------------------------------------------------------
class Wilson(ActivityModel):
    """Wilson, SIS Eqs. 9.5-11 and 9.5-12a; multicomponent form from Appendix A9.2.

    Binary, as the book prints it:

        G^ex/RT = -x1 ln(x1 + x2 L12) - x2 ln(x2 + x1 L21)

    Multicomponent, which reduces to it exactly (checked in `check_gibbs_duhem`):

        G^ex/RT   = -SUM_i x_i ln(SUM_j x_j L_ij)
        ln gamma_i = 1 - ln(SUM_j x_j L_ij) - SUM_k x_k L_ki / (SUM_j x_j L_kj)

    Parameters
    ----------
    Lambda : the 2 constants (L12, L21) for a binary, or a full n x n matrix with
        L_ii = 1. Dimensionless.

    **Wilson cannot predict liquid-liquid immiscibility** for any parameter values
    -- the reason Chapter 11 reaches for NRTL or UNIQUAC instead. The model is not
    broken; the limitation is structural, and Sec. 9.11 says so.

    **One instance is one temperature** -- see the module docstring. Wilson is the
    one guarded model with no `rescaled_to`: Lambda_ij bundles the volume ratio
    V_j/V_i in with the exponential, so the energy cannot be recovered from Lambda
    alone. Build it from lambda_ij and the molar volumes at the temperature you want.
    """

    _T_LAW = "Lambda_ij = (V_j/V_i) exp[-(lambda_ij - lambda_ii)/RT] (Eq. 9.5-11)"

    def __init__(self, Lambda, L21=None):
        if L21 is not None:
            Lambda = [[1.0, float(Lambda)], [float(L21), 1.0]]
        self.Lambda = np.asarray(Lambda, dtype=float)
        if self.Lambda.ndim != 2 or self.Lambda.shape[0] != self.Lambda.shape[1]:
            raise ValueError("Lambda must be square, or pass (L12, L21)")
        # Lambda_ij = (v_j/v_i) exp[-(lambda_ij - lambda_ii)/RT] is a ratio of volumes
        # times an exponential, so it is strictly positive. A negative value makes
        # ln(SUM_j x_j Lambda_ij) complex somewhere in the composition range -- which
        # shows up as a bare `invalid value encountered in log`, not as an error.
        if np.any(self.Lambda <= 0):
            raise ValueError("Wilson Lambda_ij must be positive; got "
                             f"{self.Lambda.tolist()}")
        self.n = self.Lambda.shape[0]

    def _params(self):
        return {"Lambda": self.Lambda.tolist()}

    def lngamma(self, x, T):
        x = self._x(x)
        S = self.Lambda @ x                      # S_i = SUM_j x_j L_ij
        return 1.0 - np.log(S) - (self.Lambda.T @ (x / S))

    def gex_over_RT(self, x, T):
        x = self._x(x)
        return float(-(x @ np.log(self.Lambda @ x)))

    def infinite_dilution(self, printed=False):
        """The two infinite-dilution values of **ln gamma**, SIS Eq. 9.5-12b.

            ln gamma_1^inf = -ln L12 + 1 - L21
            ln gamma_2^inf = -ln L21 + 1 - L12

        **ERRATUM, SIS Eq. 9.5-12b.** The book prints `-ln L12 + 1 - L12` and
        `-ln L21 + 1 - L21` -- the *same* subscript twice in each expression. Taking
        x1 -> 0 in the book's own Eq. 9.5-12a,

            ln gamma_1 = -ln(x1 + x2 L12) + x2 [L12/(x1 + x2 L12) - L21/(x1 L21 + x2)]

        gives -ln L12 + 1 - L21: the first term's L12 cancels and what survives from
        the second is L21. The printed form is also not the standard Wilson result.
        For Lambda = (1.2594, 0.3783) -- a least-squares fit to the benzene /
        2,2,4-trimethyl pentane data -- the printed expression gives ln gamma_1^inf =
        -0.490 while Eq. 9.5-12a itself gives +0.391, i.e. gamma_1^inf = 0.61 against
        1.48. The sign of the deviation from ideality is wrong, not just its size.

        Set `printed=True` to reproduce the equation as printed.
        """
        if self.n != 2:
            raise ValueError("Eq. 9.5-12b is the binary result")
        L12, L21 = self.Lambda[0, 1], self.Lambda[1, 0]
        if printed:
            return np.array([-np.log(L12) + 1 - L12, -np.log(L21) + 1 - L21])
        return np.array([-np.log(L12) + 1 - L21, -np.log(L21) + 1 - L12])

    @classmethod
    def fit(cls, x1, gamma1, gamma2, p0=(1.0, 1.0)):
        """Least squares on ln gamma, with Lambda_ij constrained positive."""
        return fit_binary(cls, x1, gamma1, gamma2, T=None, p0=p0,
                          bounds=(1e-8, np.inf))


class NRTL(ActivityModel):
    """NRTL, SIS Eqs. 9.5-13 and 9.5-14a; multicomponent form from Appendix A9.2.

    Three parameters for a binary -- (alpha, tau12, tau21) -- with
    ln G_ij = -alpha tau_ij, and

        G^ex/RT = x1 x2 [ tau21 G21/(x1 + x2 G21) + tau12 G12/(x2 + x1 G12) ]

    Multicomponent, with G_ii = 1 and tau_ii = 0:

        G^ex/RT    = SUM_i x_i [SUM_j tau_ji G_ji x_j] / [SUM_k G_ki x_k]
        ln gamma_i = SUM_j tau_ji G_ji x_j / SUM_k G_ki x_k
                     + SUM_j (x_j G_ij / SUM_k G_kj x_k)
                       (tau_ij - SUM_m x_m tau_mj G_mj / SUM_k G_kj x_k)

    `alpha` is the nonrandomness parameter; 0.2-0.47 in practice, and Renon and
    Prausnitz's own recommendation is to fix it rather than fit it, which is why it
    is a separate argument here.

    **One instance is one temperature** -- tau_ij is an energy over RT, so it
    scales as 1/T. See the module docstring; `rescaled_to` moves it.
    """

    _T_LAW = "tau_ij = (g_ij - g_jj)/RT (Eq. 9.5-13)"

    def __init__(self, tau12=None, tau21=None, alpha=0.3, tau=None, alpha_matrix=None):
        if tau is not None:
            self.tau = np.asarray(tau, dtype=float)
            if self.tau.ndim != 2 or self.tau.shape[0] != self.tau.shape[1]:
                raise ValueError("tau must be square")
        else:
            if tau12 is None or tau21 is None:
                raise ValueError("give (tau12, tau21) or a full tau matrix")
            self.tau = np.array([[0.0, float(tau12)], [float(tau21), 0.0]])
        self.n = self.tau.shape[0]
        if alpha_matrix is not None:
            self.alpha = np.asarray(alpha_matrix, dtype=float)
        else:
            self.alpha = np.full((self.n, self.n), float(alpha))
            np.fill_diagonal(self.alpha, 0.0)
        self.G = np.exp(-self.alpha * self.tau)

    def _params(self):
        return {"tau": self.tau.tolist(), "alpha": float(self.alpha[0, 1])}

    def lngamma(self, x, T):
        x = self._x(x)
        G, tau = self.G, self.tau
        den = G.T @ x                              # den_j = SUM_k G_kj x_k
        num = (G * tau).T @ x                      # num_j = SUM_k tau_kj G_kj x_k
        first = num / den                          # index i via j->i
        second = (G * (tau - (num / den))) @ (x / den)
        return first + second

    def gex_over_RT(self, x, T):
        x = self._x(x)
        den = self.G.T @ x
        num = (self.G * self.tau).T @ x
        return float(x @ (num / den))

    def infinite_dilution(self):
        """SIS Eq. 9.5-14b: ln gamma_1^inf = tau21 + tau12 exp(-alpha tau12)."""
        if self.n != 2:
            raise ValueError("Eq. 9.5-14b is the binary result")
        t12, t21 = self.tau[0, 1], self.tau[1, 0]
        a = self.alpha[0, 1]
        return np.array([t21 + t12 * np.exp(-a * t12),
                         t12 + t21 * np.exp(-a * t21)])

    def rescaled_to(self, T):
        """tau ~ 1/T, so tau(T) = tau(T_ref) T_ref/T. Eq. 9.5-13."""
        T_ref = self.__dict__.get("_T_ref")
        if T_ref is None:
            raise ValueError(
                "this NRTL does not know which temperature its tau belongs to; "
                "call `.at_T(T_ref)` first, or build it from the energies")
        return type(self)(tau=self.tau * (T_ref / float(T)),
                          alpha_matrix=self.alpha).at_T(T)

    @classmethod
    def fit(cls, x1, gamma1, gamma2, alpha=0.3, p0=(0.5, 0.5)):
        """Fit tau12 and tau21 at fixed `alpha` -- Renon and Prausnitz's advice."""
        def make(p):
            return cls(tau12=p[0], tau21=p[1], alpha=alpha)
        return _fit_with(make, x1, gamma1, gamma2, T=None, p0=p0)


class UNIQUAC(ActivityModel):
    """UNIQUAC, SIS Eqs. 9.5-19 to 9.5-23b (Abrams and Prausnitz, 1975).

    Two adjustable parameters per binary pair -- tau12 and tau21 -- on top of the
    pure-species size parameters r_i and q_i, which come from molecular structure by
    group contribution (Table 9.5-2, and `UNIQUAC.from_groups`).

        ln gamma_i(combinatorial) = ln(phi_i/x_i) - (z/2) q_i ln(phi_i/theta_i)
                                    + l_i - (phi_i/x_i) SUM_j x_j l_j
        ln gamma_i(residual) = q_i [ 1 - ln(SUM_j theta_j tau_ji)
                                     - SUM_j theta_j tau_ij / SUM_k theta_k tau_kj ]

    with theta_i = x_i q_i / SUM_j x_j q_j the surface area fraction, phi_i =
    x_i r_i / SUM_j x_j r_j the segment (volume) fraction, l_i = (z/2)(r_i - q_i)
    - (r_i - 1), and z = 10 the coordination number.

    **Watch the subscript order in the residual term.** The tau inside the
    logarithm is tau_ji; the one in the numerator of the sum is tau_ij. tau is not
    symmetric, so swapping them changes every activity coefficient -- and the
    pure-component limit does *not* catch it, because gamma_i -> 1 either way. The
    error only shows up in a mixture. `check_gibbs_duhem` does catch it, because
    Eq. 9.5-21 is written with tau_ji and would then disagree.

    Parameters
    ----------
    r, q : per-species volume and surface area parameters.
    tau : n x n matrix with tau_ii = 1, or pass (tau12, tau21) for a binary.

    **One instance is one temperature.** para. 1303 defines ln tau_ij =
    -(u_ij - u_jj)/RT, so tau_ij = exp(-a_ij/T) and the constant is a_ij, not tau.
    See the module docstring; `rescaled_to` moves it. This is the model and the
    mistake behind Fig. 11.2-6.
    """

    Z = 10.0
    _T_LAW = "ln tau_ij = -(u_ij - u_jj)/RT (para. 1303)"

    def __init__(self, r, q, tau=None, tau12=None, tau21=None):
        self.r = np.asarray(r, dtype=float).ravel()
        self.q = np.asarray(q, dtype=float).ravel()
        if self.r.size != self.q.size:
            raise ValueError("r and q must have the same length")
        self.n = self.r.size
        if tau is None:
            if tau12 is None or tau21 is None:
                raise ValueError("give a tau matrix, or (tau12, tau21) for a binary")
            if self.n != 2:
                raise ValueError("(tau12, tau21) is the binary shorthand")
            tau = [[1.0, float(tau12)], [float(tau21), 1.0]]
        self.tau = np.asarray(tau, dtype=float)
        if self.tau.shape != (self.n, self.n):
            raise ValueError(f"tau must be {self.n}x{self.n}")
        self.l = (self.Z / 2.0) * (self.r - self.q) - (self.r - 1.0)

    def _params(self):
        return {"r": self.r.tolist(), "q": self.q.tolist(), "tau": self.tau.tolist()}

    @classmethod
    def from_groups(cls, groups, tau=None, tau12=None, tau21=None, R=None, Q=None):
        """Build r and q by group contribution from Table 9.5-2.

        `groups` is one dict per species, {subgroup_no: count}, exactly as
        `thermo.UNIFAC.gamma` takes them -- so a notebook can hand the same molecular
        description to both models. R and Q default to the tabulated values.

        This is the calculation of Illustration 9.5-2: benzene as six ACH groups and
        2,2,4-trimethyl pentane as 5 CH3 + CH2 + CH + C.
        """
        if R is None or Q is None:
            from .data import load_unifac_subgroups
            subs = load_unifac_subgroups()
            R = dict(zip(subs.subgroup_no, subs.R))
            Q = dict(zip(subs.subgroup_no, subs.Q))
        r = [sum(n * R[k] for k, n in g.items()) for g in groups]
        q = [sum(n * Q[k] for k, n in g.items()) for g in groups]
        return cls(r, q, tau=tau, tau12=tau12, tau21=tau21)

    # -- the size fractions, which the illustrations ask for on their own ---
    def theta(self, x):
        """Surface area fractions, SIS Sec. 9.5. Printed as theta_i."""
        x = self._x(x)
        return x * self.q / float(x @ self.q)

    def phi(self, x):
        """Segment (volume) fractions, SIS Sec. 9.5. Printed as phi_i."""
        x = self._x(x)
        return x * self.r / float(x @ self.r)

    def lngamma(self, x, T):
        x = self._x(x)
        th, ph = self.theta(x), self.phi(x)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(x > 0, ph / np.where(x > 0, x, 1.0), 0.0)
        comb = (np.log(ratio) - (self.Z / 2) * self.q * np.log(ph / th)
                + self.l - ratio * float(x @ self.l))
        S = self.tau.T @ th                      # S_j = SUM_k theta_k tau_kj
        res = self.q * (1.0 - np.log(S) - (self.tau @ (th / S)))
        return comb + res

    def gex_over_RT(self, x, T):
        """SIS Eqs. 9.5-20 + 9.5-21, which are an independent statement of the model.

        Agreement with SUM_i x_i ln gamma_i is the check described in the module
        docstring, and it is what pins down the tau_ij / tau_ji order.
        """
        x = self._x(x)
        th, ph = self.theta(x), self.phi(x)
        nz = x > 0
        comb = float(np.sum(x[nz] * np.log(ph[nz] / x[nz]))
                     + (self.Z / 2) * np.sum(x[nz] * self.q[nz]
                                             * np.log(th[nz] / ph[nz])))
        res = -float(np.sum(self.q * x * np.log(self.tau.T @ th)))
        return comb + res

    def rescaled_to(self, T):
        """tau = exp(-a/T), so tau(T) = tau(T_ref)**(T_ref/T). para. 1303.

        tau_ii = 1 is a fixed point of the exponentiation, so the diagonal is safe.
        """
        T_ref = self.__dict__.get("_T_ref")
        if T_ref is None:
            raise ValueError(
                "this UNIQUAC does not know which temperature its tau belongs to; "
                "call `.at_T(T_ref)` first, or build it from the energies")
        return type(self)(self.r, self.q,
                          tau=self.tau ** (T_ref / float(T))).at_T(T)

    @classmethod
    def fit(cls, r, q, x1, gamma1, gamma2, p0=(1.0, 1.0)):
        """Fit tau12 and tau21 with r and q held at their structural values."""
        def make(p):
            return cls(r, q, tau12=p[0], tau21=p[1])
        return _fit_with(make, x1, gamma1, gamma2, T=None, p0=p0,
                         bounds=(1e-8, np.inf))


# ---------------------------------------------------------------------------
# Sec. 9.6 -- the predictive model that needs no mixture data
# ---------------------------------------------------------------------------
class RegularSolution(ActivityModel):
    """Regular solution theory, SIS Eqs. 9.6-8 to 9.6-11 (Scatchard-Hildebrand).

        RT ln gamma_i = V_i (delta_i - deltabar)^2      deltabar = SUM_j Phi_j delta_j
        Phi_i = x_i V_i / SUM_j x_j V_j

    which for a binary is Eq. 9.6-10, RT ln gamma_1 = V_1 Phi_2^2 (delta_1 - delta_2)^2.

    **This is the only model in the module that needs no mixture data at all** --
    both parameters are pure-component properties. That is its whole appeal, and
    Illustration 9.6-1 is the book's own warning about the price: for benzene /
    2,2,4-trimethyl pentane the prediction is qualitatively right and quantitatively
    poor. Sec. 9.6 also notes it can only give gamma_i >= 1, so it cannot describe a
    negative deviation from ideality at all.

    Parameters (SI)
    ---------------
    V : molar liquid volumes, m^3/mol.
    delta : solubility parameters, Pa^(1/2).

    Use `from_table_9_6_1` to enter the book's cc/mol and (cal/cc)^(1/2) as printed.
    """

    def __init__(self, V, delta):
        self.V = np.asarray(V, dtype=float).ravel()
        self.delta = np.asarray(delta, dtype=float).ravel()
        if self.V.size != self.delta.size:
            raise ValueError("V and delta must have the same length")
        self.n = self.V.size

    def _params(self):
        return {"V": self.V.tolist(), "delta": self.delta.tolist()}

    @classmethod
    def from_table_9_6_1(cls, V_cc_mol, delta_cal_cc_half):
        """Construct from the traditional units Table 9.6-1 prints.

            >>> rs = RegularSolution.from_table_9_6_1([89, 165], [9.2, 6.93])

        is benzene / 2,2,4-trimethyl pentane as Illustration 9.6-1 sets it up.
        """
        return cls(np.asarray(V_cc_mol, float) * CC_PER_MOL,
                   np.asarray(delta_cal_cc_half, float) * CAL_CC_HALF)

    @classmethod
    def from_book_names(cls, names):
        """Look species up in `TABLE_9_6_1` by the name the table prints."""
        rows = []
        for nm in names:
            key = nm.strip().lower()
            if key not in TABLE_9_6_1:
                raise KeyError(f"{nm!r} is not in Table 9.6-1; "
                               f"pass V and delta explicitly")
            rows.append(TABLE_9_6_1[key])
        V, d = zip(*rows)
        return cls.from_table_9_6_1(V, d)

    def volume_fraction(self, x):
        """Phi_i, SIS Eq. 9.6-8. Printed with a *capital* Phi, unlike the
        lowercase phi_i of UNIQUAC -- a different measure of size (molar volume
        rather than the group parameter r_i), and the chapter keeps both symbols."""
        x = self._x(x)
        return x * self.V / float(x @ self.V)

    def delta_bar(self, x):
        return float(self.volume_fraction(x) @ self.delta)

    def lngamma(self, x, T):
        x = self._x(x)
        return self.V * (self.delta - self.delta_bar(x)) ** 2 / (R * T)

    @staticmethod
    def delta_from_vaporization(dHvap, V, T):
        """delta = [(dHvap - RT)/V]^(1/2), the estimate used in Illustration 9.6-1
        for 2,2,4-trimethyl pentane, whose delta Table 9.6-1 does not give.

        dHvap in J/mol, V in m^3/mol, T in K; returns Pa^(1/2). The internal energy
        of vaporization is approximated as dHvap - RT (ideal vapor, liquid volume
        neglected).
        """
        return np.sqrt((float(dHvap) - R * float(T)) / float(V))


# ---------------------------------------------------------------------------
# fitting
# ---------------------------------------------------------------------------
def _fit_with(make, x1, gamma1, gamma2, T, p0, bounds=(-np.inf, np.inf)):
    """Least squares on ln gamma for both species, given a parameter->model factory."""
    x1 = np.asarray(x1, float)
    obs = np.concatenate([np.log(np.asarray(gamma1, float)),
                          np.log(np.asarray(gamma2, float))])

    def resid(p):
        m = make(p)
        ln = np.array([m.lngamma([a, 1.0 - a], T if T is not None else 300.0)
                       for a in x1])
        return np.concatenate([ln[:, 0], ln[:, 1]]) - obs

    sol = optimize.least_squares(resid, np.asarray(p0, float), bounds=bounds)
    if not sol.success:
        raise RuntimeError(f"activity-model fit did not converge: {sol.message}")
    return make(sol.x)


def fit_binary(model_cls, x1, gamma1, gamma2, T=None, p0=(1.0, 1.0), **kw):
    """Fit a two-parameter binary model to activity coefficient data.

    Works for `VanLaar` and `Wilson`, whose two constructor arguments *are* the two
    parameters. `NRTL.fit` and `UNIQUAC.fit` have their own signatures because they
    carry a third parameter (alpha) or structural constants (r, q) that are not
    fitted.

    T is accepted for interface symmetry and is only used by models whose activity
    coefficients depend on it -- which, among the two-parameter binaries, is none of
    them.
    """
    return _fit_with(lambda p: model_cls(p[0], p[1]), x1, gamma1, gamma2, T, p0, **kw)


# ---------------------------------------------------------------------------
# the book's own tables
# ---------------------------------------------------------------------------
#
# Table 9.5-1, the van Laar constants (SIS p. 451). Keyed on the printed
# "Component 1-Component 2" pair, lowercased. Entries with two temperature ranges in
# the book appear twice, distinguished by the range. The source is an adaptation of
# Perry's Chemical Engineers' Handbook, 4th ed. (1963), p. 13-7.
#
# Only the alpha/beta pairs are stored; the temperature range is the second element
# and is text, because it is what tells a reader whether the constants apply.
TABLE_9_5_1_VAN_LAAR = {
    "acetaldehyde-water":                    (1.59, 1.80, "19.8-100"),
    "acetone-benzene":                       (0.405, 0.405, "56.1-80.1"),
    "acetone-methanol":                      (0.58, 0.56, "56.1-64.6"),
    "acetone-water":                         (1.89, 1.66, "25"),
    "acetone-water@56.1-100":                (2.05, 1.50, "56.1-100"),
    "benzene-isopropanol":                   (1.36, 1.95, "71.9-82.3"),
    "carbon disulfide-acetone":              (1.28, 1.79, "39.5-56.1"),
    "carbon disulfide-carbon tetrachloride": (0.23, 0.16, "46.3-76.7"),
    "carbon tetrachloride-benzene":          (0.12, 0.11, "76.4-80.2"),
    "ethanol-benzene":                       (1.946, 1.610, "67.0-80.1"),
    "ethanol-cyclohexane":                   (2.102, 1.729, "66.3-80.8"),
    "ethanol-toluene":                       (1.757, 1.757, "76.4-110.7"),
    "ethanol-water":                         (1.54, 0.97, "25"),
    "ethyl acetate-benzene":                 (1.15, 0.92, "71.1-80.2"),
    "ethyl acetate-ethanol":                 (0.896, 0.896, "71.7-78.3"),
    "ethyl acetate-toluene":                 (0.09, 0.58, "77.2-110.7"),
    "ethyl ether-acetone":                   (0.741, 0.741, "34.6-56.1"),
    "ethyl ether-ethanol":                   (0.97, 1.27, "34.6-78.3"),
    "n-hexane-ethanol":                      (1.57, 2.58, "59.3-78.3"),
    "isobutane-furfural":                    (2.62, 3.02, "37.8"),
    "isobutane-furfural@51.7":               (2.51, 2.83, "51.7"),
    "isopropanol-water":                     (2.40, 1.13, "82.3-100"),
    "methanol-benzene":                      (0.56, 0.56, "55.5-64.6"),
    "methanol-ethyl acetate":                (1.16, 1.16, "62.1-77.1"),
    "methanol-water":                        (0.58, 0.46, "25"),
    "methanol-water@64.6-100":               (0.83, 0.51, "64.6-100"),
    "methyl acetate-methanol":               (1.06, 1.06, "53.7-64.6"),
    "methyl acetate-water":                  (2.99, 1.89, "57.0-100"),
    "n-propanol-water":                      (2.53, 1.13, "88.0-100"),
    "water-phenol":                          (0.83, 3.22, "100-181"),
}

#
# Table 9.6-1, molar liquid volumes and solubility parameters of some nonpolar
# liquids (SIS p. 464), in the units it prints: (V^L in cc/mol, delta in
# (cal/cc)^(1/2)). Source: J. M. Prausnitz, *Molecular Thermodynamics of Fluid-Phase
# Equilibria*, Prentice-Hall (1969).
#
# The table is in two blocks measured at different conditions -- liquefied gases
# at 90 K and liquid solvents at 25 C -- and the book does not say the parameters
# transfer between them. Sec. 9.6 does note that (delta_1 - delta_2) is often nearly
# temperature independent, so the table "may be used at temperatures other than the
# one at which they were obtained." That license is the book's; the block a species
# came from is recorded here so a notebook can say which it is relying on.
TABLE_9_6_1 = {
    # liquefied gases at 90 K
    "nitrogen":              (38.1, 5.3),
    "carbon monoxide":       (37.1, 5.7),
    "argon":                 (29.0, 6.8),
    "oxygen":                (28.0, 7.2),
    "methane":               (35.3, 7.4),
    "carbon tetrafluoride":  (46.0, 8.3),
    "ethane":                (45.7, 9.5),
    # liquid solvents at 25 C
    "perfluoro-n-heptane":   (226.0, 6.0),
    "neopentane":            (122.0, 6.2),
    "isopentane":            (117.0, 6.8),
    "n-pentane":             (116.0, 7.1),
    "n-hexane":              (132.0, 7.3),
    "1-hexene":              (126.0, 7.3),
    "n-octane":              (164.0, 7.5),
    "n-hexadecane":          (294.0, 8.0),
    "cyclohexane":           (109.0, 8.2),
    "carbon tetrachloride":  (97.0, 8.6),
    "ethyl benzene":         (123.0, 8.8),
    "toluene":               (107.0, 8.9),
    "benzene":               (89.0, 9.2),
    "styrene":               (116.0, 9.3),
    "tetrachloroethylene":   (103.0, 9.3),
    "carbon disulfide":      (61.0, 10.0),
    "bromine":               (51.0, 11.5),
}

TABLE_9_6_1_BLOCK = {k: ("liquefied gas at 90 K" if i < 7 else "liquid solvent at 25 C")
                     for i, k in enumerate(TABLE_9_6_1)}
