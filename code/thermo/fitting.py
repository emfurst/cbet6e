"""Redlich-Kister correlation of mixing data, and the partial molar properties
that follow from it.

This is the machinery of SIS Sec. 8.6 -- *The Experimental Determination of the
Partial Molar Volume and Enthalpy* -- and the Python replacement for the 4th-edition
`PRTLMOLR` worksheet.

The physical content is one idea. A property change on mixing must vanish at both
pure-component limits, so it is fitted with a polynomial that vanishes there by
construction (Eq. 8.6-5a),

    dmix(theta) = x1 x2 SUM_i a_i (x1 - x2)^i

and the partial molar properties then follow from the *slope* of that curve
(Eqs. 8.6-6a,b) rather than from a tangent line drawn by hand on graph paper:

    thetabar_1 - theta_1 = x2^2 SUM_i a_i [(x1-x2)^i + 2 i x1 (x1-x2)^(i-1)]
    thetabar_2 - theta_2 = x1^2 SUM_i a_i [(x1-x2)^i - 2 i x2 (x1-x2)^(i-1)]

Two things about that are worth saying out loud, because they are the whole reason
this module exists rather than a call to `numpy.polyfit`.

**The derivative is much more sensitive to the order than the fit is.** Adding a
term barely moves the curve and can move the infinite-dilution partial molar
properties -- which are the endpoints, where the extrapolation is worst -- by ten
percent. `loo_rms` and `infinite_dilution` are here so that a notebook can *show*
that instead of hiding it behind a default. See `RedlichKister.scan_order`.

**The summation identity holds exactly, by construction.** Whatever order is
chosen, x1*(thetabar_1 - theta_1) + x2*(thetabar_2 - theta_2) reproduces the fitted
dmix(theta) to machine precision, because Eqs. 8.6-6a,b are an exact rearrangement.
So agreement there is *not* evidence the fit is good -- `residuals` is. A
hand-drawn tangent construction, by contrast, does not satisfy it: the 5e's own
Table 8.6-4 misses it by up to 16 J/mol (see `code/ch8/validation/`).

Units are whatever the data are in: this module is unit-agnostic, and the caller
keeps SI (m^3/mol, J/mol) as everywhere else in the package.

Eric M. Furst
August 2026
"""
import numpy as np

__all__ = ["RedlichKister", "tangent_intercepts", "LogPolynomialActivity",
           "ILLUSTRATION_9_7_1_HEMOGLOBIN", "SIS_9_7_1_CONSTANTS"]


class RedlichKister:
    """A Redlich-Kister expansion of a property change on mixing (SIS Eq. 8.6-5a).

    Parameters
    ----------
    a : sequence of coefficients a_0 ... a_n, in the units of the property.

    Construct from data with `RedlichKister.fit`, or from published constants
    directly -- the book prints a set for the water-methanol volume data on p.389:

        >>> rk = RedlichKister([-4.0034e-6, -0.17756e-6, 0.54139e-6, 0.60481e-6])
        >>> d1, d2 = rk.partial_molar_excess(0.5266)
    """

    def __init__(self, a):
        self.a = np.asarray(a, dtype=float)
        if self.a.ndim != 1 or self.a.size == 0:
            raise ValueError("a must be a non-empty 1-D sequence of coefficients")

    # -- construction ------------------------------------------------------

    @classmethod
    def fit(cls, x1, dmix, order=3):
        """Least-squares fit of order `order` (so order+1 coefficients).

        The model is linear in the coefficients, so this is one `lstsq` call and
        has no starting guess and no convergence to worry about.
        """
        x1 = np.asarray(x1, dtype=float)
        dmix = np.asarray(dmix, dtype=float)
        if x1.shape != dmix.shape:
            raise ValueError("x1 and dmix must have the same shape")
        if x1.size < order + 1:
            raise ValueError(f"order {order} needs at least {order + 1} points, "
                             f"got {x1.size}")
        a, *_ = np.linalg.lstsq(cls._design(x1, order), dmix, rcond=None)
        return cls(a)

    @staticmethod
    def _design(x1, order):
        """Columns x1*x2*(x1-x2)^i, i = 0..order."""
        x1 = np.asarray(x1, dtype=float)
        x2 = 1.0 - x1
        u = x1 - x2
        return np.column_stack([x1 * x2 * u ** i for i in range(order + 1)])

    @property
    def order(self):
        return self.a.size - 1

    # -- the fitted curve and its slope ------------------------------------

    def __call__(self, x1):
        """dmix(theta) at composition x1 (Eq. 8.6-5a)."""
        x1 = np.asarray(x1, dtype=float)
        x2 = 1.0 - x1
        u = x1 - x2
        return x1 * x2 * sum(ai * u ** i for i, ai in enumerate(self.a))

    def derivative(self, x1):
        """d[dmix(theta)]/dx1 at constant T, P (Eq. 8.6-5b)."""
        x1 = np.asarray(x1, dtype=float)
        u = 2.0 * x1 - 1.0
        first = -sum(ai * u ** (i + 1) for i, ai in enumerate(self.a))
        second = 2.0 * x1 * (1.0 - x1) * sum(
            i * ai * u ** (i - 1) for i, ai in enumerate(self.a) if i)
        return first + second

    # -- the partial molar properties --------------------------------------

    def partial_molar_excess(self, x1):
        """(thetabar_1 - theta_1, thetabar_2 - theta_2) at x1 (Eqs. 8.6-6a,b).

        Returned as differences from the pure-component molar values, which is
        the form Tables 8.6-2 and 8.6-4 print; add the pure values for thetabar
        itself, or use `partial_molar`.
        """
        x1 = np.asarray(x1, dtype=float)
        x2 = 1.0 - x1
        u = x1 - x2
        d1 = x2 ** 2 * sum(
            ai * (u ** i + (2 * i * x1 * u ** (i - 1) if i else 0.0))
            for i, ai in enumerate(self.a))
        d2 = x1 ** 2 * sum(
            ai * (u ** i - (2 * i * x2 * u ** (i - 1) if i else 0.0))
            for i, ai in enumerate(self.a))
        return d1, d2

    def partial_molar(self, x1, theta1_pure, theta2_pure):
        """(thetabar_1, thetabar_2) at x1, given the pure-component molar values."""
        d1, d2 = self.partial_molar_excess(x1)
        return theta1_pure + d1, theta2_pure + d2

    def infinite_dilution(self):
        """(thetabar_1 - theta_1 as x1->0, thetabar_2 - theta_2 as x1->1).

        Both limits are the expansion summed at an endpoint -- SUM a_i (-1)^i and
        SUM a_i -- which is where an extra term does the most damage. Compare
        across orders before quoting either (`scan_order`).
        """
        a = self.a
        return (float(sum(ai * (-1.0) ** i for i, ai in enumerate(a))),
                float(a.sum()))

    # -- how good is it ----------------------------------------------------

    def residuals(self, x1, dmix):
        """Fitted minus observed."""
        return self(np.asarray(x1, dtype=float)) - np.asarray(dmix, dtype=float)

    def rms(self, x1, dmix):
        return float(np.sqrt(np.mean(self.residuals(x1, dmix) ** 2)))

    @classmethod
    def loo_rms(cls, x1, dmix, order):
        """Leave-one-out cross-validated rms error at this order.

        `rms` always falls as terms are added; this does not, so it is the one
        that can say a term is not paying for itself.
        """
        x1 = np.asarray(x1, dtype=float)
        dmix = np.asarray(dmix, dtype=float)
        err = np.empty(x1.size)
        for k in range(x1.size):
            keep = np.ones(x1.size, dtype=bool)
            keep[k] = False
            fit = cls.fit(x1[keep], dmix[keep], order)
            err[k] = fit(x1[k]) - dmix[k]
        return float(np.sqrt(np.mean(err ** 2)))

    @classmethod
    def scan_order(cls, x1, dmix, orders=range(2, 8)):
        """Fit at several orders and report what changes.

        Returns a list of dicts with keys `order`, `rms`, `loo_rms`,
        `inf_dilution_1`, `inf_dilution_2`, `fit`. The point of collecting them
        together is the comparison the notebooks make: `rms` improves smoothly
        while the two infinite-dilution numbers wander, and choosing an order on
        `rms` alone quietly picks an endpoint value.
        """
        out = []
        for n in orders:
            fit = cls.fit(x1, dmix, n)
            i1, i2 = fit.infinite_dilution()
            out.append({"order": n,
                        "rms": fit.rms(x1, dmix),
                        "loo_rms": cls.loo_rms(x1, dmix, n),
                        "inf_dilution_1": i1,
                        "inf_dilution_2": i2,
                        "fit": fit})
        return out

    def __repr__(self):
        return f"RedlichKister(order={self.order}, a={np.array2string(self.a, precision=5)})"


def tangent_intercepts(x1_star, dmix_star, slope):
    """The graphical construction of Sec. 8.6, as arithmetic.

    Draw the tangent to the dmix(theta)-versus-x1 curve at `x1_star`; its
    intercepts with the ordinates at x1 = 0 (point A in Fig. 8.6-1) and x1 = 1
    (point B) are (thetabar_2 - theta_2) and (thetabar_1 - theta_1) at that
    composition -- Eqs. 8.6-4a and 8.6-4b.

    Returns (A, B) = (thetabar_2 - theta_2, thetabar_1 - theta_1).

    This is the same answer `RedlichKister.partial_molar_excess` gives, by a
    different route, and the notebooks use it as an independent check: it needs
    only a point and a slope, so it will agree with *any* correlation that passes
    through the data with that slope.
    """
    A = dmix_star - x1_star * slope
    B = dmix_star + (1.0 - x1_star) * slope
    return A, B


# ---------------------------------------------------------------------------
# Illustration 9.7-1 -- the activity coefficient of a protein in solution.
#
# A protein is a *nonsimple* mixture in the sense of SIS Sec. 9.7: it does not
# exist as a pure liquid, so its activity coefficient is defined on a Henry's law
# basis and referred to concentration rather than mole fraction,
#
#     gamma(c) = a(c) / c        with gamma -> 1 as c -> 0
#
# and correlated with an expansion in c whose constant term is zero, so that the
# limit is satisfied by construction:
#
#     ln gamma = A c + B c^2 + C c^3                            (SIS Eq. 9.7-21)
#
# The order matters less than the OBJECTIVE. Fitting ln gamma and fitting
# gamma weight the fifteen points very differently -- gamma spans 1.16 to 36.8,
# so a least-squares fit on gamma is dominated by the last three points while a
# fit on ln gamma is not. The two give coefficients that differ by ~5% in B.
# `fit` therefore takes `objective` explicitly and has no default that hides the
# choice, and `LogPolynomialActivity.report` prints which one was used.
# ---------------------------------------------------------------------------

#: Illustration 9.7-1, as printed: concentration and activity of hemoglobin in
#: aqueous solution, both in g/L. The solution's gamma column is *not* stored --
#: it is gamma = a/c exactly, so it is a check on the transcription, not data.
ILLUSTRATION_9_7_1_HEMOGLOBIN = np.array([
    # c (g/L),  a (g/L)
    [ 20.0,      23.1],
    [ 40.0,      53.5],
    [ 60.0,      94.6],
    [ 80.0,     150.0],
    [100.0,     226.0],
    [120.0,     330.0],
    [140.0,     473.0],
    [160.0,     679.0],
    [180.0,     973.0],
    [200.0,    1410.0],
    [220.0,    2060.0],
    [240.0,    3040.0],
    [260.0,    4580.0],
    [280.0,    7040.0],
    [300.0,   11050.0],
])

#: The constants printed in Illustration 9.7-1's solution. Read out of the
#: MathType record in the Word source; the exponents are not printed in a form
#: any text extraction recovers, and were pinned by testing against the data.
SIS_9_7_1_CONSTANTS = (7.139e-3, 6.940e-6, 3.116e-8)


class LogPolynomialActivity:
    """ln gamma = A c + B c^2 + ... , the correlation of SIS Eq. 9.7-21.

    Parameters
    ----------
    coeffs : sequence A, B, C, ... of coefficients of c, c^2, c^3, ...
             There is deliberately no constant term: gamma -> 1 as c -> 0.
    objective : the string recorded by `fit`, or None when the coefficients were
             taken from a publication rather than fitted here.
    """

    def __init__(self, coeffs, objective=None):
        self.coeffs = np.asarray(coeffs, float)
        self.objective = objective

    @property
    def order(self):
        return len(self.coeffs)

    def ln_gamma(self, c):
        c = np.asarray(c, float)
        return sum(a * c ** (k + 1) for k, a in enumerate(self.coeffs))

    def gamma(self, c):
        return np.exp(self.ln_gamma(c))

    def residuals(self, c, gamma):
        """Relative residuals in gamma, (fitted - observed)/observed."""
        return (self.gamma(c) - np.asarray(gamma, float)) / np.asarray(gamma, float)

    def rms(self, c, gamma):
        return float(np.sqrt((self.residuals(c, gamma) ** 2).mean()))

    def max_error(self, c, gamma):
        return float(np.abs(self.residuals(c, gamma)).max())

    @classmethod
    def fit(cls, c, gamma, order=3, objective="ln gamma"):
        """Least-squares fit of the expansion to gamma(c).

        objective : "ln gamma" fits the logarithm, weighting every point equally
                    in ln gamma. "gamma" fits gamma itself, which weights the
                    high-concentration points far more heavily. The choice is
                    required rather than defaulted-away because it moves the
                    coefficients by more than the experimental scatter does.
        """
        c = np.asarray(c, float)
        gamma = np.asarray(gamma, float)
        A = np.vstack([c ** (k + 1) for k in range(order)]).T
        if objective == "ln gamma":
            coeffs, *_ = np.linalg.lstsq(A, np.log(gamma), rcond=None)
        elif objective == "gamma":
            # weight the log-space normal equations by gamma, which is the
            # first-order equivalent of least squares on gamma itself
            coeffs, *_ = np.linalg.lstsq(A * gamma[:, None], np.log(gamma) * gamma,
                                         rcond=None)
        else:
            raise ValueError("objective must be 'ln gamma' or 'gamma', "
                             f"not {objective!r}")
        return cls(coeffs, objective=objective)

    def report(self, c, gamma, label=""):
        """One line: the coefficients, the objective, and how well it does."""
        names = "ABCDEFG"
        terms = "  ".join(f"{names[k]}={a:10.4e}" for k, a in enumerate(self.coeffs))
        obj = self.objective or "as printed"
        return (f"{label:26s} {terms}   [{obj:9s}]  "
                f"max|d|={100*self.max_error(c, gamma):5.2f}%  "
                f"rms={100*self.rms(c, gamma):5.2f}%")

    def __repr__(self):
        return (f"LogPolynomialActivity(order={self.order}, "
                f"objective={self.objective!r})")
