"""Combined equation-of-state + excess Gibbs energy model -- SIS Section 9.9.

The **Wong-Sandler mixing rules**: a way of putting an activity coefficient model
*inside* a cubic equation of state, so that one model describes a highly nonideal
mixture in both phases and over a wide range of temperature and pressure.

⭐ **This is the one place the book teaches a method its own author invented**
(Wong and Sandler, *AIChE J.* **38**, 671 (1992)), and in the 5e it is prose only:
Sec. 9.9 has no illustration, no figure, and no worked number anywhere in the
chapter. Everything here therefore rests on the section's own derivation, and the
module is written so that the two boundary conditions the derivation is *built from*
can be checked numerically rather than taken on faith -- see `check_boundary`.

## The idea, in the section's own order

A cubic equation of state has two constants, `a` and `b`, so two conditions can be
imposed on them.

1. **At low density**, the second virial coefficient from the EOS,
   `B = b - a/RT` (Eq. 9.9-2), must be quadratic in composition as statistical
   mechanics requires (Eq. 9.9-1). That gives

       b - a/RT = SUM_i SUM_j x_i x_j (b_ij - a_ij/RT)              (Eq. 9.9-3)

   -- one equation, fixing the *combination*, not `a` and `b` separately.

2. **At liquid densities**, the excess Helmholtz energy from the EOS must match the
   one from an activity coefficient model (Eq. 9.9-4). Helmholtz and not Gibbs
   because `G^ex = A^ex + P V^ex` (Eq. 9.9-5) diverges as P -> infinity while `A^ex`
   is nearly pressure independent, so `A^ex(T, P->inf) ~ A^ex(T, 1 bar) =
   G^ex(T, 1 bar)` (Eqs. 9.9-6, 9.9-7). For a cubic,

       A^ex_EOS = C* [ a/b - SUM_i x_i a_i/b_i ]                    (Eq. 9.9-8)

   with `C*` a pure number set by the equation of state.

Solving the two together gives the mixing rules (Eq. 9.9-9):

    a/RT = Q D / (1 - D)        b = Q / (1 - D)
    Q = SUM_i SUM_j x_i x_j (b_ij - a_ij/RT)
    D = SUM_i x_i a_i/(b_i RT) + G^ex(T, x) / (C* RT)

## Why this is not just another mixing rule

The van der Waals one-fluid rules make `b` **linear** and `a` **quadratic** in mole
fraction, because that is the simplest thing to write. Here neither is: both come out
of `Q/(1-D)`, and `D` carries the whole activity coefficient model. The composition
dependence is inherited from the G^ex model rather than assumed -- which is why the
same equation of state can then describe acetone/water, and why Sec. 10.3 can fit at
one temperature and predict 100-200 degrees away.

⚠️ **Setting k_ij = 0 is NOT the ideal-solution case here.** With the van der Waals
rules it very nearly is; with these, `D` still carries G^ex, so a zero binary
parameter can still produce extreme nonideality. Sec. 10.3 makes exactly this point
about acetone/water, and `code/ch9/` demonstrates it.

## Units and conventions

SI throughout (T in K, P in Pa, a in Pa m^6/mol^2, b in m^3/mol, G^ex in J/mol), as
everywhere else in the package. `kij` is dimensionless and symmetric.

Eric M. Furst
August 2026
"""
import numpy as np
from scipy import constants

from .peng_robinson import PengRobinson
from .phi_phi import PhiPhiVLE

R = constants.R
_SQRT2 = np.sqrt(2.0)

__all__ = ["WongSandler", "GexFromUNIFAC", "C_STAR_PR", "C_STAR_VDW"]

# C* of Eq. 9.9-8, the P -> infinity limit of the EOS excess Helmholtz energy. The
# book gives both values in the text following the equation: -1 for van der Waals,
# and ln(sqrt(2) - 1)/sqrt(2) for Peng-Robinson, which it prints as -0.62323.
C_STAR_PR = float(np.log(_SQRT2 - 1.0) / _SQRT2)      # -0.6232252401...
C_STAR_VDW = -1.0


class GexFromUNIFAC:
    """Adapter: make a `UNIFAC` instance look like an `activity_models` model.

    `UNIFAC.gamma` takes the group assignments as its first argument, so it does not
    have the `gamma(x, T)` / `gex(x, T)` signature `WongSandler` wants. Binding the
    groups once here is cleaner than teaching the mixing rule about two shapes of
    activity model.

    >>> from thermo import UNIFAC
    >>> from thermo.data import unifac_groups
    >>> g = GexFromUNIFAC(UNIFAC(), [{18: 1, 1: 1}, unifac_groups("water")])

    ⭐ This is Sec. 10.3's fourth route -- the mixing rule fed by a *predictive* model,
    so a mixture with no data at all can be carried to high pressure.
    """

    def __init__(self, unifac, groups):
        self.unifac = unifac
        self.groups = list(groups)

    def gamma(self, x, T):
        return np.asarray(self.unifac.gamma(self.groups, np.asarray(x, float), T),
                          dtype=float)

    def gex(self, x, T):
        """G^ex = RT SUM_i x_i ln gamma_i (Eq. 9.3-12), J/mol."""
        x = np.asarray(x, dtype=float)
        return float(R * T * (x @ np.log(self.gamma(x, T))))

    def __repr__(self):
        return f"GexFromUNIFAC({self.unifac!r}, {self.groups!r})"


class WongSandler(PhiPhiVLE):
    """Peng-Robinson with the Wong-Sandler mixing rules of SIS Eq. 9.9-9.

    Parameters
    ----------
    components : list of PengRobinson
        The pure fluids. ⚠️ Sec. 10.3 says to use **PRSV** here (`kappa1=`, Eqs.
        7.5-1/7.5-2 with kappa1 fitted to each pure vapor pressure) rather than the
        generalized kappa of Eq. 6.7-4, "to ensure that the pure component vapor
        pressures are correct." Nothing in this class requires it, and it does not
        check -- but a mixing rule cannot repair a pure-component vapor pressure.
    gex : activity coefficient model
        Anything exposing `gex(x, T)` in J/mol and `gamma(x, T)`; every model in
        `thermo.activity_models` qualifies, and `GexFromUNIFAC` adapts UNIFAC.
    kij : float, (n, n) array, or None
        The single binary parameter of Eq. 9.9-10. `None` means zero -- which, see
        above, is *not* the ideal-solution case.
    combining : '9.9-10b' (default) or '9.9-10a'
        Which combining rule supplies the cross term. The book offers both and notes
        that 10b "has the advantage of being more similar to the van der Waals
        one-fluid mixing rules."
    C_star : float, optional
        Defaults to the Peng-Robinson value.
    """

    def __init__(self, components, gex, kij=None, combining="9.9-10b", C_star=None):
        self.components = list(components)
        n = len(self.components)
        if n < 2:
            raise ValueError("a mixture needs at least two components")
        if combining not in ("9.9-10a", "9.9-10b"):
            raise ValueError("combining must be '9.9-10a' or '9.9-10b'")
        self.gex_model = gex
        self.combining = combining
        self.C_star = C_STAR_PR if C_star is None else float(C_star)
        if kij is None:
            self.kij = np.zeros((n, n))
        else:
            self.kij = (np.asarray(kij, dtype=float) if np.ndim(kij)
                        else np.full((n, n), float(kij)) - np.eye(n) * float(kij))
            if self.kij.shape != (n, n):
                raise ValueError(f"kij must be {n}x{n}")
        self.b_pure = np.array([c.b for c in self.components])

    @classmethod
    def from_database(cls, keys, gex, **kw):
        """Build the pure fluids from `pure_property.csv` by name."""
        return cls([PengRobinson.from_database(k) for k in keys], gex, **kw)

    @property
    def n(self):
        return len(self.components)

    def a_pure(self, T):
        return np.array([c.a(T) for c in self.components])

    # -- the cross term, Eq. 9.9-10 ---------------------------------------
    def cross_matrix(self, T):
        """The matrix of (b - a/RT)_ij, Eq. 9.9-10a or 9.9-10b.

        Its diagonal is the pure-component second virial coefficient b_i - a_i/RT,
        under either rule -- which is the check that the combining rule is written
        correctly, and `check_boundary` runs it.
        """
        a, b = self.a_pure(T), self.b_pure
        Bi = b - a / (R * T)                       # pure second virial coefficients
        if self.combining == "9.9-10a":
            # ⚠️ THE PRINTED FORM IS SIGN-AMBIGUOUS, and taking it literally is
            # catastrophic. Below the Boyle temperature -- which is everywhere these
            # models are used -- every B_i is NEGATIVE. The product of two negatives
            # is positive, so the principal square root of Eq. 9.9-10a is POSITIVE,
            # giving a cross term of the opposite sign to both diagonals. Q then
            # passes through zero somewhere in the middle of the composition range
            # and b = Q/(1-D) collapses with it.
            #
            # The only reading consistent with the rule reducing to B_i on the
            # diagonal is to carry the sign of the pures through the root, i.e.
            #     (b - a/RT)_ij = -sqrt[(a_ii/RT - b_ii)(a_jj/RT - b_jj)](1 - k_ij)
            # for the usual case. That is what is done here.
            if np.any(np.sign(Bi) != np.sign(Bi[0])) or np.any(Bi == 0):
                raise ValueError(
                    "Eq. 9.9-10a takes the square root of a product of the pure "
                    f"second virial coefficients {Bi}, which do not share a sign at "
                    f"T = {T:g} K -- the geometric mean is undefined. Use "
                    "combining='9.9-10b'.")
            M = np.sign(Bi[0]) * np.sqrt(np.outer(Bi, Bi)) * (1.0 - self.kij)
        else:
            M = (0.5 * np.add.outer(b, b)
                 - np.sqrt(np.outer(a, a)) / (R * T) * (1.0 - self.kij))
        np.fill_diagonal(M, Bi)                    # both rules give B_i exactly at i=j
        return M

    # -- the mixing rule, Eq. 9.9-9 ---------------------------------------
    def Q(self, x, T):
        """Q = SUM_i SUM_j x_i x_j (b_ij - a_ij/RT), Eq. 9.9-9b."""
        x = self._x(x)
        return float(x @ self.cross_matrix(T) @ x)

    def D(self, x, T):
        """D = SUM_i x_i a_i/(b_i RT) + G^ex(T,x)/(C* RT), Eq. 9.9-9c."""
        x = self._x(x)
        a, b = self.a_pure(T), self.b_pure
        return float(x @ (a / (b * R * T))
                     + self.gex_model.gex(x, T) / (self.C_star * R * T))

    def b_mix(self, x, T):
        """b = Q/(1 - D), Eq. 9.9-9a. ⚠️ Temperature dependent, unlike the van der
        Waals one-fluid b -- Q carries a/RT and D carries G^ex(T).

        ⚠️ **D > 1 is the normal case, not an error.** Reading Eq. 9.9-9a for the
        first time it is tempting to guard on `D < 1`, since `1 - D` is a
        denominator. But D is dominated by SUM_i x_i a_i/(b_i RT), which is roughly
        a/(bRT) -- about 12 for liquid water at 25 C, and greater than 1 for every
        fluid below its Boyle temperature. Q is negative there by the same token, and
        the two negatives divide to a positive b. The pure-component limit is the
        proof: b = (b_i - a_i/RT)/(1 - a_i/(b_i RT)) = b_i identically, and that
        identity runs through 1 - D < 0.

        What is actually fatal is D = 1 exactly, and a b that comes out non-positive.
        """
        D = self.D(x, T)
        if abs(1.0 - D) < 1e-12:
            raise ValueError(
                f"D = {D:.12g}: Eq. 9.9-9a is singular at D = 1 exactly.")
        b = self.Q(x, T) / (1.0 - D)
        if b <= 0.0:
            raise ValueError(
                f"the mixing rule returned b = {b:.4g} <= 0 (Q = {self.Q(x, T):.4g}, "
                f"D = {D:.4g}). Q and 1 - D must share a sign; if you are using "
                f"combining='9.9-10a', see the sign note in `cross_matrix`.")
        return b

    def a_mix(self, x, T):
        """a = RT Q D/(1 - D), Eq. 9.9-9a."""
        return R * T * self.b_mix(x, T) * self.D(x, T)

    # -- the two boundary conditions the rule was built from ---------------
    def second_virial(self, x, T):
        """B_mix = b - a/RT from the mixture parameters (Eq. 9.9-2)."""
        return self.b_mix(x, T) - self.a_mix(x, T) / (R * T)

    def excess_helmholtz(self, x, T):
        """A^ex_EOS at infinite pressure, Eq. 9.9-8: C*[a/b - SUM x_i a_i/b_i]."""
        x = self._x(x)
        a, b = self.a_pure(T), self.b_pure
        return float(self.C_star * (self.a_mix(x, T) / self.b_mix(x, T)
                                    - x @ (a / b)))

    def check_boundary(self, x, T):
        """Both conditions of Sec. 9.9, as relative residuals. Should be ~1e-16.

        Returns a dict:

        ``virial``   Eq. 9.9-3 -- does b - a/RT come out quadratic in composition?
        ``helmholtz`` Eqs. 9.9-4/9.9-7/9.9-8 -- does A^ex_EOS(P->inf) equal the
                      activity coefficient model's G^ex at low pressure?
        ``diagonal``  does the combining rule reduce to the pure B_i at i = j?

        ⭐ **This is the whole content of the section, run as arithmetic.** The mixing
        rules were *derived* by imposing these two, so their holding is not evidence
        that the derivation is right -- it is evidence that this implementation of it
        is. With no printed number anywhere in Sec. 9.9 to check against, it is the
        strongest internal test available.
        """
        x = self._x(x)
        M = self.cross_matrix(T)
        quad = float(x @ M @ x)
        lhs = self.second_virial(x, T)
        gex = self.gex_model.gex(x, T)
        aex = self.excess_helmholtz(x, T)
        a, b = self.a_pure(T), self.b_pure
        scale = lambda v, s: abs(v) / max(abs(s), 1e-300)
        return {
            "virial": scale(lhs - quad, quad),
            "helmholtz": scale(aex - gex, gex) if gex != 0 else abs(aex),
            "diagonal": float(np.max(np.abs(np.diag(M) - (b - a / (R * T)))
                                     / np.abs(b - a / (R * T)))),
        }

    # -- the equation of state --------------------------------------------
    def _AB(self, x, T, P):
        return (self.a_mix(x, T) * P / (R * T) ** 2,
                self.b_mix(x, T) * P / (R * T))

    def compressibility(self, x, T, P):
        """All real roots Z of the PR cubic, ascending."""
        A, B = self._AB(x, T, P)
        coeffs = [1.0, -1 + B, A - 3 * B ** 2 - 2 * B, -A * B + B ** 2 + B ** 3]
        roots = np.roots(coeffs)
        return np.sort(roots.real[np.abs(roots.imag) < 1e-9])

    def Z(self, x, T, P, phase="vapor"):
        zs = self.compressibility(x, T, P)
        zs = zs[zs > self._AB(x, T, P)[1]]        # Z > B, or the log terms are complex
        if len(zs) == 0:
            raise ValueError("no physical root (Z > B) at these conditions")
        return float(zs.max() if phase == "vapor" else zs.min())

    def molar_volume(self, x, T, P, phase="vapor"):
        return self.Z(x, T, P, phase) * R * T / P

    # -- composition derivatives, Eqs. 9.9-12 and 9.9-13 -------------------
    def dND_dNi(self, x, T):
        """d(ND)/dN_i = a_i/(b_i RT) + ln gamma_i / C*, Eq. 9.9-13.

        ⭐ The second term is the neat part of the whole method: differentiating
        N G^ex with respect to N_i gives RT ln gamma_i by definition, so the activity
        coefficient model enters the fugacity coefficient *as itself* -- no numerical
        differentiation, and any model with a gamma can be dropped in.
        """
        x = self._x(x)
        a, b = self.a_pure(T), self.b_pure
        return (a / (b * R * T)
                + np.log(self.gex_model.gamma(x, T)) / self.C_star)

    def dN2Q_dNi(self, x, T):
        """(1/N) d(N^2 Q)/dN_i = 2 SUM_j x_j (b - a/RT)_ij, Eq. 9.9-12 (third line)."""
        x = self._x(x)
        return 2.0 * (self.cross_matrix(T) @ x)

    def dNb_dNi(self, x, T):
        """d(Nb)/dN_i, Eq. 9.9-12 (first line).

        Nb = (N^2 Q)/(N - ND), so the quotient rule gives

            d(Nb)/dN_i = [d(N^2Q)/dN_i]/(1-D) - Q[1 - d(ND)/dN_i]/(1-D)^2
        """
        D = self.D(x, T)
        return (self.dN2Q_dNi(x, T) / (1.0 - D)
                - self.Q(x, T) * (1.0 - self.dND_dNi(x, T)) / (1.0 - D) ** 2)

    def dN2a_dNi(self, x, T):
        """(1/N) d(N^2 a)/dN_i = RT[D d(Nb)/dN_i + b d(ND)/dN_i], Eq. 9.9-12 (second).

        Follows from a = RT b D, i.e. N^2 a = RT (Nb)(ND) -- a product rule, which is
        why the same two derivatives appear again.
        """
        return R * T * (self.D(x, T) * self.dNb_dNi(x, T)
                        + self.b_mix(x, T) * self.dND_dNi(x, T))

    def ln_phi(self, x, T, P, phase="vapor"):
        """ln of the fugacity coefficient of each species, Eq. 9.9-11.

        The general cubic-mixture expression, written with composition *derivatives*
        rather than with b_i and 2 SUM_j x_j a_ij: for the van der Waals one-fluid
        rules those are the same thing, and `code/ch9/` checks this against
        `PRMixture` on that case.
        """
        x = self._x(x)
        a, b = self.a_mix(x, T), self.b_mix(x, T)
        B = b * P / (R * T)
        Z = self.Z(x, T, P, phase)
        nb, n2a = self.dNb_dNi(x, T), self.dN2a_dNi(x, T)
        log_term = np.log((Z + (1 + _SQRT2) * B) / (Z + (1 - _SQRT2) * B))
        return (nb / b * (Z - 1) - np.log(Z - B)
                - a / (2 * _SQRT2 * b * R * T) * (n2a / a - nb / b) * log_term)

    def phi(self, x, T, P, phase="vapor"):
        return np.exp(self.ln_phi(x, T, P, phase))

    def fugacity(self, x, T, P, phase="vapor"):
        """Component fugacities f_i = x_i phi_i P (Pa)."""
        return self._x(x) * self.phi(x, T, P, phase) * P

    # -- helpers -----------------------------------------------------------
    def _x(self, x):
        x = np.asarray(x, dtype=float)
        if x.size != self.n:
            raise ValueError(f"expected {self.n} mole fractions, got {x.size}")
        if not np.isclose(x.sum(), 1.0, atol=1e-8):
            raise ValueError(f"mole fractions must sum to 1, got {x.sum():.6g}")
        return x

    def __repr__(self):
        names = ", ".join(getattr(c, "name", "?") for c in self.components)
        return (f"WongSandler([{names}], {type(self.gex_model).__name__}, "
                f"combining={self.combining!r})")
