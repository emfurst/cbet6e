"""Pure-fluid Redlich-Kwong equation of state -- the 1949 original, not Soave's.

    P = RT/(V - b) - a/(sqrt(T) V (V + b))
    a = 0.42748 R^2 Tc^2.5 / Pc        b = 0.08664 R Tc / Pc

THIS IS A DIFFERENT EQUATION FROM `SoaveRedlichKwong`, and confusing the two is
easy because they share a name, a `b`, and the same 0.42748. The difference is the
temperature dependence of the attraction:

    Redlich-Kwong        a/sqrt(T),  fixed once Tc and Pc are known
    Soave-Redlich-Kwong  a(T) = a_c alpha(T),  alpha fitted through the acentric
                         factor so that the equation reproduces vapor pressures

So Redlich-Kwong needs NO acentric factor -- it is a two-constant equation in the
corresponding-states sense -- and it is correspondingly poor at vapor pressure,
which is why Soave modified it. Chapter 6 uses the original in
Problems 6.42 through 6.47, and the book's Table 6.4-3 carries its cubic form.

THE SECOND VIRIAL COEFFICIENT IS ELEMENTARY HERE, which is Problem 6.45:

    B(T) = b - a/(R T^1.5)

and setting it to zero gives the Boyle temperature, Problem 6.43:

    T_Boyle = (a/(R b))^(2/3) = (0.42748/0.08664)^(2/3) Tc = 2.898 Tc

Both are exact consequences of the equation, so `second_virial` and
`boyle_temperature` are checks on the constants that need no data at all.

SI units throughout: T in K, P in Pa, V in m^3/mol, energies in J/mol, entropy in
J/(mol K).

Eric M. Furst
September 2026
"""
import numpy as np
from numpy.polynomial import Polynomial
from scipy import constants

from .cubic import CubicEOS, real_roots
from .data import get_compound

R = constants.R


class RedlichKwong(CubicEOS):
    """Pure-component Redlich-Kwong EOS (1949).

    Parameters
    ----------
    Tc, Pc : critical temperature (K) and pressure (Pa)
    name   : optional label
    cp     : optional ideal-gas Cp polynomial coefficients (a, b, c, d) for
             Cp* = a + b T + c T^2 + d T^3 in J/(mol K) (Appendix A.II form)

    There is deliberately no `omega`: the equation does not use one.
    """

    def __init__(self, Tc, Pc, name=None, cp=None):
        self.Tc = float(Tc)
        self.Pc = float(Pc)
        self.name = name
        self.cp = tuple(cp) if cp is not None else None
        self.b = 0.08664 * R * self.Tc / self.Pc
        # a carries sqrt(K) -- the Tc^2.5 is not a typo for Tc^2.
        self.a_rk = 0.42748 * R ** 2 * self.Tc ** 2.5 / self.Pc

    def __repr__(self):
        return (f"<RedlichKwong {self.name or '?'}: Tc={self.Tc} K, "
                f"Pc={self.Pc/1e5:.4g} bar>")

    @classmethod
    def from_database(cls, key, cp=None):
        """Build from `pure_property.csv` (Pc there is in bar -> converted to Pa)."""
        c = get_compound(key)
        if cp is None:
            cp = (float(c.CpA), float(c.CpB), float(c.CpC), float(c.CpD))
        return cls(Tc=float(c.Tc), Pc=float(c.Pc) * 1e5, name=str(c.Name), cp=cp)

    # --- EOS parameters --------------------------------------------------
    def a(self, T):
        """The attraction term's coefficient AS IT ENTERS THE CUBIC, a_rk/sqrt(T).

        Named to match `PengRobinson.a(T)` so that anything written against the
        one works against the other; `a_rk` is the temperature-independent
        constant underneath it.
        """
        return self.a_rk / np.sqrt(T)

    def dadT(self, T):
        return -0.5 * self.a_rk / T ** 1.5

    def pressure(self, V, T):
        """Pressure (Pa) from molar volume V (m^3/mol) and T."""
        return R * T / (V - self.b) - self.a(T) / (V * (V + self.b))

    def _AB(self, T, P):
        return self.a(T) * P / (R * T) ** 2, self.b * P / (R * T)

    # --- roots -----------------------------------------------------------
    def compressibility(self, T, P):
        """All real roots Z of the RK cubic, ascending (SIS Table 6.4-3).

            Z^3 - Z^2 + (A - B - B^2) Z - A B = 0
        """
        A, B = self._AB(T, P)
        return real_roots(Polynomial([-A * B, A - B - B ** 2, -1.0, 1.0]).roots())

    # Z, molar_volume, fugacity, spinodal_bounds and vapor_pressure come from
    # CubicEOS -- they are the same for every pure-fluid cubic.

    # --- exact consequences of the equation ------------------------------
    def second_virial(self, T):
        """B(T) = b - a/(R T^1.5), m^3/mol. Problem 6.45."""
        return self.b - self.a_rk / (R * T ** 1.5)

    def boyle_temperature(self):
        """The temperature at which B(T) = 0, K. Problem 6.43.

        (a/(R b))^(2/3) = (0.42748/0.08664)^(2/3) Tc = 2.8980 Tc, so it is a
        fixed multiple of the critical temperature for every fluid -- which is
        the corresponding-states content of a two-constant equation.
        """
        return (self.a_rk / (R * self.b)) ** (2.0 / 3.0)

    # --- fugacity and departures -----------------------------------------
    def _log_term(self, Z, B):
        return np.log(1 + B / Z)

    def ln_phi(self, T, P, phase="vapor"):
        """ln of the fugacity coefficient."""
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return Z - 1 - np.log(Z - B) - A / B * self._log_term(Z, B)

    def departure_H(self, T, P, phase="vapor"):
        """(H - H_ideal-gas) at (T, P), J/mol.

        With a(T) = a_rk/sqrt(T), T da/dT - a = -1.5 a(T), so the bracket that
        Soave-Redlich-Kwong carries as a general expression collapses here.
        """
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return (R * T * (Z - 1)
                - 1.5 * self.a(T) / self.b * self._log_term(Z, B))

    def departure_S(self, T, P, phase="vapor"):
        """(S - S_ideal-gas) at (T, P), J/(mol K)."""
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return R * np.log(Z - B) + self.dadT(T) / self.b * self._log_term(Z, B)
