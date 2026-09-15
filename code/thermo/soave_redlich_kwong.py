"""Pure-fluid Soave-Redlich-Kwong equation of state.

The book introduces it beside Peng-Robinson as the other generalized cubic (SIS
Sec. 6.7, Table 6.4-3), and chapter 6 sets eight problems on it -- 6.35 through
6.41 and 6.48 -- while nothing in this package could evaluate it.

Same shape as `PengRobinson`, one term different. The attraction is divided by
V(V + b) rather than by V(V + b) + b(V - b), which makes the departure integrals
elementary and their logarithm ln(1 + B/Z) rather than the Peng-Robinson ratio:

    P = RT/(V - b) - a(T)/(V(V + b))
    a(T) = 0.42748 R^2 Tc^2 / Pc * alpha(T)
    alpha(T) = [1 + m (1 - sqrt(T/Tc))]^2,  m = 0.480 + 1.574 w - 0.176 w^2
    b = 0.08664 R Tc / Pc

CAUTION: `m` HERE IS `kappa` THERE, and the two correlations are different functions of
the acentric factor -- 0.480 + 1.574w - 0.176w^2 against Peng-Robinson's
0.37464 + 1.54226w - 0.26992w^2. For oxygen (w = 0.021) they are 0.5130 and
0.4069, which is not a small difference and is not interchangeable.

VALIDATED AGAINST A NUMBER THIS PACKAGE DID NOT PRODUCE. The 5e solution manual
solved Problem 6.37 -- Illustration 6.4-1 repeated with this equation -- inside a
Mathcad worksheet, and the worksheet prints its result for oxygen at 173.15 K and
1 bar: Z = 0.9952, V = 0.014327 m^3/mol, H = -3602.73 J/mol, S = -15.62 J/(mol K),
b = 2.207e-5 m^3/mol. This class reproduces all five to every printed digit when
fed that worksheet's own ideal-gas Cp and its R = 8.31451. The check is in
`solution-manual/notebooks/ch06/6.48.ipynb`.

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


class SoaveRedlichKwong(CubicEOS):
    """Pure-component Soave-Redlich-Kwong EOS.

    Parameters
    ----------
    Tc, Pc : critical temperature (K) and pressure (Pa)
    omega  : Pitzer acentric factor
    name   : optional label
    cp     : optional ideal-gas Cp polynomial coefficients (a, b, c, d) for
             Cp* = a + b T + c T^2 + d T^3 in J/(mol K) (Appendix A.II form)
    """

    def __init__(self, Tc, Pc, omega, name=None, cp=None):
        self.Tc = float(Tc)
        self.Pc = float(Pc)
        self.omega = float(omega)
        self.name = name
        self.cp = tuple(cp) if cp is not None else None
        # SIS Sec. 6.7. Called `m` in the Soave literature and in the 5e's own
        # worksheets; the Peng-Robinson analogue is `kappa` and is a DIFFERENT
        # function of omega.
        self.m = 0.480 + 1.574 * self.omega - 0.176 * self.omega ** 2
        self.b = 0.08664 * R * self.Tc / self.Pc
        self._a_c = 0.42748 * R ** 2 * self.Tc ** 2 / self.Pc

    def __repr__(self):
        return (f"<SoaveRedlichKwong {self.name or '?'}: Tc={self.Tc} K, "
                f"Pc={self.Pc/1e5:.4g} bar, omega={self.omega}>")

    @classmethod
    def from_database(cls, key, cp=None):
        """Build from `pure_property.csv` (Pc there is in bar -> converted to Pa).

        `cp` overrides the ideal-gas heat capacity; pass
        `thermo.data.APPENDIX_A2_CP[key]` for the book's Appendix A.II set, which
        is what the printed tables are computed with.
        """
        c = get_compound(key)
        if cp is None:
            cp = (float(c.CpA), float(c.CpB), float(c.CpC), float(c.CpD))
        return cls(Tc=float(c.Tc), Pc=float(c.Pc) * 1e5, omega=float(c.Omega),
                   name=str(c.Name), cp=cp)

    # --- EOS parameters --------------------------------------------------
    def _sqrt_alpha(self, T):
        return 1 + self.m * (1 - np.sqrt(T / self.Tc))

    def a(self, T):
        return self._a_c * self._sqrt_alpha(T) ** 2

    def dadT(self, T):
        # d/dT [sqrt_alpha^2] = 2 sqrt_alpha * d(sqrt_alpha)/dT, and
        # d(sqrt_alpha)/dT = -m / (2 sqrt(T Tc))
        return -self._a_c * self.m * self._sqrt_alpha(T) / np.sqrt(T * self.Tc)

    def pressure(self, V, T):
        """Pressure (Pa) from molar volume V (m^3/mol) and T."""
        return R * T / (V - self.b) - self.a(T) / (V * (V + self.b))

    def _AB(self, T, P):
        return self.a(T) * P / (R * T) ** 2, self.b * P / (R * T)

    # --- roots -----------------------------------------------------------
    def compressibility(self, T, P):
        """All real roots Z of the SRK cubic, ascending.

            Z^3 - Z^2 + (A - B - B^2) Z - A B = 0
        """
        A, B = self._AB(T, P)
        return real_roots(Polynomial([-A * B, A - B - B ** 2, -1.0, 1.0]).roots())

    # Z, molar_volume, fugacity, spinodal_bounds and vapor_pressure come from
    # CubicEOS -- they are the same for every pure-fluid cubic.

    # --- fugacity and departures -----------------------------------------
    def _log_term(self, Z, B):
        """ln(1 + B/Z). The SRK counterpart of Peng-Robinson's root-2 ratio."""
        return np.log(1 + B / Z)

    def ln_phi(self, T, P, phase="vapor"):
        """ln of the fugacity coefficient."""
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return Z - 1 - np.log(Z - B) - A / B * self._log_term(Z, B)

    def departure_H(self, T, P, phase="vapor"):
        """(H - H_ideal-gas) at (T, P), J/mol."""
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return (R * T * (Z - 1)
                + (T * self.dadT(T) - self.a(T)) / self.b * self._log_term(Z, B))

    def departure_S(self, T, P, phase="vapor"):
        """(S - S_ideal-gas) at (T, P), J/(mol K)."""
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return R * np.log(Z - B) + self.dadT(T) / self.b * self._log_term(Z, B)
