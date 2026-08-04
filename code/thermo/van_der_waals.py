"""Pure-fluid van der Waals equation of state.

The book's first cubic equation of state (SIS Eq. 6.4-1 and Sec. 7.3), refactored
from `code/ch7/vapor_pressure_n_butane.ipynb` -- the notebook behind Figure 7.5-2,
where the van der Waals vapor pressure is the control case that shows what the
Peng-Robinson temperature dependence alpha(T) buys.

The API mirrors `PengRobinson` so that one calculation can be handed either
equation. In particular `a(T)` is a method here too, even though the van der Waals
a is a constant, and `dadT(T)` returns zero -- which is why the entropy departure
below is simply R ln(Z - B) and the enthalpy departure carries no dadT term.

SI units throughout: T in K, P in Pa, V in m^3/mol.

Eric M. Furst
August 2026
"""
import numpy as np
from numpy.polynomial import Polynomial
from scipy import constants

from .cubic import CubicEOS, real_roots
from .data import get_compound

R = constants.R


class VanDerWaals(CubicEOS):
    """Pure-component van der Waals EOS.

    Parameters
    ----------
    Tc, Pc : critical temperature (K) and pressure (Pa)
    omega  : accepted and stored for API symmetry with `PengRobinson`; the van
             der Waals equation does not use the acentric factor, which is the
             whole of its trouble with vapor pressure (Fig. 7.5-2, curve a)
    name   : optional label
    cp     : optional ideal-gas Cp coefficients (a, b, c, d), J/(mol K)
    """

    def __init__(self, Tc, Pc, omega=None, name=None, cp=None):
        self.Tc = float(Tc)
        self.Pc = float(Pc)
        self.omega = None if omega is None else float(omega)
        self.name = name
        self.cp = tuple(cp) if cp is not None else None
        # from the critical-point constraints (dP/dV)_Tc = (d2P/dV2)_Tc = 0
        self._a = 27 * R ** 2 * self.Tc ** 2 / (64 * self.Pc)
        self.b = R * self.Tc / (8 * self.Pc)

    def __repr__(self):
        return (f"<VanDerWaals {self.name or '?'}: Tc={self.Tc} K, "
                f"Pc={self.Pc/1e5:.4g} bar>")

    @classmethod
    def from_database(cls, key):
        """Build from `pure_property.csv` (Pc there is in bar -> converted to Pa)."""
        c = get_compound(key)
        cp = (float(c.CpA), float(c.CpB), float(c.CpC), float(c.CpD))
        return cls(Tc=float(c.Tc), Pc=float(c.Pc) * 1e5, omega=float(c.Omega),
                   name=str(c.Name), cp=cp)

    # --- EOS parameters --------------------------------------------------
    def a(self, T=None):
        """The van der Waals a, which does not depend on temperature."""
        return self._a

    def dadT(self, T=None):
        return 0.0

    def pressure(self, V, T):
        """Pressure (Pa) from molar volume V (m^3/mol) and T."""
        return R * T / (V - self.b) - self._a / V ** 2

    def _AB(self, T, P):
        return self._a * P / (R * T) ** 2, self.b * P / (R * T)

    # --- roots -----------------------------------------------------------
    def compressibility(self, T, P):
        """All real roots Z of Z^3 - (1 + B) Z^2 + A Z - A B = 0, ascending."""
        A, B = self._AB(T, P)
        return real_roots(Polynomial([-A * B, A, -(1 + B), 1]).roots())

    # --- fugacity and departures ----------------------------------------
    def ln_phi(self, T, P, phase="vapor"):
        """ln of the fugacity coefficient."""
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return Z - 1 - np.log(Z - B) - A / Z

    def departure_H(self, T, P, phase="vapor"):
        """(H - H_ideal-gas) at (T, P), J/mol.

        RT(Z - 1) - a/V. There is no dadT term: for van der Waals, a is constant.
        """
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return R * T * (Z - 1) - self._a * P / (Z * R * T)

    def departure_S(self, T, P, phase="vapor"):
        """(S - S_ideal-gas) at (T, P), J/(mol K).

        R ln(Z - B) = R ln[(V - b) P / RT], again with no dadT term.
        """
        A, B = self._AB(T, P)
        return R * np.log(self.Z(T, P, phase) - B)
