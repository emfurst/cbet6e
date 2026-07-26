"""Pure-fluid Peng-Robinson equation of state.

Refactored from the validated Chapter 6/7 notebooks and the book (SIS Eqs. 6.4-2,
6.4-29, 6.4-30, and Table 6.4-3). An Aspen-optional substitute for pure-fluid
property calculations.

SI units throughout: temperature in K, pressure in Pa, molar volume in m^3/mol,
energies in J/mol, entropy in J/(mol K). R = 8.314 J/(mol K).

Eric M. Furst
July 2026
"""
import numpy as np
from numpy.polynomial import Polynomial
from scipy import constants, optimize

from .data import get_compound

R = constants.R
_SQRT2 = np.sqrt(2.0)


class PengRobinson:
    """Pure-component Peng-Robinson EOS.

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
        self.kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega ** 2
        self.b = 0.07780 * R * self.Tc / self.Pc

    @classmethod
    def from_database(cls, key):
        """Build from `pure_property.csv` (Pc there is in bar -> converted to Pa)."""
        c = get_compound(key)
        cp = (float(c.CpA), float(c.CpB), float(c.CpC), float(c.CpD))
        return cls(Tc=float(c.Tc), Pc=float(c.Pc) * 1e5, omega=float(c.Omega),
                   name=str(c.Name), cp=cp)

    # --- EOS parameters --------------------------------------------------
    def _sqrt_alpha(self, T):
        return 1 + self.kappa * (1 - np.sqrt(T / self.Tc))

    def a(self, T):
        return 0.45724 * R ** 2 * self.Tc ** 2 / self.Pc * self._sqrt_alpha(T) ** 2

    def dadT(self, T):
        return (-0.45724 * R ** 2 * self.Tc ** 2 / self.Pc
                * self.kappa * self._sqrt_alpha(T) / np.sqrt(T * self.Tc))

    def pressure(self, V, T):
        """Pressure (Pa) from molar volume V (m^3/mol) and T (SIS Eq. 6.4-2)."""
        a, b = self.a(T), self.b
        return R * T / (V - b) - a / (V * (V + b) + b * (V - b))

    def _AB(self, T, P):
        return self.a(T) * P / (R * T) ** 2, self.b * P / (R * T)

    # --- roots -----------------------------------------------------------
    def compressibility(self, T, P):
        """All real roots Z of the PR cubic (SIS Table 6.4-3), ascending."""
        A, B = self._AB(T, P)
        alpha = -1 + B
        beta = A - 3 * B ** 2 - 2 * B
        gamma = -A * B + B ** 2 + B ** 3
        roots = Polynomial([gamma, beta, alpha, 1]).roots()
        return np.sort(roots.real[np.abs(roots.imag) < 1e-9])

    def Z(self, T, P, phase="vapor"):
        """Vapor (largest) or liquid (smallest) real compressibility root."""
        zs = self.compressibility(T, P)
        if len(zs) == 0:
            raise ValueError("no real root")
        return zs.max() if phase == "vapor" else zs.min()

    def molar_volume(self, T, P, phase="vapor"):
        return self.Z(T, P, phase) * R * T / P

    # --- fugacity and departures ----------------------------------------
    def _log_term(self, Z, B):
        return np.log((Z + (1 + _SQRT2) * B) / (Z + (1 - _SQRT2) * B))

    def ln_phi(self, T, P, phase="vapor"):
        """ln of the fugacity coefficient (standard PR)."""
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return Z - 1 - np.log(Z - B) - A / (2 * _SQRT2 * B) * self._log_term(Z, B)

    def fugacity(self, T, P, phase="vapor"):
        return P * np.exp(self.ln_phi(T, P, phase))

    def departure_H(self, T, P, phase="vapor"):
        """(H - H_ideal-gas) at (T, P), J/mol (SIS Eq. 6.4-29)."""
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return (R * T * (Z - 1)
                + (T * self.dadT(T) - self.a(T)) / (2 * _SQRT2 * self.b)
                * self._log_term(Z, B))

    def departure_S(self, T, P, phase="vapor"):
        """(S - S_ideal-gas) at (T, P), J/(mol K) (SIS Eq. 6.4-30).

        Note: uses R*ln(Z - B); the ch6 throttle notebook's calc_depS instead
        wrote R*(Z - B), which appears to be a typo (flagged for the editor).
        """
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return R * np.log(Z - B) + self.dadT(T) / (2 * _SQRT2 * self.b) * self._log_term(Z, B)

    # --- saturation ------------------------------------------------------
    def vapor_pressure(self, T, P_guess=None):
        """Saturation pressure (Pa) at T by equating liquid and vapor fugacities."""
        if T >= self.Tc:
            raise ValueError("T must be below Tc")
        if P_guess is None:  # Pitzer/Lee-Kesler-style initial estimate
            Tr = T / self.Tc
            P_guess = self.Pc * np.exp(5.373 * (1 + self.omega) * (1 - 1 / Tr))

        def resid(lnP):
            P = np.exp(lnP)
            return self.ln_phi(T, P, "liquid") - self.ln_phi(T, P, "vapor")

        lnP = optimize.newton(resid, np.log(P_guess), tol=1e-10, maxiter=200)
        return float(np.exp(lnP))
