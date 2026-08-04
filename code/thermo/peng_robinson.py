"""Pure-fluid Peng-Robinson equation of state, with the PRSV option.

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
from scipy import constants

from .cubic import CubicEOS, real_roots
from .data import get_compound

R = constants.R
_SQRT2 = np.sqrt(2.0)


class PengRobinson(CubicEOS):
    """Pure-component Peng-Robinson EOS, optionally in its PRSV form.

    Parameters
    ----------
    Tc, Pc : critical temperature (K) and pressure (Pa)
    omega  : Pitzer acentric factor
    name   : optional label
    cp     : optional ideal-gas Cp polynomial coefficients (a, b, c, d) for
             Cp* = a + b T + c T^2 + d T^3 in J/(mol K) (Appendix A.II form)
    kappa1 : optional PRSV parameter (SIS Eq. 7.5-1). Leave it None for the
             standard Peng-Robinson equation, where kappa comes from the
             generalized correlation in omega and is temperature independent.
             Supply it (0.0 counts) to switch to PRSV, where

                 kappa(T) = kappa0 + kappa1 (1 + sqrt(Tr)) (0.7 - Tr)

             with kappa0 the Stryjek-Vera cubic in omega (SIS Eq. 7.5-2).
             PRSV with kappa1 = 0 is NOT the standard Peng-Robinson equation:
             kappa0 and the PR kappa differ (for water, by 0.16%).

    Notes
    -----
    Values of kappa1 are substance specific and fitted to vapor-pressure data;
    the book gives -0.0665 for water (Illustration 7.5-3).
    """

    def __init__(self, Tc, Pc, omega, name=None, cp=None, kappa1=None):
        self.Tc = float(Tc)
        self.Pc = float(Pc)
        self.omega = float(omega)
        self.name = name
        self.cp = tuple(cp) if cp is not None else None
        self.kappa1 = None if kappa1 is None else float(kappa1)
        if self.kappa1 is None:
            # standard PR: kappa is a constant (SIS Eq. 6.7-4)
            self.kappa0 = (0.37464 + 1.54226 * self.omega
                           - 0.26992 * self.omega ** 2)
        else:
            # PRSV: Stryjek and Vera, Can. J. Chem. Eng. 64, 323 (1986).
            # NOTE THE MINUS on the omega^2 term. The 5e prints a plus in
            # Eq. 7.5-2, which is an erratum -- it does not reproduce the book's
            # own PRSV column (see revision_notes/c07_manuscript_edits_7.5-3.md).
            self.kappa0 = (0.378893 + 1.4897153 * self.omega
                           - 0.17131848 * self.omega ** 2
                           + 0.0196554 * self.omega ** 3)
        # `kappa` is the temperature-independent part. For standard PR it IS
        # kappa; for PRSV it is kappa0. Kept as an attribute because the ch6
        # `_thermo` notebooks print it.
        self.kappa = self.kappa0
        self.b = 0.07780 * R * self.Tc / self.Pc

    @property
    def is_prsv(self):
        return self.kappa1 is not None

    def __repr__(self):
        form = f"PRSV, kappa1={self.kappa1}" if self.is_prsv else "PR"
        return (f"<{type(self).__name__} {self.name or '?'}: Tc={self.Tc} K, "
                f"Pc={self.Pc/1e5:.4g} bar, omega={self.omega} ({form})>")

    @classmethod
    def from_database(cls, key, kappa1=None, cp=None):
        """Build from `pure_property.csv` (Pc there is in bar -> converted to Pa).

        `cp` overrides the ideal-gas heat capacity. The database carries the
        Reid-Prausnitz-Poling coefficients; pass `thermo.data.APPENDIX_A2_CP[key]`
        to use the book's Appendix A.II set instead, which is what the printed
        tables of Illustrations 6.4-1, 7.5-1 and 7.5-2 are computed with.
        """
        c = get_compound(key)
        if cp is None:
            cp = (float(c.CpA), float(c.CpB), float(c.CpC), float(c.CpD))
        return cls(Tc=float(c.Tc), Pc=float(c.Pc) * 1e5, omega=float(c.Omega),
                   name=str(c.Name), cp=cp, kappa1=kappa1)

    # --- EOS parameters --------------------------------------------------
    def kappa_T(self, T):
        """kappa at temperature T. Constant unless this is a PRSV fluid."""
        if self.kappa1 is None:
            return self.kappa0
        Tr = T / self.Tc
        return self.kappa0 + self.kappa1 * (1 + np.sqrt(Tr)) * (0.7 - Tr)

    def _dkappa_dT(self, T):
        if self.kappa1 is None:
            return 0.0
        Tr = T / self.Tc
        return self.kappa1 * ((0.7 - Tr) / (2 * np.sqrt(T * self.Tc))
                              - (1 + np.sqrt(Tr)) / self.Tc)

    def _sqrt_alpha(self, T):
        return 1 + self.kappa_T(T) * (1 - np.sqrt(T / self.Tc))

    def a(self, T):
        return 0.45724 * R ** 2 * self.Tc ** 2 / self.Pc * self._sqrt_alpha(T) ** 2

    def dadT(self, T):
        if self.kappa1 is None:
            return (-0.45724 * R ** 2 * self.Tc ** 2 / self.Pc
                    * self.kappa0 * self._sqrt_alpha(T) / np.sqrt(T * self.Tc))
        # PRSV: kappa itself depends on T, so d(sqrt_alpha)/dT picks up a term.
        # sqrt_alpha = 1 + kappa(T) (1 - sqrt(Tr))
        #   d/dT = kappa'(T) (1 - sqrt(Tr)) - kappa(T) / (2 sqrt(T Tc))
        d_sqrt_alpha = (self._dkappa_dT(T) * (1 - np.sqrt(T / self.Tc))
                        - self.kappa_T(T) / (2 * np.sqrt(T * self.Tc)))
        return (2 * 0.45724 * R ** 2 * self.Tc ** 2 / self.Pc
                * self._sqrt_alpha(T) * d_sqrt_alpha)

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
        return real_roots(Polynomial([gamma, beta, alpha, 1]).roots())

    # Z, molar_volume, fugacity, spinodal_bounds and vapor_pressure come from
    # CubicEOS -- they are the same for every pure-fluid cubic.

    # --- fugacity and departures ----------------------------------------
    def _log_term(self, Z, B):
        return np.log((Z + (1 + _SQRT2) * B) / (Z + (1 - _SQRT2) * B))

    def ln_phi(self, T, P, phase="vapor"):
        """ln of the fugacity coefficient (standard PR)."""
        A, B = self._AB(T, P)
        Z = self.Z(T, P, phase)
        return Z - 1 - np.log(Z - B) - A / (2 * _SQRT2 * B) * self._log_term(Z, B)

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
