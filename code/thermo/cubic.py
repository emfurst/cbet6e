"""Shared machinery for the pure-fluid cubic equations of state.

A subclass supplies the equation itself -- the attributes `Tc`, `Pc`, `b`, and the
methods `pressure(V, T)`, `compressibility(T, P)`, and `ln_phi(T, P, phase)` -- and
inherits everything here: root selection, molar volume, fugacity, the spinodal
pressures, and the saturation solver.

`vapor_pressure` is the algorithm of SIS Figure 7.5-1, in the form the figure is
drawn for the 6th edition: the equal-fugacity root is **bracketed between the two
spinodal pressures** (the local maximum and minimum of the sub-critical isotherm)
rather than reached from a guess. Any pressure strictly inside that interval gives
two distinct roots of the cubic, so the trivial solution Sec. 7.5 warns about --
guessing outside the van der Waals loop, getting one root twice, and accepting the
guess as the answer -- cannot occur. It is also the only bracket that always works:
the fugacity residual changes sign exactly once across the loop.

SI units throughout: T in K, P in Pa, V in m^3/mol.

Eric M. Furst
August 2026
"""
import numpy as np
from scipy import constants, optimize

R = constants.R


def real_roots(roots):
    """The real roots of a cubic, ascending.

    The tolerance is scaled to the root magnitude: near the critical point the
    three real roots crowd together and pick up tiny numerical imaginary parts
    that a fixed tolerance would discard.
    """
    keep = np.abs(roots.imag) < 1e-9 * np.maximum(1.0, np.abs(roots.real))
    return np.sort(roots.real[keep])


class CubicEOS:
    """Base class: what every pure-fluid cubic equation of state has in common."""

    # --- roots -----------------------------------------------------------
    def _B(self, T, P):
        return self.b * P / (R * T)

    def physical_Z(self, T, P):
        """Real compressibility roots with Z > B, ascending.

        A root at or below B corresponds to a molar volume smaller than the
        excluded volume b; it is a root of the cubic but not a state of the fluid.
        """
        zs = self.compressibility(T, P)
        return zs[zs > self._B(T, P)]

    def Z(self, T, P, phase="vapor"):
        """Vapor (largest) or liquid (smallest) physical compressibility root."""
        zs = self.physical_Z(T, P)
        if len(zs) == 0:
            raise ValueError(f"no physical root at T = {T} K, P = {P} Pa")
        return zs.max() if phase == "vapor" else zs.min()

    def molar_volume(self, T, P, phase="vapor"):
        return self.Z(T, P, phase) * R * T / P

    def fugacity(self, T, P, phase="vapor"):
        return P * np.exp(self.ln_phi(T, P, phase))

    # --- saturation ------------------------------------------------------
    def spinodal_bounds(self, T, n=20000):
        """The pressures at the two turning points of the isotherm at T.

        These are the limits of mechanical stability, (dP/dV)_T = 0, and they
        bracket the vapor pressure. Returns None once the loop has vanished
        (T >= Tc), which is the "does a two-phase solution exist?" test.
        """
        V = np.logspace(np.log10(1.0000001 * self.b), 3.0, n)
        P = self.pressure(V, T)
        turns = np.where(np.diff(np.sign(np.diff(P))) != 0)[0]
        if len(turns) < 2:
            return None
        return P[turns].min(), P[turns].max()

    def ln_fugacity_ratio(self, T, P):
        """ln(f_L / f_V) at T, P. NaN where the cubic gives no two distinct roots.

        This is the quantity the Fig. 7.5-1 convergence test is written on: it is
        zero at the vapor pressure, positive below it and negative above.
        """
        if len(self.physical_Z(T, P)) < 2:
            return np.nan
        return self.ln_phi(T, P, "liquid") - self.ln_phi(T, P, "vapor")

    def vapor_pressure(self, T, P_guess=None):
        """Saturation pressure (Pa) at T by equal fugacity (SIS Fig. 7.5-1).

        The root is bracketed between the spinodal pressures and found by Brent's
        method, so no initial guess is needed; `P_guess` is accepted and ignored,
        and remains only so that older call sites keep working.
        """
        if T >= self.Tc:
            raise ValueError(f"T = {T} K is at or above Tc = {self.Tc} K")
        bounds = self.spinodal_bounds(T)
        if bounds is None:
            raise ValueError(f"no van der Waals loop at T = {T} K: "
                             "the equation of state has no two-phase solution here")
        P_lo, P_hi = bounds
        P_lo = max(P_lo, 1e-12 * self.Pc)   # the vapor spinodal pressure can be negative

        # Walk the ends inward until both give a usable residual of opposite sign.
        # The spinodal pressures themselves are the degenerate limits, where two of
        # the three roots merge, so the residual is evaluated just inside them.
        for _ in range(200):
            r_lo = self.ln_fugacity_ratio(T, P_lo)
            r_hi = self.ln_fugacity_ratio(T, P_hi)
            if np.isfinite(r_lo) and np.isfinite(r_hi) and r_lo * r_hi < 0:
                break
            if not np.isfinite(r_lo):
                P_lo *= 1.5
            if not np.isfinite(r_hi):
                P_hi *= 0.999
            if P_hi <= P_lo:
                raise ValueError(f"could not bracket the vapor pressure at T = {T} K")
        else:
            raise ValueError(f"could not bracket the vapor pressure at T = {T} K")

        return float(optimize.brentq(lambda P: self.ln_fugacity_ratio(T, P),
                                     P_lo, P_hi, xtol=1e-12, rtol=1e-12))

    def saturation_volumes(self, T, P=None):
        """Saturated liquid and vapor molar volumes (m^3/mol) at T."""
        if P is None:
            P = self.vapor_pressure(T)
        return self.molar_volume(T, P, "liquid"), self.molar_volume(T, P, "vapor")
