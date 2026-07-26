"""Peng-Robinson equation of state for mixtures.

Extends the pure-fluid `PengRobinson` (peng_robinson.py) to multicomponent
mixtures using the **van der Waals one-fluid mixing rules** with binary
interaction parameters k_ij (SIS Chapter 9):

    a_ij = (1 - k_ij) sqrt(a_i a_j)      a_mix = sum_i sum_j x_i x_j a_ij
                                         b_mix = sum_i x_i b_i

and the fugacity coefficient of a species in the mixture (SIS Eq. 9.4-9):

    ln phi_i = (b_i/b) (Z - 1) - ln(Z - B)
               - A / (2 sqrt2 B) [ 2 sum_j x_j a_ij / a - b_i / b ] * ln(...)

On top of that it provides the bubble/dew/flash drivers of Chapter 10
(the high-pressure phi-phi VLE method, SIS Section 10.3).

The Aspen-optional substitute for mixture PR property and VLE calculations.
SI units throughout: T in K, P in Pa, molar volume in m^3/mol, R = 8.314.

Eric M. Furst
July 2026
"""
import numpy as np
from numpy.polynomial import Polynomial
from scipy import constants, optimize

from .peng_robinson import PengRobinson

R = constants.R
_SQRT2 = np.sqrt(2.0)


class PRMixture:
    """Peng-Robinson EOS for a mixture with van der Waals one-fluid mixing.

    Parameters
    ----------
    components : sequence of PengRobinson
        The pure-component EOS objects, in a fixed order.
    kij : 2-D array, optional
        Symmetric matrix of binary interaction parameters (default all zeros).
        k_ii is ignored (a_ii = a_i). Book k_ij values are tabulated per system
        (SIS Table 9.4-1); pass them here.
    """

    def __init__(self, components, kij=None):
        self.components = list(components)
        n = len(self.components)
        if n < 2:
            raise ValueError("a mixture needs at least two components")
        if kij is None:
            self.kij = np.zeros((n, n))
        else:
            self.kij = np.asarray(kij, dtype=float)
            if self.kij.shape != (n, n):
                raise ValueError(f"kij must be {n}x{n}")
        self.b = np.array([c.b for c in self.components])  # pure b_i (T-independent)

    @classmethod
    def from_database(cls, keys, kij=None):
        """Build from `pure_property.csv` keys, e.g. ["benzene", "toluene"]."""
        return cls([PengRobinson.from_database(k) for k in keys], kij=kij)

    @property
    def n(self):
        return len(self.components)

    # --- mixing rules ----------------------------------------------------
    def a_pure(self, T):
        """Vector of pure-component a_i(T)."""
        return np.array([c.a(T) for c in self.components])

    def a_matrix(self, T):
        """The a_ij = (1 - k_ij) sqrt(a_i a_j) matrix (SIS Eq. 9.4-8)."""
        a = self.a_pure(T)
        return (1.0 - self.kij) * np.sqrt(np.outer(a, a))

    def a_mix(self, x, T):
        """Mixture a = sum_i sum_j x_i x_j a_ij (SIS Eq. 9.4-6)."""
        x = np.asarray(x, dtype=float)
        return float(x @ self.a_matrix(T) @ x)

    def b_mix(self, x):
        """Mixture b = sum_i x_i b_i (SIS Eq. 9.4-7)."""
        return float(np.asarray(x, dtype=float) @ self.b)

    def _AB(self, x, T, P):
        return self.a_mix(x, T) * P / (R * T) ** 2, self.b_mix(x) * P / (R * T)

    # --- roots -----------------------------------------------------------
    def compressibility(self, x, T, P):
        """All real roots Z of the PR cubic (same cubic as the pure fluid),
        ascending. The mixture enters only through A and B."""
        A, B = self._AB(x, T, P)
        alpha = -1 + B
        beta = A - 3 * B ** 2 - 2 * B
        gamma = -A * B + B ** 2 + B ** 3
        roots = Polynomial([gamma, beta, alpha, 1]).roots()
        return np.sort(roots.real[np.abs(roots.imag) < 1e-9])

    def Z(self, x, T, P, phase="vapor"):
        """Vapor (largest) or liquid (smallest) real root."""
        zs = self.compressibility(x, T, P)
        if len(zs) == 0:
            raise ValueError("no real root")
        return zs.max() if phase == "vapor" else zs.min()

    def molar_volume(self, x, T, P, phase="vapor"):
        return self.Z(x, T, P, phase) * R * T / P

    # --- component fugacities -------------------------------------------
    def _log_term(self, Z, B):
        return np.log((Z + (1 + _SQRT2) * B) / (Z + (1 - _SQRT2) * B))

    def ln_phi(self, x, T, P, phase="vapor"):
        """Vector of ln(fugacity coefficient) of each species in the mixture
        (SIS Eq. 9.4-9)."""
        x = np.asarray(x, dtype=float)
        a_ij = self.a_matrix(T)
        a = float(x @ a_ij @ x)
        b = self.b_mix(x)
        A, B = a * P / (R * T) ** 2, b * P / (R * T)
        Z = self.Z(x, T, P, phase)
        # 2 sum_j x_j a_ij for each i
        sum_j = 2.0 * (a_ij @ x)
        return (self.b / b * (Z - 1) - np.log(Z - B)
                - A / (2 * _SQRT2 * B) * (sum_j / a - self.b / b)
                * self._log_term(Z, B))

    def phi(self, x, T, P, phase="vapor"):
        return np.exp(self.ln_phi(x, T, P, phase))

    def fugacity(self, x, T, P, phase="vapor"):
        """Vector of component fugacities f_i = x_i phi_i P (Pa)."""
        return np.asarray(x, dtype=float) * self.phi(x, T, P, phase) * P

    # --- K-value initial guess ------------------------------------------
    def _wilson_K(self, T, P):
        """Wilson-correlation K_i = phi_i^L / phi_i^V estimate for initializing
        VLE iterations (SIS Eq. 10.2-... style ideal K-value seed)."""
        Tc = np.array([c.Tc for c in self.components])
        Pc = np.array([c.Pc for c in self.components])
        w = np.array([c.omega for c in self.components])
        return Pc / P * np.exp(5.373 * (1 + w) * (1 - Tc / T))

    # --- bubble / dew points --------------------------------------------
    def _equilibrate_y(self, x, T, P, y0, max_iter=200, tol=1e-10):
        """Inner successive-substitution: at fixed (x, T, P) find the incipient
        vapor y in equilibrium with liquid x. Returns (y_unnormalized_sum, y)."""
        y = np.array(y0, dtype=float)
        y /= y.sum()
        for _ in range(max_iter):
            phiL = self.phi(x, T, P, "liquid")
            phiV = self.phi(y, T, P, "vapor")
            K = phiL / phiV
            y_new = K * np.asarray(x, dtype=float)
            s = y_new.sum()
            y_new = y_new / s
            if np.max(np.abs(y_new - y)) < tol:
                return s, y_new
            y = y_new
        return s, y

    def bubble_pressure(self, x, T, P_guess=None, max_iter=100, tol=1e-9):
        """Bubble-point pressure and incipient vapor composition at (x, T).

        Returns (P, y). Liquid x is given; the vapor is in equilibrium.
        """
        x = np.asarray(x, dtype=float)
        K = self._wilson_K(T, 1e5)  # rough K to seed P and y
        if P_guess is None:
            P_guess = float(1.0 / np.sum(x / K)) if np.all(K > 0) else 1e5
        y = K * x
        y /= y.sum()
        P = P_guess
        for _ in range(max_iter):
            s, y = self._equilibrate_y(x, T, P, y)
            # sum(K x) = s must equal 1 at the bubble point; Newton on ln P
            f = s - 1.0
            if abs(f) < tol:
                return float(P), y
            dP = P * 1e-6
            s2, _ = self._equilibrate_y(x, T, P + dP, y)
            dfdP = (s2 - s) / dP
            P = P - f / dfdP
            if P <= 0:
                P = P_guess / 2
                P_guess = P
        return float(P), y

    def dew_pressure(self, y, T, P_guess=None, max_iter=100, tol=1e-9):
        """Dew-point pressure and incipient liquid composition at (y, T).

        Returns (P, x). Vapor y is given; the liquid is in equilibrium.
        """
        y = np.asarray(y, dtype=float)
        K = self._wilson_K(T, 1e5)
        if P_guess is None:
            P_guess = float(np.sum(y * K)) * 1e5 if np.all(K > 0) else 1e5
        x = y / K
        x /= x.sum()
        P = P_guess

        def sum_x(P, x):
            for _ in range(200):
                phiV = self.phi(y, T, P, "vapor")
                phiL = self.phi(x, T, P, "liquid")
                K = phiL / phiV
                x_new = y / K
                s = x_new.sum()
                x_new = x_new / s
                if np.max(np.abs(x_new - x)) < 1e-10:
                    return s, x_new
                x = x_new
            return s, x

        for _ in range(max_iter):
            s, x = sum_x(P, x)
            f = s - 1.0
            if abs(f) < tol:
                return float(P), x
            dP = P * 1e-6
            s2, _ = sum_x(P + dP, x)
            dfdP = (s2 - s) / dP
            P = P - f / dfdP
            if P <= 0:
                P = P_guess / 2
                P_guess = P
        return float(P), x

    def bubble_temperature(self, x, P, T_guess=None, max_iter=100, tol=1e-9):
        """Bubble-point temperature and vapor composition at (x, P).

        Returns (T, y).
        """
        x = np.asarray(x, dtype=float)
        if T_guess is None:
            Tb = np.array([c.Tc * 0.7 for c in self.components])  # crude
            T_guess = float(x @ Tb)
        T = T_guess
        K = self._wilson_K(T, P)
        y = K * x
        y /= y.sum()
        for _ in range(max_iter):
            s, y = self._equilibrate_y(x, T, P, y)
            f = s - 1.0
            if abs(f) < tol:
                return float(T), y
            dT = T * 1e-6
            s2, _ = self._equilibrate_y(x, T + dT, P, y)
            dfdT = (s2 - s) / dT
            T = T - f / dfdT
        return float(T), y

    def dew_temperature(self, y, P, T_guess=None, max_iter=100, tol=1e-9):
        """Dew-point temperature and liquid composition at (y, P).

        Returns (T, x).
        """
        y = np.asarray(y, dtype=float)
        if T_guess is None:
            Tb = np.array([c.Tc * 0.7 for c in self.components])
            T_guess = float(y @ Tb)
        T = T_guess

        def sum_x(T, x):
            for _ in range(200):
                phiV = self.phi(y, T, P, "vapor")
                phiL = self.phi(x, T, P, "liquid")
                K = phiL / phiV
                x_new = y / K
                s = x_new.sum()
                x_new = x_new / s
                if np.max(np.abs(x_new - x)) < 1e-10:
                    return s, x_new
                x = x_new
            return s, x

        K = self._wilson_K(T, P)
        x = y / K
        x /= x.sum()
        for _ in range(max_iter):
            s, x = sum_x(T, x)
            f = s - 1.0
            if abs(f) < tol:
                return float(T), x
            dT = T * 1e-6
            s2, _ = sum_x(T + dT, x)
            dfdT = (s2 - s) / dT
            T = T - f / dfdT
        return float(T), x

    # --- isothermal flash ------------------------------------------------
    def flash(self, z, T, P, max_iter=200, tol=1e-10):
        """Isothermal (T, P) flash of a feed of composition z.

        Returns (beta, x, y) where beta is the vapor molar fraction and x, y are
        the liquid and vapor compositions. beta == 0 -> subcooled liquid,
        beta == 1 -> superheated vapor (feed is single-phase at (T, P)).
        """
        z = np.asarray(z, dtype=float)
        K = self._wilson_K(T, P)
        beta = 0.5
        for _ in range(max_iter):
            beta = self._rachford_rice(z, K, beta)
            x = z / (1 + beta * (K - 1))
            y = K * x
            x = x / x.sum()
            y = y / y.sum()
            phiL = self.phi(x, T, P, "liquid")
            phiV = self.phi(y, T, P, "vapor")
            K_new = phiL / phiV
            if np.max(np.abs(K_new / K - 1)) < tol:
                K = K_new
                break
            K = K_new
        beta = self._rachford_rice(z, K, beta)
        beta = min(max(beta, 0.0), 1.0)
        x = z / (1 + beta * (K - 1))
        y = K * x
        return beta, x / x.sum(), y / y.sum()

    @staticmethod
    def _rachford_rice(z, K, beta0=0.5):
        """Solve sum_i z_i (K_i - 1) / (1 + beta (K_i - 1)) = 0 for beta in [0,1].
        Falls back to the single-phase edges when the feed does not split."""
        z, K = np.asarray(z, dtype=float), np.asarray(K, dtype=float)

        def g(beta):
            return np.sum(z * (K - 1) / (1 + beta * (K - 1)))

        # single-phase checks (SIS: bubble if g(0)<0, dew if g(1)>0)
        if g(0.0) <= 0:
            return 0.0
        if g(1.0) >= 0:
            return 1.0
        lo = 1.0 / (1.0 - K.max()) + 1e-10
        hi = 1.0 / (1.0 - K.min()) - 1e-10
        return float(optimize.brentq(g, max(lo, 0.0), min(hi, 1.0), xtol=1e-12))
