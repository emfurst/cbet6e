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
from .phi_phi import PhiPhiVLE

R = constants.R
_SQRT2 = np.sqrt(2.0)


class PRMixture(PhiPhiVLE):
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
    def from_database(cls, keys, kij=None, warn_missing=True):
        """Build from `pure_property.csv` keys, e.g. ["benzene", "toluene"].

        `kij` may be an explicit matrix, `None` for all zeros, or the string
        ``"table"`` to look the pairs up in **Table 9.4-1** (`code/data/pr_kij.csv`):

            >>> m = PRMixture.from_database(["ethane", "n-butane"], kij="table")

        **A pair Table 9.4-1 does not list is set to zero and warned about**, not
        filled in silently. The table is 65% blank and its own footnote says to
        substitute an estimate from a similar mixture -- a judgment the caller has to
        make, so it is surfaced rather than buried. Pass `warn_missing=False` once
        you have decided that zero is what you want.
        """
        if isinstance(kij, str):
            if kij != "table":
                raise ValueError("kij must be a matrix, None, or 'table'")
            from .data import pr_kij_matrix
            kij, missing = pr_kij_matrix(keys)
            if missing and warn_missing:
                import warnings
                pairs = ", ".join(f"{a}/{b}" for a, b in missing)
                warnings.warn(
                    f"Table 9.4-1 gives no k_ij for {pairs}; using 0. Its footnote "
                    f"says to estimate from a similar mixture instead.",
                    stacklevel=2)
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

    # --- the VLE drivers are in phi_phi.PhiPhiVLE ------------------------
    # `bubble_pressure`, `dew_pressure`, `bubble_temperature`, `dew_temperature`
    # and `flash` are inherited. They were written here, but they touch this class
    # only through `components` and `phi`, so Chapter 10 shares them with the
    # Wong-Sandler mixing rule rather than keeping two copies. See phi_phi.py.
