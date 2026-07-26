"""UNIFAC activity-coefficient model (group contribution).

Two variants, selectable:
  - "modified"  : temperature-dependent (Dortmund) UNIFAC, Psi = exp[-(a/T + b + c T)],
                  combinatorial with r^(3/4). Fully parameterized from the legacy tables.
  - "original"  : classic UNIFAC, single parameter, Psi = exp(-a/T), combinatorial r.
                  NOT yet usable: the original subgroup R/Q constants were not in the
                  legacy .mat (only the a_mn matrix was). See code/data/README.md.

Faithful translation of the legacy MATLAB `calc_coeff` routine. An Aspen-optional
substitute for activity-coefficient calculations.

Eric M. Furst
July 2026
"""
import numpy as np

from .data import load_unifac_subgroups, load_unifac_interactions


class UNIFAC:
    def __init__(self, kind="modified"):
        if kind not in ("modified", "original"):
            raise ValueError("kind must be 'modified' or 'original'")
        self.kind = kind

        subs = load_unifac_subgroups()
        self.R = dict(zip(subs.subgroup_no, subs.R))
        self.Q = dict(zip(subs.subgroup_no, subs.Q))
        self.main = dict(zip(subs.subgroup_no, subs.main_group_no))  # subgroup -> main group

        inter = load_unifac_interactions(kind)
        self.A = {(int(i), int(j)): a for i, j, a in
                  zip(inter.main_i, inter.main_j, inter.a_ij)}
        if kind == "modified":
            self.B = {(int(i), int(j)): b for i, j, b in
                      zip(inter.main_i, inter.main_j, inter.b_ij)}
            self.C = {(int(i), int(j)): c for i, j, c in
                      zip(inter.main_i, inter.main_j, inter.c_ij)}
        else:
            self.B = self.C = {}

    # --- interaction matrix over the subgroups present -------------------
    def _psi(self, mains, T):
        n = len(mains)
        Psi = np.empty((n, n))
        for i, mi in enumerate(mains):
            for j, mj in enumerate(mains):
                a = self.A.get((mi, mj), 0.0) - self.A.get((mj, mj), 0.0)
                b = self.B.get((mi, mj), 0.0) - self.B.get((mj, mj), 0.0)
                c = self.C.get((mi, mj), 0.0) - self.C.get((mj, mj), 0.0)
                Psi[i, j] = np.exp(-(a / T + b + c * T))
        return Psi

    @staticmethod
    def _ln_Gamma(Qk, theta, Psi):
        """Group residual ln Gamma_k for group-fraction weights theta."""
        s1 = theta @ Psi                       # sum_m theta_m Psi[m,k]  (over k)
        s2 = Psi @ (theta / s1)                # sum_m theta_m Psi[k,m]/s1[m]  (over k)
        return Qk * (1 - np.log(s1) - s2)

    def gamma(self, groups, x, T):
        """Activity coefficients.

        groups : list (one per species) of {subgroup_no: count}
        x      : mole fractions (same order as `groups`)
        T      : temperature, K
        """
        if self.kind == "original":
            raise NotImplementedError(
                "original-UNIFAC R/Q subgroup constants are not yet loaded "
                "(only a_mn is available); see code/data/README.md. Use kind='modified'.")

        x = np.asarray(x, float)
        subs = sorted({k for g in groups for k in g})
        Rk = np.array([self.R[k] for k in subs])
        Qk = np.array([self.Q[k] for k in subs])
        mains = [int(self.main[k]) for k in subs]
        nu = np.array([[g.get(k, 0) for k in subs] for g in groups], float)  # species x subgroup

        # combinatorial term
        r = nu @ Rk
        q = nu @ Qk
        f = r ** 0.75 if self.kind == "modified" else r
        phi = r / (x @ r)
        phip = f / (x @ f)
        theta = q / (x @ q)
        lng_C = 1 - phip + np.log(phip) - 5 * q * (1 + np.log(phi / theta) - phi / theta)

        # residual term
        Psi = self._psi(mains, T)
        Xmix = x @ nu
        theta_mix = Qk * Xmix / np.sum(Qk * Xmix)
        lnG_mix = self._ln_Gamma(Qk, theta_mix, Psi)

        lng_R = np.empty(len(groups))
        for i in range(len(groups)):
            Xi = nu[i] / nu[i].sum()
            theta_i = Qk * Xi / np.sum(Qk * Xi)
            lnG_i = self._ln_Gamma(Qk, theta_i, Psi)
            lng_R[i] = np.sum(nu[i] * (lnG_mix - lnG_i))

        return np.exp(lng_C + lng_R)
