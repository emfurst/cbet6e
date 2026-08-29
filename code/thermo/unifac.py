"""UNIFAC activity-coefficient model (group contribution).

**Two parameter sets, and both are yours to choose** -- the same choice a process
simulator gives you, and for the same reason: UNIFAC is one idea regressed more than
once, and which regression you are standing on changes the answer.

  - "modified"  : temperature-dependent (Dortmund) UNIFAC, Psi = exp[-(a/T + b + c T)],
                  combinatorial with r^(3/4). **The default, and the book's model**
                  -- Table 9.5-2 is the Dortmund parameter set and Eq. 9.6-12a's
                  r^(3/4) is the Dortmund modification. Aspen Plus calls it UNIF-DMD.
  - "original"  : classic UNIFAC (Fredenslund, Jones and Prausnitz, 1975), single
                  temperature-independent a_mn, Psi = exp(-a/T), combinatorial in r.
                  What earlier editions of this book computed with, and what Aspen
                  Plus calls plain UNIFAC.

Everything downstream takes the model as an argument, so switching sets is one word:

    UNIFAC()              # Dortmund -- the book's default everywhere
    UNIFAC("original")    # the classic model, for comparison

Figure 10.2-8 is the illustration of why anyone would: on n-pentane/propionaldehyde
the two sets differ by 2-3 % in pressure and by almost nothing in vapor composition.

Both sets are complete: 92 subgroups over 47 main groups for Dortmund, 85 over 44 for
original. Neither is a stub.

The R/Q and the a_mn must come from the same parameter set, which is why `kind`
selects both together and neither table is meant to be loaded on its own.

**And the two sets do not agree on what a subgroup NUMBER means.** 27, 37-39 and
78-85 name different groups in each -- 27 is `FCH2O` in the original set and
`cy-CH2 OCH2` in Dortmund; 78-80 are `SiH3`/`SiH2`/`SiH` against `cy-CH2`/`cy-CH`/`cy-C`;
81-82 are `Si`/`SiH2O` against `OH(s)`/`OH(t)`; 85 is `NMP` against `CNH2`.

**Dortmund also splits the hydroxyl group and the original set does not**, which is the
divergence inside the otherwise-shared low numbers: 14 is a plain `OH` in the original
set and specifically `OH(p)` in Dortmund, where a secondary or tertiary alcohol takes 81
or 82 instead. So 2-propanol decomposes differently in the two sets even though every
number involved exists in both.

Both numbers always exist in both tables, so handing one set's group assignment to the
other model does not raise. Use `data.unifac_groups(name, kind)`, which reads the
matching table, rather than reusing a dict across kinds.

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

        subs = load_unifac_subgroups(kind)
        self.R = dict(zip(subs.subgroup_no, subs.R))
        self.Q = dict(zip(subs.subgroup_no, subs.Q))
        self.main = dict(zip(subs.subgroup_no, subs.main_group_no))  # subgroup -> main group
        # Kept so that coverage can be ASKED about rather than discovered by crashing.
        # The two sets do not cover the same subgroups, and a reader switching kinds
        # deserves to see that before the KeyError, not after it.
        self.subgroup_names = dict(zip(subs.subgroup_no, subs.subgroup_name))

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
        x = np.asarray(x, float)
        subs = sorted({k for g in groups for k in g})
        missing = [k for k in subs if k not in self.R]
        if missing:
            raise KeyError(
                f"subgroup(s) {missing} are not in the {self.kind!r} R/Q table "
                f"({len(self.R)} subgroups). "
                + ("The original set has no cyclic (cy-CH2), OH(s)/OH(t) or cy-CON "
                   "subgroups -- those are Dortmund's, and Dortmund numbers them 78-89, "
                   "which the original set uses for its silicon groups. Assign the "
                   "species with unifac_groups(name, 'original'), or use "
                   "kind='modified'."
                   if self.kind == "original" else
                   "Check the subgroup number against Table 9.5-2."))
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
