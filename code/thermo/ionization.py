"""ionization -- ionization, acidity and the charge on biomolecules, SIS Secs. 13.5-13.6.

`reaction.py` deliberately stopped at Sec. 13.4. Sections 13.5 to 13.6 are a different
calculation: no formation properties, no extent of reaction, no equation of state. What
they need instead is a proton balance, an apparent (concentration) equilibrium constant,
and -- when the solution is not dilute -- a Debye-Huckel activity coefficient. That is
this module.

    from thermo.ionization import WeakAcid, AminoAcid, PolyproticAcid

    WeakAcid(4.76).pH(0.1)                      # Illustration 13.5-3, 2.883
    WeakAcid(4.76).pH_with_base(0.1, 0.07)      # Illustration 13.5-7, 5.128
    AminoAcid(2.34, 9.6).pI                     # Illustration 13.6-2, 5.97
    PolyproticAcid([2.95, 5.41]).charge(7.0)    # Fig. 13.6-2

## Where each printed equation is

| SIS | what | entry point |
|-----|------|-------------|
| Eq. 13.5-2b, Table 13.5-1 | Ka,W and pKw against temperature | `Ka_water`, `pKw` |
| Eq. 13.5-3, 13.5-4, 13.5-5 | pK, pH, pOH | `pK`, `pH_from_activity` |
| Eqs. 13.5-6, 13.5-7a, 13.5-7b | strong acid, water ionization kept | `strong_acid` |
| Eqs. 13.5-8a, 13.5-8b | strong base, the same | `strong_base` |
| Eq. 13.5-9 | the two apparent constants K_HA and K_W | `WeakAcid` |
| Eqs. 13.5-10, 13.5-11 | weak acid, dissociation neglected in M_HA | `WeakAcid.pH(mode="simple")` |
| Eqs. 13.5-12a, 13.5-12b | weak acid, M_HA corrected | `WeakAcid.pH(mode="quadratic")` |
| Eq. 13.5-13 | the full cubic, water kept | `WeakAcid.pH(mode="cubic")` |
| Illustration 13.5-6 | dissociation with added salt, Debye-Huckel | `WeakAcid.dissociation` |
| Eqs. 13.5-14, 13.5-15, 13.5-16 | weak acid + strong base | `WeakAcid.alpha`, `.pH_with_base` |
| Eqs. 13.5-17, 13.5-18 | a buffered solution, exactly | `WeakAcid.pH_buffer` |
| Eq. 13.5-19 | Henderson-Hasselbalch | `WeakAcid.henderson_hasselbalch` |
| Illustration 13.5-2 | strong acid titrated with strong base | `titration_strong_acid` |
| Illustration 13.5-8 | weak acid titrated with strong base | `titration_weak_acid` |
| Fig. 13.5-1 | fraction deprotonated against pH and I | `fraction_deprotonated` |
| Eqs. 13.6-1 to 13.6-4c | a dibasic acid, and any polyprotic one | `PolyproticAcid` |
| Eqs. 13.6-5c, 13.6-6b, 13.6-7a | an amino acid's three states | `AminoAcid` |
| Eq. 13.6-7b | the isoelectric point | `AminoAcid.pI` |
| Eqs. 13.6-9, 13.6-10 | the zwitterion square, K1 K3 = K2 K4 | `zwitterion_constants` |
| Fig. 13.6-7 | net charge on a protein | `protein_charge`, `LYSOZYME` |
| Illustration 13.6-3 | an amino acid titrated, as a buffer | `amino_acid_titration` |
| Sec. 13.6, after Fig. 13.6-9 | the apparent constant kappa(pH) | `apparent_K` |
| Sec. 13.6, last page | pK(T) from the heat of ionization | `pK_at_temperature` |

## Molality, molarity, and why the module does not police the difference

Footnote 10 of Sec. 13.5 says it plainly: the section works in molality M, but for dilute
aqueous solutions one liter of solution holds one kilogram of water, "and therefore in
some of the calculations that follow, especially the titration calculations in this
section, we may ignore the distinction". Every concentration here is a number in mol/kg
that the illustrations also read as mol/L. The `(M = 1)` factors that make the printed
equilibrium constants dimensionless are carried implicitly, exactly as the illustrations
carry them.

## Apparent constants, not thermodynamic ones

Sec. 13.5 draws the distinction and then uses the apparent constant K_HA nearly
everywhere, because that is what is measured. So does this module: `WeakAcid` holds an
apparent pK and every routine that takes activity coefficients takes them explicitly.
Nothing here silently converts between K and K_a.

**Sec. 13.5 states that K_HA "has units of M^-1".** It has units of M: K_HA =
M_A- M_H+ / M_HA is (molality)^2 / molality, and it is K_HA/(M = 1) that the same
paragraph makes dimensionless. Filed as an erratum; the code uses the correct power.

## The sign switches, and why they exist

Three printed results in Secs. 13.5-13.6 do not agree with the figure they are attached
to, and one printed equation does not agree with the equation two lines above it. Rather
than silently correct them, the routines that touch them take a `convention` (or
`printed=`) argument so that a notebook can draw both curves and let them disagree
[[paired-notebooks-contrast]]. The three are:

* `fraction_deprotonated` -- Sec. 13.5's F(pH, I) for benzoyl tyrosine, Fig. 13.5-1.
  The general form and the worked BT form are printed with **different** signs on pH,
  and neither reproduces the figure; the thermodynamic result disagrees with the figure
  in the sign of the ionic-strength term as well.
* `amino_acid_titration` -- Illustration 13.6-3's X > 50 mL branch, which omits the
  5 mmol of hydroxide that has already neutralized the glycine.
* `apparent_K` -- the kappa(pH) of the benzoyl tyrosine + glycinamide reaction, whose
  printed form peaks four orders of magnitude below Fig. 13.6-10's plateau.

Each switch defaults to the form that is **thermodynamically correct**, and each names
the printed alternative explicitly. The validation notebook
`code/ch13/validation/ionization_module_validation.ipynb` is the audit trail.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

__all__ = [
    "TABLE_13_5_1", "Ka_water", "pKw", "pK", "pH_from_activity",
    "DEBYE_HUCKEL_LN", "DEBYE_HUCKEL_LOG10",
    "gamma_pm", "ionic_strength_1_1",
    "strong_acid", "strong_base",
    "WeakAcid", "titration_strong_acid", "titration_weak_acid",
    "fraction_deprotonated", "fraction_protonated_amine",
    "PolyproticAcid", "AminoAcid", "zwitterion_constants",
    "protein_charge", "IonizableGroup", "LYSOZYME", "isoelectric_point",
    "amino_acid_titration", "water_titration_pH", "apparent_K", "pK_at_temperature",
    "GAS_CONSTANT",
]

GAS_CONSTANT = 8.314                       # J/(mol K), the value Sec. 13.6 uses

# ---------------------------------------------------------------------------
# Table 13.5-1 -- the ionization of water
# ---------------------------------------------------------------------------
# T (deg C) : Ka,W.  The table also prints pKW, and every one of its eight rows is
# -log10 of the Ka,W beside it to the two decimals printed, so pKW is derived here
# rather than stored twice.
TABLE_13_5_1 = {
    0: 1.15e-15, 10: 2.88e-15, 20: 6.76e-15, 25: 1.00e-14,
    30: 1.48e-14, 40: 2.88e-14, 50: 5.50e-14, 60: 9.55e-14,
}


def Ka_water(T_C=25.0):
    """Ka,W at `T_C` degrees Celsius, Table 13.5-1.

    Off the tabulated temperatures the interpolation is **linear in pKw against
    1/T**, not linear in Ka,W: pKw is what the van 't Hoff equation makes nearly
    linear in 1/T, and over 0-60 C the table's own eight rows are straight in those
    coordinates to better than 0.01 in pKw. Interpolating Ka,W itself would be
    interpolating a quantity that changes by two orders of magnitude across the
    table.
    """
    T_C = np.asarray(T_C, dtype=float)
    T_tab = np.array(sorted(TABLE_13_5_1), dtype=float)
    if np.any(T_C < T_tab[0] - 1e-9) or np.any(T_C > T_tab[-1] + 1e-9):
        raise ValueError(
            f"Table 13.5-1 covers {T_tab[0]:g} to {T_tab[-1]:g} C; "
            f"{np.min(T_C):g} to {np.max(T_C):g} C was asked for. Extrapolating the "
            f"ionization of water past the table is not something the section supports")
    pK_tab = np.array([-np.log10(TABLE_13_5_1[int(t)]) for t in T_tab])
    inv_tab = 1.0 / (T_tab + 273.15)
    inv = 1.0 / (T_C + 273.15)
    # np.interp needs an increasing abscissa; 1/T decreases as T rises.
    return 10.0 ** (-np.interp(inv, inv_tab[::-1], pK_tab[::-1]))


def pKw(T_C=25.0):
    """pKw = -log10 Ka,W, Eq. 13.5-3 applied to Table 13.5-1."""
    return -np.log10(Ka_water(T_C))


def pK(K):
    """pK = -log10 K, Eq. 13.5-3. The minus sign is in the definition so that the
    value comes out positive."""
    return -np.log10(np.asarray(K, dtype=float))


def pH_from_activity(M_H, gamma_pm_=1.0):
    """pH = -log10 a_H+ = -log10(M_H+ gamma_pm), Eq. 13.5-4.

    Sec. 13.5 is not consistent about this. Illustration 13.5-1b and Illustration
    13.5-6 both include gamma_pm; Illustration 13.5-4's nonideal part reports
    -log10(M_H+) and calls it pH. Eq. 13.5-4 settles it: the activity.
    """
    return -np.log10(np.asarray(M_H, dtype=float) * np.asarray(gamma_pm_, dtype=float))


# ---------------------------------------------------------------------------
# Debye-Huckel, in the two forms Sec. 13.5 writes
# ---------------------------------------------------------------------------
# The section prints the limiting law twice on the same line, once as
# gamma_pm = exp(-1.176 sqrt(I)) and once as 10^(-0.5116 sqrt(I)). Those are not the
# same law: 1.176/ln(10) = 0.5107, and 0.5116*ln(10) = 1.1780.
#
# Which is the book's own? **Chapter 9's Table 9.10-1 gives alpha = 1.175 for water at
# 25 C**, i.e. a log10 slope of 0.5103 -- so Sec. 13.5's 1.176 agrees with Chapter 9 and
# 0.5116 is the outlier. But Sec. 13.5's WORKED numbers all use 0.5116: it is the only
# one of the three that gives Illustration 13.5-1's printed gamma_pm = 0.768 at I = 0.05
# (1.175 and 1.176 both give 0.769), and 1.178 is what the extended form after
# Fig. 13.5-1 and Illustration 13.5-6 print.
#
# So the value below is 1.178, because this module's job is to reproduce what Sec. 13.5
# prints. The disagreement with Chapter 9 is 0.26 % and moves nothing the book prints --
# it is a cross-chapter consistency question for the author, not a numerical one.
# `thermo.electrolytes.water_parameters` carries Chapter 9's temperature-dependent alpha
# for anything that should agree with Sec. 9.10 instead.
DEBYE_HUCKEL_LN = 1.178        # ln gamma_pm = -1.178 |z+ z-| sqrt(I) / (1 + A sqrt(I))
DEBYE_HUCKEL_LOG10 = DEBYE_HUCKEL_LN / np.log(10.0)     # 0.5116, to four places


def gamma_pm(I, z_product=1.0, extended=False, A=1.0, k_s=0.0):
    """Mean ionic activity coefficient at ionic strength `I` (mol/kg).

    `extended=False` (default) is the limiting law Sec. 13.5 uses for its 1:1
    electrolytes; `extended=True` is the form printed after Fig. 13.5-1, with the
    1 + A sqrt(I) denominator and the optional salting term k_s I. Illustration
    13.5-6 uses `extended=True, A=1`.

    This is Sec. 13.5's own constant, not Sec. 9.10's alpha(T). `thermo.electrolytes`
    carries the temperature-dependent version and the 9.10-17 / 9.10-18 forms; use it
    when the temperature is not 25 C. The two agree at 25 C to the digits Sec. 13.5
    prints.
    """
    I = np.asarray(I, dtype=float)
    if np.any(I < 0):
        raise ValueError("ionic strength cannot be negative")
    root = np.sqrt(I)
    ln_g = -DEBYE_HUCKEL_LN * z_product * root
    if extended:
        ln_g = ln_g / (1.0 + A * root) + k_s * I
    return np.exp(ln_g)


def ionic_strength_1_1(*molalities):
    """I = 1/2 sum z_i^2 M_i for singly charged ions only, which is what Sec. 13.5's
    illustrations need: pass the molality of each 1:1 salt or of each ion pair.

    Illustration 13.5-6 writes it out -- I = 1/2 (M_H+ + M_Ac- + M_Na+ + M_Cl-) =
    alpha M_HAc,o + M_NaCl -- and that is exactly `ionic_strength_1_1(alpha*M, M_NaCl)`.
    For any other charge type use `thermo.electrolytes.ionic_strength`.
    """
    return float(np.sum([float(m) for m in molalities]))


# ---------------------------------------------------------------------------
# strong acids and strong bases -- Eqs. 13.5-6 to 13.5-8b
# ---------------------------------------------------------------------------
def strong_acid(M_HA, gamma=1.0, Kw=1.0e-14):
    """The hydrogen ion molality of a strong, completely ionized acid, Eq. 13.5-7b.

    Returns ``(M_H, M_H_from_water)``: the total, and the part contributed by the
    ionization of water, Eq. 13.5-7a.

    The quadratic behind both is Eq. 13.5-6,

        Ka,W (M=1)^2 = (M_HA + [M_H+]_W) [M_H+]_W gamma_pm^2

    Eq. 13.5-7b is printed as ``M_H+ = M_HA + [M_H+]_HA``. The right-hand side is
    right; the middle label is not -- the term added to M_HA is ``[M_H+]_W``, the
    water's contribution of Eq. 13.5-7a, not the acid's own. Filed as an erratum.

    At M_HA = 0 this reduces to sqrt(Ka,W)/gamma, i.e. pH 7 at 25 C, which is the
    property the Comment on Illustration 13.5-2 claims for it -- and the reason the
    simpler M_H+ = M_HA is singular at the neutral point while this is not.
    """
    M_HA = np.asarray(M_HA, dtype=float)
    disc = M_HA**2 + 4.0 * Kw / gamma**2
    M_H_water = 0.5 * (-M_HA + np.sqrt(disc))
    return M_HA + M_H_water, M_H_water


def strong_base(M_BOH, gamma=1.0, Kw=1.0e-14):
    """Hydroxyl and hydrogen ion molalities for a strong base, Eqs. 13.5-8a and 8b.

    Returns ``(M_OH, M_H)``. M_OH is Eq. 13.5-8a, the mirror of Eq. 13.5-7b; M_H then
    follows from the water equilibrium, which is Eq. 13.5-8b written exactly rather
    than in its M_OH- ~ M_BOH approximation.
    """
    M_BOH = np.asarray(M_BOH, dtype=float)
    M_OH = 0.5 * (M_BOH + np.sqrt(M_BOH**2 + 4.0 * Kw / gamma**2))
    return M_OH, Kw / (M_OH * gamma**2)


# ---------------------------------------------------------------------------
# the weak acid -- Eqs. 13.5-9 through 13.5-19
# ---------------------------------------------------------------------------
class WeakAcid:
    """A weak acid HA <-> H+ + A-, held by its **apparent** pK.

    Parameters
    ----------
    pKa : float
        pK_HA = -log10 K_HA, with K_HA = M_H+ M_A- / M_HA the apparent
        (concentration) constant of Eq. 13.5-9.
    Kw : float
        The apparent water constant, default 1e-14 (25 C, Table 13.5-1).

    Sec. 13.5 offers three levels of approximation for the pH of a solution of this
    acid alone, and they are three modes of `pH`:

    ``"simple"``     Eq. 13.5-11.  M_HA held at its initial value.
    ``"quadratic"``  Eq. 13.5-12b. M_HA reduced by what dissociates; water neglected.
    ``"cubic"``      Eq. 13.5-13.  Nothing neglected. Rarely needed, as the text says.
    """

    def __init__(self, pKa, Kw=1.0e-14):
        self.pKa = float(pKa)
        self.Kw = float(Kw)

    @property
    def K(self):
        """The apparent dissociation constant K_HA, in molality units."""
        return 10.0 ** (-self.pKa)

    def __repr__(self):
        return f"WeakAcid(pKa={self.pKa:g}, Kw={self.Kw:g})"

    # -- the acid on its own ------------------------------------------------
    def M_H(self, M0, mode="quadratic"):
        """Hydrogen ion molality of a solution `M0` molal in the acid."""
        M0 = np.asarray(M0, dtype=float)
        K = self.K
        if mode == "simple":                                     # Eq. 13.5-10
            return np.sqrt(K * M0)
        if mode == "quadratic":                                  # Eq. 13.5-12a
            return 0.5 * K * (np.sqrt(1.0 + 4.0 * M0 / K) - 1.0)
        if mode == "cubic":                                      # Eq. 13.5-13
            return self._cubic(M0)
        raise ValueError(f"mode must be simple, quadratic or cubic, not {mode!r}")

    def _cubic(self, M0):
        """Eq. 13.5-13, solved as a bracketed root rather than by a cubic formula.

        M_H^3 + K M_H^2 - M_H (K M_HA,o + Kw) - K Kw = 0

        The bracket is physical and closed on both ends: the hydrogen ion molality
        cannot be below what pure water gives (sqrt(Kw), since the acid can only add
        protons) and cannot be above complete dissociation (M_HA,o + sqrt(Kw)). The
        left-hand side is monotone increasing across it, so the root is unique --
        which matters, because the cubic has two other roots that converge perfectly
        well and are negative [[spurious-root-is-the-default]].
        """
        K, Kw = self.K, self.Kw
        def f(x):
            return x**3 + K * x**2 - x * (K * M0 + Kw) - K * Kw
        lo, hi = np.sqrt(Kw), float(M0) + np.sqrt(Kw)
        if f(lo) > 0:                     # only if M0 is 0 to within rounding
            return lo
        return brentq(f, lo, hi, xtol=1e-300, rtol=8.9e-16)

    def pH(self, M0, mode="quadratic", gamma=1.0):
        """pH of a solution `M0` molal in the acid. See `M_H` for the modes."""
        if np.ndim(M0) and mode == "cubic":
            return np.array([self.pH(m, mode, gamma) for m in np.asarray(M0)])
        return pH_from_activity(self.M_H(M0, mode), gamma)

    # -- with added salt, Illustration 13.5-6 -------------------------------
    def dissociation(self, M0, M_salt=0.0, ideal=False, A=1.0):
        """Fractional dissociation alpha and pH with an added 1:1 salt.

        Illustration 13.5-6 exactly: the thermodynamic relation

            K_a,HA = alpha^2 M_HA,o gamma_pm^2 / [(1 - alpha)(1 M)]

        with I = alpha M_HA,o + M_salt and the extended Debye-Huckel law. `ideal=True`
        sets gamma_pm to 1, which is the illustration's part (a) -- and in that case
        the salt drops out entirely, as the illustration notes.

        Returns ``(alpha, pH, gamma_pm)``. The pH is -log10(alpha M_HA,o gamma_pm),
        i.e. an activity, which is what the illustration prints.

        This is the calculation that settles the direction of the salt effect for
        the whole section: alpha **rises** with added salt, so the acid is apparently
        stronger and its apparent pK falls. Fig. 13.5-1 shows the opposite.
        """
        M0 = float(M0)
        Ka = self.K

        def residual(alpha):
            g = 1.0 if ideal else gamma_pm(alpha * M0 + M_salt, extended=True, A=A)
            return alpha**2 * M0 / (1.0 - alpha) * g**2 - Ka

        alpha = brentq(residual, 1e-14, 1.0 - 1e-14, xtol=1e-15, rtol=8.9e-16)
        g = 1.0 if ideal else gamma_pm(alpha * M0 + M_salt, extended=True, A=A)
        return alpha, -np.log10(alpha * M0 * g), g

    # -- with a strong base, Eqs. 13.5-14 to 13.5-16 ------------------------
    def alpha(self, M0, M_base):
        """Fraction of the acid dissociated in the presence of `M_base` strong base,
        Eq. 13.5-15.

        The quadratic is  alpha^2 M_HA,o + alpha (K_HA - M_BOH) - K_HA = 0, whose
        two roots straddle zero (the product of the roots is -K_HA/M_HA,o < 0), so
        the physical root is unambiguously the ``+`` one.
        """
        M0, B, K = np.asarray(M0, float), np.asarray(M_base, float), self.K
        return (-(K - B) + np.sqrt((K - B) ** 2 + 4.0 * K * M0)) / (2.0 * M0)

    def M_H_with_base(self, M0, M_base):
        """Eq. 13.5-16, M_H+ = alpha M_HA,o - M_BOH written out.

        The asymmetry is real and is **not** a typo: the leading bracket carries
        ``K_HA + M_BOH`` while the one under the radical carries ``K_HA - M_BOH``.
        Substituting Eq. 13.5-15 into M_H+ = alpha M_HA,o - M_BOH produces exactly
        that. Illustration 13.5-7 prints ``(K_HA + M_BOH)^2`` inside the radical,
        which is the typo -- and its own printed answer, 7.445e-6, comes from the
        form here, not from the form it prints.
        """
        M0, B, K = np.asarray(M0, float), np.asarray(M_base, float), self.K
        return (-(K + B) + np.sqrt((K - B) ** 2 + 4.0 * K * M0)) / 2.0

    def pH_with_base(self, M0, M_base, gamma=1.0):
        """pH of a weak acid + strong base mixture, Eq. 13.5-16."""
        return pH_from_activity(self.M_H_with_base(M0, M_base), gamma)

    # -- buffered, Eqs. 13.5-18 and 13.5-19 ---------------------------------
    def pH_buffer(self, M_acid, M_salt):
        """pH of a buffered solution, Eq. 13.5-18, solved rather than iterated.

        Eq. 13.5-18 has pH on both sides. Rather than substitute repeatedly, the
        equation it came from is solved directly for M_H+: electroneutrality plus
        the two equilibria give

            K_HA (a M_H+ + Kw - M_H+^2) = M_H+ (b M_H+ + M_H+^2 - Kw)

        which is the cubic in the middle of p. 816 cleared of fractions. The root is
        bracketed between the pure-water value and the unbuffered acid's value: adding
        the conjugate salt can only *raise* the pH, and can never raise it past
        neutrality for an acid.
        """
        a, b, K, Kw = float(M_acid), float(M_salt), self.K, self.Kw

        def f(h):
            return K * (a * h + Kw - h * h) - h * (b * h + h * h - Kw)

        hi = float(self.M_H(a, "quadratic")) if a > 0 else np.sqrt(Kw)
        lo = np.sqrt(Kw) * 1e-8
        return -np.log10(brentq(f, lo, max(hi, np.sqrt(Kw)) * (1 + 1e-12),
                                xtol=1e-300, rtol=8.9e-16))

    def henderson_hasselbalch(self, M_acid, M_salt):
        """pH = pK_HA + log10(M_salt/M_acid), Eq. 13.5-19.

        The reduction of Eq. 13.5-18 valid for 4 <= pH <= 10, where both 10^-pH and
        10^(pH-pKw) are negligible beside a and b.
        """
        return self.pKa + np.log10(float(M_salt) / float(M_acid))


# ---------------------------------------------------------------------------
# titration curves -- Illustrations 13.5-2 and 13.5-8
# ---------------------------------------------------------------------------
def titration_strong_acid(X, V_acid=10.0, M_acid=0.20, M_base=0.20, Kw=1.0e-14):
    """pH of `V_acid` mL of a strong acid titrated with `X` mL of a strong base.

    Illustration 13.5-2, whose numbers are V_acid = 10 mL and both solutions 0.20 M.
    Below the equivalence point the excess acid sets the pH; above it the excess base
    does, through the water equilibrium; at the equivalence point the answer is 7.

    The equivalence point is handled by a tolerance rather than by an exact test on a
    float, and the returned pH there is exactly ``pKw/2``.
    """
    X = np.asarray(X, dtype=float)
    V = V_acid + X
    n_acid = V_acid * M_acid                      # mmol
    n_base = X * M_base
    excess = n_acid - n_base                      # mmol of H+ (>0) or OH- (<0)
    with np.errstate(divide="ignore", invalid="ignore"):
        M_H = np.where(excess > 0, excess / V, np.nan)
        M_OH = np.where(excess < 0, -excess / V, np.nan)
        pH = np.where(excess > 0, -np.log10(M_H),
                      np.where(excess < 0, 14.0 + np.log10(M_OH), 0.5 * -np.log10(Kw)))
    neutral = np.isclose(excess, 0.0, atol=1e-12 * max(n_acid, 1.0))
    return np.where(neutral, 0.5 * -np.log10(Kw), pH)


def titration_weak_acid(X, acid, V_acid=10.0, M_acid=0.20, M_base=0.20,
                        Kw=1.0e-14, printed=False):
    """pH of `V_acid` mL of a weak acid titrated with `X` mL of a strong base.

    Illustration 13.5-8: 10 mL of 0.20 M acetic acid, 0.20 M NaOH. Below the
    equivalence point Eq. 13.5-16 is used on the **diluted** concentrations

        M_HA,o = V_acid M_acid / (V_acid + X)     M_NaOH = M_base X / (V_acid + X)

    and above it the excess hydroxide sets the pH through the water equilibrium.

    `printed=True` reproduces the illustration's printed X > equivalence formula,

        pH = -log10[ (10 - X) 1e-14 / (2 X - 2) ]

    which is **undefined** past the equivalence point: its argument is negative for
    every X > 10. Two errors in one line -- the numerator's ``(10 - X)`` should be
    ``(10 + X)``, and the denominator has been multiplied by ten without the
    numerator. The line two above it in the same illustration,
    M_H+ = (10 + X) 1e-14 / (0.2 X - 2), is right, and is what `printed=False` uses.
    """
    X = np.asarray(X, dtype=float)
    V = V_acid + X
    n_acid, n_base = V_acid * M_acid, X * M_base
    equiv = n_acid / M_base                        # mL of base at the equivalence point

    if printed:
        # the illustration's own printed expression, reproduced without repair
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(
                X > equiv,
                -np.log10((V_acid - X) * Kw / (2.0 * X - 2.0)),
                _weak_below(X, acid, V_acid, M_acid, M_base, equiv, Kw))

    with np.errstate(divide="ignore", invalid="ignore"):
        above = 14.0 + np.log10((n_base - n_acid) / V)
    return np.where(X > equiv, above,
                    _weak_below(X, acid, V_acid, M_acid, M_base, equiv, Kw))


def _weak_below(X, acid, V_acid, M_acid, M_base, equiv, Kw):
    """The X <= equivalence branch: Eq. 13.5-16 on the diluted concentrations."""
    V = V_acid + X
    M_HA_o = V_acid * M_acid / V
    M_BOH = M_base * X / V
    with np.errstate(divide="ignore", invalid="ignore"):
        M_H = acid.M_H_with_base(M_HA_o, M_BOH)
        pH = -np.log10(np.where(M_H > 0, M_H, np.nan))
    return np.where(np.isclose(X, equiv, atol=1e-12 * max(equiv, 1.0)),
                    0.5 * -np.log10(Kw), pH)


# ---------------------------------------------------------------------------
# the fraction deprotonated -- Sec. 13.5 and Fig. 13.5-1
# ---------------------------------------------------------------------------
def fraction_deprotonated(pH, pKa, I=0.0, convention="thermodynamic"):
    """Fraction of a carboxylic acid RCOOH that is deprotonated, at `pH` and ionic
    strength `I`.

    Sec. 13.5 develops this twice -- ideally, then with the Debye-Huckel correction --
    and Fig. 13.5-1 plots it for benzoyl tyrosine, pK = 3.7, at I = 0, 1, 5 and 10.
    **The two printed forms disagree with each other and neither reproduces the
    figure.** All three readings are available here so that a notebook can plot them
    together:

    ``"ideal"``
        F = 10^(pH - pK)/(1 + 10^(pH - pK)); `I` is ignored. This one is right and
        is printed correctly, half-deprotonation at pH = pK.

    ``"thermodynamic"`` (default)
        F built from the section's own K1 = [RCOO-] gamma_pm 10^-pH / [RCOOH],
        i.e. the ratio [RCOO-]/[RCOOH] = K1/(gamma_pm 10^-pH), giving an exponent
        ``pH - pK + 0.5116 sqrt(I)``. Half-deprotonation at pH = pK - 0.5116 sqrt(I):
        adding salt makes a neutral acid apparently **stronger**, so the curve moves
        to **lower** pH. That is the direction Illustration 13.5-6 and Table 13.5-2
        of this same section compute for acetic acid.

    ``"figure"``
        Exponent ``pH - pK - 0.5116 sqrt(I)``. Half-deprotonation at
        pH = pK + 0.5116 sqrt(I) -- 3.70, 4.21, 4.84, 5.32 for the figure's four
        ionic strengths, which is what Fig. 13.5-1 plots and what the sentence under
        it describes ("a slight shift ... to higher pH with increasing ionic
        strength"). The gamma_pm sits on the wrong side.

    ``"printed_general"``
        Exponent ``-pH - pK - 0.5116 sqrt(I)``, the general F printed on p. 803. It
        is a decreasing function of pH: it says a carboxylic acid *protonates* as the
        solution is made more basic.

    ``"printed_bt"``
        Exponent ``pH + pK - 0.5116 sqrt(I)``, the benzoyl tyrosine F printed four
        lines below the general one. Half-deprotonation at pH = -pK, off the left of
        the figure; F is above 0.999 everywhere on the figure's axes.
    """
    pH = np.asarray(pH, dtype=float)
    s = DEBYE_HUCKEL_LOG10 * np.sqrt(np.asarray(I, dtype=float))
    exponents = {
        "ideal":           pH - pKa,
        "thermodynamic":   pH - pKa + s,
        "figure":          pH - pKa - s,
        "printed_general": -pH - pKa - s,
        "printed_bt":      pH + pKa - s,
    }
    try:
        e = exponents[convention]
    except KeyError:
        raise ValueError(
            f"convention must be one of {sorted(exponents)}, not {convention!r}"
        ) from None
    r = 10.0 ** e
    return r / (1.0 + r)


def fraction_protonated_amine(pH, pKa):
    """Fraction of an ammonium residue R-NH3+ that is still protonated (charged),
    10^(pK-pH)/(1 + 10^(pK-pH)).

    Sec. 13.5's last equation before Illustration 13.5-2. The section calls this "the
    fraction of the ammonium residue ionized", which is the protonated, charged form
    -- the opposite sense from `fraction_deprotonated`, and the printed expression is
    correct in that sense.
    """
    r = 10.0 ** (pKa - np.asarray(pH, dtype=float))
    return r / (1.0 + r)


# ---------------------------------------------------------------------------
# Sec. 13.6 -- polyprotic acids
# ---------------------------------------------------------------------------
class PolyproticAcid:
    """An acid with n ionizable protons, held by its n successive apparent pKs.

    Eqs. 13.6-1 to 13.6-4c do this for two protons (phthalic acid, pK1 = 2.95 and
    pK2 = 5.41); the algebra generalizes without change, and Problem 13.77 asks for
    six. The fractions are the terms of

        f_j = 10^(-[j pH + sum_{i<=j} pK_i]) / sum over all states

    normalized, with state 0 the fully protonated acid.

    ``charge`` is the average charge, Illustration 13.6-1's

        Charge = 0 f_H2A - 1 f_HA- - 2 f_A2-

    which need not be an integer even though no single molecule can carry a
    non-integer charge.
    """

    def __init__(self, pKs, charge_0=0):
        self.pKs = [float(p) for p in pKs]
        if any(b < a for a, b in zip(self.pKs, self.pKs[1:])):
            raise ValueError(
                f"successive pKs must not decrease -- got {self.pKs}. A dibasic acid "
                f"loses its first proton more easily than its second")
        self.charge_0 = int(charge_0)

    @property
    def n(self):
        """Number of ionizable protons."""
        return len(self.pKs)

    def __repr__(self):
        return f"PolyproticAcid(pKs={self.pKs}, charge_0={self.charge_0})"

    def fractions(self, pH):
        """Array of shape ``(n+1, ...)``: the fraction in each ionization state,
        most protonated first. Eqs. 13.6-4a, 4b, 4c for the dibasic case.

        Built in log space and shifted by the maximum before exponentiating, so that
        a pH ten units from every pK does not overflow -- the printed form,
        10^-2pH over a sum of 10^-2pH terms, is 10^-30 over 10^-30 at pH 15.
        """
        pH = np.asarray(pH, dtype=float)
        logs = [np.zeros_like(pH)]
        run = 0.0
        for j, p in enumerate(self.pKs, start=1):
            run += p
            logs.append(-(run - j * pH))
        L = np.stack(logs)                    # log10 of each unnormalized term
        L = L - L.max(axis=0, keepdims=True)
        w = 10.0 ** L
        return w / w.sum(axis=0)

    def charge(self, pH):
        """Average net charge, Illustration 13.6-1."""
        f = self.fractions(pH)
        z = self.charge_0 - np.arange(self.n + 1)      # 0, -1, -2, ... for an acid
        return np.tensordot(z, f, axes=(0, 0))

    def pH_of_max(self, j):
        """The pH at which the fraction in state `j` is largest.

        For 0 < j < n this is (pK_j + pK_{j+1})/2 -- Illustration 13.6-1 derives it
        for j = 1 by differentiating the denominator, and gets 4.18 for phthalic acid.
        The end states have no interior maximum and this raises for them.
        """
        if not 0 < j < self.n:
            raise ValueError(
                f"state {j} is an end state of a {self.n}-proton acid; its fraction "
                f"is monotone in pH and has no interior maximum")
        return 0.5 * (self.pKs[j - 1] + self.pKs[j])


# ---------------------------------------------------------------------------
# Sec. 13.6 -- amino acids
# ---------------------------------------------------------------------------
class AminoAcid:
    """An amino acid in the three-state description of Fig. 13.6-4: A (+1, both
    groups protonated), B (neutral), C (-1, both deprotonated).

    Sec. 13.6 solves the two reactions **separately** rather than simultaneously,
    on the argument that "the equilibrium constants are generally so different in
    value for amino acids that one reaction will have gone to completion before the
    other reaction has occurred to an appreciable extent". Eqs. 13.6-5c and 13.6-6b
    are the two halves of that, and `fractions` reproduces them.

    `exact=True` drops the approximation and uses the coupled result -- which is
    `PolyproticAcid([pK1, pK2])` -- so that a notebook can measure what the
    approximation costs. For glycine (pK1 = 2.34, pK2 = 9.6) it costs nothing
    visible; for an amino acid whose pKs are closer together it does not.
    """

    def __init__(self, pK1, pK2, name=None):
        self.pK1, self.pK2 = float(pK1), float(pK2)
        if self.pK2 < self.pK1:
            raise ValueError(
                f"pK1 ({pK1}) is the carboxyl group and pK2 ({pK2}) the ammonium; "
                f"pK2 must be the larger")
        self.name = name

    def __repr__(self):
        tag = f", name={self.name!r}" if self.name else ""
        return f"AminoAcid(pK1={self.pK1:g}, pK2={self.pK2:g}{tag})"

    def fractions(self, pH, exact=False):
        """``(f_A, f_B, f_C)`` -- positive, neutral and negative forms.

        `exact=False` is Eqs. 13.6-5c and 13.6-6b with f_B taken by difference, which
        is what Illustration 13.6-2 does. `exact=True` is the coupled two-proton
        result.
        """
        if exact:
            return PolyproticAcid([self.pK1, self.pK2]).fractions(pH)
        pH = np.asarray(pH, dtype=float)
        rA = 10.0 ** (self.pK1 - pH)                     # Eq. 13.6-5c
        f_A = rA / (1.0 + rA)
        f_C = 1.0 / (1.0 + 10.0 ** (self.pK2 - pH))      # Eq. 13.6-6b
        return np.stack([f_A, 1.0 - f_A - f_C, f_C])

    def charge(self, pH, exact=False):
        """Average net charge, Eq. 13.6-7a: (+1) f_A + (-1) f_C."""
        f = self.fractions(pH, exact=exact)
        return f[0] - f[2]

    @property
    def pI(self):
        """Isoelectric point, Eq. 13.6-7b: (pK1 + pK2)/2.

        Setting Eq. 13.6-7a to zero gives 10^(pK1 - pH) 10^(pK2 - pH) = 1 exactly,
        with no approximation -- so 13.6-7b is not the "left to the reader" estimate
        the text implies but an identity of the three-state model.
        """
        return 0.5 * (self.pK1 + self.pK2)

    def at_temperature(self, T_C, dH1=0.0, dH2=-44770.0, T1_C=25.0):
        """A new `AminoAcid` with both pKs shifted to `T_C` by the van 't Hoff
        relation, Sec. 13.6's last page.

        The defaults are its glycine example: the ammonium group's heat of ionization
        is about -44.77 kJ/mol and the carboxylic group's is approximately zero.
        """
        return AminoAcid(pK_at_temperature(self.pK1, dH1, T1_C, T_C),
                         pK_at_temperature(self.pK2, dH2, T1_C, T_C),
                         name=self.name)


def zwitterion_constants(K1, K2, K3, K4, rtol=1e-9):
    """The two measurable constants of the four-state zwitterion square, Eq. 13.6-8.

    Returns ``(K1_bar, K2_bar)`` where, by Eqs. 13.6-10a and 13.6-10b,

        K1_bar = K1 + K2            K2_bar = K3 K4 / (K3 + K4)

    GH and GH' cannot be told apart in the laboratory, so what is measured is the sum
    in the first case and the parallel combination in the second.

    The four constants are not independent: Eqs. 13.6-9a and 13.6-9b both evaluate to
    M_G- (M_H+)^2 / M_GH2+ along the two different paths, so K1 K3 = K2 K4. That is
    checked here rather than assumed, because a set of four measured constants that
    violates it is a data error and not a new state.
    """
    if not np.isclose(K1 * K3, K2 * K4, rtol=rtol):
        raise ValueError(
            f"K1 K3 = {K1 * K3:.6g} but K2 K4 = {K2 * K4:.6g}; Eqs. 13.6-9a and 9b "
            f"require them to be equal, since both are M_G- (M_H+)^2 / M_GH2+")
    return K1 + K2, K3 * K4 / (K3 + K4)


# ---------------------------------------------------------------------------
# Sec. 13.6 -- proteins
# ---------------------------------------------------------------------------
class IonizableGroup:
    """One kind of ionizable group on a protein: how many, what pK, and whether it
    is charged when protonated (cationic) or when deprotonated (anionic)."""

    __slots__ = ("name", "count", "pK", "cationic")

    def __init__(self, name, count, pK, cationic):
        self.name, self.count, self.pK = name, int(count), float(pK)
        self.cationic = bool(cationic)

    def charge(self, pH):
        """This group's contribution to the net charge at `pH`."""
        pH = np.asarray(pH, dtype=float)
        if self.cationic:                      # +1 while protonated
            return self.count / (1.0 + 10.0 ** (pH - self.pK))
        return -self.count / (1.0 + 10.0 ** (self.pK - pH))     # -1 once deprotonated

    def __repr__(self):
        kind = "cationic" if self.cationic else "anionic"
        return f"IonizableGroup({self.name!r}, {self.count}, pK={self.pK:g}, {kind})"


# Hen egg white lysozyme, the protein of Fig. 13.6-7.
#
# Sec. 13.6 gives three facts about it and no pK values: 32 ionizable groups, "19 basic
# groups and 13 acid groups", and a charge running from +19 in very acidic solutions to
# -13 in very basic ones. The composition below is the one that satisfies all three --
# the amino acid counts of hen egg white lysozyme (UniProt P00698, residues 19-147)
# plus the two chain termini -- and the pKs are the standard model values for free
# residues, since the section's own treatment ignores the shifts a folded structure
# imposes.  The figure is therefore NOT strictly reproducible: the book prints no
# pKs. What is reproducible is the shape, the two end values, and the statement that
# the isoelectric point is near 11.2.
LYSOZYME = [
    IonizableGroup("alpha-amino (N terminus)",  1,  7.7,  True),
    IonizableGroup("histidine",                 1,  6.0,  True),
    IonizableGroup("lysine",                    6, 10.5,  True),
    IonizableGroup("arginine",                 11, 12.5,  True),
    IonizableGroup("alpha-carboxyl (C term.)",  1,  3.1,  False),
    IonizableGroup("aspartic acid",             7,  3.9,  False),
    IonizableGroup("glutamic acid",             2,  4.3,  False),
    IonizableGroup("tyrosine",                  3, 10.1,  False),
]


def protein_charge(pH, groups=None):
    """Average net charge on a protein at `pH`, summed over its ionizable groups.

    Fig. 13.6-7. Each group contributes independently -- the model has no interaction
    between sites, which is the same approximation Eq. 13.6-3 makes for a dibasic acid
    and the reason Problem 13.77 can be answered by writing six terms.
    """
    groups = LYSOZYME if groups is None else groups
    return sum(g.charge(pH) for g in groups)


def isoelectric_point(groups=None, bracket=(0.0, 14.0)):
    """The pH at which `protein_charge` vanishes, by bisection on the bracket.

    The charge is strictly decreasing in pH -- every term is -- so the root is unique
    and bisection cannot land on a spurious one.
    """
    groups = LYSOZYME if groups is None else groups
    lo, hi = bracket
    if protein_charge(lo, groups) < 0 or protein_charge(hi, groups) > 0:
        raise ValueError(
            f"the net charge does not change sign on pH {lo} to {hi}: "
            f"{protein_charge(lo, groups):+.2f} to {protein_charge(hi, groups):+.2f}")
    return brentq(lambda p: protein_charge(p, groups), lo, hi, xtol=1e-12)


# ---------------------------------------------------------------------------
# Illustration 13.6-3 -- an amino acid titrated with a strong base
# ---------------------------------------------------------------------------
def _amino_acid_titration_exact(X, amino, V0, M0, M_base, Kw):
    """The proton condition for a diprotic amino acid titrated with a strong base,
    solved as one equation over the whole titration -- `amino_acid_titration`'s
    `water=True` branch.

    The solution holds the fully protonated amino acid (GH2+, as its hydrochloride)
    with its counter-ion, plus the Na+ that comes in with the base. Charge balance is

        [Na+] + [H+] + [GH2+] = [Cl-] + [OH-] + [G-]

    and with [Cl-] = C_T and f_GH2+ + f_GH + f_G- = 1 that is

        C_B + [H+] - C_T (f_GH + 2 f_G-) - Kw/[H+] = 0.

    The left side is strictly decreasing in pH -- every term is -- so the root is
    unique and bisection cannot land on a spurious one [[spurious-root-is-the-default]].
    """
    X = np.atleast_1d(np.asarray(X, dtype=float))
    K1, K2 = 10.0 ** (-amino.pK1), 10.0 ** (-amino.pK2)
    n0 = V0 * M0

    pH = np.empty_like(X)
    frac = np.empty((3,) + X.shape)

    for i, x in enumerate(X):
        V = V0 + x
        C_T, C_B = n0 / V, M_base * x / V

        def residual(p):
            h = 10.0 ** (-p)
            D = h * h + K1 * h + K1 * K2
            return C_B + h - C_T * (K1 * h + 2.0 * K1 * K2) / D - Kw / h

        p = brentq(residual, -2.0, 16.0, xtol=1e-13)
        h = 10.0 ** (-p)
        D = h * h + K1 * h + K1 * K2
        pH[i] = p
        frac[:, i] = [h * h / D, K1 * h / D, K1 * K2 / D]

    return pH, frac


def amino_acid_titration(X, amino, V0=25.0, M0=0.1, M_base=0.1, Kw=1.0e-14,
                         printed=False, water=False):
    """Illustration 13.6-3: `V0` mL of `M0` molal amino acid titrated with `X` mL of
    a `M_base` molal strong base.

    Returns ``(pH, fractions)`` with `fractions` of shape ``(3, ...)`` holding
    f_GH2+, f_GH and f_G-.

    The illustration's two-regime approximation is used verbatim -- the first
    ionization is taken to complete before the second begins, so that below the first
    equivalence point only

        K1 = alpha (alpha - M_base X) / [(n0 - alpha)(V0 + X)]

    is solved, and above it only the corresponding equation in beta. Amounts are in
    mmol and volumes in mL throughout, so that a molality is a ratio of the two.

    `printed=False` (default) subtracts, past the **second** equivalence point, the
    hydroxide that has already neutralized the amino acid:

        M_OH- = (M_base X - 2 n0) / (V0 + X)

    The illustration prints ``M_OH- = 0.1 X/(25 + X)`` with no such subtraction, even
    though the sentence introducing it says "all the hydrogen ions resulting from the
    ionization of glycine have been neutralized". `printed=True` reproduces that. The
    difference is large and visible: at X = 50.5 mL the printed form gives pH 12.83
    and the corrected one 10.82, and Fig. 13.6-8 plots a curve that leaves the second
    plateau near 11 -- so the figure was drawn from the corrected form.

    A second and separate disagreement is not repaired here and cannot be: over
    25 < X < 50 mL these equations put the second buffer plateau at pH = pK2 = 9.60,
    rising from 8.2 to 11.0 across the interval, while Fig. 13.6-8 plots a nearly flat
    8.3 to 8.7. The equations are self-consistent and reproduce every number the
    illustration states in words (81 % GH2+, 19 % GH and pH 1.72 at X = 0); what drew
    that branch of the figure is not established.

    At the two equivalence points themselves the returned pH is **NaN**. These
    equations neglect the ionization of water, so at X = 25 and X = 50 mL they give
    M_H+ = 0 and no logarithm -- the same singularity the Comment on Illustration
    13.5-2 identifies for the simpler strong-acid equations.

    `water=True` solves the full proton condition instead, keeping the ionization of
    water and both glycine equilibria at once,

        C_B + [H+] - C_T (f_GH + 2 f_G-) - Kw/[H+] = 0

    which is one equation valid over the whole titration. Use it for anything that is
    PLOTTED. The two-regime form above is not merely singular AT the equivalence
    points, it is wrong in a neighborhood of them, and the error has the wrong sign:
    driving M_H+ to zero sends the pH up, so the curve rises to 13.0 at X = 49.99 mL,
    breaks, and restarts at 9.1 -- a titration curve that runs backwards as base is
    added. `water=True` passes smoothly through pH 5.99 at the first equivalence point
    (the isoelectric point, 5.97) and 11.05 at the second.

    Nothing is given up in fidelity. The exact form reproduces every number the
    illustration states in words -- pH 1.716 at X = 0, 80.8 % GH2+ and 19.2 % GH,
    pH 9.598 = pK2 at the half-point. It also makes `printed` moot: the hydroxide
    already spent on the amino acid is subtracted by construction, so `water=True`
    ignores `printed` rather than silently combining with it.

    Where the two forms differ is NOT symmetric about the two equivalence points, and
    it is wider than a first look suggests. Measured for glycine on a 0.05 mL grid,
    they differ by more than 0.02 pH units only at X = 25.0 itself and over roughly
    46 to 54 mL -- a band about 4 mL either side of the SECOND equivalence point,
    where the missing 5 mmol of hydroxide compounds the divergence. Five milliliters
    away from either point they agree to 0.013.

    It does NOT reconcile Fig. 13.6-8. The exact curve still rises 8.2 to 11.0
    across 25 < X < 50 mL, on pK2 = 9.60, where the figure plots a flat 8.3 to 8.7.
    That disagreement is the separate one above and is untouched by this switch.
    """
    if water:
        return _amino_acid_titration_exact(X, amino, V0, M0, M_base, Kw)
    X = np.atleast_1d(np.asarray(X, dtype=float))
    n0 = V0 * M0                                   # mmol of amino acid
    K1, K2 = 10.0 ** (-amino.pK1), 10.0 ** (-amino.pK2)
    eq1, eq2 = n0 / M_base, 2.0 * n0 / M_base      # mL at the two equivalence points

    pH = np.empty_like(X)
    frac = np.empty((3,) + X.shape)

    for i, x in enumerate(X):
        V = V0 + x
        added = M_base * x                         # mmol of OH- added
        if x < eq1 * (1 - 1e-12):
            a = brentq(lambda al: al * (al - added) / ((n0 - al) * V) - K1,
                       added + 1e-16, n0 - 1e-13, xtol=1e-300, rtol=8.9e-16)
            b = 0.0
        elif x < eq2 * (1 - 1e-12):
            a = n0
            lo = max(added - n0, 0.0)
            b = brentq(lambda be: be * (n0 + be - added) / ((n0 - be) * V) - K2,
                       lo + 1e-16, n0 - 1e-14, xtol=1e-300, rtol=8.9e-16)
        else:
            a, b = n0, n0
        M_H = (a + b - added) / V
        if x > eq2:
            excess = added if printed else added - 2.0 * n0
            M_H = Kw * V / excess if excess > 0 else np.nan
        # At an equivalence point every proton released has been neutralized, so
        # these equations give M_H+ = 0 and no pH at all. That is not a solver
        # failure -- it is the same singularity the Comment on Illustration 13.5-2
        # identifies for the simpler strong-acid equations, and it has the same
        # cause: the ionization of water has been dropped. The honest answer is
        # NaN, which leaves a gap in the plotted curve. Returning pKw/2 = 7 instead
        # would put a downward spike into the figure at the one X where the
        # solution is in fact strongly basic.
        pH[i] = -np.log10(M_H) if M_H > 0 else np.nan
        frac[:, i] = [(n0 - a) / n0, (a - b) / n0, b / n0]

    return pH, frac


def water_titration_pH(X, V0=25.0, M_base=0.1, Kw=1.0e-14, water=False):
    """pH of `V0` mL of pure water on adding `X` mL of `M_base` strong base --
    the dashed "without glycine" curve of Fig. 13.6-8,

        pH = 14 + log10(M_base X / (V0 + X))

    printed correctly in the illustration.

    That printed form counts only the hydroxide added and none of the water's own,
    so it has the same defect as the two-regime amino acid equations at the low end:
    below about pH 7.5 it returns values BELOW 7 for a solution to which base has
    been added, and it diverges to -infinity as X goes to zero. The `X > 0` guard
    below patches the single point X = 0 to 7 and leaves the neighborhood wrong, so
    on a plotting grid the dashed curve steps discontinuously out of 7.

    `water=True` solves the charge balance with the water equilibrium kept,

        [H+] + C_B = Kw/[H+]   ->   h = (-C_B + sqrt(C_B^2 + 4 Kw)) / 2

    which is continuous from pH 7 up and never falls below it. Use it for anything
    that is PLOTTED; the two agree to three decimals above pH 8, which is over
    essentially the whole figure.
    """
    X = np.asarray(X, dtype=float)
    if water:
        C_B = M_base * X / (V0 + X)
        h = 0.5 * (-C_B + np.sqrt(C_B * C_B + 4.0 * Kw))
        return -np.log10(h)
    with np.errstate(divide="ignore"):
        return np.where(X > 0, 14.0 + np.log10(M_base * X / (V0 + X)), 7.0)


# ---------------------------------------------------------------------------
# the apparent equilibrium constant of a reaction between ionizable species
# ---------------------------------------------------------------------------
def apparent_K(pH, K, pK_acid, pK_base, convention="printed"):
    """kappa(pH), the apparent constant for A + B -> AB when A is a weak acid and B a
    weak base, Sec. 13.6's benzoyl tyrosine + glycinamide example.

    ``kappa = [AB] / {([AH] + [A-])([B] + [BH+])}``, so kappa = K f_AH f_B, a product
    of two fractions each at most one. Hence **kappa <= K always**, with equality only
    where both reacting forms are simultaneously dominant.

    ``"printed"`` (default)
        The book's own algebra, kappa = K/[(1 + 10^(pH - pK_acid))(1 + 10^(pK_base -
        pH))]. It follows correctly from the printed definitions: the acid must be
        protonated (pH below pK_acid = 3.7) and the amine free (pH above pK_base =
        7.93) at the same time. Those windows do not overlap, so kappa never exceeds
        2.8e-5 -- about K/17000.

    ``"figure"``
        Both exponents reversed, kappa = K/[(1 + 10^(pK_acid - pH))(1 + 10^(pH -
        pK_base))], which does have a plateau at K between the two pKs. This is the
        only form that produces Fig. 13.6-10's shape, and it reproduces the figure's
        rise at pH 3.7; the figure's fall is near pH 9.5 rather than at 7.93, so even
        this does not reproduce the figure completely. What drew Fig. 13.6-10 is not
        established, and no story is offered for it here.

    Sec. 13.6 gives K = 0.49 "M" for this reaction. For A + B -> AB the units are
    M^-1, not M -- the same slip as K_HA's in Sec. 13.5.
    """
    pH = np.asarray(pH, dtype=float)
    if convention == "printed":
        d = (1.0 + 10.0 ** (pH - pK_acid)) * (1.0 + 10.0 ** (pK_base - pH))
    elif convention == "figure":
        d = (1.0 + 10.0 ** (pK_acid - pH)) * (1.0 + 10.0 ** (pH - pK_base))
    else:
        raise ValueError(f"convention must be 'printed' or 'figure', not {convention!r}")
    return K / d


def pK_at_temperature(pK_1, dH_rxn, T1_C=25.0, T2_C=37.0, R=GAS_CONSTANT):
    """pK at `T2_C` from its value at `T1_C` and the heat of ionization, Sec. 13.6:

        pK(T2) = pK(T1) + [dH_rxn/(2.303 R)] (1/T2 - 1/T1)

    with the temperatures in kelvin. `dH_rxn` is in J/mol.

    For glycine's ammonium group (dH = -44.77 kJ/mol, pK = 9.60 at 25 C) this gives
    9.90 at 37 C and 8.88 at 0 C, as the section prints. An exothermic ionization has
    a pK that **rises** with temperature.
    """
    T1, T2 = float(T1_C) + 273.15, np.asarray(T2_C, dtype=float) + 273.15
    return pK_1 + dH_rxn / (np.log(10.0) * R) * (1.0 / T2 - 1.0 / T1)
