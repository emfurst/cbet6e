"""partition -- the distribution of a solute between two liquid phases, SIS Sec. 11.4.

Section 11.4 is Sec. 11.2's equilibrium condition with one thing thrown away. The
solute still satisfies the equality of species fugacities, and because both phases are
liquids the pure-component fugacity still cancels, leaving SIS Eq. 11.4-4,

    x_1^I gamma_1(x^I, T) = x_1^II gamma_1(x^II, T)

What is thrown away is the *solvents'* equilibrium: the section assumes the mutual
solubility of the two solvents is unchanged by the solute, either because the solvents
are nearly immiscible or because so little solute was added. That assumption is the
only difference between this module and `thermo.lle`, and it buys a great deal --
Eq. 11.4-4 rearranges to Eq. 11.4-5,

    K_x = x_1^I / x_1^II = gamma_1^II / gamma_1^I

so a measured distribution coefficient *is* a ratio of activity coefficients. No solver
appears anywhere in this module. Everything here is that rearrangement, plus the unit
conversions that stand between a laboratory concentration and a mole fraction.

    from thermo.partition import x_from_concentration, kow_from_gamma

    x_from_concentration(2.5, 153.84, 1595.0, 159.83, 3119.0)   # Ill. 11.4-1
    kow_from_gamma(3.76e8)                                      # Ill. 11.4-2

## What Section 11.4 asks for, and where each piece is

| SIS | what | entry point |
|-----|------|-------------|
| Eq. 11.4-5 | gamma ratio from a distribution coefficient | `gamma_ratio_from_Kx` |
| Ill. 11.4-1, Eq. (1) | mole fraction from a molar concentration | `x_from_concentration` |
| Ill. 11.4-1 | Kc -> Kx for a solute in two solvents | `Kx_from_Kc` |
| Eqs. 11.4-12, 11.4-13 | K_OW from gamma_i^(W,infinity) | `kow_from_gamma` |
| Eq. 11.4-13, inverted | gamma_i^(W,infinity) from a measured K_OW | `gamma_from_kow` |
| Ill. 11.4-3 | how much solute ends up in each phase | `solute_split` |
| Ill. 11.4-4 | gamma ratio from two saturation solubilities | `gamma_ratio_from_solubility` |
| Eq. 11.4-17 | Gibbs energy of transfer between solvents | `gibbs_energy_of_transfer` |

## The one trap in the arithmetic, and it is a unit trap

A distribution coefficient reported on **concentrations** is not the distribution
coefficient that equals a ratio of activity coefficients; only the **mole fraction**
one is (Eq. 11.4-5). Converting between them needs the molar volumes of both solvents
*and* of the solute, because the volume of a solution of known molar concentration
depends on what is dissolved in it. Illustration 11.4-1 is a worked example of exactly
that conversion, and its two coefficients -- 0.096 45 and 0.045 21 for the carbon
tetrachloride phase, 0.018 and -0.033 24 for the aqueous phase -- come out of
`x_from_concentration`'s two arguments rather than being typed in.

Note the *sign* of the second one. The denominator of Eq. (1) is
1 + C(m_S/rho_S - m_B/rho_B), and whether that correction raises or lowers the mole
fraction depends on which of the solute and the solvent has the larger molar volume.
Bromine is denser and heavier than water but has the smaller molar volume, so the
aqueous coefficient is negative where the organic one is positive. A conversion written
with an assumed sign is wrong for one of the two phases in this illustration.

## Units

Concentrations and molar volumes are unit-agnostic and must simply agree: if C is in
kmol/m^3 then molecular weights are in kg/kmol and densities in kg/m^3, which is how
Illustration 11.4-1 prints them. `gibbs_energy_of_transfer` is the one function with
fixed units, J/mol, because it takes a temperature.

Eric Furst
August 2026
"""
import numpy as np
from scipy import constants

R = constants.R

__all__ = ["gamma_ratio_from_Kx", "x_from_concentration", "Kx_from_Kc",
           "kow_from_gamma", "gamma_from_kow", "solute_split",
           "gamma_ratio_from_solubility", "x_from_molarity",
           "gibbs_energy_of_transfer", "OCTANOL_WATER",
           "air_water_partition", "compartment_partition",
           "compartment_concentrations", "COMPARTMENTS", "EQ_12_5_8_CONSTANT"]

# n-octanol and water as Sec. 11.4 gives them, in the paragraph that derives
# Eq. 11.4-12: density in g/cc and molecular weight in g/mol, from which the ratio
# of total molar concentrations C^O/C^W = 0.114 follows.
OCTANOL_WATER = {"rho_octanol": 0.827, "mw_octanol": 130.22,
                 "rho_water": 1.0, "mw_water": 18.0,
                 "x_water_in_octanol": 0.26, "gamma_octanol_inf": 5.0}


def gamma_ratio_from_Kx(Kx):
    """gamma_1^II / gamma_1^I from the mole-fraction distribution coefficient.

    SIS Eq. 11.4-5, which is the identity function. It is here so that a notebook
    naming this quantity does not have to explain in a comment that the two are the
    same number; the equation number is the explanation.
    """
    return np.asarray(Kx, dtype=float)


def x_from_concentration(C, mw_solvent, rho_solvent, mw_solute, rho_solute):
    """Mole fraction of a solute from its molar concentration, no volume change on mixing.

    SIS Eq. (1) of Illustration 11.4-1,

        x_B = C_B (m_S/rho_S) / [1 + C_B (m_S/rho_S - m_B/rho_B)]

    which follows from V^mix = x_B V_B + (1 - x_B) V_S with V_i = m_i/rho_i, and from
    the statement that C_B V^mix is the number of moles of solute in one mole of
    solution. Units are the caller's, and must agree: C in kmol/m^3 with molecular
    weights in kg/kmol and densities in kg/m^3 is the illustration's own choice.

    The one modeling assumption is Delta_mix V = 0. It is not incidental -- without it
    the molar volume of the solution is not the mole-fraction average and there is no
    closed form at all.
    """
    C = np.asarray(C, dtype=float)
    v_s = mw_solvent / rho_solvent
    v_b = mw_solute / rho_solute
    return C * v_s / (1.0 + C * (v_s - v_b))


def Kx_from_Kc(C_I, Kc, solvent_I, solvent_II, solute):
    """Convert a concentration distribution coefficient to a mole-fraction one.

    Illustration 11.4-1's whole calculation. `C_I` is the solute concentration in
    phase I -- the entries of the printed table -- and `Kc = C_I/C_II`, so the
    concentration in phase II is C_I/Kc. Each of `solvent_I`, `solvent_II` and
    `solute` is a `(molecular weight, density)` pair in units consistent with `C_I`.

    Returns `(x_I, x_II, Kx)`. By Eq. 11.4-5 the third of those is also
    gamma_1^II/gamma_1^I, which is what the illustration is after.

    Both phases are converted with the *same* function and the same solute, which is
    the reason this is not two lines of arithmetic in a notebook: the aqueous
    conversion differs from the organic one only in which solvent is named, and
    writing it twice is how a sign or a molecular weight gets transposed.
    """
    x_I = x_from_concentration(C_I, solvent_I[0], solvent_I[1], solute[0], solute[1])
    x_II = x_from_concentration(np.asarray(C_I, dtype=float) / np.asarray(Kc, float),
                                solvent_II[0], solvent_II[1], solute[0], solute[1])
    return x_I, x_II, x_I / x_II


def kow_from_gamma(gamma_water_inf, equation="11.4-13"):
    """Octanol-water partition coefficient from the aqueous infinite-dilution gamma.

    Two correlations, and they are not the same kind of statement:

    `"11.4-12"`  K_OW = 0.0228 gamma_i^(W,inf), equivalently
                 log10 K_OW = -1.642 + log10 gamma_i^(W,inf).
                 Derived, not fitted: 0.0228 is (C^O/C^W)/gamma_i^(O,inf) = 0.114/5,
                 with C^O/C^W from the densities and molecular weights of the two
                 solvents and gamma_i^(O,inf) taken as "of order 5" because the
                 octanol-rich phase is itself organic. Its slope in log-log is
                 therefore exactly 1 by construction.

    `"11.4-13"`  log10 K_OW = -0.486 + 0.806 log10 gamma_i^(W,inf), which SIS gives as
                 the correlation of measured data. The slope is *not* 1, so the two
                 agree only near log10 gamma^(W,inf) = 6 and diverge on either side:
                 11.4-12 is high for a very hydrophobic species, which is the case
                 Illustration 11.4-2 works.

    Returns K_OW itself, not its logarithm.
    """
    g = np.asarray(gamma_water_inf, dtype=float)
    if equation in ("11.4-12", "1", 1):
        return 0.0228 * g
    if equation in ("11.4-13", "2", 2):
        return 10.0 ** (-0.486 + 0.806 * np.log10(g))
    raise ValueError("equation must be '11.4-12' or '11.4-13'")


def gamma_from_kow(Kow, equation="11.4-13"):
    """Invert `kow_from_gamma`: gamma_i^(W,inf) from a measured K_OW.

    The paragraph after Eq. 11.4-13 is explicit that this is the point of the
    correlation -- "knowing any one among the infinite-dilution activity coefficient,
    octanol-water partition coefficient, and saturation solubility in water, the other
    two can be estimated." The third leg of that triangle is `lle.gamma_from_solubility`.
    """
    K = np.asarray(Kow, dtype=float)
    if equation in ("11.4-12", "1", 1):
        return K / 0.0228
    if equation in ("11.4-13", "2", 2):
        return 10.0 ** ((np.log10(K) + 0.486) / 0.806)
    raise ValueError("equation must be '11.4-12' or '11.4-13'")


def solute_split(n_solute, V_I, V_II, K):
    """Distribute a fixed amount of solute between two liquid phases of known volume.

    Illustration 11.4-3. With K = C_I/C_II the concentration distribution coefficient
    and the solute conserved,

        n = C_II V_II + C_I V_I = C_II (V_II + K V_I)

    so `(C_I, C_II)` follows in one line. Units are the caller's; the illustration
    works in mol and mL.

    This is the calculation that decides whether a purification works, and the reason
    it is worth writing down rather than doing in one's head is that the answer depends
    on K and on the *volume ratio* together. A K of 65 still leaves 1.5 % of the solute
    behind when the phases have equal volume, which is why Problem 11.4-2 asks for the
    same extraction in several small batches instead of one large one.
    """
    C_II = np.asarray(n_solute, dtype=float) / (V_II + np.asarray(K, float) * V_I)
    return C_II * np.asarray(K, dtype=float), C_II


def x_from_molarity(molarity, moles_solvent_per_kg):
    """Mole fraction of a dilute solute from its molarity.

    Illustration 11.4-4's conversion: `moles_solvent_per_kg` is 55.51 for water and
    21.707 for ethanol, and the mole fraction is M/(M + that). Strictly this treats a
    molarity (per liter of solution) as a molality (per kg of solvent), which for the
    0.186 M and 2.3e-5 M solutions of that illustration is well inside the precision
    of the solubility data.
    """
    M = np.asarray(molarity, dtype=float)
    return M / (M + moles_solvent_per_kg)


def gamma_ratio_from_solubility(x_sat_1, x_sat_2):
    """gamma^inf in solvent 2 over gamma^inf in solvent 1 from saturation solubilities.

    Illustration 11.4-4. When the *same* solid is in equilibrium with two saturated
    liquid solutions, its fugacity is the same in both, so x_1 gamma_1 = x_2 gamma_2
    and the ratio of activity coefficients is the inverse ratio of the solubilities.
    The solid's own properties -- its fugacity, its melting point, its enthalpy of
    fusion -- cancel and never have to be known.

        gamma_2^inf / gamma_1^inf = x_1^sat / x_2^sat

    Both mole fractions must be small enough that the saturated solutions can stand in
    for infinitely dilute ones; the illustration says so out loud, at 3.34e-3 and
    1.06e-6.
    """
    return np.asarray(x_sat_1, dtype=float) / np.asarray(x_sat_2, dtype=float)


def gibbs_energy_of_transfer(gamma_inf_ratio, T):
    """Gibbs energy of transfer of one mole of solute between solvents, J/mol.

    SIS Eq. 11.4-17,

        Delta_tfr G = RT ln[ gamma_S^inf(solvent 2) / gamma_S^inf(solvent 1) ]

    for transfer *from* solvent 1 *to* solvent 2, at infinite dilution in both. A
    positive value means the transfer is unfavorable -- work must be supplied -- and a
    negative one means Gibbs energy is released.

    Only the ratio enters, which is why `gamma_ratio_from_solubility` is enough to
    evaluate it and why neither activity coefficient has to be known separately.
    """
    return R * T * np.log(np.asarray(gamma_inf_ratio, dtype=float))


# --- Sec. 12.5: the environmental compartments -------------------------------
#
# Section 12.5 is Sec. 11.4's argument pointed at a river instead of a separatory
# funnel. Every one of Eqs. 12.5-9 to 12.5-11 is "K_OW times the fraction of this
# compartment that behaves like octanol," and the whole model is that one idea:
#
#   biota     the lipid fraction IS octanol                     w_B  K_OW
#   soil      organic carbon is 40 % as good as octanol   0.4 * w_S  K_OW
#   sediment  the same, with a different organic fraction 0.4 * w_D  K_OW
#
# The 0.4 is empirical and Sec. 12.5 says so. It is not a thermodynamic result, and
# it is the single largest soft number in the section.

#: Sec. 12.5's own average organic/lipid fractions, and the densities Ill. 12.5-3 gives.
COMPARTMENTS = {"biota": dict(w=0.05, rho=1000.0, factor=1.0),
                "soil": dict(w=0.02, rho=1500.0, factor=0.4),
                "sediment": dict(w=0.05, rho=1420.0, factor=0.4)}

#: The constant of Eq. 12.5-8, K_AW = C_CONSTANT * gamma^inf P^vap / T, for P^vap in
#: bar and T in K. It is not fitted: 1.218e4 / 5.556e4 / 1.013, where 1.218e4 comes
#: from the ideal gas law (Eq. 12.5-6), 5.556e4 from 1e6/18 (Eq. 12.5-7), and 1.013 is
#: the atmospheric pressure at which Eq. 12.5-5 is evaluated.
EQ_12_5_8_CONSTANT = (1.218e4 / 5.556e4) / 1.013


def air_water_partition(gamma_inf, p_vap, T, constant=EQ_12_5_8_CONSTANT):
    """K_AW from an infinite-dilution activity coefficient -- SIS Eq. 12.5-8.

        K_AW,i = 0.2164 gamma_i^inf P_i^vap / T

    for P^vap in **bar** and T in K, giving a ratio of g/m^3 to g/m^3. The molecular
    weight cancels between Eqs. 12.5-6 and 12.5-7, which is why it does not appear.

    Watch the pressure unit. Illustration 12.5-1 states benzo[a]pyrene's vapor
    pressure in **pascals** (2.13e-5 Pa) and converts inside the equation; hand it bar.
    """
    return constant * np.asarray(gamma_inf, float) * np.asarray(p_vap, float) / np.asarray(T, float)


def compartment_partition(Kow, compartment):
    """K_BW, K_SW or K_DW from K_OW -- SIS Eqs. 12.5-9, 12.5-10 and 12.5-11.

        biota     K_BW = w_B K_OW
        soil      K_SW = 0.4 w_S K_OW
        sediment  K_DW = 0.4 w_D K_OW

    `compartment` is "biota", "soil" or "sediment"; the weight fractions and the 0.4
    are Sec. 12.5's own averages, in `COMPARTMENTS`.
    """
    if compartment not in COMPARTMENTS:
        raise KeyError(f"{compartment!r} is not one of {sorted(COMPARTMENTS)}")
    c = COMPARTMENTS[compartment]
    return c["factor"] * c["w"] * np.asarray(Kow, float)


def compartment_concentrations(C_water, Kow, K_aw, compartments=None):
    """Concentration of a pollutant in every compartment, given its concentration in water.

    Illustration 12.5-3 in one call. Returns a dict keyed by compartment name, each
    entry holding

        K          the partition coefficient against water
        per_mass   concentration per GRAM of the compartment, in whatever mass unit
                   C_water carried. With C_water in ng/m^3 this is ng/g, so it is ppb
                   by weight and ppm by weight is `per_mass / 1000`.
        per_volume concentration per CUBIC METER of the compartment, in the same unit
                   as C_water -- per_mass times the density in g/m^3

    **The two are not interchangeable and the illustration reports both.** A soil
    concentration of 0.247 ppm by weight is 3.71e8 ng/m^3, and the factor between them
    is the soil density -- 1500 kg/m^3, not 1000. Quoting one where the other is meant
    is a factor of 1.5 for soil and 1.42 for sediment on top of the unit conversion.
    Water is the reference, so 1 m^3 of water is 1e3 kg by definition here, which is the
    identity Illustration 12.5-3's closing parenthesis points out.

    Checked against Illustration 12.5-3's four printed answers: air 1.66, soil 3.71e8,
    sediment 8.78e8 and biota 1.55e9 ng/m^3, and its 0.247, 0.618 and 1.55 ppm.

    Air is not a partition against octanol and so is not in `compartments`; it comes
    from `K_aw` directly.
    """
    names = sorted(COMPARTMENTS) if compartments is None else list(compartments)
    C_water = np.asarray(C_water, float)
    out = {"air": dict(K=float(K_aw), per_mass=None,
                       per_volume=float(K_aw) * C_water)}
    for name in names:
        K = float(compartment_partition(Kow, name))
        # K is (mass of i per 1e6 g of compartment) / (mass of i per m^3 of water), so
        # dividing by 1e6 puts per_mass on a per-GRAM basis. Going from there to a
        # per-m^3 basis needs the density in g/m^3, which is rho[kg/m^3] x 1e3 -- and
        # forgetting that 1e3 is a factor-of-1000 error that still looks plausible.
        per_mass = K * C_water / 1e6            # per gram of compartment
        out[name] = dict(K=K, per_mass=per_mass,
                         per_volume=per_mass * COMPARTMENTS[name]["rho"] * 1e3)
    return out
