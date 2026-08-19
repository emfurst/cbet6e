"""osmotic -- osmotic equilibrium and osmotic pressure, SIS Sec. 11.5.

A membrane that passes the solvent and not the solute removes one equilibrium
condition and adds one unknown. The solvent still satisfies equality of fugacities,
Eq. 11.5-1; the solute does not, because it is held back mechanically rather than
thermodynamically. What is left after the Poynting correction (Eq. 11.5-3) is
SIS Eq. 11.5-4,

    Pi = P^II - P^I = -(RT / V_solvent) ln(x_solvent gamma_solvent)

and the rest of the section is that equation in three limits: the dilute limit that
turns it into a molecular weight (Eqs. 11.5-6, 11.5-7), the osmotic virial expansion
that turns the departure from that limit into a second virial coefficient
(Eq. 11.5-8), and the electrolyte case where complete dissociation doubles the solute
mole fraction (Illustration 11.5-4).

    from thermo.osmotic import pressure_from_height, virial_fit

    pressure_from_height(980.0, 0.0085)        # Ill. 11.5-1: 81.72 Pa
    virial_fit(C, Pi, 298.15).m_S              # Ill. 11.5-3: 29 500 g/mol

## Why the logarithm matters

Eq. 11.5-4 is the logarithm of a number very close to one, multiplied by RT/V, which
for water at 25 C is 1377 bar. A solution that is 99.44 mol % water -- physiological
saline -- has an osmotic pressure above 7 bar, and the fifth decimal place of the
solvent activity coefficient is worth 0.05 bar in the answer. That is the reason the
dilute forms below are written as separate functions rather than as one function with
an `ideal=True` flag: they are not small corrections to each other, and a notebook
that reports Pi to three figures has to say which one it used.

## Units

Every function that takes a solute mass concentration takes it in the units osmometry
is reported in, which are the units SIS prints: `C_S` in **g/L**, osmotic pressure in
Pa, molecular weight in g/mol, `B2` in L/mol and `B3` in L^2/mol^2. So R appears as
8314 Pa L/(mol K) rather than 8.314 in `osmotic_pressure_dilute`,
`molecular_weight`, `virial_fit` and `multicomponent_dilute`. `osmotic_pressure` and
`pressure_from_height` take no concentration and are SI throughout (m^3/mol, Pa,
kg/m^3, m).

Eric Furst
August 2026
"""
import numpy as np
from scipy import constants

R = constants.R
G_STANDARD = 9.81          # m/s^2, the value SIS uses in Illustration 11.5-1
MOLES_WATER_PER_KG = 55.51  # SIS uses this throughout Chapters 9, 11 and 15

__all__ = ["osmotic_pressure", "osmotic_pressure_dilute", "pressure_from_height",
           "molecular_weight", "virial_fit", "VirialFit",
           "multicomponent_dilute", "electrolyte_mole_fraction",
           "MOLES_WATER_PER_KG"]


def osmotic_pressure(x_solvent, V_solvent, T, gamma_solvent=1.0):
    """SIS Eq. 11.5-4, the osmotic pressure with no approximation beyond Poynting.

    `V_solvent` is the molar volume of the *solvent* in m^3/mol (1.8e-5 for water),
    taken independent of pressure -- which is the incompressibility assumption behind
    Eq. 11.5-3, not an extra one. Returns Pa.

    The solvent activity coefficient defaults to 1, which is the ideal-solution
    estimate the section works first. Illustration 11.5-4 is the argument for not
    stopping there: for 0.9 wt % aqueous sodium chloride, gamma_water differs from
    unity in the fourth decimal place and moves the answer by several percent.
    """
    x = np.asarray(x_solvent, dtype=float)
    g = np.asarray(gamma_solvent, dtype=float)
    return -R * T / V_solvent * np.log(x * g)


def osmotic_pressure_dilute(C_S, m_S, T):
    """SIS Eq. 11.5-6, Pi = RT C_S / m_S -- the van't Hoff limit.

    `C_S` is the solute mass concentration in g/L, `m_S` its molecular weight in
    g/mol; returns Pa. Two approximations get here from Eq. 11.5-4: ln(1 - z) is
    replaced by -z, and the solvent activity coefficient by 1. Both fail in the same
    direction, so the agreement of this equation with Eq. 11.5-4 at moderate
    concentration is not evidence that either approximation is good.
    """
    return (R * 1e3 * T * np.asarray(C_S, dtype=float)
            / np.asarray(m_S, dtype=float))


def pressure_from_height(rho, h, g=G_STANDARD):
    """Pi = rho g h, the osmometer of SIS Fig. 11.5-1.

    `rho` in kg/m^3, `h` in m; returns Pa. Illustration 11.5-1's 0.85 cm of
    cyclohexanone (density 980 kg/m^3) is 81.72 Pa, which is what makes osmometry
    practical: a 60 000 g/mol polymer at 2 g/L lifts a measurable centimeter of
    liquid where the same solution depresses the freezing point by microkelvin.
    """
    return np.asarray(rho, dtype=float) * g * np.asarray(h, dtype=float)


def molecular_weight(C_S, Pi, T):
    """SIS Eq. 11.5-7, m_S = RT C_S / Pi, in g/mol for `C_S` in g/L and `Pi` in Pa.

    The section is explicit that one measurement is not the intended use: "for high
    accuracy, this measurement is repeated several times at varying solute
    concentrations, and the limiting value of C_solute/Pi as C_solute approaches zero
    is used." `virial_fit` is that extrapolation, and it also returns what the single
    point cannot -- how far from the limit the measurement was.
    """
    return (R * 1e3 * T * np.asarray(C_S, dtype=float)
            / np.asarray(Pi, dtype=float))


class VirialFit:
    """The result of `virial_fit`: a molecular weight, virial coefficients, and errors.

    Attributes
    ----------
    m_S : molecular weight of the solute, g/mol, from the intercept.
    B2, B3 : osmotic virial coefficients, L/mol and L^2/mol^2. `B3` is None for a
        first- or second-order fit.
    coeffs : the polynomial coefficients of Pi/C_S in C_S, lowest order first, in
        Pa L/g, Pa L^2/g^2, ...
    stderr : standard errors of `coeffs`, from the least-squares covariance.
    m_S_stderr, B2_stderr : those errors propagated. Reported because they are the
        point: the intercept of an osmometry series is usually well determined and
        its slope usually is not, and a B2 quoted without them cannot be compared
        with another laboratory's.
    order : the polynomial order fitted.
    T : temperature, K.
    """

    def __init__(self, coeffs, stderr, T, order):
        self.coeffs = np.asarray(coeffs, dtype=float)
        self.stderr = np.asarray(stderr, dtype=float)
        self.T = float(T)
        self.order = int(order)
        RT = R * 1e3 * self.T                      # Pa L / mol
        a0 = self.coeffs[0]
        self.m_S = RT / a0
        self.m_S_stderr = self.m_S * self.stderr[0] / abs(a0)
        self.B2 = self.coeffs[1] * self.m_S ** 2 / RT if order >= 1 else None
        if order >= 1:
            # dB2/da1 at fixed m_S, plus the m_S contribution through a0
            rel = np.hypot(self.stderr[1] / abs(self.coeffs[1]),
                           2 * self.stderr[0] / abs(a0))
            self.B2_stderr = abs(self.B2) * rel
        else:
            self.B2_stderr = None
        self.B3 = (self.coeffs[2] * self.m_S ** 3 / RT) if order >= 2 else None

    def pi_over_C(self, C_S):
        """The fitted Pi/C_S at concentration `C_S` (g/L), in Pa L/g."""
        return np.polyval(self.coeffs[::-1], np.asarray(C_S, dtype=float))

    def pressure(self, C_S):
        """The fitted osmotic pressure at `C_S` (g/L), in Pa."""
        C = np.asarray(C_S, dtype=float)
        return C * self.pi_over_C(C)

    def rms(self, C_S, Pi):
        """RMS residual in Pi/C_S, Pa L/g."""
        y = np.asarray(Pi, dtype=float) / np.asarray(C_S, dtype=float)
        return float(np.sqrt(np.mean((y - self.pi_over_C(C_S)) ** 2)))

    def __repr__(self):
        b2 = "None" if self.B2 is None else f"{self.B2:.1f}"
        return (f"VirialFit(order={self.order}, m_S={self.m_S:.0f} g/mol, "
                f"B2={b2} L/mol)")


def virial_fit(C_S, Pi, T, order=2):
    """Fit the osmotic virial expansion, SIS Eq. 11.5-8, to (C_S, Pi) data.

    `C_S` in g/L, `Pi` in Pa. Written as

        Pi / C_S = (RT/m_S) [1 + B2 (C_S/m_S) + B3 (C_S/m_S)^2 + ...]

    the expansion is a polynomial in C_S whose *intercept* carries the molecular
    weight and whose higher coefficients carry the virial coefficients, so the fit is
    linear least squares in Pi/C_S and needs no starting guess. Returns a `VirialFit`.

    `order=2` is what Illustration 11.5-3 fits, and it is the smallest order that
    separates m_S from B2 -- a straight line through data that curve puts the
    curvature into the slope and biases both.

    Two properties of this fit bear on any B2 quoted from it.

    The intercept is an extrapolation to a concentration at which nothing was
    measured, and the data closest to it carry the least signal, since Pi itself goes
    to zero there. The order therefore moves B2 much more than it moves m_S: over
    Illustration 11.5-3's six points, going from order 1 to 2 to 3 moves m_S by 4 %
    and B2 by a factor of three. `VirialFit.B2_stderr` is the honest report of that,
    and a notebook comparing two data sets should quote it.

    The fit is unweighted in Pi/C_S, which weights the dilute points equally with the
    concentrated ones. Fitting Pi itself instead -- also least squares, also
    defensible -- weights the concentrated end and returns a different answer. Neither
    is more correct; they answer slightly different questions, and this module fits the
    quantity SIS plots.
    """
    C = np.asarray(C_S, dtype=float)
    P = np.asarray(Pi, dtype=float)
    if C.shape != P.shape:
        raise ValueError(f"got {C.size} concentrations and {P.size} pressures")
    if C.size < order + 1:
        raise ValueError(f"order {order} needs at least {order + 1} points, "
                         f"got {C.size}")
    if np.any(C <= 0):
        raise ValueError("concentrations must be positive")
    y = P / C
    coeffs, cov = np.polyfit(C, y, order, cov=True)
    return VirialFit(coeffs[::-1], np.sqrt(np.diag(cov))[::-1], T, order)


def multicomponent_dilute(C_S, m_S, T):
    """Osmotic pressure of a dilute multi-solute solution, Pi = RT SUM_i C_i/m_i.

    The equation SIS gives after Eq. 11.5-8, for the case where the solvent activity
    coefficient is unity and every virial coefficient is neglected. `C_S` and `m_S`
    are sequences in g/L and g/mol; returns Pa.

    Each solute contributes in proportion to its *molar* concentration, so a
    milligram of salt outweighs a milligram of protein by three orders of magnitude.
    That is why an isotonic fluid is specified by its salt content and not by its
    total dissolved solids.
    """
    C = np.asarray(C_S, dtype=float)
    m = np.asarray(m_S, dtype=float)
    return R * 1e3 * T * np.sum(C / m)


def electrolyte_mole_fraction(molality, nu, moles_solvent_per_kg=MOLES_WATER_PER_KG):
    """Solvent mole fraction in a completely dissociated electrolyte solution.

        x_S = n_S / (n_S + nu * M)

    with `nu` the number of ions per formula unit -- 2 for NaCl, 3 for CaCl2 -- as
    Illustration 11.5-4 and Appendix A9.3 both spell out. Getting `nu` wrong is the
    largest single error available in an osmotic-pressure calculation: for saline it
    is a factor of two in the answer, and it is a *counting* error rather than a
    thermodynamic one.
    """
    M = np.asarray(molality, dtype=float)
    return moles_solvent_per_kg / (moles_solvent_per_kg + nu * M)
