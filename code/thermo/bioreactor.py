"""bioreactor -- the fermenter atom balances, energy balance and second-law
constraint of SIS Sec. 15.7.

**Eleven illustrations, one object.** Illustrations 15.7-1 through 15.7-11 are the
same calculation asked eleven ways: write C, H, N and O balances for a reactor whose
streams are named only by their atom ratios, close them with an energy balance built
from heats of *combustion* rather than heats of formation, and bound the answer with
the second law. Written once, that object serves the whole section -- which is 48
percent of Chapter 15.

    from thermo.bioreactor import CMole, Fermentation, TABLE_15_7_2

    glucose = CMole("CH2O", name="glucose")
    biomass = CMole("CH1.8O0.5N0.2", name="average biomass")
    f = Fermentation(substrate=glucose, biomass=biomass,
                     nitrogen=CMole("NH3"), product=CMole("CH3O0.5"))
    f.solve(Y_B_S=0.14)          # Illustration 15.7-3

## What Sec. 15.7 asks for, and where each piece is

| SIS | what | entry point |
|-----|------|-------------|
| Eq. 15.7-9 | the generalized degree of reduction xi | `CMole.xi` |
| Eqs. 15.7-1 .. 15.7-4 | the four atom balances in yield-factor form | `Fermentation.residuals` |
| Eq. 15.7-12 | the 4C + H - 2O combination, i.e. the oxygen balance as xi | `Fermentation.Y_O2_from_xi` |
| Eq. 15.7-8b | the energy balance from heats of combustion | `Fermentation.heat_load` |
| Eq. 15.7-10 | the energy regularity principles | `regularity_G`, `regularity_H` |
| Eq. 15.7-11 | the energy balance with regularity for every species | `Fermentation.heat_load_regularity` |
| Eq. 15.7-19 | the second-law constraint on Gibbs energies of combustion | `Fermentation.second_law` |
| Eq. 15.7-20 | the same constraint under the regularity approximation | `Fermentation.second_law_regularity` |
| Eq. 15.7-26a | the entropy generated per C-mole of substrate | `Fermentation.entropy_generated` |
| Table 15.7-2 | xi, dcG, dcH for 52 compounds | `TABLE_15_7_2` |

## THE C-MOLE IS THE UNIT, AND IT IS NOT A MOLE

Everything in this module is per **C-mole**: one mole of the species divided by the
number of carbon atoms in it. Glucose C6H12O6 is CH2O, ethanol C2H5OH is CH3O0.5,
xylose C5H10O5 is CH2O. Species with no carbon -- oxygen, ammonia, water, molecular
nitrogen -- keep ordinary moles, and `CMole` carries a flag saying which it is. Mixing
the two is the arithmetic mistake this section invites: Illustration 15.7-6 as printed
calls ethanol's 684.5 kJ/C-mole a "kJ/mol", and Illustration 15.7-9's yield limit of
0.727 is C-moles per C-mole, which is 0.545 g/g only after both molar masses are put
back in.

## THE SECOND LAW IS AN INEQUALITY, SO IT CANNOT BE "SOLVED"

`second_law` returns a `Constraint`, not a boolean -- the two sides, their difference,
and the efficiency ratio -- because every use of Eq. 15.7-19 in Sec. 15.7 wants a
different one of those. Illustration 15.7-7 wants the ratio (85.1 percent),
Illustration 15.7-9 wants the equality case (the maximum yield), and Illustration
15.7-10 wants the constraint line drawn on a figure. A function that returned
`True` would have thrown away all three.

## WHY THE ATOM BALANCES ARE SOLVED AS A LINEAR SYSTEM AND NOT SUBSTITUTED

Illustrations 15.7-3, 15.7-4 and 15.7-11 each specify a *different* subset of the
yield factors and ask for the rest. Sec. 15.7 walks the substitutions by hand each
time, in a different order each time -- which is where the printed slips live
(Illustration 15.7-4's oxygen balance is set as `1 + 2*Y_W/S + Y_W/S`, and
Illustration 15.7-11's nitrogen balance as `Y_N/S = 0.15 * Y_N/S`). Here the four
balances are assembled as `A y = b` over whichever factors are unknown and solved
once, so the order never matters and an over- or under-specified problem **raises**
instead of quietly returning one arbitrary solution of many.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

GAS_CONSTANT = 8.314        # J/(mol K)

# Eq. 15.7-10, the "energy regularity principles" -- the averages of the last two
# columns of Table 15.7-2, over the organic compounds only.
REGULARITY_G = 112.0        # kJ per C-mole per unit xi
REGULARITY_H = 110.9        # kJ per C-mole per unit xi

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


def parse_formula(formula):
    """{element: count} from a formula string. Fractional counts are the point.

    `CH1.8O0.5N0.2` is the Roels average biomass and is not a molecule -- it is an
    elemental analysis normalized to one carbon. A formula parser that insists on
    integers cannot read this chapter.
    """
    counts = {}
    pos = 0
    for m in _TOKEN.finditer(formula):
        if m.start() != pos:
            raise ValueError(f"cannot parse {formula!r} at position {pos}")
        pos = m.end()
        el, n = m.group(1), m.group(2)
        counts[el] = counts.get(el, 0.0) + (float(n) if n else 1.0)
    if pos != len(formula):
        raise ValueError(f"cannot parse {formula!r} at position {pos}")
    return counts


def degree_of_reduction(C=0.0, H=0.0, O=0.0):
    """Eq. 15.7-9. xi = (4C + H - 2O)/C with carbon, xi = H - 2O without.

    **Nitrogen does not appear.** Sec. 15.7 says so in the parenthesis under
    Eq. 15.7-9 and it is easy to miss: this is one of several definitions of the
    generalized degree of reduction in the literature, and in *this* one xi(N2) = 0
    and xi(NH3) = 3. Using a definition in which ammonia's nitrogen is reduced --
    xi(NH3) = 0, with the nitrogen valence taken as -3 -- changes every second-law
    line in the section.
    """
    C, H, O = float(C), float(H), float(O)
    if C:
        return (4.0 * C + H - 2.0 * O) / C
    return H - 2.0 * O


def regularity_G(xi):
    """Eq. 15.7-10, Gibbs energy of combustion: 112 xi kJ per C-mole."""
    return REGULARITY_G * np.asarray(xi, float)


def regularity_H(xi):
    """Eq. 15.7-10, heat of combustion: 110.9 xi kJ per C-mole."""
    return REGULARITY_H * np.asarray(xi, float)


@dataclass(frozen=True)
class CMole:
    """One species of a fermentation, on the C-mole basis of Sec. 15.7.

    `formula` is the elemental analysis per C-mole for a carbon-containing species
    (`CH1.8O0.5N0.2`) and the ordinary molecular formula for one without (`NH3`,
    `H2O`, `O2`, `N2`). `dcG` and `dcH` are the measured Gibbs energy and heat of
    combustion in kJ per C-mole (or per mole, for a carbon-free species); leave them
    None to fall back on the energy regularity estimate, which is what Sec. 15.7 does
    for biomass.
    """
    formula: str
    name: str = ""
    dcG: float | None = None
    dcH: float | None = None
    _atoms: dict = field(init=False, repr=False, default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "_atoms", parse_formula(self.formula))
        if self.C not in (0.0, 1.0):
            raise ValueError(
                f"{self.formula!r} carries {self.C:g} carbon atoms. Sec. 15.7 works "
                f"per C-MOLE, so a carbon-containing species must be normalized to "
                f"exactly one carbon -- glucose C6H12O6 is CH2O, ethanol C2H5OH is "
                f"CH3O0.5. See the module docstring.")

    @property
    def C(self):
        return self._atoms.get("C", 0.0)

    @property
    def H(self):
        return self._atoms.get("H", 0.0)

    @property
    def O(self):
        return self._atoms.get("O", 0.0)

    @property
    def N(self):
        return self._atoms.get("N", 0.0)

    @property
    def xi(self):
        """Eq. 15.7-9."""
        return degree_of_reduction(self.C, self.H, self.O)

    def combustion_G(self, regularity=False):
        """dcG in kJ/C-mole. Falls back on Eq. 15.7-10 when no datum is carried."""
        if regularity or self.dcG is None:
            return float(regularity_G(self.xi))
        return float(self.dcG)

    def combustion_H(self, regularity=False):
        """dcH in kJ/C-mole. Falls back on Eq. 15.7-10 when no datum is carried."""
        if regularity or self.dcH is None:
            return float(regularity_H(self.xi))
        return float(self.dcH)

    def __repr__(self):
        tag = f" ({self.name})" if self.name else ""
        return f"CMole({self.formula}{tag}, xi = {self.xi:g})"


# --- Table 15.7-2 ----------------------------------------------------------------
#
# xi, dcG and dcH per C-mole, kJ. The two right-hand columns of the printed table
# (dcG/xi and dcH/xi) are DERIVED, so they are not stored -- they are regenerated in
# the notebook, which is what makes every row self-checking. Transcribing a derived
# column is how Table 14.6-1's four bad cells stayed invisible.
#
# The printed compound column reads `Trytophane`; the amino acid is tryptophan.
# The key here is the correct spelling and the printed form is recorded beside it.
TABLE_15_7_2 = {
    # name:               (xi,   dcG,    dcH)
    "acetic acid":        (4.00, 447.0, 438.0),
    "propionic acid":     (4.67, 511.0, 509.7),
    "butyric acid":       (5.00, 543.3, 548.5),
    "valeric acid":       (5.20, 562.6, 568.2),
    "palmitic acid":      (5.75, 612.5, 624.3),
    "lactic acid":        (4.00, 459.0, 456.3),
    "oxalic acid":        (1.00, 163.5, 123.0),
    "succinic acid":      (3.50, 399.8, 373.3),
    "fumaric acid":       (3.00, 362.0, 334.3),
    "malic acid":         (3.00, 361.0, 332.3),
    "citric acid":        (3.00, 357.8, 327.2),
    "glucose":            (4.00, 478.7, 467.8),
    "ethanol":            (6.00, 659.5, 684.5),
    "i-propanol":         (6.00, 648.7, 663.0),
    "n-butanol":          (6.00, 648.0, 670.0),
    "ethylene glycol":    (5.00, 585.0, 590.5),
    "glycerol":           (4.67, 547.7, 554.3),
    "glucitol":           (4.33, 514.0, 508.2),
    "acetone":            (5.33, 578.0, 597.7),
    "acetaldehyde":       (5.00, 561.5, 584.0),
    "alanine":            (5.00, 547.3, 569.0),
    "arginine":           (5.67, 631.0, 624.0),
    "asparagine":         (4.50, 499.8, 484.0),
    "glutamic acid":      (4.20, 463.0, 450.0),
    "aspartic acid":      (3.75, 421.5, 402.0),
    "glutamine":          (4.80, 525.6, 514.0),
    "glycine":            (4.50, 505.5, 487.0),
    "leucine":            (5.50, 594.2, 598.0),
    "isoleucine":         (5.50, 594.0, 598.0),
    "phenylalanine":      (4.78, 516.3, 517.0),
    "serine":             (4.33, 500.7, 485.0),
    "threonine":          (4.75, 532.5, 526.0),
    "tryptophan":         (4.73, 513.6, 512.0),      # printed `Trytophane`
    "tyrosine":           (4.56, 498.1, 493.0),
    "valine":             (5.40, 584.0, 584.0),
    "guanine":            (4.60, 522.4, 500.0),
    "n-hexane":           (6.33, 670.5, 693.9),
    "n-heptane":          (6.29, 665.7, 688.1),
    "n-octane":           (6.25, 662.2, 683.9),
    "n-decane":           (6.20, 657.1, 677.7),
    "n-dodecane":         (6.17, 653.7, 673.5),
    "n-hexadecane":       (6.13, 649.5, 667.9),
    "n-eicosane":         (6.10, 646.8, 663.9),
    "toluene":            (5.14, 547.7, 563.0),
    "cyclohexane":        (6.00, 636.1, 652.0),
    "ethyl acetate":      (5.00, 552.0, 563.6),
    # "Other compounds" -- these are per MOLE, not per C-mole, and none of them is
    # an organic, so none is in the average that produced Eq. 15.7-10's 112 and 110.9.
    "ammonia":            (3.00, 391.9, 348.1),      # l, 0.01 M
    "hydrogen":           (2.00, 238.0, 286.0),
    "carbon monoxide":    (2.00, 257.0, 283.0),
    "nitric acid":       (-5.00,   7.3, -30.0),
    "hydrazine":          (4.00, 602.4, 622.0),
    "hydrogen sulfide":   (2.00, 323.0, 247.0),
    "sulfuric acid":     (-6.00, -507.4, -602.0),
}

# The 46 organics above the rule in the printed table -- the rows Eq. 15.7-10's two
# constants are averaged over. The seven "other compounds" are excluded, which is what
# the printed table's `Average` line does and what makes 112.0/110.9 reproducible.
TABLE_15_7_2_ORGANIC = tuple(list(TABLE_15_7_2)[:46])


# The MOLECULAR formula of each row, which the printed table does not carry. It is what
# makes the xi column self-checking: Eq. 15.7-9 regenerates every one of the 53 values
# from the formula alone, so a transcription error in xi cannot survive.
TABLE_15_7_2_FORMULA = {
    "acetic acid": "C2H4O2", "propionic acid": "C3H6O2", "butyric acid": "C4H8O2",
    "valeric acid": "C5H10O2", "palmitic acid": "C16H32O2", "lactic acid": "C3H6O3",
    "oxalic acid": "C2H2O4", "succinic acid": "C4H6O4", "fumaric acid": "C4H4O4",
    "malic acid": "C4H6O5", "citric acid": "C6H8O7", "glucose": "C6H12O6",
    "ethanol": "C2H6O", "i-propanol": "C3H8O", "n-butanol": "C4H10O",
    "ethylene glycol": "C2H6O2", "glycerol": "C3H8O3", "glucitol": "C6H14O6",
    "acetone": "C3H6O", "acetaldehyde": "C2H4O", "alanine": "C3H7NO2",
    "arginine": "C6H14N4O2", "asparagine": "C4H8N2O3", "glutamic acid": "C5H9NO4",
    "aspartic acid": "C4H7NO4", "glutamine": "C5H10N2O3", "glycine": "C2H5NO2",
    "leucine": "C6H13NO2", "isoleucine": "C6H13NO2", "phenylalanine": "C9H11NO2",
    "serine": "C3H7NO3", "threonine": "C4H9NO3", "tryptophan": "C11H12N2O2",
    "tyrosine": "C9H11NO3", "valine": "C5H11NO2", "guanine": "C5H5N5O",
    "n-hexane": "C6H14", "n-heptane": "C7H16", "n-octane": "C8H18",
    "n-decane": "C10H22", "n-dodecane": "C12H26", "n-hexadecane": "C16H34",
    "n-eicosane": "C20H42", "toluene": "C7H8", "cyclohexane": "C6H12",
    "ethyl acetate": "C4H8O2", "ammonia": "NH3", "hydrogen": "H2",
    "carbon monoxide": "CO", "nitric acid": "HNO3", "hydrazine": "N2H4",
    "hydrogen sulfide": "H2S", "sulfuric acid": "H2SO4",
}


def c_mole_formula(formula):
    """Divide a molecular formula through by its carbon count. C6H12O6 -> CH2O.

    A species with no carbon is returned unchanged, which is Sec. 15.7's convention:
    its properties stay on a per-mole basis.
    """
    atoms = parse_formula(formula)
    for el in atoms:
        if el not in ("C", "H", "N", "O", "S"):
            raise ValueError(f"{formula}: element {el} is outside C/H/N/O/S")
    C = atoms.get("C", 0.0)
    if not C:
        return formula
    out = ""
    for el in ("C", "H", "N", "O", "S"):
        if el in atoms:
            n = atoms[el] / C
            out += el + ("" if abs(n - 1.0) < 1e-12 else f"{n:.6g}")
    return out


def from_table(name):
    """A `CMole` on Table 15.7-2's measured dcG and dcH, on the C-mole basis.

    **Name the species exactly.** There is no fuzzy match here, deliberately:
    `from_database('hydrogen')` returning hydrogen bromide drew two ch13 figures
    before anyone noticed. A miss raises and lists the near neighbors.

    The row's own xi is checked against Eq. 15.7-9 applied to the formula, and a
    disagreement **raises** rather than picking one.
    """
    key = name.strip().lower()
    if key not in TABLE_15_7_2:
        near = [k for k in TABLE_15_7_2 if key[:4] in k or k[:4] in key]
        raise KeyError(f"{name!r} is not in Table 15.7-2."
                       + (f" Did you mean one of {near}?" if near else ""))
    xi, dcG, dcH = TABLE_15_7_2[key]
    c = CMole(c_mole_formula(TABLE_15_7_2_FORMULA[key]), name=key, dcG=dcG, dcH=dcH)
    if abs(c.xi - xi) > 6e-3:
        raise ValueError(
            f"{name}: {TABLE_15_7_2_FORMULA[key]} gives xi = {c.xi:.4g}, Table 15.7-2 "
            f"prints {xi:g}. One of the two is wrong -- do not average them.")
    return c


@dataclass
class Constraint:
    """One side of Eq. 15.7-19 against the other, with what each caller wants."""
    available: float        # dcG of substrate + nitrogen source consumed
    required: float         # dcG of biomass + product produced
    label: str = ""

    @property
    def satisfied(self):
        return self.available >= self.required

    @property
    def slack(self):
        return self.available - self.required

    @property
    def efficiency(self):
        """The fraction of the feed's Gibbs energy that appears in the products."""
        return self.required / self.available

    def __repr__(self):
        rel = ">=" if self.satisfied else "<"
        return (f"Constraint({self.available:.4g} {rel} {self.required:.4g} kJ, "
                f"efficiency {100 * self.efficiency:.1f} %"
                + (f", {self.label}" if self.label else "") + ")")


# The yield factors the balances are written over, in a fixed order so that a caller
# can read a solution vector without guessing.
FACTORS = ("Y_B_S", "Y_P_S", "Y_C_S", "Y_N_S", "Y_W_S", "Y_O2_S")


@dataclass
class Fermentation:
    """Substrate + nitrogen source + O2 + water -> biomass + product + CO2.

    Every species is a `CMole`. `product` may be None (Illustration 15.7-10 produces
    only biomass) and `nitrogen` may be None (Illustration 15.7-9 produces no biomass
    and therefore needs no nitrogen source).
    """
    substrate: CMole
    biomass: CMole | None = None
    nitrogen: CMole | None = None
    product: CMole | None = None
    T: float = 298.15

    # --- the four atom balances, Eqs. 15.7-1 .. 15.7-4 ---------------------------
    def residuals(self, y):
        """C, H, N and O balance residuals at the yield factors `y` (a dict).

        Written so that all four are zero at the solution, in the printed order.
        `y` uses the `FACTORS` names.
        """
        S, B, P, N = self.substrate, self.biomass, self.product, self.nitrogen
        zero = CMole("H2")          # a stand-in with every count zero but H
        B = B or CMole("C")
        P = P or CMole("C")
        N = N or CMole("H2")
        YB, YP, YC = y["Y_B_S"], y["Y_P_S"], y["Y_C_S"]
        YN, YW, YO = y["Y_N_S"], y["Y_W_S"], y["Y_O2_S"]
        if self.biomass is None:
            YB = 0.0
        if self.product is None:
            YP = 0.0
        if self.nitrogen is None:
            YN = 0.0
        del zero
        return {
            # Eq. 15.7-1  (C_S = 1 by construction)
            "C": 1.0 + YN * N.C - YB * B.C - YP * P.C - YC,
            # Eq. 15.7-2
            "H": S.H + YN * N.H + 2.0 * YW - YB * B.H - YP * P.H,
            # Eq. 15.7-3
            "N": S.N + YN * N.N - YB * B.N - YP * P.N,
            # Eq. 15.7-4
            "O": S.O + YN * N.O + 2.0 * YO + YW - 2.0 * YC - YB * B.O - YP * P.O,
        }

    def solve(self, **known):
        """Solve the four balances for the yield factors not given in `known`.

        **Raises rather than guessing.** Four balances determine four unknowns; a
        problem that leaves five open is under-specified and one that leaves three is
        over-specified, and in both cases `numpy.linalg.solve` on a non-square system
        would either throw an unreadable error or -- worse, with `lstsq` -- return a
        least-squares answer that satisfies none of the balances exactly. Sec. 15.7
        never says out loud how many factors each illustration fixes, so the check
        has to be here.
        """
        for k in known:
            if k not in FACTORS:
                raise KeyError(f"{k!r} is not a yield factor; expected {FACTORS}")
        # Factors that are structurally zero because the species is absent.
        fixed = dict(known)
        if self.biomass is None:
            fixed.setdefault("Y_B_S", 0.0)
        if self.product is None:
            fixed.setdefault("Y_P_S", 0.0)
        if self.nitrogen is None:
            fixed.setdefault("Y_N_S", 0.0)
        unknown = [f for f in FACTORS if f not in fixed]
        if len(unknown) != 4:
            raise ValueError(
                f"four balances determine four unknowns; {len(unknown)} are open "
                f"({unknown}). Fix or free one and try again.")

        base = {f: fixed.get(f, 0.0) for f in FACTORS}
        r0 = self.residuals(base)
        order = ("C", "H", "N", "O")
        b = -np.array([r0[k] for k in order])
        A = np.empty((4, 4))
        for j, f in enumerate(unknown):
            probe = dict(base)
            probe[f] = base[f] + 1.0
            r1 = self.residuals(probe)
            A[:, j] = [r1[k] - r0[k] for k in order]
        if abs(np.linalg.det(A)) < 1e-12:
            raise ValueError(
                f"the balances are singular in {unknown} -- this set of unknowns is "
                f"not determined by C, H, N and O. (Illustration 15.7-2's substrate "
                f"has unknown H and O, which is exactly why Sec. 15.7 says only two "
                f"of its yield factors can be found.)")
        sol = np.linalg.solve(A, b)
        out = dict(base)
        out.update(dict(zip(unknown, sol)))
        res = self.residuals(out)
        worst = max(abs(v) for v in res.values())
        if worst > 1e-9:
            raise AssertionError(f"balances not closed: {res}")
        return out

    # --- Eq. 15.7-12, the 4C + H - 2O combination --------------------------------
    def Y_O2_from_xi(self, Y_B_S=0.0, Y_P_S=0.0, Y_N_S=0.0):
        """Eq. 15.7-12: Y_O2/S = (xi_S + xi_N Y_N/S - xi_B Y_B/S - xi_P Y_P/S)/4.

        The same number the oxygen balance gives, obtained without knowing the
        substrate's H and O separately. Sec. 15.7 derives it as 4 x (carbon balance)
        + (hydrogen balance) - 2 x (oxygen balance), and the multipliers are the
        valences -- the same numbers that define xi.
        """
        xi_N = self.nitrogen.xi if self.nitrogen else 0.0
        xi_B = self.biomass.xi if self.biomass else 0.0
        xi_P = self.product.xi if self.product else 0.0
        return (self.substrate.xi + xi_N * Y_N_S
                - xi_B * Y_B_S - xi_P * Y_P_S) / 4.0

    # --- Eq. 15.7-8b, the energy balance -----------------------------------------
    def heat_load(self, Y_B_S=0.0, Y_P_S=0.0, Y_N_S=0.0, regularity=False):
        """Y_Q/S in kJ per C-mole of substrate consumed -- Eq. 15.7-8b.

        Negative means heat must be REMOVED. Sec. 15.7 proves from the second law
        that it always is (Eq. 15.7-24 onward), so a positive return is a signal that
        the yield factors handed in are not thermodynamically possible.
        """
        h = lambda s: s.combustion_H(regularity) if s else 0.0
        return (Y_B_S * h(self.biomass) + Y_P_S * h(self.product)
                - h(self.substrate) - Y_N_S * h(self.nitrogen))

    def heat_load_regularity(self, Y_B_S=0.0, Y_P_S=0.0, Y_N_S=0.0):
        """Eq. 15.7-11 -- the energy balance with Eq. 15.7-10 used for EVERY species.

        Kept separate from `heat_load(regularity=True)` reading the same way on
        purpose: Illustration 15.7-6's Note compares the two, and the difference
        (-51.0 against -34.8 kJ) is almost entirely the regularity estimate's error
        on glucose, not on the biomass it was introduced for.
        """
        xi_B = self.biomass.xi if self.biomass else 0.0
        xi_P = self.product.xi if self.product else 0.0
        xi_N = self.nitrogen.xi if self.nitrogen else 0.0
        return REGULARITY_H * (Y_B_S * xi_B + Y_P_S * xi_P
                               - self.substrate.xi - Y_N_S * xi_N)

    # --- Eqs. 15.7-19 and 15.7-20, the second law --------------------------------
    def second_law(self, Y_B_S=0.0, Y_P_S=0.0, Y_N_S=0.0, regularity=False):
        """Eq. 15.7-19 as a `Constraint`."""
        g = lambda s: s.combustion_G(regularity) if s else 0.0
        return Constraint(
            available=g(self.substrate) + Y_N_S * g(self.nitrogen),
            required=Y_B_S * g(self.biomass) + Y_P_S * g(self.product),
            label="Eq. 15.7-19" + (" (energy regularity)" if regularity else ""))

    def second_law_regularity(self, Y_B_S=0.0, Y_P_S=0.0, Y_N_S=0.0):
        """Eq. 15.7-20, the constraint written on degrees of reduction alone."""
        xi_B = self.biomass.xi if self.biomass else 0.0
        xi_P = self.product.xi if self.product else 0.0
        xi_N = self.nitrogen.xi if self.nitrogen else 0.0
        return Constraint(
            available=self.substrate.xi + xi_N * Y_N_S,
            required=xi_B * Y_B_S + xi_P * Y_P_S,
            label="Eq. 15.7-20")

    def max_product_yield(self):
        """The equality case of Eq. 15.7-20 with no biomass: Y_P/S <= xi_S/xi_P.

        Illustration 15.7-9(b). Returns None when there is no product.
        """
        if self.product is None:
            return None
        return self.substrate.xi / self.product.xi

    def max_biomass_yield(self, regularity=False):
        """The equality case of the second law when biomass is the only product.

        Illustration 15.7-10(c). With the nitrogen balance Y_N/S = N_B * Y_B/S
        folded in, Eq. 15.7-19 becomes linear in Y_B/S and the maximum is closed
        form. It is NOT xi_S/xi_B: the nitrogen source carries Gibbs energy into
        the reactor too, which is what makes the printed coefficient 0.2439 rather
        than 1/4.8 = 0.2083.
        """
        if self.biomass is None:
            raise ValueError("no biomass to bound")
        gB = self.biomass.combustion_G(regularity)
        gN = self.nitrogen.combustion_G(regularity) if self.nitrogen else 0.0
        gS = self.substrate.combustion_G(regularity)
        denom = gB - self.biomass.N * gN
        if denom <= 0:
            raise ValueError(
                "the nitrogen source carries more Gibbs energy into the reactor than "
                "the biomass carries out; the second law places no upper bound here")
        return gS / denom

    # --- Eq. 15.7-26a, the entropy balance ---------------------------------------
    def entropy_generated(self, Y_B_S=0.0, Y_P_S=0.0, Y_N_S=0.0, T=None,
                          regularity=False):
        """S_gen per C-mole of substrate, kJ/(C-mol K) -- Eq. 15.7-26a."""
        c = self.second_law(Y_B_S, Y_P_S, Y_N_S, regularity)
        return c.slack / (self.T if T is None else T)


def metabolic_rate(oxygen_m3_per_day):
    """Sec. 15.7's mammalian correlation: 21 800 kJ per m3 of oxygen per day.

    The animal-kingdom analogue of Eq. 15.7-25's -443.6 kJ per mole of O2 in a
    fermenter, and the reason Fig. 15.7-5 belongs in this section rather than in a
    physiology text.
    """
    return 21_800.0 * np.asarray(oxygen_m3_per_day, float)
