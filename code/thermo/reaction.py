"""reaction -- chemical reaction equilibrium, SIS Chapter 13.

⭐ **This module is the replacement for CHEMEQ** and for "the chemical equilibrium
constant calculation programs of Appendix B.I or B.II", the Visual Basic executables the
5e shipped on its website. Their whole numerical content was two screens of arithmetic:
sum the formation properties of Appendix A.IV against the stoichiometric coefficients,
sum the Appendix A.II heat-capacity coefficients the same way, and integrate the van 't
Hoff equation analytically. `Reaction.table` below is that program, and the 5e source was
read to establish it -- not guessed from the printed output.

    from thermo.reaction import Reaction

    rxn = Reaction.parse("0.5 N2 + 1.5 H2 = NH3")      # Illustrations 13.1-4, 13.1-8
    rxn.Ka(450.0)                                       # 1.218, as the book prints
    rxn.table(np.arange(500, 900, 100))                 # what CHEMEQ printed

    from thermo.reaction import equilibrium_extent
    equilibrium_extent(rxn, 450.0, {"N2": 0.5, "H2": 1.5}, P=4.0)

## What Chapter 13 asks for, and where each piece is

| SIS | what | entry point |
|-----|------|-------------|
| Eq. 13.1-5 | N_i = N_i,0 + nu_i X, and the mole fractions | `Reaction.moles`, `Reaction.balance_table` |
| Eq. 13.1-18 | Ka = exp(-dG_rxn(T)/RT) | `Reaction.Ka`, `Reaction.ln_Ka` |
| Eq. 13.1-19 | Ka = prod a_i^nu_i, the relation that is solved | `equilibrium_extent` |
| Sec. 13.1, after Eq. 13.1-19 | dG_rxn from the Appendix A.IV formation data | `Reaction.delta_G` |
| Eq. 13.1-20b | the van 't Hoff equation itself, d ln Ka/dT | `Reaction.dlnKa_dT` |
| Eq. 13.1-21 | dH_rxn(T) as an integral of dCp | `Reaction.delta_H` |
| Eq. 13.1-22a | the integrated van 't Hoff, as an integral | `Reaction.ln_Ka(mode="quadrature")` |
| Eq. 13.1-22b | the same, with dH_rxn constant | `Reaction.ln_Ka(mode="constant_H")` |
| the second 13.1-22a (see below) | dH_rxn(T), Cp polynomial carried through | `Reaction.delta_H` |
| the second 13.1-22b (see below) | ⭐ ln Ka(T), fully integrated -- **this is CHEMEQ** | `Reaction.ln_Ka` |
| Eqs. 13.1-22c, 13.1-22d | the same with dCp constant | `Reaction.ln_Ka(mode="constant_cp")` |
| Eqs. 13.1-23a to 23c | the ratios Kc, Kx, Ky, Kp | `Reaction.K_ratio` |
| Eq. 13.1-23d | Knu and Kgamma, the nonideality corrections | `K_nu_from_eos`, `K_gamma_from_model` |
| Table 13.1-1 | the mass balance table the chapter prints | `Reaction.balance_table` |
| Table 13.1-2 | species activity by choice of standard state | `activity` |
| Table 13.1-3 | Ka to Kc/Kx/Ky/Kp | `Reaction.K_ratio` |
| Sec. 13.1 at constant T and V | Illustration 13.1-4c, P floating with mole number | `equilibrium_extent(basis="TV")` |
| Sec. 13.2, Eq. 13.2-11 | dG_rxn(T) for heterogeneous reactions, the Ellingham line | `ellingham` |
| Sec. 13.3, Eqs. 13.3-4, 13.3-5 | several reactions at once, as coupled extents | `multireaction_extents` |
| Sec. 13.3 / Fig. 13.1-1 | the same state as a minimum of G | `gibbs_minimization`, `gibbs_curve` |

## The sign convention, and why it does not matter

Stoichiometric coefficients are **negative for reactants and positive for products** --
CHEMEQ's convention, and the book's. Sec. 13.1 makes the point that the choice is
arbitrary: replacing every nu_i by -nu_i interchanges reactants and products and leads to
the same equilibrium state, and multiplying the equilibrium relation by any constant
changes nothing. `Reaction.reversed` and `Reaction.scaled` exist so that a notebook can
demonstrate that, rather than assert it.

## ⛔ The trap: the equilibrium relation has roots outside the physical range

`equilibrium_extent` never takes a bare initial guess. The extent is bounded by
exhaustion -- no species may end up with a negative mole number -- which gives a closed
bracket

    X_lo = max over products  of (-N_i,0/nu_i)      X_hi = min over reactants of (N_i,0/-nu_i)

and the solve is a bracketed Brent on that interval. Inside it, for an ideal mixture at
fixed T and P, prod a_i^nu_i rises monotonically with X, so the root is unique; the
solver **checks that monotonicity** and says so in the result rather than trusting it.
Outside the bracket the same algebra has further roots that converge perfectly well and
mean nothing -- Illustration 13.1-5's cubic in X is the chapter's own example, where the
physical root is X = 0.0173 and the approximate analysis gives X proportional to
P^(-1/3). ⛔ On failure the solver **raises**; it does not return the last iterate.

## ⚠️ Two equations in Sec. 13.1 are printed with duplicate numbers

The labels **(13.1-22a)** and **(13.1-22b)** each appear twice in Sec. 13.1: once on the
general integrated van 't Hoff relation, and again about a page later on the pair that
carries the heat-capacity polynomial through the integration. The prose then refers to
the second pair as *"Eqs. 13.1-23a and b"* (twice -- in Illustration 13.1-3 and again in
Illustration 13.2-2), while 13.1-23a through 23e are the Kc/Kx/Ky/Knu/Kgamma group. One
of the two references has to move, and which one is the author's call, so the table above
names the second pair by position rather than by number. This is inherited from the 5e --
both labels are duplicated in the printed book, so it is not an artifact of the .docx.

## ⛔ Illustration 13.1-8's printed Ka table is not reproducible from Appendix A

Illustration 13.1-8 prints Ka for the ammonia reaction as 0.258, 0.02946, 0.005754 and
0.001595 at 500, 600, 700 and 800 K, and attributes them to "the programs in Appendix B.I
or B.II" -- that is, to CHEMEQ. They do not come out of CHEMEQ. This module reproduces
Illustration 13.1-3 to every printed digit and Illustration 13.1-4's Ka(450 K) = 1.218
exactly, from the same data and the same algorithm, and then gives 0.324, 0.0425, 0.00955
and 0.00304 for this table -- high by factors of 1.26 rising to 1.91. The printed values
imply a dH_rxn of -54 to -60 kJ/mol over that range where Appendix A.II gives -50 to -54,
and they are inconsistent with the book's own Ka(450 K): extrapolated back, the table
gives 1.06, not 1.218. The species data are not the explanation -- CHEMEQ's own database
was read and it is identical to `code/data/react_property.csv`, all 99 rows and all seven
columns. `code/ch13/validation/` records the comparison; the source of the printed table
is **not** established, and nothing here is tuned to match it.

## Units

J/mol and K throughout, as Appendix A.IV and A.II print them once scaled (`thermo.data`
does the scaling). Pressures are in **bar**, because the standard state of Chapter 13 is
1 bar and the activity of a low-pressure gas is then just y_i P -- the form every
illustration in the chapter uses. `K_nu_from_eos` is the one place a pressure crosses
into an equation of state, and it converts to Pa itself.

Eric Furst
August 2026
"""
import re
from collections import defaultdict

import numpy as np
from scipy import constants, integrate, optimize

from .data import get_reaction_species, reaction_cp

R = constants.R
T_REF = 298.15                  # the standard state of Appendix A.IV: 25 C, 1 bar
P_REF = 1.0                     # bar

__all__ = ["Reaction", "Extent", "equilibrium_extent", "multireaction_extents",
           "gibbs_minimization", "gibbs_curve", "ellingham", "activity",
           "K_nu_from_eos", "K_gamma_from_model", "elements", "molecular_weight",
           "formation_gibbs_T", "ATOMIC_MASS", "ELEMENT_REFERENCE_STATE",
           "PHASE_TAGS", "FORMULA_ALIASES"]


# --- formulas and elements ---------------------------------------------------
#
# ⚠️ A formula parser is a liability, and `electrolytes.py` deliberately does without
# one. It is unavoidable here: the element balance is what checks a reaction's
# stoichiometry, and Gibbs minimization is posed subject to element conservation. So the
# parser is written to be *checkable* rather than clever -- every one of the 99 names in
# react_property.csv is parsed and the result weighed against the molecular weight in the
# Reid-Prausnitz-Poling table, which is an independent column from an independent source.
# `code/ch13/validation/` runs that check. The names below are the ones no rule gets
# right, and they are listed rather than pattern-matched.

#: Tags that name a phase or a crystal form, not a group of atoms. ⚠️ `PbO(red)` and
#: `SiO2qrtz` are why this cannot be "strip anything in parentheses" -- and why
#: `Ca(OH)2` cannot be handled by stripping either.
PHASE_TAGS = ("(g)", "(l)", "(s)", "(aq)", "(red)", "(yellow)", "qrtz")

#: Names whose formula no rule recovers. `EtC6H6` is ethylbenzene, so `Et` is C2H5 and
#: the ring keeps its own hydrogens; the rest are isomer or position prefixes that carry
#: no atoms at all. Verified against the RPP molecular weights.
FORMULA_ALIASES = {
    "EtC6H6": "C8H10",          # ethylbenzene, not "Et" + benzene
    "CaOSiO2": "CaSiO3",        # calcium metasilicate, written as a double oxide
    "FeOSiO2": "FeSiO3",
    "CH4O": "CH4O",             # methanol, already atomic
    "(CH3)2CO": "C3H6O",        # acetone -- the leading group would trip the parser
}

#: The reference state of each element: the species in which its formation properties are
#: zero. Needed to turn a species into its formation reaction, which is how
#: `formation_gibbs_T` gets dG_f(T) for the Gibbs-minimization route.
ELEMENT_REFERENCE_STATE = {
    "Al": "Al", "Br": "Br2(g)", "C": "C", "Ca": "Ca", "Cl": "Cl2", "Cu": "Cu",
    "F": "F2", "Fe": "Fe", "H": "H2", "Hg": "Hg(l)", "I": "I2(g)", "Mg": "Mg",
    "N": "N2", "Na": "Na", "O": "O2", "Pb": "Pb", "S": "S", "Si": "Si",
}

#: Atomic masses (IUPAC 2021), for the molecular-weight check on the parser only.
ATOMIC_MASS = {
    "Al": 26.9815, "Br": 79.904, "C": 12.011, "Ca": 40.078, "Cl": 35.45,
    "Cu": 63.546, "F": 18.9984, "Fe": 55.845, "H": 1.008, "Hg": 200.592,
    "I": 126.904, "Mg": 24.305, "N": 14.007, "Na": 22.9898, "O": 15.999,
    "Pb": 207.2, "S": 32.06, "Si": 28.085,
}

_TOKEN = re.compile(r"([A-Z][a-z]?)([0-9]*)|(\()|(\))([0-9]*)")
_PREFIX = re.compile(r"^(?:[a-z]+[0-9]*-|[cs]s?[0-9]*-|[ni](?=[A-Z]))")


def _strip_decoration(name):
    """A species name reduced to a bare formula: phase tags and isomer prefixes gone."""
    formula = str(name)
    for tag in PHASE_TAGS:
        formula = formula.replace(tag, "")
    # "i-C4H10", "n-C4H10", "m-C8H10", "cs2-C4H8", "trs2-C4H8", "nC6H14", "iC3H8O".
    # ⚠️ The lookahead on the bare "n"/"i" form is what keeps "Na" intact.
    formula = _PREFIX.sub("", formula)
    return formula


def elements(name):
    """The atom counts of a species, as a dict -- 'C2H5OH' -> {'C': 2, 'H': 6, 'O': 1}.

    Phase tags and isomer prefixes are discarded, parenthesized groups are expanded, and
    the handful of names no rule recovers come from `FORMULA_ALIASES`. The counts are
    checked against the RPP molecular weights in `code/ch13/validation/`.
    """
    formula = FORMULA_ALIASES.get(str(name)) or _strip_decoration(name)
    counts, stack = defaultdict(float), []
    pos, current = 0, defaultdict(float)
    while pos < len(formula):
        m = _TOKEN.match(formula, pos)
        if m is None or m.end() == pos:
            raise ValueError(f"cannot parse the formula of species {name!r} "
                             f"(stopped at {formula[pos:]!r})")
        pos = m.end()
        symbol, count, opening, closing, multiplier = m.groups()
        if opening:
            stack.append(current)
            current = defaultdict(float)
        elif closing:
            if not stack:
                raise ValueError(f"unbalanced parentheses in {name!r}")
            factor = float(multiplier) if multiplier else 1.0
            outer = stack.pop()
            for k, v in current.items():
                outer[k] += v * factor
            current = outer
        else:
            if symbol not in ATOMIC_MASS:
                raise ValueError(f"{name!r} names element {symbol!r}, which is not one "
                                 f"of the elements in Appendix A.IV")
            current[symbol] += float(count) if count else 1.0
    if stack:
        raise ValueError(f"unbalanced parentheses in {name!r}")
    counts.update(current)
    return {k: v for k, v in counts.items() if v}


def molecular_weight(name):
    """Molecular weight from the parsed formula -- the check on `elements`, not a datum."""
    return sum(ATOMIC_MASS[el] * n for el, n in elements(name).items())


# --- the reaction -----------------------------------------------------------

class Reaction:
    """A single reaction, its standard-state properties, and its equilibrium constant.

    Parameters
    ----------
    stoichiometry : dict
        Species name -> stoichiometric coefficient, **negative for reactants**. Names
        are the `Name` column of `code/data/react_property.csv`, which is Appendix A.IV
        joined to Appendix A.II -- formulas, e.g. `'N2O4'`, `'H2O(g)'`, `'nC6H14(l)'`.
    name : str, optional
        A label for printing. Built from the stoichiometry when omitted.

    Notes
    -----
    Standard state: each species pure, at the temperature of interest and 1 bar, in the
    state of aggregation Appendix A.IV lists -- the choice Eq. 13.1-15 makes.
    """

    def __init__(self, stoichiometry, name=None):
        if not stoichiometry:
            raise ValueError("a reaction needs at least one species")
        self.nu = {str(k): float(v) for k, v in stoichiometry.items()}
        if any(v == 0.0 for v in self.nu.values()):
            raise ValueError("a stoichiometric coefficient of zero means the species "
                             "does not take part; leave it out")
        if all(v < 0 for v in self.nu.values()) or all(v > 0 for v in self.nu.values()):
            raise ValueError("every coefficient has the same sign, so this reaction has "
                             "no products (or no reactants); reactants are negative")
        self.species = list(self.nu)
        for s in self.species:                      # fail now, not inside a solver
            get_reaction_species(s)
        self.name = name or self.equation

    # -- construction --------------------------------------------------------

    @classmethod
    def parse(cls, equation, name=None):
        """A reaction from the way the book writes it: `'N2 + 3 H2 = 2 NH3'`.

        `=`, `==`, `->` and `=` with an arrow all separate the two sides. Coefficients
        may be written `3 H2`, `3H2` or omitted. ⚠️ `3H2` is ambiguous in principle -- no
        species name in Appendix A.IV begins with a digit, which is what makes it safe.
        """
        parts = re.split(r"<?=+>?|->|→|⇌|⇄", str(equation))
        if len(parts) != 2:
            raise ValueError(f"cannot find one reaction arrow in {equation!r}; write it "
                             f"as 'N2 + 3 H2 = 2 NH3'")
        stoich = {}
        for side, sign in zip(parts, (-1.0, +1.0)):
            for term in side.split("+"):
                term = term.strip()
                if not term:
                    continue
                m = re.match(r"^([0-9]*\.?[0-9]*)\s*(.+)$", term)
                coefficient, species = m.group(1), m.group(2).strip()
                factor = float(coefficient) if coefficient else 1.0
                stoich[species] = stoich.get(species, 0.0) + sign * factor
        return cls({k: v for k, v in stoich.items() if v != 0.0}, name=name)

    @property
    def equation(self):
        """The reaction as a string, reactants first."""
        def side(items):
            out = []
            for s, v in items:
                mag = abs(v)
                shown = "" if mag == 1.0 else (f"{mag:g} ")
                out.append(f"{shown}{s}")
            return " + ".join(out)
        left = side([(s, v) for s, v in self.nu.items() if v < 0])
        right = side([(s, v) for s, v in self.nu.items() if v > 0])
        return f"{left} = {right}"

    def reversed(self):
        """The same reaction written backwards -- and the same equilibrium state.

        Sec. 13.1's point that the choice of reactants and products is arbitrary. Use it
        to show that `Ka` is inverted while `equilibrium_extent` lands on the same
        composition, rather than taking the book's word for it.
        """
        return Reaction({s: -v for s, v in self.nu.items()})

    def scaled(self, factor):
        """The reaction multiplied through by a constant. Ka is raised to that power."""
        if factor == 0:
            raise ValueError("scaling a reaction by zero deletes it")
        return Reaction({s: v * float(factor) for s, v in self.nu.items()})

    # -- stoichiometry -------------------------------------------------------

    @property
    def delta_nu(self):
        """Sum of the stoichiometric coefficients -- the mole number change on reaction.

        Zero means pressure has no direct effect on the extent (Illustration 13.1-7);
        nonzero is the primary pressure effect of Illustration 13.1-4.
        """
        return float(sum(self.nu.values()))

    def element_balance(self):
        """Net atoms created by the reaction as written -- all zero if it is balanced.

        ⭐ This is the reaction's own check on itself, and it is the only thing standing
        between a mistyped coefficient and a converged, meaningless answer.
        """
        net = defaultdict(float)
        for s, v in self.nu.items():
            for el, n in elements(s).items():
                net[el] += v * n
        return {el: n for el, n in net.items() if abs(n) > 1e-9}

    def check_balance(self):
        """Raise unless the reaction conserves every element. Returns self, so it chains."""
        net = self.element_balance()
        if net:
            raise ValueError(f"{self.equation!r} does not balance: net "
                             + ", ".join(f"{n:+g} {el}" for el, n in sorted(net.items())))
        return self

    def moles(self, initial, X):
        """N_i = N_i,0 + nu_i X -- SIS Eq. 13.1-5. `initial` is a dict of mole numbers.

        Species absent from `initial` start at zero. Species in `initial` that take no
        part in the reaction are carried through unchanged -- that is how the diluent of
        Illustration 13.1-4b and Illustration 13.1-5 enters.
        """
        X = float(X)
        out = {s: float(n) for s, n in initial.items()}
        for s, v in self.nu.items():
            out[s] = out.get(s, 0.0) + v * X
        return out

    def mole_fractions(self, initial, X):
        """The mole fractions at extent X -- the last column of the balance table."""
        N = self.moles(initial, X)
        total = sum(N.values())
        if total <= 0:
            raise ValueError(f"no moles left at X = {X:g}")
        return {s: n / total for s, n in N.items()}

    def extent_bracket(self, initial, inset=1e-9, in_excess=()):
        """The physically possible range of X: no species may go negative.

        Returns `(X_lo, X_hi)`, pulled in by `inset` so that a log of zero cannot happen
        at an endpoint. ⛔ Everything about `equilibrium_extent`'s reliability rests on
        this bracket -- see the module docstring.

        `in_excess` names species that do not bound the extent because there is assumed
        to be plenty of them -- the condition Illustration 13.3-1 states outright for its
        solid carbon, *"as long as there is sufficient solid carbon present to ensure
        equilibrium."* A species is put in excess by marking it a solid phase and giving
        it no initial mole number.
        """
        in_excess = set(in_excess)
        lo, hi = -np.inf, np.inf
        for s, v in self.nu.items():
            if s in in_excess:
                continue
            N0 = float(initial.get(s, 0.0))
            if v < 0:
                hi = min(hi, N0 / -v)               # reactant runs out
            else:
                lo = max(lo, -N0 / v)               # product runs out going backwards
        if not np.isfinite(lo) or not np.isfinite(hi):
            raise ValueError("the initial state bounds the extent on one side only; "
                             "give a mole number for every species in the reaction, or "
                             "mark the unbounded one as a solid present in excess")
        if hi - lo <= 0:
            raise ValueError(f"the initial state leaves no room to react "
                             f"(bracket [{lo:g}, {hi:g}])")
        span = hi - lo
        return lo + inset * span, hi - inset * span

    def balance_table(self, initial, X=None):
        """The mass balance table Chapter 13 prints -- Table 13.1-1 and its siblings.

        Returns a DataFrame indexed by species with the initial mole number, the final
        mole number as a function of X, and the mole fraction. With `X` given, the two
        numeric columns are evaluated there as well; without it, the table is the
        symbolic one the book prints.
        """
        import pandas as pd

        rows, total_terms = {}, []
        for s in list(dict.fromkeys(list(initial) + self.species)):
            N0, v = float(initial.get(s, 0.0)), self.nu.get(s, 0.0)
            if v == 0:
                final = f"{N0:g}"
            elif N0 == 0:
                final = f"{v:g} X" if v != 1 else "X"
            else:
                final = f"{N0:g} {'+' if v > 0 else '-'} {abs(v):g} X"
                if abs(v) == 1:
                    final = f"{N0:g} {'+' if v > 0 else '-'} X"
            rows[s] = {"initial": N0, "final": final}
            total_terms.append((N0, v))
        N_total0 = sum(n for n, _ in total_terms)
        dnu = self.delta_nu
        total = f"{N_total0:g}" if dnu == 0 else \
                f"{N_total0:g} {'+' if dnu > 0 else '-'} {abs(dnu):g} X"
        for s, row in rows.items():
            row["mole fraction"] = f"({row['final']})/({total})"
        table = pd.DataFrame(rows).T
        table.loc["Total"] = {"initial": N_total0, "final": total, "mole fraction": "1"}
        if X is not None:
            N, y = self.moles(initial, X), self.mole_fractions(initial, X)
            table["N at X"] = [N.get(s, 0.0) for s in table.index[:-1]] + [sum(N.values())]
            table["y at X"] = [y.get(s, 0.0) for s in table.index[:-1]] + [1.0]
        return table

    # -- standard-state properties -------------------------------------------

    @property
    def delta_G_ref(self):
        """dG_rxn at 25 C, from the Appendix A.IV Gibbs energies of formation, in J/mol."""
        return float(sum(v * float(get_reaction_species(s)["DG"]) * 1e3
                         for s, v in self.nu.items()))

    @property
    def delta_H_ref(self):
        """dH_rxn at 25 C, from the Appendix A.IV enthalpies of formation, in J/mol."""
        return float(sum(v * float(get_reaction_species(s)["DH"]) * 1e3
                         for s, v in self.nu.items()))

    @property
    def delta_cp(self):
        """`(da, db, dc, dd, de)`: the Appendix A.II Cp coefficients summed with nu.

        Already scaled to SI by `thermo.data.reaction_cp`, so

            dCp = da + db T + dc T^2 + dd T^3 + de/T^2        J/(mol K)

        ⚠️ The A.II coefficients are *printed* pre-scaled by powers of ten; using the raw
        columns is wrong by three orders of magnitude and silently so. `thermo.data`
        is the only place that factor lives.
        """
        return sum(v * reaction_cp(s) for s, v in self.nu.items())

    def delta_cp_at(self, T):
        """dCp evaluated at T, in J/(mol K)."""
        a, b, c, d, e = self.delta_cp
        T = np.asarray(T, dtype=float)
        return a + b * T + c * T**2 + d * T**3 + e / T**2

    @property
    def _enthalpy_constant(self):
        """The constant in dH_rxn(T) = const + da T + db T^2/2 + ... -- Eq. 13.1-21.

        ⭐ Worth its own name because the chapter prints it twice for the same reaction
        and in two different guises: Illustration 13.1-3 gives it as 56 189 J/mol, and
        gives the coefficient of (1/T - 1/T1) in ln Ka as -6758.4, which is this constant
        over -R. Reproducing both from one expression is the check that this is right.
        """
        a, b, c, d, e = self.delta_cp
        T1 = T_REF
        return self.delta_H_ref - (a * T1 + b * T1**2 / 2 + c * T1**3 / 3
                                   + d * T1**4 / 4 - e / T1)

    def delta_H(self, T):
        """dH_rxn(T) in J/mol -- SIS Eq. 13.1-21, integrated.

            dH(T) = dH(T1) + da(T - T1) + db/2 (T^2 - T1^2) + dc/3 (T^3 - T1^3)
                            + dd/4 (T^4 - T1^4) - de (1/T - 1/T1)

        Illustration 13.1-3 prints this for N2O4 = 2 NO2 as
        56 189 + 12.80 T - 3.62e-2 T^2 + 1.434e-5 T^3 + 3.933e-9 T^4, which this
        reproduces to every digit shown.
        """
        a, b, c, d, e = self.delta_cp
        T = np.asarray(T, dtype=float)
        return (self._enthalpy_constant + a * T + b * T**2 / 2 + c * T**3 / 3
                + d * T**4 / 4 - e / T)

    def dlnKa_dT(self, T):
        """The van 't Hoff equation itself, d ln Ka/dT = dH_rxn(T)/(R T^2) -- Eq. 13.1-20b.

        Its sign is the chapter's rule of thumb: an exothermic reaction (dH < 0) has an
        equilibrium constant that falls with temperature.
        """
        T = np.asarray(T, dtype=float)
        return self.delta_H(T) / (R * T**2)

    def ln_Ka(self, T, mode="full"):
        """ln Ka at temperature T. Four routes, all from Eq. 13.1-20b.

        Parameters
        ----------
        T : float or array
        mode : {'full', 'constant_cp', 'constant_H', 'quadrature'}
            `'full'` -- the heat-capacity polynomial carried analytically through the
            integration. ⭐ **This is what CHEMEQ computed**, and it is the chapter's
            "General equation for the variation of the equilibrium constant with
            temperature" (printed as Eq. 13.1-22b, a duplicated label -- see the module
            docstring).

            `'constant_cp'` -- Eq. 13.1-22d, dCp frozen at its 25 C value.

            `'constant_H'` -- Eq. 13.1-22b (the first one of that number), dH_rxn frozen
            at its 25 C value. The straight line on a plot of ln Ka against 1/T.

            `'quadrature'` -- Eq. 13.1-22a integrated numerically instead of
            analytically. Returns the same numbers as `'full'`; it is here so a notebook
            can show that it does, which is the only honest way to claim the analytic
            form was integrated correctly.
        """
        T = np.asarray(T, dtype=float)
        if np.any(T <= 0):
            raise ValueError("temperature must be positive and in K")
        T1, lnKa_ref = T_REF, -self.delta_G_ref / (R * T_REF)
        a, b, c, d, e = self.delta_cp

        if mode == "full":
            poly = (a * np.log(T / T1) + b * (T - T1) / 2
                    + c * (T**2 - T1**2) / 6 + d * (T**3 - T1**3) / 12
                    + e * (1 / T**2 - 1 / T1**2) / 2)
            return lnKa_ref + (poly - self._enthalpy_constant * (1 / T - 1 / T1)) / R
        if mode == "constant_H":
            return lnKa_ref - self.delta_H_ref / R * (1 / T - 1 / T1)
        if mode == "constant_cp":
            dCp = float(self.delta_cp_at(T1))
            return (lnKa_ref - self.delta_H_ref / R * (1 / T - 1 / T1)
                    + dCp / R * np.log(T / T1) + dCp / R * (T1 / T - 1))
        if mode == "quadrature":
            scalar = np.ndim(T) == 0
            out = np.array([lnKa_ref + integrate.quad(self.dlnKa_dT, T1, float(t))[0]
                            for t in np.atleast_1d(T)])
            return float(out[0]) if scalar else out
        raise ValueError(f"mode must be 'full', 'constant_cp', 'constant_H' or "
                         f"'quadrature', got {mode!r}")

    def Ka(self, T, mode="full"):
        """The equilibrium constant -- SIS Eq. 13.1-18. See `ln_Ka` for the modes.

        ⚠️ Returns `inf` where ln Ka overflows. CHEMEQ printed ">10^38" in that case
        rather than a number, and a reaction that far to one side is better read from
        `ln_Ka` anyway.
        """
        return np.exp(self.ln_Ka(T, mode=mode))

    def delta_G(self, T, mode="full"):
        """dG_rxn(T) = -R T ln Ka(T), in J/mol.

        At 25 C this is the Appendix A.IV sum exactly; away from it, the temperature
        dependence comes through the same integration as Ka. Sec. 13.2's Eq. 13.2-11
        needs this function, not the 25 C number -- see `ellingham`.
        """
        T = np.asarray(T, dtype=float)
        return -R * T * self.ln_Ka(T, mode=mode)

    def delta_S(self, T, mode="full"):
        """dS_rxn(T) = (dH_rxn - dG_rxn)/T, in J/(mol K)."""
        T = np.asarray(T, dtype=float)
        return (self.delta_H(T) - self.delta_G(T, mode=mode)) / T

    def table(self, T, mode="full"):
        """⭐ What CHEMEQ printed: T, ln Ka, log10 Ka, Ka, dG_rxn, dH_rxn.

        The 5e program stepped 21 temperatures from a start and a step size; pass whatever
        array you want instead. dG and dH are in kJ/mol, as CHEMEQ printed them and as
        Appendix A.IV tabulates them.

            Reaction.parse("N2O4 = 2 NO2").table([300, 350, 400])
        """
        import pandas as pd

        T = np.atleast_1d(np.asarray(T, dtype=float))
        lnKa = np.atleast_1d(self.ln_Ka(T, mode=mode))
        return pd.DataFrame({"T (K)": T,
                             "ln Ka": lnKa,
                             "log10 Ka": lnKa / np.log(10.0),
                             "Ka": np.exp(lnKa),
                             "dG_rxn (kJ/mol)": self.delta_G(T, mode=mode) / 1e3,
                             "dH_rxn (kJ/mol)": np.atleast_1d(self.delta_H(T)) / 1e3})

    # -- the equilibrium ratios ----------------------------------------------

    def K_ratio(self, T, kind, P=None, C=None, K_nu=1.0, K_gamma=1.0, mode="full"):
        """Ka converted to one of the measurable ratios -- Table 13.1-3, Eqs. 13.1-23.

        Parameters
        ----------
        kind : {'Ky', 'Kp', 'Kx', 'Kc'}
            `Ky` and `Kx` are mole fraction ratios (vapor, liquid), `Kp` a
            partial-pressure ratio with units of pressure^dnu, `Kc` a concentration ratio
            with units of concentration^dnu.
        P : float
            Total pressure in bar. Needed for `Ky` and `Kp`.
        C : float
            Total molar concentration, for `Kc`, in whatever units the answer should
            carry.
        K_nu, K_gamma : float
            The nonideality corrections of Eq. 13.1-23d. Leave at 1 for an ideal mixture.

        Notes
        -----
        ⚠️ Ka depends only on temperature and on the choice of standard state. These
        ratios do not: they move with pressure, with concentration and with the mixture
        nonidealities, which is why the chapter warns that a literature Kc "has meaning
        only in the situation in which it was obtained."
        """
        Ka, dnu = self.Ka(T, mode=mode), self.delta_nu
        if kind == "Ky":
            if P is None:
                raise ValueError("Ky needs the total pressure P, in bar")
            return Ka / K_nu * (P / P_REF) ** -dnu
        if kind == "Kp":
            if P is None:
                raise ValueError("Kp needs the total pressure P, in bar")
            return Ka / K_nu * P_REF ** dnu
        if kind == "Kx":
            return Ka / K_gamma
        if kind == "Kc":
            if C is None:
                raise ValueError("Kc needs the total molar concentration C")
            # ⚠️ Watch this exponent. Table 13.1-3 prints the liquid relation as
            #     Ka = C^(-sum nu) Kc Kgamma
            # so Kc = Ka C^(+dnu)/Kgamma, with C RAISED to +dnu, not lowered. This had
            # the sign inverted until it was checked against the printed table: x_i =
            # C_i/C makes Kx = Kc C^(-dnu), and it is easy to carry that minus sign one
            # step too far. Illustration 13.1-6 exercises only Ka, Kp and Ky, so nothing
            # in the chapter's worked numbers catches it -- the identity in the
            # validation notebook does.
            return Ka / K_gamma * C ** dnu
        raise ValueError(f"kind must be 'Ky', 'Kp', 'Kx' or 'Kc', got {kind!r}")


# --- activities -------------------------------------------------------------

def activity(y, P, phase="gas", phi=1.0, gamma=1.0, molality=None):
    """The activity of one species -- Table 13.1-2, the low-to-moderate-pressure column.

    gas    : a = phi y P / 1 bar, and a = y P / 1 bar when phi = 1
    liquid : a = gamma x
    solid  : a = 1 for a pure solid or pure liquid (Sec. 13.2's whole simplification)
    molal  : a = gamma M / 1 molal

    Kept as a function rather than folded into the solver because Chapter 13's standard
    states are the part a reader gets wrong, and a notebook should be able to show the
    four side by side.
    """
    if phase == "gas":
        return phi * y * P / P_REF
    if phase == "liquid":
        return gamma * y
    if phase in ("solid", "pure"):
        return 1.0
    if phase == "molal":
        if molality is None:
            raise ValueError("the 1-molal standard state needs a molality")
        return gamma * molality
    raise ValueError(f"phase must be 'gas', 'liquid', 'solid' or 'molal', got {phase!r}")


def K_nu_from_eos(reaction, T, P, y, eos, order, phase="vapor"):
    """Knu = prod phi_i^nu_i from an equation of state -- Eq. 13.1-23d.

    Parameters
    ----------
    eos : PRMixture
    order : sequence of str
        The reaction's species names, in the order the EOS components were built in.
        ⚠️ Required, and not inferred: `react_property.csv` names species by formula
        (`'H2O(g)'`) and `pure_property.csv` by name (`'water'`), so there is no reliable
        way to line the two lists up automatically. Getting the order wrong gives a
        converged, wrong Knu, so the caller states it.
    y : dict
        Mole fractions by reaction species name.
    P : float
        Pressure in bar; converted to Pa here, since the EOS works in SI.

    Notes
    -----
    This is the quantity that makes Illustrations 13.1-7b and 13.1-8 iterative: phi_i
    needs the composition, and the composition is what the equilibrium relation is being
    solved for. `equilibrium_extent` does that outer loop; this is one evaluation of it.
    """
    order = list(order)
    if len(order) != getattr(eos, "n", len(order)):
        raise ValueError(f"the EOS has {eos.n} components but `order` names "
                         f"{len(order)}")
    missing = [s for s in reaction.nu if s not in order]
    if missing:
        raise KeyError(f"`order` does not name {missing}, which take part in the reaction")
    x = np.array([y[s] for s in order], dtype=float)
    phi = np.asarray(eos.phi(x / x.sum(), T, P * 1e5, phase=phase), dtype=float)
    lookup = dict(zip(order, phi))
    return float(np.prod([lookup[s] ** v for s, v in reaction.nu.items()]))


def K_gamma_from_model(reaction, x, model, order=None):
    """Kgamma = prod gamma_i^nu_i from an activity coefficient model -- Eq. 13.1-23d.

    `model` is any of the Chapter 9 models in `thermo.activity_models`; it is called as
    `model.gamma(x_array, T)` if that signature works and `model.gamma(x_array)`
    otherwise, because the models in the package differ on whether T is an argument.
    """
    order = list(order or x)
    xa = np.array([x[s] for s in order], dtype=float)
    xa = xa / xa.sum()
    try:
        gam = np.asarray(model.gamma(xa), dtype=float)
    except TypeError:
        gam = np.asarray(model.gamma(xa, None), dtype=float)
    lookup = dict(zip(order, gam))
    missing = [s for s in reaction.nu if s not in lookup]
    if missing:
        raise KeyError(f"the activity model does not cover {missing}")
    return float(np.prod([lookup[s] ** v for s, v in reaction.nu.items()]))


# --- the extent of reaction -------------------------------------------------

class Extent:
    """The result of an equilibrium-extent solve: the answer, and how it was reached.

    Attributes
    ----------
    X : float
        The equilibrium molar extent of reaction.
    y : dict
        Mole fractions at equilibrium.
    moles : dict
        Mole numbers at equilibrium.
    P : float
        Pressure at equilibrium, in bar. Equal to the pressure given for a constant-T,P
        solve, and floating for a constant-T,V one (Illustration 13.1-4c).
    Ka : float
        The equilibrium constant used.
    bracket : tuple
        The physical bracket the root was found in.
    monotone : bool
        Whether prod a_i^nu_i rose monotonically across the bracket. ⛔ `False` means the
        root may not be unique and the answer needs looking at, not trusting.
    history : list
        The outer iterates on Knu (or Kgamma). One entry for a solve with no nonideality.
    K_nu, K_gamma : float
        The converged nonideality corrections.
    """

    def __init__(self, X, y, moles, P, Ka, bracket, monotone, history,
                 K_nu=1.0, K_gamma=1.0):
        self.X = float(X)
        self.y = dict(y)
        self.moles = dict(moles)
        self.P = float(P)
        self.Ka = float(Ka)
        self.bracket = tuple(bracket)
        self.monotone = bool(monotone)
        self.history = list(map(float, history))
        self.K_nu = float(K_nu)
        self.K_gamma = float(K_gamma)

    @property
    def iterations(self):
        return len(self.history) - 1

    @property
    def direction(self):
        """'rising', 'falling' or 'non-monotone' for the outer iteration on X."""
        if len(self.history) < 2:
            return "single pass"
        d = np.diff(self.history)
        return ("rising" if np.all(d > 0) else
                "falling" if np.all(d < 0) else "non-monotone")

    def __float__(self):
        return self.X

    def __repr__(self):
        return (f"Extent(X={self.X:.6g}, P={self.P:.6g} bar, Ka={self.Ka:.6g}, "
                f"K_nu={self.K_nu:.6g}, iterations={self.iterations}, "
                f"monotone={self.monotone})")


def equilibrium_extent(reaction, T, initial, P=None, basis="TP", Ka=None,
                       K_nu=1.0, K_gamma=1.0, phases=None, mode="full",
                       tol=1e-10, max_iter=50):
    """Solve Ka = prod a_i^nu_i for the extent of reaction -- SIS Eq. 13.1-19.

    Parameters
    ----------
    reaction : Reaction
    T : float
        Temperature in K.
    initial : dict
        Initial mole numbers, including any inert diluent.
    P : float
        Pressure in bar. With `basis='TP'` it is the reaction pressure; with
        `basis='TV'` it is the *initial* pressure, and the answer's `P` is the final one.
    basis : {'TP', 'TV'}
        Constant temperature and pressure, or constant temperature and volume. The
        chapter's Illustration 13.1-4 does the same reaction both ways and gets different
        answers -- parts (a) and (c) -- which is the point of having the option.
    Ka : float, optional
        Use this equilibrium constant instead of computing one from Appendix A. Needed
        when the book quotes an experimental value, as in Illustrations 13.1-5 and
        13.1-7.
    K_nu : float or callable
        Gas-phase nonideality, Eq. 13.1-23d. A callable is treated as
        `K_nu(y, P, T) -> float` and the solve becomes the iteration Illustrations 13.1-7b
        and 13.1-8 describe: assume Knu, solve for the composition, recompute Knu, repeat.
    K_gamma : float or callable
        Liquid-phase nonideality, same protocol as `K_nu`, called as `K_gamma(x, T)`.
    phases : dict, optional
        Species name -> `'gas'`, `'liquid'`, `'solid'` or `'molal'`. Anything marked
        `'solid'` or `'pure'` has unit activity and drops out of the product, which is
        Sec. 13.2's simplification for heterogeneous reactions. Default: all gas.
    mode : str
        Passed to `Reaction.ln_Ka`.

    Returns
    -------
    Extent

    Raises
    ------
    ValueError
        If the bracket holds no root -- which means the reaction is driven to one
        boundary and the answer is exhaustion, not equilibrium -- or if the outer
        iteration on Knu does not converge. ⛔ It does not return the last iterate.
    """
    if basis not in ("TP", "TV"):
        raise ValueError(f"basis must be 'TP' or 'TV', got {basis!r}")
    if P is None:
        raise ValueError("give the pressure P in bar")
    Ka_value = float(reaction.Ka(T, mode=mode)) if Ka is None else float(Ka)
    if not np.isfinite(Ka_value):
        raise ValueError("Ka overflows at this temperature; the reaction goes to "
                         "completion and the extent is set by exhaustion, not by "
                         "equilibrium (see Sec. 13.1 on reactions that go to completion)")
    phases = dict(phases or {})
    # ⭐ A condensed pure phase is not part of the fluid mixture: it has unit activity, it
    # does not dilute anything, and it must stay out of the mole fraction denominator.
    # Illustration 13.3-1 says so for its solid carbon, and Aspen getting this wrong is
    # the reason that illustration's own note tells the reader to add a separator.
    condensed = {s for s, kind in phases.items() if kind in ("solid", "pure")}
    in_excess = {s for s in condensed if s not in initial}
    fluid = [s for s in dict.fromkeys(list(initial) + reaction.species)
             if s not in condensed]
    N_fluid0 = sum(float(initial.get(s, 0.0)) for s in fluid)
    if N_fluid0 <= 0:
        raise ValueError("the fluid phase is empty at the start, so there are no mole "
                         "fractions to solve for; give a feed of at least one fluid "
                         "species (an inert diluent will do)")
    lo, hi = reaction.extent_bracket(initial, in_excess=in_excess)

    def state(X):
        """Fluid-phase mole fractions and pressure at extent X."""
        N = reaction.moles(initial, X)
        total = sum(N[s] for s in fluid if s in N)
        if total <= 0:
            raise ValueError(f"no fluid moles left at X = {X:g}")
        y = {s: N[s] / total for s in fluid if s in N}
        if basis == "TV":
            return y, P * total / N_fluid0, N          # ideal gas at fixed T and V
        return y, P, N

    def log_product(X, knu, kgam):
        """ln prod a_i^nu_i at extent X, with the nonidealities held fixed."""
        y, P_x, _ = state(X)
        total = 0.0
        for s, v in reaction.nu.items():
            kind = phases.get(s, "gas")
            if kind in ("solid", "pure"):
                continue                               # unit activity
            if kind == "gas":
                a = y[s] * P_x / P_REF
            else:
                a = y[s]                               # gamma carried in kgam
            if a <= 0:
                return -np.inf if v > 0 else np.inf
            total += v * np.log(a)
        return total + np.log(knu) + np.log(kgam)

    history, knu, kgam = [], 1.0, 1.0
    if not callable(K_nu):
        knu = float(K_nu)
    if not callable(K_gamma):
        kgam = float(K_gamma)

    X = None
    for iteration in range(max_iter):
        def residual(X_try):
            return log_product(X_try, knu, kgam) - np.log(Ka_value)

        f_lo, f_hi = residual(lo), residual(hi)
        if not (np.isfinite(f_lo) or np.isfinite(f_hi)):
            raise ValueError("the equilibrium relation is undefined across the whole "
                             "bracket; check the initial mole numbers")
        if np.sign(f_lo) == np.sign(f_hi):
            raise ValueError(
                f"no root in the physical bracket [{lo:.6g}, {hi:.6g}]: the residual is "
                f"{f_lo:+.4g} at one end and {f_hi:+.4g} at the other. Ka = {Ka_value:.6g} "
                f"drives this reaction to a boundary, so the answer is exhaustion of a "
                f"species rather than an interior equilibrium.")
        X = optimize.brentq(residual, lo, hi, xtol=tol, rtol=1e-14)
        history.append(X)
        y, P_x, _ = state(X)
        if not (callable(K_nu) or callable(K_gamma)):
            break
        knu_new = float(K_nu(y, P_x, T)) if callable(K_nu) else knu
        kgam_new = float(K_gamma(y, T)) if callable(K_gamma) else kgam
        if (abs(knu_new - knu) <= 1e-10 * max(1.0, abs(knu))
                and abs(kgam_new - kgam) <= 1e-10 * max(1.0, abs(kgam))
                and iteration > 0):
            knu, kgam = knu_new, kgam_new
            break
        knu, kgam = knu_new, kgam_new
    else:
        raise ValueError(f"the nonideality iteration did not converge in {max_iter} "
                         f"passes; the extents so far are {history}. A fixed point that "
                         f"has not converged is not an answer.")

    # ⛔ Uniqueness is not assumed. For an ideal mixture at fixed T and P the product of
    # activities rises monotonically with X, so a sign change means exactly one root;
    # record whether that actually held rather than relying on it.
    probe = np.linspace(lo, hi, 65)
    values = np.array([log_product(x, knu, kgam) for x in probe])
    finite = np.isfinite(values)
    monotone = bool(np.all(np.diff(values[finite]) > 0))

    y, P_x, N = state(X)
    return Extent(X, y, N, P_x, Ka_value, (lo, hi), monotone, history, knu, kgam)


# --- several reactions at once ----------------------------------------------

def multireaction_extents(reactions, T, initial, P=None, basis="TP", Ka=None,
                          phases=None, mode="full", guess=None, tol=1e-12):
    """Sec. 13.3: R reactions in one phase, solved as R coupled extents.

    Each reaction contributes one equilibrium relation, Ka_j = prod_i a_i^nu_ij, and the
    mole numbers are N_i = N_i,0 + sum_j nu_ij X_j (Eq. 13.3-1). Returns a dict with the
    extents, the composition and the residuals.

    ⚠️ The coupled set is where a bare Newton solve earns its bad reputation, on two
    counts. The residuals are logarithms, so an iterate that pushes any mole number
    negative makes them undefined rather than merely wrong -- so the feasible region
    `N_i,0 + sum_j nu_ij X_j >= 0` is imposed as a constraint, not hoped for. And the
    obvious starting point, all extents zero, is **degenerate**: with no product present
    yet every log is minus infinity. So the default start is not zero but a strictly
    interior point of the feasible polytope, found by a small linear program that
    maximizes the smallest mole number.

    ⛔ If the residuals do not come down to `tol` it raises, and it also raises on a
    converged root that violates feasibility -- see the comment at the end of this
    function for the book's own statement of that trap. Compare `gibbs_minimization`,
    which is posed as a constrained minimization and is the more robust route.
    """
    reactions = list(reactions)
    if not reactions:
        raise ValueError("give at least one reaction")
    phases = dict(phases or {})
    if P is None:
        raise ValueError("give the pressure P in bar")
    Ka_values = np.array([float(r.Ka(T, mode=mode)) for r in reactions] if Ka is None
                         else [float(k) for k in Ka], dtype=float)
    if len(Ka_values) != len(reactions):
        raise ValueError("give one Ka per reaction")
    species = list(dict.fromkeys(list(initial)
                                 + [s for r in reactions for s in r.species]))
    condensed = {s for s, kind in phases.items() if kind in ("solid", "pure")}
    fluid_mask = np.array([s not in condensed for s in species])
    N0 = np.array([float(initial.get(s, 0.0)) for s in species])
    nu = np.array([[r.nu.get(s, 0.0) for s in species] for r in reactions])
    N_fluid0 = N0[fluid_mask].sum()
    if N_fluid0 <= 0:
        raise ValueError("the fluid phase is empty at the start")

    def composition(X):
        N = N0 + nu.T @ np.asarray(X, dtype=float)
        total = N[fluid_mask].clip(0.0, None).sum()
        safe = np.clip(N, 1e-300, None)            # keep the logs defined
        P_x = P * total / N_fluid0 if basis == "TV" else P
        return N, safe / max(total, 1e-300), P_x

    def residuals(X):
        _, y, P_x = composition(X)
        out = []
        for r, Ka_j in zip(reactions, Ka_values):
            acc = 0.0
            for s, v in r.nu.items():
                kind = phases.get(s, "gas")
                if kind in ("solid", "pure"):
                    continue
                a = y[species.index(s)] * (P_x / P_REF if kind == "gas" else 1.0)
                acc += v * np.log(max(a, 1e-300))
            out.append(acc - np.log(Ka_j))
        return np.array(out)

    nu_fluid = nu[:, fluid_mask]                   # R x (fluid species)
    N0_fluid = N0[fluid_mask]

    def interior_start():
        """A strictly interior point of `N0 + nu^T X >= 0`, by maximizing the slack.

        Maximize t subject to N0_i + sum_j nu_ij X_j >= t for every fluid species, with t
        capped so the program stays bounded. ⭐ This is what replaces the degenerate
        all-zero guess: at X = 0 any species with no feed is at exactly zero, and the log
        of that is what makes the residuals infinite before the solver has taken a step.
        """
        R_count = len(reactions)
        cap = 0.05 * max(N_fluid0, 1e-12)
        # variables [X_1..X_R, t];  -nu^T X + t <= N0
        A_ub = np.hstack([-nu_fluid.T, np.ones((nu_fluid.shape[1], 1))])
        result = optimize.linprog(
            c=np.concatenate([np.zeros(R_count), [-1.0]]),
            A_ub=A_ub, b_ub=N0_fluid,
            bounds=[(None, None)] * R_count + [(0.0, cap)])
        if result.success and result.x[-1] > 0:
            return result.x[:R_count]
        return np.full(R_count, 1e-3 * max(N_fluid0, 1e-12))

    def feasible(X):
        """Every fluid mole number non-negative -- the book's acceptability conditions."""
        return bool(np.all(N0_fluid + nu_fluid.T @ np.asarray(X, float) >= -1e-9))

    # ⭐ Multi-start Newton, then a constrained fallback. A Newton root-find reaches
    # machine precision on these equations when it starts anywhere sensible, so it is
    # tried from several points rather than from one: the interior point of the feasible
    # polytope, and a few fractions of the way across it. Any start that converges to a
    # *feasible* root is accepted.
    #
    # ⛔ For an ideal mixture G is strictly convex in the extents, so there is only one
    # feasible root -- the equilibrium state. Distinct starts landing on distinct
    # feasible roots would contradict that, so it is checked rather than assumed, and it
    # raises instead of quietly returning whichever one came first.
    starts = []
    if guess is not None:
        starts.append(np.asarray(guess, dtype=float))
    interior = interior_start()
    starts.append(interior)
    starts += [interior * fraction for fraction in (0.5, 0.25, 0.1, 0.05)]
    starts += [np.full(len(reactions), f * max(N_fluid0, 1e-12))
               for f in (0.1, 0.05, 0.01)]

    roots, best = [], None
    for start in starts:
        attempt = optimize.root(residuals, start, method="hybr", tol=1e-14)
        worst = float(np.max(np.abs(residuals(attempt.x))))
        if best is None or worst < best[0]:
            best = (worst, np.asarray(attempt.x, dtype=float))
        # ⚠️ `attempt.success` is deliberately **not** consulted. hybr reports failure
        # ("not making good progress") on roots it has already driven to 1e-15, because
        # its own step test cannot improve on machine precision. The residual is the
        # criterion, so the residual is what is measured.
        if worst <= max(tol, 1e-8) and feasible(attempt.x):
            if not any(np.allclose(attempt.x, seen, rtol=1e-6, atol=1e-9)
                       for seen in roots):
                roots.append(np.asarray(attempt.x, dtype=float))
    if len(roots) > 1:
        raise ValueError(
            f"different starting points converged on {len(roots)} distinct feasible sets "
            f"of extents: {roots}. For an ideal mixture the equilibrium state is unique, "
            f"so this needs looking at rather than picking one.")

    if roots:
        X_solved = roots[0]
    else:
        # Nothing converged cleanly. Fall back to a least squares that cannot leave the
        # feasible polytope, then try to polish it.
        constraint = optimize.LinearConstraint(nu_fluid.T, -N0_fluid, np.inf)
        fallback = optimize.minimize(lambda X: float(np.sum(residuals(X) ** 2)),
                                     best[1], method="SLSQP", constraints=[constraint],
                                     options={"maxiter": 800, "ftol": 1e-16})
        X_solved = np.asarray(fallback.x, dtype=float)
        polished = optimize.root(residuals, X_solved, method="hybr", tol=1e-14)
        if (polished.success
                and np.max(np.abs(residuals(polished.x))) < np.max(np.abs(residuals(X_solved)))
                and feasible(polished.x)):
            X_solved = np.asarray(polished.x, dtype=float)
    res = residuals(X_solved)
    sol = type("Solution", (), {"x": X_solved})()
    if np.max(np.abs(res)) > max(tol, 1e-8):
        raise ValueError(f"the coupled equilibrium relations did not converge: residuals "
                         f"{res}. Try `gibbs_minimization`, which cannot walk out of the "
                         f"feasible region the way this formulation can.")
    N, y, P_x = composition(sol.x)
    # ⛔ Illustration 13.3-1 states the trap in the book's own words: "Since these
    # equations are nonlinear, there will be more than one set of solutions for the molar
    # extents of reaction." Its two acceptability conditions -- no more steam consumed
    # than was supplied, no more hydrogen consumed than was produced -- are both just
    # "every mole number is non-negative", so that is what is checked, and a converged
    # root that fails it is rejected rather than returned.
    negative = {s: float(n) for s, n, is_fluid in zip(species, N, fluid_mask)
                if is_fluid and n < -1e-9}
    if negative:
        raise ValueError(
            f"the solver converged on a root outside the feasible region: "
            f"{negative} would be negative. The coupled relations have more than one "
            f"solution and this is one of the unphysical ones -- restart from a different "
            f"`guess`, or use `gibbs_minimization`, which is constrained.")
    return {"extents": sol.x,
            "moles": dict(zip(species, N)),
            "y": dict(zip(species, y)),
            "P": P_x,
            "Ka": dict(zip([r.name for r in reactions], Ka_values)),
            "residuals": res,
            "feasible": not negative}


def formation_gibbs_T(species, T, mode="full"):
    """dG_f(T) for one species, in J/mol -- its formation reaction from the elements.

    Appendix A.IV tabulates the Gibbs energy of formation at 25 C only. Gibbs
    minimization needs it at the reaction temperature, so the species is written as its
    formation reaction from the element reference states, and that reaction is carried to
    T by exactly the same integration as any other -- `Reaction.delta_G`.

    ⭐ This is checkable, and it is worth checking: for any balanced reaction,
    sum_i nu_i dG_f,i(T) must equal that reaction's own `delta_G(T)`. The two paths share
    no arithmetic beyond the Cp table -- the element terms have to cancel for them to
    agree -- so it tests the element bookkeeping in `elements` as well.
    `code/ch13/validation/` runs it.
    """
    species = str(species)
    counts = elements(species)
    stoich = {species: 1.0}
    for el, n in counts.items():
        if el not in ELEMENT_REFERENCE_STATE:
            raise KeyError(f"no reference state recorded for element {el!r}")
        ref = ELEMENT_REFERENCE_STATE[el]
        per_mole = elements(ref)[el]               # H2 carries two H, O2 two O
        stoich[ref] = stoich.get(ref, 0.0) - n / per_mole
    stoich = {k: v for k, v in stoich.items() if abs(v) > 1e-12}
    # An element in its own reference state cancels to nothing, by definition: its
    # formation reaction is empty and its formation Gibbs energy is zero at every
    # temperature, which is the convention Appendix A.IV is built on.
    if len(stoich) <= 1:
        T = np.asarray(T, dtype=float)
        return 0.0 if T.ndim == 0 else np.zeros_like(T)
    return Reaction(stoich).delta_G(T, mode=mode)


def gibbs_curve(reaction, T, initial, P=None, X=None, mode="full", phases=None):
    """G of the reacting mixture against the extent of reaction -- SIS Fig. 13.1-1.

    Returns `(X, G)` with G in J, measured relative to the elements in their standard
    states at 25 C, so the curve is directly the one the chapter plots for
    CO2 + H2 = CO + H2O:

        G(X) = sum_i N_i(X) [ dG_f,i(T) + R T ln (y_i P / 1 bar) ]

    ⭐ The figure's whole point is that the mixing term -- the logarithms -- is what puts
    the minimum at an interior X, and that dropping it (the dashed line the chapter draws)
    leaves a straight line with no minimum at all. `ideal` in the return is that dashed
    line, so a notebook can draw both.
    """
    if P is None:
        raise ValueError("give the pressure P in bar")
    phases = dict(phases or {})
    condensed = {s for s, kind in phases.items() if kind in ("solid", "pure")}
    in_excess = {s for s in condensed if s not in initial}
    lo, hi = reaction.extent_bracket(initial, inset=1e-6, in_excess=in_excess)
    X = np.linspace(lo, hi, 400) if X is None else np.atleast_1d(np.asarray(X, float))
    species = list(dict.fromkeys(list(initial) + reaction.species))
    gf = {s: float(formation_gibbs_T(s, T, mode=mode)) for s in species}
    G, G_ideal = [], []
    for x in X:
        N = reaction.moles(initial, x)
        total = sum(n for s, n in N.items() if s not in condensed)   # fluid phase only
        pure = sum(n * gf[s] for s, n in N.items())
        mixing = sum(n * R * T * np.log(max(n / total, 1e-300) * P / P_REF)
                     for s, n in N.items() if n > 0 and s not in condensed)
        G.append(pure + mixing)
        G_ideal.append(pure)
    return {"X": X, "G": np.array(G), "no_mixing": np.array(G_ideal)}


def gibbs_minimization(species, T, initial, P=None, mode="full", phases=None):
    """The equilibrium state as the minimum of G, subject to element conservation.

    The alternative to writing reactions at all: given a list of species that may be
    present, minimize

        G/RT = sum_i N_i [ dG_f,i(T)/RT + ln (y_i P / 1 bar) ]

    over the mole numbers, holding every element's total fixed at what the feed supplies
    and every mole number non-negative. No stoichiometry is chosen and no extent appears,
    so nothing has to be decided about which reactions are independent -- Sec. 13.3's
    difficulty. Returns the same shape of answer as `multireaction_extents`.

    ⚠️ The species list is a modeling choice and the answer depends on it entirely: a
    species left out cannot form, and one put in will form if it lowers G. That is a
    feature when it is deliberate and a trap when it is not.
    """
    if P is None:
        raise ValueError("give the pressure P in bar")
    species = list(dict.fromkeys(list(species) + list(initial)))
    phases = dict(phases or {})
    gf = np.array([float(formation_gibbs_T(s, T, mode=mode)) for s in species])
    N0 = np.array([float(initial.get(s, 0.0)) for s in species])
    if N0.sum() <= 0:
        raise ValueError("the feed is empty")

    els = sorted({el for s in species for el in elements(s)})
    A = np.array([[elements(s).get(el, 0.0) for s in species] for el in els])
    b = A @ N0

    # ⛔ The mole fraction denominator is the **fluid phase only**. A pure condensed phase
    # is its own phase: it has unit activity and it does not dilute the gas. Letting a
    # solid into the total is exactly the error Illustration 13.3-1 warns about in its own
    # note -- "otherwise Aspen Plus includes the solid in the mole fraction calculation
    # giving incorrect results" -- and it moves the answer by ~0.1 in mole fraction, which
    # is how it was caught here: the extents route and this one stopped agreeing.
    condensed = {s for s, kind in phases.items() if kind in ("solid", "pure")}
    fluid_index = [i for i, s in enumerate(species) if s not in condensed]
    scale = N0.sum()

    def objective(N):
        total = sum(max(N[i], 0.0) for i in fluid_index)
        out = 0.0
        for i, s in enumerate(species):
            if N[i] <= 0:
                continue
            kind = phases.get(s, "gas")
            if kind in ("solid", "pure"):
                a = 1.0
            elif total <= 0:
                continue
            elif kind == "gas":
                a = (N[i] / total) * (P / P_REF)
            else:
                a = N[i] / total
            out += N[i] * (gf[i] / (R * T) + np.log(max(a, 1e-300)))
        return out / scale

    start = np.where(N0 > 0, N0, scale * 1e-4)
    sol = optimize.minimize(objective, start, method="SLSQP",
                            bounds=[(0.0, None)] * len(species),
                            constraints=[{"type": "eq",
                                          "fun": lambda N: A @ N - b,
                                          "jac": lambda N: A}],
                            options={"maxiter": 500, "ftol": 1e-14})
    if not sol.success:
        raise ValueError(f"the Gibbs minimization did not converge: {sol.message}")
    N = np.clip(sol.x, 0.0, None)
    if np.max(np.abs(A @ N - b)) > 1e-7 * max(1.0, float(np.max(np.abs(b)))):
        raise ValueError("the converged answer does not conserve the elements; "
                         "the constraint was not satisfied")
    total = N.sum()
    return {"moles": dict(zip(species, N)),
            "y": dict(zip(species, N / total)),
            "P": P,
            "G_over_RT": float(objective(N) * scale),
            "elements": dict(zip(els, b))}


# --- Sec. 13.2: heterogeneous reactions and the Ellingham diagram -----------

def ellingham(reactions, T, mode="full", per_mole_O2=False):
    """dG_rxn(T) against temperature for a set of reactions -- SIS Sec. 13.2, Fig. 13.2-3.

    The Ellingham diagram plots the standard Gibbs energy change of metal oxidation
    against temperature, one straight-ish line per metal, and reads off which metal
    reduces which oxide: the lower line wins. Returns a DataFrame of dG_rxn in kJ/mol,
    one column per reaction, indexed by temperature.

    Parameters
    ----------
    reactions : sequence of Reaction, or dict of label -> Reaction
    per_mole_O2 : bool
        Scale each reaction so it consumes exactly one mole of O2, which is what makes
        the lines comparable -- the convention the diagram is drawn in. Raises for a
        reaction that has no O2.

    Notes
    -----
    ⓘ Figure 13.2-3 in the book is third-party art, redrawn from Lupis, and its caption
    says it is also on the website as an enlargeable PDF. This function is the numerical
    construction, not a reproduction of that figure: the lines it returns come from
    Appendix A.IV and A.II, so they carry the appendices' own accuracy and their own
    temperature limits.

    ⚠️ Appendix A.II's Cp correlations were fitted over roughly 273-1800 K. The classical
    Ellingham diagram runs to 2000 K and beyond, and a line drawn out there is an
    extrapolation. It is also drawn straight through melting and boiling points, where
    the real dG has a kink that these correlations know nothing about, since Appendix
    A.IV carries one state of aggregation per species.
    """
    import pandas as pd

    if isinstance(reactions, dict):
        labels, rxns = list(reactions), list(reactions.values())
    else:
        rxns = list(reactions)
        labels = [r.name for r in rxns]
    T = np.atleast_1d(np.asarray(T, dtype=float))
    columns = {}
    for label, rxn in zip(labels, rxns):
        if per_mole_O2:
            nu_O2 = rxn.nu.get("O2")
            if not nu_O2:
                raise ValueError(f"{label!r} does not contain O2, so it cannot be put "
                                 f"on a per-mole-of-O2 basis")
            rxn = rxn.scaled(1.0 / abs(nu_O2))
        columns[label] = np.atleast_1d(rxn.delta_G(T, mode=mode)) / 1e3
    return pd.DataFrame(columns, index=pd.Index(T, name="T (K)"))
