# `thermo` — Aspen-optional Python substitutes

A small, dependency-light package (numpy + scipy + pandas) providing the pure-fluid
**Peng–Robinson EOS** and **UNIFAC** activity coefficients, for use when Aspen Plus is
not available. Parameters and constants come from the tables in `code/data/`.

```python
import sys; sys.path.append("..")      # so `import thermo` works from a chapter folder
from thermo import PengRobinson, UNIFAC
```

## Peng–Robinson (pure fluids)

```python
pr = PengRobinson.from_database("CO2")          # or PengRobinson(Tc, Pc, omega, cp=...)
pr.molar_volume(300, 5e6, phase="vapor")        # m^3/mol
pr.fugacity(300, 5e6, "vapor")                  # Pa
pr.departure_H(300, 5e6, "vapor")               # J/mol   (SIS Eq. 6.4-29)
pr.departure_S(300, 5e6, "vapor")               # J/(mol K)(SIS Eq. 6.4-30)
pr.vapor_pressure(250)                          # Pa, by equal fugacity
```

SI units throughout (T in K, P in Pa, V in m³/mol). `from_database` reads
`pure_property.csv` (converting its Pc from bar to Pa) and carries the Appendix A.II
ideal-gas Cp coefficients. Verified: `vapor_pressure(Tb)` returns ≈ 1 atm for benzene,
n-hexane, ammonia, and propane (within 1%).

> Note: `departure_S` uses the correct `R·ln(Z−B)` (SIS 6.4-30); the ch6 throttle
> notebook's `calc_depS` wrote `R·(Z−B)`, which looks like a typo to fix in that notebook.

## Peng–Robinson (mixtures)

```python
from thermo import PRMixture
mix = PRMixture.from_database(["benzene", "toluene"])   # or PRMixture([pr1, pr2], kij=...)
mix.fugacity([0.4, 0.6], 373.15, 1.1e5, "liquid")       # component fugacities f_i (Pa)
P, y = mix.bubble_pressure([0.4, 0.6], 373.15)          # bubble P (Pa) + incipient vapor
P, x = mix.dew_pressure([0.6, 0.4], 373.15)             # dew P (Pa) + incipient liquid
T, y = mix.bubble_temperature([0.4, 0.6], 1.013e5)      # bubble T (K)
beta, x, y = mix.flash([0.4, 0.6], 373.15, 1.064e5)     # isothermal flash: vapor frac + phases
```

van der Waals one-fluid mixing rules with binary interaction parameters `kij`
(default 0), the species-in-mixture fugacity coefficient (SIS Eq. 9.4-9), and the
Chapter-10 φ–φ VLE drivers (bubble/dew P & T, Rachford–Rice flash). Verified on the
near-ideal benzene/toluene pair: bubble P within 0.3% of Raoult's law, equal
component fugacities in both phases to ~1e-14, dew inverts the bubble exactly, and
the flash mass balance closes with single-phase edges (β→0/1) correct.

## UNIFAC (activity coefficients)

```python
u = UNIFAC("modified")                          # temperature-dependent (Dortmund)
ethanol = {1: 1, 2: 1, 14: 1}                   # CH3 + CH2 + OH(p), by subgroup number
hexane  = {1: 2, 2: 4}                           # 2 CH3 + 4 CH2
u.gamma([ethanol, hexane], x=[0.3, 0.7], T=298.15)   # -> array of activity coefficients
```

Subgroup numbers are in `code/data/unifac_subgroups.csv`. Verified against the pure
limit (γ=1), a near-ideal alkane pair (γ≈1), and a strong positive-deviation system
(γ_EtOH^∞ ≈ 53 in hexane).

- **`kind="modified"`** — fully working (Dortmund parameters, `Ψ = exp[-(a/T + b + cT)]`,
  combinatorial with `r^0.75`).
- **`kind="original"`** — not yet usable: the classic-UNIFAC `a_mn` matrix was recovered
  (`unifac_interactions_original.csv`), but the original subgroup R/Q constants and
  44-group names were not in the legacy `.mat`. They still need to be sourced (the 5e
  chapter/appendix PDFs in `legacy-5e/Sandler_CBET_5e files/`), after which the same class handles
  both variants.

## Provenance

Refactored from the validated ch6/ch7 notebooks (PR) and the legacy MATLAB `unifac.m`
(a faithful translation of its `calc_coeff`). See Appendix B and `code/data/README.md`.
