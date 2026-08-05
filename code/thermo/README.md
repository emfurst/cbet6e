# `thermo` — Aspen-optional Python substitutes

A small, dependency-light package (numpy + scipy + pandas) providing the pure-fluid
**cubic equations of state** (Peng–Robinson, PRSV, van der Waals) and **UNIFAC**
activity coefficients, for use when Aspen Plus is not available. Parameters and
constants come from the tables in `code/data/`.

```python
import sys; sys.path.append("..")      # so `import thermo` works from a chapter folder
from thermo import PengRobinson, VanDerWaals, UNIFAC
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
`pure_property.csv` (converting its Pc from bar to Pa). Verified:
`vapor_pressure(Tb)` returns ≈ 1 atm for benzene, n-hexane, ammonia, and propane
(within 1%).

> **⚠️ The Cp coefficients in `pure_property.csv` are the Reid–Prausnitz–Poling set,
> not the book's Appendix A.II set.** This README previously said Appendix A.II; that
> was wrong, and the same mislabel is in two ch6 notebooks. It matters below room
> temperature: for oxygen the two sets differ by 152 J/mol in H̲ and 1.23 J/(mol·K) in
> S̲ at 73 K, which is the difference between reproducing printed Tables 7.5-1 and
> 7.5-2 and not. Pass `cp=` explicitly when you need the book's numbers. See
> `revision_notes/c07.md`, *ch6 heat-capacity inconsistency* — **decided 2026-08-03**: ch6
> uses the book's Appendix A.II set, via `cp=APPENDIX_A2_CP[...]`, so ch6 and ch7 sit on
> one basis. `pure_property.csv` is left as the RPP table it is.

> **⚠️ Below 273 K, use `APPENDIX_A2_CP_CRYO` — not `APPENDIX_A2_CP`.** Appendix A.II's
> familiar row is valid **273–1800 K**, and Illustrations 6.4-1, 7.5-1 and 7.5-2 all run
> below that. The extrapolation is *not* benign: for oxygen at 173 K it is 1.235 J/(mol·K)
> (4.2%) low, worth **67 J/mol in H̲ and 0.32 J/(mol·K) in S̲** — a hundred times the error
> from the gas constant or from the 5e's own arithmetic. `APPENDIX_A2_CP_CRYO` is the same
> cubic form refitted to NIST-JANAF over 100–700 K (oxygen and nitrogen only), which cuts
> those to 2.8 J/mol and 0.011 J/(mol·K). Derivation and accuracy table are in
> `data.py`. **Appendix A.II must gain these rows in print** — until it does, the printed
> book does not contain the constants these notebooks use.

> **⚠️ For oxygen, `from_database` is not just a different Cp — it is a different EOS.**
> `pure_property.csv` carries ω = 0.025 and Pc = 50.4 bar against the book's Table 6.6-1
> values ω = **0.021** (κ = 0.4069, the number Illustration 6.4-1 prints) and
> Pc = **5.046 MPa**. Build the object explicitly when reproducing a printed table.

> Note: `departure_S` uses the correct `R·ln(Z−B)` (SIS 6.4-30); the ch6 throttle
> notebook's `calc_depS` wrote `R·(Z−B)`, which looks like a typo to fix in that notebook.

### PRSV

Pass `kappa1` to switch the same class to the Stryjek–Vera modification (SIS
Eqs. 7.5-1 and 7.5-2), where κ becomes temperature dependent:

```python
w = PengRobinson.from_database("water", kappa1=-0.0665)   # Illustration 7.5-3
w.is_prsv        # True
w.kappa_T(300)   # kappa at 300 K; kappa0 alone for standard PR
```

κ₀ is coded with the **minus** sign on ω², following Stryjek and Vera, *Can. J. Chem.
Eng.* **64**, 323 (1986). The 5e prints a plus in Eq. 7.5-2, which is an erratum — see
`revision_notes/c07_manuscript_edits_7.5-3.md`. Verified against Illustration 7.5-3:
with the book's own Tc and Pc (SIS Table 6.6-1), κ₀ = 0.87188 and the printed PRSV
column reproduces to four figures (0.6092 vs 0.6094 kPa at 273.15 K; 1550.1 vs 1550.0
at 473.15 K). Note that PRSV with `kappa1=0` is **not** standard Peng–Robinson: κ₀ and
the PR κ differ by 0.16% for water.

## van der Waals (pure fluids)

```python
from thermo import VanDerWaals
vdw = VanDerWaals.from_database("n-butane")
vdw.vapor_pressure(300) / 1e5      # 7.93 bar — the book's Fig. 7.5-2, curve a
```

Same API as `PengRobinson`, so one calculation can be handed either equation; `a(T)`
is a method here too, and `dadT` returns zero. Verified against
`ch7/vapor_pressure_n_butane.ipynb`, which computes Fig. 7.5-2 from hand-coded
equations: the vapor pressures agree at every tabulated temperature (7.930 bar at
300 K, 29.62 bar at 400 K).

## What the two share: `CubicEOS`

Both classes inherit root selection, molar volume, fugacity, the spinodal pressures
and the saturation solver from `thermo.CubicEOS` (`cubic.py`). A subclass supplies
only `Tc`, `Pc`, `b`, `pressure(V, T)`, `compressibility(T, P)` and `ln_phi(T, P, phase)`.

`vapor_pressure` is **the algorithm of SIS Figure 7.5-1 as redrawn for the 6e**: the
equal-fugacity root is bracketed between the two turning points of the isotherm and
found by Brent's method, so no initial guess is needed and the trivial solution cannot
occur. `spinodal_bounds(T)` returns that bracket, or `None` at and above Tc.

Verified: identical to the previous Newton-from-a-Pitzer-guess solver to 2×10⁻¹³
relative over five fluids × 18 temperatures, and it converges for all of n-butane,
water, oxygen, n-decane, methanol and ethanol at every Tr from 0.30 to 0.98, for both
equations of state.

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
