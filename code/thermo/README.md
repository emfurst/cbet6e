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

## Redlich–Kister correlation of mixing data — `fitting`

```python
from thermo import RedlichKister, tangent_intercepts, load_mixing_data

d = load_mixing_data("water-methanol-volume")           # SIS Table 8.6-1
rk = RedlichKister.fit(d.x1, d.dmixV_m3_mol, order=3)   # Eq. 8.6-5a
d1, d2 = rk.partial_molar_excess(0.5266)                # Eqs. 8.6-6a,b
rk.infinite_dilution()                                  # the two endpoint limits
RedlichKister.scan_order(d.x1, d.dmixV_m3_mol)          # rms, LOO-rms, endpoints vs order
```

The machinery of **SIS Sec. 8.6**, and the Python replacement for the 4e `PRTLMOLR`
worksheet. A property change on mixing must vanish at both pure limits, so it is fitted
with a polynomial that vanishes there by construction; the partial molar properties are
then the *slope* of that curve rather than a tangent line drawn by hand.
`tangent_intercepts` is the graphical route (Eqs. 8.6-4a,b) as arithmetic — a point and a
slope — and the ch8 notebooks run both to show they agree.

**Unit-agnostic**: it never mentions volume or enthalpy, which is the point of Eq. 8.6-10.
The caller keeps SI.

> ⚠️ **The derivative is far more sensitive to the fit order than the fit is.** Adding a
> term barely moves the curve and can move the infinite-dilution partial molar properties
> by ten percent or more — those are the endpoints, the one place the correlation
> extrapolates. `loo_rms` and `scan_order` exist so a notebook can *show* that instead of
> hiding it behind a default.
>
> ⚠️ **The summation identity is exact by construction**, so agreement there is not
> evidence the fit is good — `residuals` is. It is still worth checking, because a
> *graphical* construction does not satisfy it: the 5e's Table 8.6-4 misses it by up to
> 16 J/mol (`code/ch8/validation/`).

Verified against the book's own printed constants: with $a = (-4.0034, -0.17756, 0.54139,
0.60481)\times10^{-6}$ m³/mol, Eqs. 8.6-6a,b reproduce **eleven of the twelve rows of
Table 8.6-2 exactly**, and both infinite-dilution values. The twelfth disagreement is a
**sign error in the printed table** (row $x_1 = 0.9489$), confirmed against the book's own
$\bar V_1$ column.

## Reaction and mixing data — `data`

```python
from thermo.data import reaction_cp, formation_enthalpy, formation_gibbs
from thermo import load_mixing_data, get_reaction_species

a, b, c, d, e = reaction_cp("NO2")        # SCALED Appendix A.II coefficients, SI
formation_enthalpy("N2O4")                # J/mol  (the CSV stores kJ/mol)
```

> ⚠️ **`react_property.csv` stores the heat-capacity coefficients as Appendix A.II
> *prints* them** — $b$, $c$, $d$ multiplied by $10^{2}$, $10^{5}$, $10^{9}$. Evaluating
> the raw columns is wrong by three orders of magnitude and nothing warns you. Always go
> through `reaction_cp`. Verified 2026-08-09: the scaled columns reproduce Illustration
> 8.5-2's own combined coefficients and all five of its printed heats of reaction exactly.

## Property charts — `charts`, `ph_chart`, `steam_chart`

Added 2026-08-08. The machinery behind the book's property charts, promoted out of the
chapter-3 notebooks so that a figure's notebook can live in the figure's own chapter
(author, 2026-08-08) without copying ~1,000 lines three ways.

| module | what it holds | draws |
|---|---|---|
| `charts` | the drawing craft: line weights, chart-paper grids, the label-placement layer, `use_book_style()` | — |
| `ph_chart` | `ChartFluid` (a per-kilogram PR wrapper on the chart's datum) + the dome / isotherm / isentrope / isochore / quality families | Figs. **3.3-2**, **3.3-3**; **5.1-3**; **`c06uf002`** |
| `steam_chart` | `SteamTables` — Appendix A.III digitized — and the isobar / isotherm / isenthalp / quality families built from it | Figs. **3.3-1a**, **3.3-1b**; **`c05uf001`** |

```python
from thermo.charts import use_book_style
from thermo.ph_chart import ChartFluid, ph_chart

use_book_style()                                    # Computer Modern, TeX if present
n2 = ChartFluid("nitrogen", M=28.014, T_triple=63.15)
fig, ax = plt.subplots(figsize=(7.0, 4.7))
ph_chart(ax, n2, H_lim=(0, 900), isotherms=range(80, 401, 20))
```

**Two different sources, deliberately.** `ph_chart` computes from Peng–Robinson;
`steam_chart` reads Appendix A.III. For water the book *tabulates* the properties, so a
cubic would be a step backwards from data the reader already has — and for methane and
nitrogen it tabulates nothing, which is why those come from the equation of state the
reader is about to be taught.

⚠️ **PR puts saturated-liquid density about 12 % high**, the textbook failing of a cubic.
It lands on exactly one family — the constant-volume curves in the compressed-liquid
region — which should not be read quantitatively there. The dome, isotherms and
isentropes are unaffected.

**These modules are loaded lazily** (`thermo/__init__.py`, PEP 562), so
`from thermo import PengRobinson` does not pull matplotlib. None of the three imports
it at module level either; it arrives only when a drawing function runs.

**Validation.** Every curve family was checked against the chapter-3 notebooks it was
promoted from, and reproduces them **bit-for-bit** — `max|diff| = 0` on the dome,
isotherms, isentropes, isochores, quality lines, isobars, liquid branches, merged
supercritical isobars and lines of constant enthalpy, for methane, nitrogen and water.
The full charts render identically too: Fig. 3.3-2 gives 85 paths and 57 labels, and
Fig. 3.3-1(b) 99 paths and 23 labels, with identical geometry, text, position and
rotation in both cases.

## Provenance

Refactored from the validated ch6/ch7 notebooks (PR), the legacy MATLAB `unifac.m`
(a faithful translation of its `calc_coeff`), and the ch3 chart notebooks (`charts`,
`ph_chart`, `steam_chart`). See Appendix B and `code/data/README.md`.
