# `thermo` — the book's models, in Python

A small, dependency-light package (numpy + scipy + pandas) providing the thermodynamic
**models** the book teaches — cubic equations of state for pure fluids and mixtures, the
activity coefficient models, UNIFAC, the property-chart machinery and the book's own data
tables. Parameters and constants come from `code/data/`.

⚠️ **This is not an "Aspen substitute," which is how this file used to describe it.**
Python and Aspen are complementary parts of the toolbox — general-purpose scientific
computing versus specialized process modeling, not a learning tool versus a real one. This
package is open, inspectable and hackable and runs with nothing but a Python install;
Aspen owns flowsheets and plant-scale design. `revision_notes/bapp02.md` is the print-facing
version of that argument.

```python
import sys; sys.path.append("..")      # so `import thermo` works from a chapter folder
from thermo import PengRobinson, VanDerWaals, UNIFAC
```

The package holds the *method*; the per-chapter notebooks are the *applications* that call
it. That split, and everything still to be built, is [`ROADMAP.md`](ROADMAP.md).

| module | what it holds | § |
|---|---|---|
| [`cubic`](#what-the-two-share-cubiceos) | `CubicEOS` — roots, fugacity, spinodals, the $P^{vap}$ solver | 6, 7 |
| [`peng_robinson`](#pengrobinson-pure-fluids) | pure-fluid PR and PRSV | 6, 7 |
| [`van_der_waals`](#van-der-waals-pure-fluids) | pure-fluid vdW | 6, 7 |
| [`pr_mixture`](#pengrobinson-mixtures) | vdW one-fluid mixing rules, species fugacity, bubble/dew/flash | 9, 10 |
| [`activity_models`](#activity-coefficient-models--activity_models) | nine correlative and predictive γ models on one interface | 9 |
| [`unifac`](#unifac-activity-coefficients) | modified (Dortmund) UNIFAC | 9 |
| [`electrolytes`](#electrolyte-solutions--electrolytes) | Debye–Hückel and its extensions; ionic strength, $\gamma_\pm$ | 9, 15 |
| [`wong_sandler`](#combined-eos--gex-the-wongsandler-mixing-rule--wong_sandler) | an activity coefficient model inside a cubic EOS | 9, 10 |
| [`fitting`](#redlichkister-correlation-of-mixing-data--fitting) | Redlich–Kister correlation, partial molar properties | 8 |
| [`data`](#the-data-tables--data) | every table above, plus the property and reaction databases | all |
| [`charts`, `ph_chart`, `steam_chart`](#property-charts--charts-ph_chart-steam_chart) | the book's property charts | 3, 5, 6 |

⬜ **Not built yet:** chemical equilibrium and electrochemistry, for ch13–15. Chapters
6–9 are complete; ch10's VLE drivers are next. See the roadmap.

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
mix.ln_phi([0.4, 0.6], 373.15, 1.1e5, "liquid")         # ln phi-bar_i  (SIS Eq. 9.4-9)
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

### $k_{ij}$ from the book's own table

`kij="table"` looks the pairs up in **Table 9.4-1** (`code/data/pr_kij.csv`, 127 pairs
over 20 species, digitized 2026-08-12) instead of leaving them at zero:

```python
mix = PRMixture.from_database(["ethane", "n-butane"], kij="table")   # k12 = 0.010
```

Species are matched on **either** the name or the formula as Table 9.4-1 prints them, so
`["methane", "CO2"]` and `["CH4", "carbon dioxide"]` both resolve. Verified against the
book's three worked values: 0.010 ethane/n-butane (Illustration 9.4-3), 0.09
methane/CO₂ (9.4-4), 0.018 n-pentane/benzene (9.4-5).

> ⚠️ **A blank in Table 9.4-1 is not a zero, and this is where that bites.** The table is
> 65% blank, and its own footnote says to estimate the missing value from a similar
> mixture. So a pair the table lacks is set to zero **and warned about** — that judgment
> belongs to the caller, not to a default. Pass `warn_missing=False` once you have decided
> zero is what you want, or call `pr_kij_matrix(keys, strict=True)` to raise instead.

## UNIFAC (activity coefficients)

```python
u = UNIFAC("modified")                          # temperature-dependent (Dortmund)
ethanol = {1: 1, 2: 1, 14: 1}                   # CH3 + CH2 + OH(p), by subgroup number
hexane  = {1: 2, 2: 4}                           # 2 CH3 + 4 CH2
u.gamma([ethanol, hexane], x=[0.3, 0.7], T=298.15)   # -> array of activity coefficients
```

Subgroup numbers are in `code/data/unifac_subgroups.csv`, which since 2026-08-12 is a
digitization of the **book's own Table 9.5-2** — 92 subgroups over 46 main groups, with
the example-assignment column. Rather than transcribing a subgroup number by hand, ask
for it:

```python
from thermo.data import unifac_groups
unifac_groups("benzene")        # -> {9: 6}, read out of Table 9.5-2's own examples
unifac_groups("cyclohexane")    # -> {78: 6}
```

Verified against Illustration 9.5-2's printed group sums (r_B = 2.2578, q_B = 2.5926,
r_TMP = 5.0600, q_TMP = 6.3675), the pure limit (γ=1), a near-ideal alkane pair (γ≈1),
and a strong positive-deviation system (γ_EtOH^∞ ≈ 53 in hexane).

- **`kind="modified"`** — the book's model (Dortmund parameters,
  `Ψ = exp[-(a/T + b + cT)]`, combinatorial with `r^0.75`).
- **`kind="original"`** — ⛔ **retired 2026-08-12.** The book teaches modified UNIFAC
  only: Table 9.5-2 *is* the Dortmund set and Eq. 9.6-12a carries the `r^0.75`
  modification (`revision_notes/c09.md` **D1**). The `original` branch stays
  unimplemented and `unifac_interactions_original.csv` is retained but unused.

⚠️ **The book and this module do not compute the same Ψ.** Eq. 9.6-13b prints
`Ψ = exp(-a/T)`, a single temperature-independent constant per pair, where modified
UNIFAC uses `a + bT + cT²` — which the parameter file carries and this module evaluates.
Where they differ, the module is right; see c09.md **T5**.

## Activity coefficient models — `activity_models`

Every correlative model of SIS Sec. 9.5 and the predictive regular solution model of
Sec. 9.6, behind **one interface**:

```python
from thermo import VanLaar, NRTL, RegularSolution, UNIQUACModel

vl = VanLaar(0.415, 0.706)
vl.gamma([0.6, 0.4], T=328.15)      # activity coefficients
vl.gex([0.6, 0.4], T=328.15)        # molar excess Gibbs energy, J/mol
```

| model | class | parameters |
|---|---|---|
| one-constant Margules (9.5-1, 9.5-5) | `OneConstantMargules` | `A` (J/mol) |
| two-constant Margules (9.5-6, 9.5-7) | `TwoConstantMargules` | `A`, `B` (J/mol) |
| Redlich–Kister G^ex (9.5-6, any order) | `RedlichKisterGex` | `a_0 … a_n` (J/mol) |
| van Laar (9.5-9) | `VanLaar` | `alpha`, `beta` |
| Wilson (9.5-11, 9.5-12a) | `Wilson` | `Lambda_ij`, any *n* |
| NRTL (9.5-13, 9.5-14a) | `NRTL` | `tau_ij`, `alpha`, any *n* |
| Flory–Huggins (9.5-17, 9.5-18) | `FloryHuggins` | `chi`, `m` |
| UNIQUAC (9.5-19 … 9.5-23b) | `UNIQUACModel` | `r`, `q`, `tau_ij`, any *n* |
| regular solution (9.6-8 … 9.6-11) | `RegularSolution` | `V`, `delta`, any *n* |

⚠️ **`UNIQUACModel` is the molecular model of Sec. 9.5**; `UNIFAC` is the
group-contribution model of Sec. 9.6. Exported under distinct names so the two cannot be
confused. `UNIQUACModel.from_groups` builds *r* and *q* from Table 9.5-2, so the same
molecular description feeds both.

`Wilson`, `NRTL`, `UNIQUACModel` and `RegularSolution` take any number of species — pass a
square parameter matrix in place of the two binary constants; the rest are binary, as the
book writes them. Fitting: `Model.fit(x1, gamma1, gamma2)` throughout, plus
`VanLaar.from_single_point` for Eq. 9.5-10 and `RedlichKisterGex.fit_gex` for the route of
Illustration 10.2-4. `fit_binary(model_cls, ...)` is the same fit as a free function, for
`VanLaar` and `Wilson`, whose two constructor arguments *are* the two fitted parameters —
useful for looping several models over one data set.

⭐ **`check_gibbs_duhem(T)` is worth running.** Where the book prints both a G^ex
expression *and* activity coefficient expressions, they must satisfy
`G^ex/RT = Σ x_i ln γ_i` exactly — so coding both and comparing them is a real test rather
than a tautology. Eight of the nine models have that independent expression to check
against (regular solution does not, and returns `0.0` saying so); **all eight now return
≤6e-15**, worst case UNIQUAC.

⛔ **Two of them did not, and both were errata in the book** — Eqs. 9.5-12b and 9.5-18,
in print through five editions and corrected in the 6e manuscript
(`revision_notes/c09.md` §12.1). The classes carry the corrected forms, and each takes a
keyword to reproduce the printed one for comparison. ⚠️ **This is the check to run first on
every model added from here on**; it is what makes the book disagree with itself.

Tables 9.5-1 (van Laar constants for 30 binaries) and 9.6-1 (molar volumes and solubility
parameters) are included as `TABLE_9_5_1_VAN_LAAR` and `TABLE_9_6_1`;
`RegularSolution.from_table_9_6_1` takes the traditional cc/mol and (cal/cc)^(1/2) as
printed and converts to SI.

## Electrolyte solutions — `electrolytes`

Added 2026-08-12. The three equations of SIS Sec. 9.10, behind one class, because the
book presents the last two as *modifications of the first* rather than as rival models:

```python
from thermo import DebyeHuckel, Electrolyte, ELECTROLYTES, ionic_strength

DebyeHuckel("NaCl")                                # Eq. 9.10-15, the limiting law
DebyeHuckel("NaCl", beta_a=1)                      # Eq. 9.10-17
m = DebyeHuckel("NaCl", beta_a=1, delta=0.137)     # Eq. 9.10-18
m.equation                                         # '9.10-18' -- which one you built
m.gamma_pm(0.5)                                    # mean ionic activity coefficient
m.fit_delta(M, gamma_exp)                          # returns a NEW model, fitted
```

Concentrations are **molalities** (mol per kg of *solvent*), $\alpha$ and $\beta$ come
from Table 9.10-1 for water at the temperature you pass, and `Electrolyte` checks
electroneutrality (Eq. 9.10-2) in its constructor — a salt that fails it is not a salt,
and every quantity downstream would return a number anyway.

**This module does not inherit `ActivityModel`, deliberately.** Everything in Sec. 9.5
is a function of mole fraction returning one coefficient per species; this is a function
of *ionic strength* returning a single mean for the salt. The difference is the physics
of the section: electroneutrality means $\gamma_+$ and $\gamma_-$ can never be varied
apart, so only $\gamma_\pm$ is measurable and only $\gamma_\pm$ is what a model can be
asked for. Sharing the interface would be a lie about what is being computed.

Verified against the chapter's own numbers: $I = M$, $3M$, $4M$ for NaCl, CaCl₂ and
CuSO₄ (the three cases Fig. 9.10-1's caption states); Illustration 9.10-1's model ranking
in the right order; and Illustration 9.10-2's fitted **δ = 0.137**, recovered as 0.1371
and stable to the fourth decimal under either weighting.

⭐ **Table 9.10-1 was audited, and it is correct.** Its two columns are not independent —
Debye–Hückel theory makes $\alpha\rho/\beta^3$ a pure constant, 33.03, at every
temperature. All fourteen rows return 33.04 ± 0.08%. Drop the solvent density and the
ratio instead drifts monotonically by 4.3% from 0 to 100 °C, which is exactly how much
water's density falls: the table reproduces a density curve it never mentions, so its
numbers were computed from theory and transcribed correctly. The two illustrations only
ever use the 25 °C row, so nothing else in the chapter would have caught an error in the
other thirteen.

⚠️ **The book's 25 °C constant is not self-consistent to four digits.** Illustration
9.10-2's Eqs. (1)–(3) use $\alpha = 1.178$ where Table 9.10-1 prints 1.175. It is a 0.26%
discrepancy that does not reach the answer — δ comes out 0.137 either way — but pass
`alpha=1.178` if you want the illustration's printed equations exactly.

## Combined EOS + $G^{ex}$: the Wong–Sandler mixing rule — `wong_sandler`

Added 2026-08-12. SIS Sec. 9.9 — an activity coefficient model placed *inside* a cubic
equation of state, so one model covers a highly nonideal mixture in both phases and over
a wide range of temperature and pressure.

```python
from thermo import WongSandler, GexFromUNIFAC, UNIFAC
from thermo.data import unifac_groups

gex = GexFromUNIFAC(UNIFAC(), [{18: 1, 1: 1}, unifac_groups("water")])   # acetone/water
ws = WongSandler.from_database(["acetone", "water"], gex, kij=0.05)

ws.a_mix(x, T), ws.b_mix(x, T)      # Eq. 9.9-9a -- both depend on T AND on G^ex
ws.ln_phi(x, T, P, "liquid")        # Eq. 9.9-11
ws.check_boundary(x, T)             # the two conditions the rule is derived from
```

`gex` is anything exposing `gex(x, T)` and `gamma(x, T)` — every model in
`activity_models`, or `GexFromUNIFAC` to use UNIFAC (the *predictive* route of §10.3,
which needs no mixture data at all). Both combining rules of Eq. 9.9-10 are available
via `combining=`.

⚠️ **`b` is temperature dependent here**, unlike the van der Waals one-fluid `b`, and
neither `a` nor `b` is a polynomial in mole fraction — both come out of $Q/(1-D)$, and
$D$ carries the whole activity coefficient model. That is the point of the method, not
a side effect.

⚠️ **$D > 1$ is normal, not an error.** $D$ is dominated by $\sum_i x_i a_i/(b_i RT)$,
about 12 for liquid water at 25 °C and greater than 1 for any fluid below its Boyle
temperature. $Q$ is negative for the same reason and the two divide to a positive `b`;
the pure-component limit runs through $1-D < 0$ and returns $b_i$ exactly. Guard on
$D = 1$, never on $D > 1$.

### Verified without a single printed number

§9.9 has **no illustration, no figure and no worked value** anywhere in the chapter, so
there is nothing to reproduce. What can be done instead is to check the derivation
against itself, which `code/ch9/` does:

| check | result |
|---|---|
| Eq. 9.9-3, the second virial condition | ~1e-16 at every composition |
| Eqs. 9.9-4/-7/-8, the excess Helmholtz condition | ~1e-15 |
| pure-component limit, $a \to a_i$ and $b \to b_i$ | 1e-10 |
| Eqs. 9.9-12, 9.9-13 vs numerical differentiation | 1e-10 |
| Euler identities on the composition derivatives | exact |
| **Eq. 9.9-11 vs `PRMixture`, fed vdW one-fluid input** | **0 to 3e-16** |

⭐ **The last one is the check that reaches outside §9.9.** Eq. 9.9-11 claims to be the
general cubic-mixture fugacity coefficient, written with composition derivatives instead
of $b_i$ and $2\sum_j x_j a_{ij}$. Feeding it the van der Waals one-fluid rules must
therefore reproduce `PRMixture`, which was verified against **Illustration 9.4-3** — and
it does, exactly.

⭐ **And Eq. 9.9-8 was confirmed by brute force.** Taking the equation of state to
$10^{11}$ bar, $\underline A^{ex}$ converges on $C^*[a/b - \sum x_i a_i/b_i]$ to seven
figures while $\underline G^{ex}$ diverges to $-5\times10^{10}$ J/mol — which verifies
the equation, the Peng–Robinson value of $C^*$ (left as Problem 9.31), and the section's
argument for using Helmholtz rather than Gibbs, all at once.

⚠️ **Eq. 9.9-10a is sign-ambiguous as printed.** Its cross term is the square root of a
product of two pure second virial coefficients, and below the Boyle temperature both are
negative — so the principal root is *positive* where both diagonals are negative, and
$Q$ collapses through zero mid-composition. `cross_matrix` carries the sign of the pures
through the root, which is the only reading consistent with the rule reducing to
$b_i - a_i/RT$ on the diagonal. Worth a footnote in print; see `revision_notes/c09.md`.

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

## The data tables — `data`

⭐ **The database is the linchpin.** These CSVs are the book's own tables — effectively a
digitized Appendix A plus the model-parameter tables of ch6–9 — in plain text, so a student
can read them, check one against the NIST WebBook, or add a compound. That is the thing a
licensed black box will not allow, and it is why the data layer is documented here rather
than treated as an implementation detail.

| what | loader | source |
|---|---|---|
| pure-component constants (Tc, Pc, ω, Cp, Antoine) | `load_pure_properties()`, `get_compound(key)` | `pure_property.csv` — the RPP set ⚠️ below |
| ideal-gas $C_p^*$, **as the book prints it** | `APPENDIX_A2_CP`, `APPENDIX_A2_CP_CRYO` | Appendix A.II (+ the cryogenic refit ⚠️ below) |
| PR constants **as the book prints them** | `TABLE_6_6_1` | Table 6.6-1 — oxygen, nitrogen, methane |
| PR binary interaction parameters | `load_pr_kij()`, `pr_kij_matrix(keys)` | Table 9.4-1 — 127 pairs / 20 species |
| UNIFAC subgroups: R, Q, main group, examples | `load_unifac_subgroups()`, `unifac_groups(name)` | Table 9.5-2 — 92 subgroups / 46 main groups |
| UNIFAC group-interaction parameters $a, b, c$ | `load_unifac_interactions()` | legacy `UNIFAC_data.mat` ⚠️ below |
| formation $\Delta_f G$, $\Delta_f H$ + $C_p^*$ | `load_reaction_properties()`, `get_reaction_species(name)`, `reaction_cp`, `formation_enthalpy`, `formation_gibbs` | Appendices A.IV + A.II — 99 species |
| property change on mixing | `load_mixing_data(key)` | Tables 8.6-1, 8.6-3 |
| steam tables (`steam_*.csv`) | `SteamTables` — in `steam_chart`, not `data`, because it interpolates as well as loads | Appendix A.III |

```python
from thermo.data import reaction_cp, formation_enthalpy, formation_gibbs
from thermo import load_mixing_data, get_reaction_species

a, b, c, d, e = reaction_cp("NO2")        # SCALED Appendix A.II coefficients, SI
formation_enthalpy("N2O4")                # J/mol  (the CSV stores kJ/mol)
```

The two UNIFAC tables and the two mixture tables are loaded through `thermo.data` rather
than re-exported at top level, because a notebook that wants them is doing data work and
the longer path says so:

```python
from thermo.data import unifac_groups, load_unifac_subgroups, pr_kij_matrix
```

⚠️ **`load_unifac_subgroups` validates on every load, not on request.** A duplicate
subgroup number or an R/Q that disagrees with its main group changes every activity
coefficient it touches, silently — so the integrity rules run each time the file is read.
They are not decorative: they caught duplicate subgroup numbers and silicon names sitting
on the Dortmund cyclic groups in the legacy `.mat` extraction (`revision_notes/c09.md`
§12.3). ⬜ **Check the other legacy-derived tables the same way as ch10–15 open.**

> ⚠️ **R and Q come from the book; the $a_{mn}$ do not.** `unifac_subgroups.csv` is the
> printed Table 9.5-2, but the *numbers* identifying each subgroup and main group are
> printed nowhere in the chapter, and `unifac_interactions_modified.csv` is keyed on the
> main-group number — so both come from the legacy extraction, which the book cannot check.
> The asymmetry is deliberate and is recorded in `code/data/README.md`.

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
`ph_chart`, `steam_chart`). `pr_mixture`, `activity_models` and `fitting` were written
from the book's own equations and verified against its printed numbers, not ported.
See Appendix B and `code/data/README.md`.

⭐ **Writing them that way is what found the errata.** Making two printed expressions for
the same quantity compute against each other — $\underline G^{ex}/RT = \sum_i x_i\ln\gamma_i$
for the activity models, $P\hat V = \hat H - \hat U$ for the steam tables, the summation
identity for Redlich–Kister — has now caught errors in Eqs. 9.5-12b and 9.5-18, in Table
8.6-2, and in Tables 9.5-2 and A.III. A model that merely runs has not been checked.
