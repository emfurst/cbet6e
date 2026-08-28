# `thermo` — the book's models, in Python

A small, dependency-light package (numpy + scipy + pandas) providing the thermodynamic
**models** the book teaches — cubic equations of state for pure fluids and mixtures, the
activity coefficient models, UNIFAC, the property-chart machinery and the book's own data
tables. Parameters and constants come from `code/data/`.

**This is not an "Aspen substitute," which is how this file used to describe it.**
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
| [`reaction`](#chemical-reaction-equilibrium--reaction) | $K_a(T)$, the extent of reaction, Ellingham, Gibbs minimization — **replaces CHEMEQ** | 13 |
| [`data`](#the-data-tables--data) | every table above, plus the property and reaction databases | all |
| [`charts`, `ph_chart`, `steam_chart`](#property-charts--charts-ph_chart-steam_chart) | the book's property charts | 3, 5, 6 |

**Not built yet:** electrochemistry, and the ionization half of ch13 (Secs. 13.5–13.7 —
$K_a$ coupled to $\gamma_\pm$ through the ionic strength, which belongs with `electrolytes`).
Reaction equilibrium itself landed 2026-08-18 as [`reaction`](#chemical-reaction-equilibrium--reaction).
See the roadmap.

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

> **The Cp coefficients in `pure_property.csv` are the Reid–Prausnitz–Poling set,
> not the book's Appendix A.II set.** This README previously said Appendix A.II; that
> was wrong, and the same mislabel is in two ch6 notebooks. It matters below room
> temperature: for oxygen the two sets differ by 152 J/mol in H̲ and 1.23 J/(mol·K) in
> S̲ at 73 K, which is the difference between reproducing printed Tables 7.5-1 and
> 7.5-2 and not. Pass `cp=` explicitly when you need the book's numbers. See
> `revision_notes/c07.md`, *ch6 heat-capacity inconsistency* — **decided 2026-08-03**: ch6
> uses the book's Appendix A.II set, via `cp=APPENDIX_A2_CP[...]`, so ch6 and ch7 sit on
> one basis. `pure_property.csv` is left as the RPP table it is.

> **Below 273 K, use `APPENDIX_A2_CP_CRYO` — not `APPENDIX_A2_CP`.** Appendix A.II's
> familiar row is valid **273–1800 K**, and Illustrations 6.4-1, 7.5-1 and 7.5-2 all run
> below that. The extrapolation is *not* benign: for oxygen at 173 K it is 1.235 J/(mol·K)
> (4.2%) low, worth **67 J/mol in H̲ and 0.32 J/(mol·K) in S̲** — a hundred times the error
> from the gas constant or from the 5e's own arithmetic. `APPENDIX_A2_CP_CRYO` is the same
> cubic form refitted to NIST-JANAF over 100–700 K (oxygen and nitrogen only), which cuts
> those to 2.8 J/mol and 0.011 J/(mol·K). Derivation and accuracy table are in
> `data.py`. **Appendix A.II must gain these rows in print** — until it does, the printed
> book does not contain the constants these notebooks use.

> **For oxygen, `from_database` is not just a different Cp — it is a different EOS.**
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

> **A blank in Table 9.4-1 is not a zero, and this is where that bites.** The table is
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
- **`kind="original"`** — **retired 2026-08-12.** The book teaches modified UNIFAC
  only: Table 9.5-2 *is* the Dortmund set and Eq. 9.6-12a carries the `r^0.75`
  modification (`revision_notes/c09.md` **D1**). The `original` branch stays
  unimplemented and `unifac_interactions_original.csv` is retained but unused.

**The book and this module do not compute the same Ψ.** Eq. 9.6-13b prints
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

**`UNIQUACModel` is the molecular model of Sec. 9.5**; `UNIFAC` is the
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

**`check_gibbs_duhem(T)` is worth running.** Where the book prints both a G^ex
expression *and* activity coefficient expressions, they must satisfy
`G^ex/RT = Σ x_i ln γ_i` exactly — so coding both and comparing them is a real test rather
than a tautology. Eight of the nine models have that independent expression to check
against (regular solution does not, and returns `0.0` saying so); **all eight now return
≤6e-15**, worst case UNIQUAC.

**Two of them did not, and both were errata in the book** — Eqs. 9.5-12b and 9.5-18,
in print through five editions and corrected in the 6e manuscript
(`revision_notes/c09.md` §12.1). The classes carry the corrected forms, and each takes a
keyword to reproduce the printed one for comparison. **This is the check to run first on
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

**Table 9.10-1 was audited, and it is correct.** Its two columns are not independent —
Debye–Hückel theory makes $\alpha\rho/\beta^3$ a pure constant, 33.03, at every
temperature. All fourteen rows return 33.04 ± 0.08%. Drop the solvent density and the
ratio instead drifts monotonically by 4.3% from 0 to 100 °C, which is exactly how much
water's density falls: the table reproduces a density curve it never mentions, so its
numbers were computed from theory and transcribed correctly. The two illustrations only
ever use the 25 °C row, so nothing else in the chapter would have caught an error in the
other thirteen.

**The book's 25 °C constant is not self-consistent to four digits.** Illustration
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

**`b` is temperature dependent here**, unlike the van der Waals one-fluid `b`, and
neither `a` nor `b` is a polynomial in mole fraction — both come out of $Q/(1-D)$, and
$D$ carries the whole activity coefficient model. That is the point of the method, not
a side effect.

**$D > 1$ is normal, not an error.** $D$ is dominated by $\sum_i x_i a_i/(b_i RT)$,
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

**The last one is the check that reaches outside §9.9.** Eq. 9.9-11 claims to be the
general cubic-mixture fugacity coefficient, written with composition derivatives instead
of $b_i$ and $2\sum_j x_j a_{ij}$. Feeding it the van der Waals one-fluid rules must
therefore reproduce `PRMixture`, which was verified against **Illustration 9.4-3** — and
it does, exactly.

**And Eq. 9.9-8 was confirmed by brute force.** Taking the equation of state to
$10^{11}$ bar, $\underline A^{ex}$ converges on $C^*[a/b - \sum x_i a_i/b_i]$ to seven
figures while $\underline G^{ex}$ diverges to $-5\times10^{10}$ J/mol — which verifies
the equation, the Peng–Robinson value of $C^*$ (left as Problem 9.31), and the section's
argument for using Helmholtz rather than Gibbs, all at once.

**Eq. 9.9-10a is sign-ambiguous as printed.** Its cross term is the square root of a
product of two pure second virial coefficients, and below the Boyle temperature both are
negative — so the principal root is *positive* where both diagonals are negative, and
$Q$ collapses through zero mid-composition. `cross_matrix` carries the sign of the pures
through the root, which is the only reading consistent with the rule reducing to
$b_i - a_i/RT$ on the diagonal. Worth a footnote in print; see `revision_notes/c09.md`.

## Phase equilibrium — `vle` (Ch. 10) and `lle` (Ch. 11)

`vle` is the γ–φ method: an activity coefficient model in the liquid, ideal gas in the
vapor, with `GammaPhi.bubble_pressure / dew_pressure / bubble_temperature /
dew_temperature / flash` and the diagram generators `pxy`, `txy`, `azeotrope`. Its
solvers carry the **same names and signatures as `PRMixture`'s**, which is what lets one
generator draw both Chapter 10's low-pressure diagrams and Sec. 10.3's high-pressure
ones. It also holds the vapor-pressure correlations (`Antoine`, `ClausiusClapeyron`,
`Wagner`, `Riedel`, `TabulatedPsat`, and `psat_from_database`).

`lle`, added 2026-08-16, is SIS Secs. 11.2 and 11.3. The equilibrium condition loses its
vapor pressure entirely — the pure-component fugacity cancels from both sides — and what
is left is Eq. 11.2-2, $x_i^{\rm I}\gamma_i^{\rm I} = x_i^{\rm II}\gamma_i^{\rm II}$:

```python
from thermo import VanLaar, binary_lle, vlle_binary
from thermo.lle import binary_lle_envelope, consolute_temperature, tie_line_split

xI, xII = binary_lle(VanLaar(2.62, 3.02), 310.95)       # Illustration 11.2-2
P, y, xI, xII = vlle_binary(gp, 310.95)                 # Illustration 11.3-1
```

| what | function | SIS |
|---|---|---|
| Gibbs energy of mixing, its curvature, stability | `gmix_over_RT`, `d2gmix_over_RT`, `is_stable` | 11.2-16, 11.2-9 |
| limits of stability; the common-tangent brackets | `spinodal`, `common_tangent` | 11.2-10, Fig. 11.2-5 |
| binary LLE, and the coexistence curve over T | `binary_lle`, `binary_lle_envelope` | 11.2-2 |
| upper / lower consolute temperature | `consolute_temperature` | 11.2-14 |
| multicomponent two-liquid flash | `lle_flash` | 11.2-2 + 11.2-24 |
| lever rule and stream mixing (mass or mole) | `tie_line_split`, `mix_streams` | 11.2-1b |
| activity coefficient from a solubility | `gamma_from_solubility`, `solubility_at_T` | 11.2-18, 11.2-22 |
| three-phase P and T, immiscible limit, steam distillation | `vlle_binary`, `vlle_temperature`, `immiscible_pressure`, `steam_distillation` | 11.3-3, 11.3-5 |
| P–x and T–x with the LLE region resolved | `pxy_lle`, `txy_lle` | Figs. 11.3-1, 11.3-2, 11.3-4 |
| **the equation-of-state route** — liquid root in both phases | `eos_binary_lle`, `eos_vlle_binary` | 11.2-5, 11.3-2 |

**Equation 11.2-2 has a trivial root.** $x^{\rm I} = x^{\rm II}$ satisfies it at every
temperature for every model, and it is what a solver handed a plain guess converges to.
`binary_lle` is therefore not seeded by a guess: it takes the **lower convex hull** of the
Gibbs energy of mixing — the common-tangent construction of Fig. 11.2-5, done numerically
— which cannot produce a trivial seed, and reports *no split* instead. After the solve it
checks that the tangent it found lies below the curve everywhere. A pair of compositions
that satisfies Eq. 11.2-2 to 1e-12 and is still not the equilibrium state is the ordinary
failure of this calculation, not a hypothetical one.

**Models whose parameters depend on temperature.** Every function that varies T
accepts either an `ActivityModel` or a **callable `T -> model`**, which is what
Illustration 11.2-6 needs (χ = 1473/T).

**Two routes, two signatures, and the difference is physical.** `eos_binary_lle` and
`eos_vlle_binary` take anything with an `ln_phi(x, T, P, phase=...)` — `PRMixture` or the
Wong-Sandler mixing rule — and they take **P as an argument**, because a cubic has no
incompressible-liquid shortcut. `eos_vlle_binary` therefore solves all four unknowns
(both liquid compositions, the vapor composition, the pressure) at once rather than
solving LLE first and taking a bubble point, which `vlle_binary` legitimately can:

```python
from thermo import eos_vlle_binary
P, y, xI, xII = eos_vlle_binary(pr_mixture, 235.65)     # Illustration 11.3-2
```

**The vapor equation has a trivial root of its own**, and it is a different one:
$y = x^{\rm II}$ satisfies it *exactly* wherever the cubic has a single real root, since
then the "vapor" root and the liquid root are the same number. That is the situation
above the three-phase pressure, so a badly seeded solve returns a converged answer with
three liquids in it. `eos_vlle_binary` rejects any solution whose vapor is
indistinguishable from either liquid. `PhiPhiVLE.bubble_pressure` has the same hole
and is **not** guarded — its default seed is far below, so the ordinary path is
unaffected; see `revision_notes/c11.md` §13.

Validated in `code/ch11/validation/lle_module_validation.ipynb` against Illustrations
11.2-1, 11.2-2, 11.2-3, 11.2-4, 11.2-6, 11.2-7, 11.2-8, 11.2-9, 11.3-1, 11.3-2 and
11.3-3, the fourteen-row $x_i\gamma_i$ table of Illustration 11.2-2, the thirty-two
printed values of Illustration 11.3-2, and — with no rounding in them — the analytic
results of Eqs. 11.2-11, 11.2-13, 11.2-14 and Problem 11.2-1(a,b). That pass found five
defects in the chapter, and a sixth in this module; they are filed in
`revision_notes/c11.md`.

**Illustration 11.3-2 needs SIS Table 6.6-1's constants, not `pure_property.csv`'s.**
The database is Reid, Prausnitz and Poling and gives ω = 0.239 for carbon dioxide where
the book prints 0.225 — a 2 % systematic error in the three-phase pressure. Build the
`PengRobinson` objects explicitly when reproducing a printed table, the same way
`APPENDIX_A2_CP` exists for heat capacities.

**Still not here:** Illustration 11.2-5's *Wong-Sandler* curve. Its UNIQUAC parameters
were "fit only to the data at 235.65 K" and are not printed, so that curve has to be
refitted rather than reproduced. The van der Waals half of the same figure runs today
with `eos_binary_lle`.

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

> **The derivative is far more sensitive to the fit order than the fit is.** Adding a
> term barely moves the curve and can move the infinite-dilution partial molar properties
> by ten percent or more — those are the endpoints, the one place the correlation
> extrapolates. `loo_rms` and `scan_order` exist so a notebook can *show* that instead of
> hiding it behind a default.
>
> **The summation identity is exact by construction**, so agreement there is not
> evidence the fit is good — `residuals` is. It is still worth checking, because a
> *graphical* construction does not satisfy it: the 5e's Table 8.6-4 misses it by up to
> 16 J/mol (`code/ch8/validation/`).

Verified against the book's own printed constants: with $a = (-4.0034, -0.17756, 0.54139,
0.60481)\times10^{-6}$ m³/mol, Eqs. 8.6-6a,b reproduce **eleven of the twelve rows of
Table 8.6-2 exactly**, and both infinite-dilution values. The twelfth disagreement is a
**sign error in the printed table** (row $x_1 = 0.9489$), confirmed against the book's own
$\bar V_1$ column.

## The data tables — `data`

**The database is the linchpin.** These CSVs are the book's own tables — effectively a
digitized Appendix A plus the model-parameter tables of ch6–9 — in plain text, so a student
can read them, check one against the NIST WebBook, or add a compound. That is the thing a
licensed black box will not allow, and it is why the data layer is documented here rather
than treated as an implementation detail.

| what | loader | source |
|---|---|---|
| pure-component constants (Tc, Pc, ω, Cp, Antoine) | `load_pure_properties()`, `get_compound(key)` | `pure_property.csv` — the RPP set below |
| ideal-gas $C_p^*$, **as the book prints it** | `APPENDIX_A2_CP`, `APPENDIX_A2_CP_CRYO` | Appendix A.II (+ the cryogenic refit below) |
| PR constants **as the book prints them** | `TABLE_6_6_1` | Table 6.6-1 — oxygen, nitrogen, methane |
| PR binary interaction parameters | `load_pr_kij()`, `pr_kij_matrix(keys)` | Table 9.4-1 — 127 pairs / 20 species |
| UNIFAC subgroups: R, Q, main group, examples | `load_unifac_subgroups()`, `unifac_groups(name)` | Table 9.5-2 — 92 subgroups / 46 main groups |
| UNIFAC group-interaction parameters $a, b, c$ | `load_unifac_interactions()` | legacy `UNIFAC_data.mat` below |
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

**`load_unifac_subgroups` validates on every load, not on request.** A duplicate
subgroup number or an R/Q that disagrees with its main group changes every activity
coefficient it touches, silently — so the integrity rules run each time the file is read.
They are not decorative: they caught duplicate subgroup numbers and silicon names sitting
on the Dortmund cyclic groups in the legacy `.mat` extraction (`revision_notes/c09.md`
§12.3). **Check the other legacy-derived tables the same way as ch10–15 open.**

> **R and Q come from the book; the $a_{mn}$ do not.** `unifac_subgroups.csv` is the
> printed Table 9.5-2, but the *numbers* identifying each subgroup and main group are
> printed nowhere in the chapter, and `unifac_interactions_modified.csv` is keyed on the
> main-group number — so both come from the legacy extraction, which the book cannot check.
> The asymmetry is deliberate and is recorded in `code/data/README.md`.

> **`react_property.csv` stores the heat-capacity coefficients as Appendix A.II
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

**PR puts saturated-liquid density about 12 % high**, the textbook failing of a cubic.
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

## Triangular composition diagrams — `ternary`

Added 2026-08-16, for Figs. **11.2-7** through **11.2-11** — five figures nothing else
in the book could draw. `charts` rules rectangular chart paper and `vle_chart` draws
binary envelopes; a ternary diagram is neither. This module is **draftsmanship only**:
the thermodynamics is in `lle` (`lle_flash`, `tie_line_split`, `mix_streams`) and
nothing here solves anything.

```python
from thermo.charts import use_book_style
from thermo import ternary as tern

use_book_style()
tern.ternary_axes(ax, top="A", left="MIK", right="W", symbol="w", percent=True)
tern.plot(ax, acetone, mik, water)                  # the binodal
tern.tie_line(ax, xI, xII)                          # one tie line
tern.lever_arm(ax, z, xI, xII)                      # Illustrations 11.2-7/-8/-9
tern.check_labels(ax)                               # before saving
```

| function | what it draws |
|---|---|
| `to_xy` / `from_xy` | the coordinate pair, `(top, left, right)` ↔ `(x, y)` |
| `ternary_axes` | triangle, three scales, corner names, edge symbols, two-tier grid |
| `plot`, `scatter` | binodal curves; measured points as open circles |
| `tie_line`, `lever_arm` | one tie line; the lever rule with the feed point marked |
| `point`, `text`, `region_label` | a labeled composition; free text; a `2L`/`3L` region tag |
| `read_construction` | **Fig. 11.2-7** — how a point is read, one line per species |
| `check_labels` | fails if any two labels overlap |

**Which edge carries which species is not arbitrary, and the obvious guess is
wrong.** A species' fraction is proportional to the perpendicular distance from the
edge **opposite** its corner — that is the geometry. But the printed **ruler** is on an
**adjacent** edge, the one ending at its own corner, so each scale runs *toward* the
thing it measures. The first draft put the apex species on the left edge; the 5e's own
Figs. 11.2-7 and 11.2-8 put it on the right, and the convention was taken from those
two figures rather than invented.

**`check_labels` is `check_print_art.py`'s blind spot, made reusable.** That gate
reads color and font family and cannot see two labels sitting on top of each other —
and a triangular diagram crowds three scales, three corner names and three edge symbols
around one small perimeter. It caught `x_C = 0.25` running through the `0.2` tick at a
0.4 px overlap, which looked deliberate at thumbnail size. Worth copying into the other
figure notebooks.

**Not `mpltern`.** The author's practice files use it, and it is a good library. It
would be a new student-facing dependency in a `pyproject.toml` that is deliberately a
short list of floors with no lock file, to keep resolving years from now — and the
house rules (pure black ink, Computer Modern, nothing below 7 pt) are easier to obey
owning ~200 lines of geometry than bending a third-party projection. Fig. 11.2-7 is not
a data plot at all.

**Validation.** The geometry is exact: the three corners land on their vertices, all
three edges are unit length, `from_xy(to_xy(c)) == c` to machine precision on 500 random
compositions, and lines of constant composition come out horizontal / 60° / 120° to
within 6e-17. The drawing was checked against the two figures it imitates, and the
lever-rule path reproduces **Illustration 11.2-8**: feed 15 % A / 75 % MIK / 10 % W, and
3.720 / 0.280 kg against the book's 3.721 / 0.279.

## Chemical reaction equilibrium — `reaction`

**This module is the replacement for CHEMEQ** and for "the chemical equilibrium constant
calculation programs of Appendix B.I or B.II" — the Visual Basic executables the 5e shipped
on its website. It was written from that source, not reverse-engineered from the printed
output: `chemeq_source.zip` was unpacked and `frmChemEq.frm` read, and the two routines
that carry the whole program (`Ka_T298` and `cmdCalKa_Trange_Click`) are the analytic van 't
Hoff integration below. CHEMEQ's own 99-row species database was read too, and it is
identical to `code/data/react_property.csv` — every row, every column.

```python
from thermo import Reaction, equilibrium_extent

rxn = Reaction.parse("0.5 N2 + 1.5 H2 = NH3").check_balance()
rxn.Ka(450.0)                                   # 1.218, as Illustration 13.1-4 prints
rxn.table([500, 600, 700, 800])                 # what CHEMEQ printed
equilibrium_extent(rxn, 450.0, {"N2": 0.5, "H2": 1.5}, P=4.0)   # X = 0.6306
```

Coefficients are **negative for reactants**, and species are named by the formulas of
Appendix A.IV (`'N2O4'`, `'H2O(g)'`, `'nC6H14(l)'`).

| what | entry point | SIS |
|---|---|---|
| $K_a(T)$ four ways — full Cp polynomial, constant $\Delta C_P$, constant $\Delta H$, quadrature | `Reaction.ln_Ka`, `Reaction.Ka` | 13.1-18, 13.1-22 |
| $\Delta_{rxn}G^\circ(T)$, $\Delta_{rxn}H^\circ(T)$, $\Delta_{rxn}S^\circ(T)$ | `Reaction.delta_G`, `.delta_H`, `.delta_S` | 13.1-21 |
| the CHEMEQ output table | `Reaction.table` | — |
| the chapter's mass-balance tables | `Reaction.balance_table`, `.moles` | Table 13.1-1, Eq. 13.1-5 |
| the extent of reaction, at fixed $T,P$ or fixed $T,V$ | `equilibrium_extent` | 13.1-19 |
| $K_c$, $K_x$, $K_y$, $K_p$ | `Reaction.K_ratio` | Table 13.1-3, 13.1-23 |
| $K_\nu$ from an EOS, $K_\gamma$ from a $\gamma$ model | `K_nu_from_eos`, `K_gamma_from_model` | 13.1-23d |
| several reactions at once | `multireaction_extents` | 13.3 |
| the same state without choosing reactions | `gibbs_minimization` | 13.3 |
| $G$ against the extent — the figure the chapter opens with | `gibbs_curve` | Fig. 13.1-1 |
| the Ellingham construction | `ellingham` | 13.2, Fig. 13.2-3 |
| the adiabatic reaction temperature | `adiabatic_reaction_temperature` | 14.3-10a, Ill. 14.3-2 |
| that Sec. 14.3's two energy-balance forms are one equation | `energy_balance_forms_agree` | 14.3-10a |

### Guards, and what they caught

**`adiabatic_reaction_temperature` brackets DIRECTIONALLY, and that is not a nicety.** The
energy-balance branch passes through zero at the inlet temperature by construction, so
`T = T_in` is *always* a root of the difference between the two branches — an undirected
bracket can converge on it, report a temperature rise of zero, and look perfectly converged.
The default bracket therefore starts one degree away from `T_in`, on the side the sign of
$\Delta_{rxn}H$ requires: an exothermic reaction heats its own effluent, an endothermic one
cools it. If the branches do not cross inside the bracket it raises, with both end-point
residuals in the message, rather than returning its last iterate.

**`energy_balance_forms_agree` asserts an identity rather than trusting it.** Eq. 14.3-10a is
usually written as the sensible heat of the *feed* plus $\Delta_{rxn}H$ at the *effluent*
temperature; the total-enthalpy form heats the *effluent* composition and takes
$\Delta_{rxn}H$ at the *inlet*. Those are two paths between the same pair of states, so they
must give the same extent — and checking it costs one call, whereas discovering later that a
sign or a limit was wrong costs a chapter. It was this check that made it safe to report
Illustration 14.3-2's adiabatic temperature as 701.2 K against two different printed values.


**`equilibrium_extent` never takes a bare initial guess.** The extent is bounded by
exhaustion — no species may end with a negative mole number — which gives a closed bracket,
and the solve is a bracketed Brent on it. Monotonicity of $\prod a_i^{\nu_i}$ across that
bracket is *checked and reported*, not assumed, and on failure the solver raises rather than
returning its last iterate. Illustration 13.1-5 is the chapter's own example of an equation
whose extra roots converge perfectly well and mean nothing.

**A pure condensed phase stays out of the mole-fraction denominator.** Mark it
`phases={"C": "solid"}` and it has unit activity, does not dilute the fluid, and does not
bound the extent. Getting this wrong moves the answer by 0.11 in mole fraction — and it is
the same error Illustration 13.3-1's own note attributes to Aspen Plus.

### Validated against Chapter 13 — 33 of 34 checks

`code/ch13/validation/reaction_module_validation.ipynb`. Illustration 13.1-3 reproduces to
**every printed digit** (including the internal identity that 56 189 and −6758.4 are the same
number twice), Illustration 13.1-4's three extents come out exactly (0.6306, 0.4072, 0.5741),
and Illustration 13.2-2's decomposition-pressure table reproduces across **25 orders of
magnitude** to 0.66 %. Ten of the checks are identities with no rounding in them at all —
the analytic integration against quadrature, $\Delta_{rxn}G^\circ(T)$ by two routes that
share no arithmetic, the minimum of $G(X)$ against the root of $K_a=\prod a_i^{\nu_i}$, and
the coupled-extent solve against Gibbs minimization.

**The one failure is a finding about the book, not the module**: Illustration 13.1-8's
printed $K_a$ table is not reproducible from Appendix A by any route the chapter offers, and
it contradicts the book's own $K_a(450\ \mathrm{K})$ for the same reaction. It is filed as
Q8 in `revision_notes/c13.md`, and nothing here is tuned to match it.

## Bioreactors — `bioreactor` (Sec. 15.7)

Eleven illustrations, one object. Sec. 15.7 balances **atoms** rather than species, because
its products are cells and cells have no molecular formula — what they have is an elemental
analysis, Roels' average biomass `CH1.8O0.5N0.2`. Everything is per **C-mole**: the species'
formula divided through by its carbon count, so glucose C₆H₁₂O₆ is `CH2O` and ethanol
C₂H₅OH is `CH3O0.5`.

```python
from thermo.bioreactor import CMole, Fermentation, from_table

f = Fermentation(substrate=from_table("glucose"),
                 biomass=CMole("CH1.8O0.5N0.2"),
                 nitrogen=from_table("ammonia"),
                 product=from_table("ethanol"))
f.solve(Y_B_S=0.14, Y_O2_S=0.0)      # Illustration 15.7-3
f.heat_load(0.14, 0.569, 0.028)      # Eq. 15.7-8b, kJ per C-mole
f.second_law(0.14, 0.569, 0.028)     # Eq. 15.7-19, as a Constraint
```

| SIS | what | entry point |
|---|---|---|
| Eq. 15.7-9 | the generalized degree of reduction ξ | `CMole.xi` |
| Eqs. 15.7-1 … 15.7-4 | the four atom balances in yield-factor form | `Fermentation.solve` |
| Eq. 15.7-12 | the 4C + H − 2O combination, i.e. the oxygen balance as ξ | `Y_O2_from_xi` |
| Eq. 15.7-8b | the energy balance, from heats of **combustion** | `heat_load` |
| Eq. 15.7-10 | the energy regularity principles, 112 ξ and 110.9 ξ | `regularity_G`, `regularity_H` |
| Eq. 15.7-19 | the second-law constraint | `second_law` |
| Eq. 15.7-26a | the entropy generated per C-mole | `entropy_generated` |
| Table 15.7-2 | ξ, Δ_cG, Δ_cH and the molecular formula, 53 compounds | `TABLE_15_7_2`, `from_table` |

**Three design decisions, each of them a guard.**

*The balances are solved as a linear system, not substituted by hand.* Each illustration
fixes a different pair of the six yield factors, and Sec. 15.7 walks the substitutions in a
different order every time — which is where its printed slips are (Illustration 15.7-4's
oxygen balance is set as `1 + 2 Y_W/S + Y_W/S`, with no oxygen term in it). `solve` assembles
`A y = b` over whichever factors are open and **raises** on an over- or under-specified
problem instead of returning one arbitrary solution of many.

*The second law returns a `Constraint`, not a boolean* — both sides, the slack and the
efficiency ratio — because the three illustrations that use Eq. 15.7-19 each want a
different one of those, and a `True` would have thrown away all three.

*Table 15.7-2 ships each compound's molecular formula, which the printed table does not
carry.* That is what makes its ξ column self-checking: Eq. 15.7-9 regenerates all 53 values
from the formulas alone, and all 53 agree. The printed table's two right-hand columns are
derived and are regenerated rather than transcribed — the Table 14.6-1 contract, written
after four of that table's cells turned out not to follow from their own data.

**What it found.** Illustration 15.7-6 quotes ammonia's heat of combustion as 383.0 where
Table 15.7-2 gives 348.1 and the arithmetic uses 348.1; Illustration 15.7-8 prints
"6 × 14.5 = 81" where the heat load one line above is 13.5, and a second-law sum of 450.2
where the ratio two clauses later uses the correct 450.5; Illustration 15.7-11's entropy
balance multiplies the ethanol yield by 684.5, which is ethanol's Δ_cH where Δ_cG = 659.5 is
needed and is what the printed answer used; Eq. 15.7-26b is missing the 112 kJ/C-mole that
makes it an entropy; and Figure 15.7-3's vertical axis runs −200 to −700 kJ/C-mole where the
equation beside it runs −443.6 to 0.

## Provenance

Refactored from the validated ch6/ch7 notebooks (PR), the legacy MATLAB `unifac.m`
(a faithful translation of its `calc_coeff`), and the ch3 chart notebooks (`charts`,
`ph_chart`, `steam_chart`). `pr_mixture`, `activity_models` and `fitting` were written
from the book's own equations and verified against its printed numbers, not ported.
See Appendix B and `code/data/README.md`.

**Writing them that way is what found the errata.** Making two printed expressions for
the same quantity compute against each other — $\underline G^{ex}/RT = \sum_i x_i\ln\gamma_i$
for the activity models, $P\hat V = \hat H - \hat U$ for the steam tables, the summation
identity for Redlich–Kister — has now caught errors in Eqs. 9.5-12b and 9.5-18, in Table
8.6-2, and in Tables 9.5-2 and A.III. A model that merely runs has not been checked.
