# `thermo` package roadmap — a toolkit for CBET 6e, ch6–15

A living plan for building out the Python `thermo` package through the whole book.
Derived from the **4th-edition Appendix B** program catalog
(`legacy-4e/Sandler_4e_app_B.pdf`), read as a specification of *what the historical
computing suite could do*, re-expressed as Python modules.

**Status 2026-08-12 — chapters 6–9 are complete.** T13 (`electrolytes.py`) and T14
(`wong_sandler.py`) both landed; thirteen modules, the property and model-parameter
databases, and notebook sets in `ch6`–`ch9`. **Next is ch10** — the VLE drivers, which
are the first application of `wong_sandler` and the first place UNIFAC gets a P–xy
driver.

## The organizing principle (from 4e Appendix B)

The 4e suite shipped the **same capabilities four times** — Windows Visual Basic
(B.I), DOS Basic (B.II), MATHCAD worksheets (B.III), and MATLAB (B.IV). Across all
four, the software was organized **by thermodynamic model/tool, not by chapter**:

- **Property** — a pure-fluid constants database (600+ compounds).
- **Peng–Robinson (pure)** — properties, fugacity, departures, vapor pressure, P–V/P–H/P–S plots, VLE envelope.
- **Peng–Robinson (mixtures)** — vdW one-fluid mixing rules + $k_{ij}$, component fugacities, bubble/dew/flash.
- **UNIFAC** — activity coefficients, $G^{ex}$, P–xy.
- **Activity-coefficient models** — van Laar, Margules, Wilson, NRTL, Redlich–Kister (+ fitting from $G^{ex}$ data).
- **SRK / RK** — cubic-EOS variants (an $a(T)$ swap on PR).
- **Chemical equilibrium (CHEMEQ)** — $K_a(T)$ from $\Delta_f G$, $\Delta_f H$, $C_p$.
- **Fitting utilities** — Antoine ($P^{vap}$) fits, partial molar properties via Redlich–Kister.

Each tool supported a **consistent operation vocabulary** (properties → fugacity →
bubble/dew/flash → plots) and was **reused across many chapters**. Chapters were
*application sites*, not separate programs.

**This is the model for the 6e Python code.** The `thermo` package is the toolkit
(the "programs"); per-chapter notebooks are the worked *applications* that call it —
exactly the ch6 pattern, where `PR_*_thermo.ipynb` notebooks apply
`thermo.PengRobinson`. See [`../ch6/README.md`](../ch6/README.md) for the
self-contained/`_thermo`-twin convention and [`README.md`](README.md) for the API.

## Positioning: general-purpose computing vs specialized process modeling

Python and Aspen are **two complementary parts of a practicing engineer's toolbox — not a
learning tool and a "real" one.** Code is itself part of engineering practice (increasingly
so with accessible languages and coding agents), not merely a step toward it. The axis that
separates them is *general-purpose vs specialized*, not *learning vs practice*.

- **`thermo` (Python + Jupyter) — general-purpose scientific computing.** Both how we
  *learn* the thermodynamics in this book *and* a broad engineering skill used constantly
  (analysis, data, custom calculations, prototyping). Open, inspectable, hackable — read
  every step, run it anywhere, extend it with AI coding tools, accessible with nothing more
  than a Python install (no license, no gatekeeper).
- **Aspen — specialized process modeling.** The right tool for flowsheets, unit operations,
  and plant-scale design; the book treats it as first-class (Appendix C ships Aspen
  illustration input files). Note the access reality: at UD it is delivered through
  **AppsAnywhere**, a real barrier — a further reason general-purpose computation, not a
  process simulator, is the everyday tool.
- **NIST Chemistry WebBook — authoritative data lookup.** Use it to *verify and extend*
  the local property database.

**Scope implication.** The Python stack aims for **completeness on the book's own
calculations** — a clear, verifiable path for its examples, illustrations, and homework
across ch6–15 (general-purpose computation). Aspen owns the genuinely process-scale work
(full-process/flowsheet design) and is introduced where the book teaches that, not as the
default for routine calculations.

**The property database is the linchpin of the foundation.** `pure_property.csv` /
`react_property.csv` are the book's own curated data — effectively the **digitized
Appendix A** (A.II ideal-gas heat capacities, A.IV enthalpies/Gibbs of formation) and the
same tables that backed the legacy Visual Basic *Property* program. Because they are plain
CSV and the models are plain Python, a student (with or without an AI assistant) can add a
compound from the NIST WebBook, cross-check a value against Aspen, or swap in a different
model — exactly what a licensed black box will not let them do. Keeping this layer
authoritative and easy to extend is what makes the Python stack a foundation students can
build on rather than a demo.

## Capability → module → status → chapters

| Capability — from the 4e suite unless marked *new* | Target module | Status | Chapter(s) |
|---|---|---|---|
| Property database (constants, $C_p^*$, $P^{vap}$, $\Delta_f G/\Delta_f H$) | `data.py` | done (`get_compound`, CSV loaders; `react_property.csv` has formation data) | all |
| Model-parameter tables digitized from the book | `data.py` | **done 2026-08-12** — Table 9.5-2 UNIFAC subgroups (92 over 46 main groups; `load_unifac_subgroups`, `unifac_groups`), PR $k_{ij}$ (127 pairs / 20 species; `load_pr_kij`, `pr_kij_matrix`), §8.6 mixing data (`load_mixing_data`) | 8, 9, 10 |
| PR pure: props, fugacity, departures, $P^{vap}$ | `peng_robinson.py` | done | 6, 7 |
| PR pure: P–H charts, dome + isotherm/isentrope/isochore/quality families | **`ph_chart.py`** (new) + **`charts.py`** (new) | **done 2026-08-08** (`ChartFluid`, `ph_chart`; reproduces Figs. 3.3-2/3.3-3 bit-for-bit) | 3, 5, 6 |
| Steam charts from Appendix A.III (Mollier, T–S) | **`steam_chart.py`** (new) | **done 2026-08-08** (`SteamTables`, `temperature_entropy`; reproduces Fig. 3.3-1b bit-for-bit) | 3, 5 |
| PR pure: P–S plots, VLE (T–V) envelope | `ph_chart.py` | to add (P–V isotherms exist as a notebook) | 6, 7 |
| **PRSV** (Stryjek–Vera $\kappa(T)$, Eqs. 7.5-1/7.5-2) | `peng_robinson.py` (`kappa1=`) | done (verified against Illustration 7.5-3) | 7 |
| **van der Waals** | **`van_der_waals.py`** (new) | done (verified against the Fig. 7.5-2 notebook) | 6, 7 |
| Shared cubic machinery: roots, fugacity, spinodals, $P^{vap}$ | **`cubic.py`** (new) | done (`CubicEOS`; bracketed Fig. 7.5-1 solver) | 6, 7 |
| SRK / RK variants (a(T) swap) | `cubic.py` — subclass as vdW and PR now do | to add | 6 |
| PR **mixtures**: mixing rules + $k_{ij}$, mixture fugacity | **`pr_mixture.py`** (new) | done (`PRMixture`, vdW one-fluid + `ln_phi`, verified) | 9, 10 |
| Bubble/dew T & P; iso-T flash | `pr_mixture.py` + phase-eq driver | bubble/dew P&T + Rachford-Rice flash done; iso-H/S flash | 10 |
| UNIFAC: $\gamma$ (modified) | `unifac.py` | works; R/Q now from the book's Table 9.5-2 | 9, 10 |
| ~~UNIFAC: original variant~~ | `unifac.py` | **retired 2026-08-12 (EMF)** — the book teaches modified (Dortmund) UNIFAC; the `original` branch stays unimplemented. See `revision_notes/c09.md` D1 | — |
| UNIFAC: $G^{ex}$–x, P–xy VLE driver | `unifac.py` (+ driver) | to add | 10 |
| Activity models: van Laar, Margules (1- & 2-constant), Redlich–Kister, Wilson, NRTL, Flory–Huggins, UNIQUAC, regular solution (+ $G^{ex}$ fitting) | **`activity_models.py`** | **done 2026-08-12** — one `gamma(x,T)`/`gex(x,T)` interface; multicomponent Wilson, NRTL, UNIQUAC and regular solution; `check_gibbs_duhem` on every model | 9, 10, 11 |
| **Combined EOS + $G^{ex}$: the Wong–Sandler mixing rule** (§9.9) — *new, not in the 4e suite* | **`wong_sandler.py`** | **done 2026-08-12 (T14)** — Eqs. 9.9-9/-10/-11/-12/-13, both combining rules, `GexFromUNIFAC`, and `check_boundary` for the two conditions the rule is derived from. **Verified without a single printed number**: the boundary conditions hold to 1e-15, Eq. 9.9-8 and $C^*$ confirmed by taking the EOS to $10^{11}$ bar, and Eq. 9.9-11 reproduces `PRMixture` exactly on vdW one-fluid input | 9, 10 |
| **Electrolytes: Debye–Hückel and its extensions** (§9.10) — ionic strength, mean ionic activity coefficients — *new, not in the 4e suite* | **`electrolytes.py`** | **done 2026-08-12 (T13)** — Eqs. 9.10-15/17/18 on one `DebyeHuckel` class, `Electrolyte` with an electroneutrality check, Table 9.10-1, mixed-electrolyte `ionic_strength`. Recovers Illustration 9.10-2's fitted δ = 0.137. **Its audit proved Table 9.10-1 correct** (α ρ/β³ = 33.03 on all fourteen rows) and found the printed art for `c09uf003` never existed. Reused by ch15 | 9, 15 |
| VLLE / LLE / gas solubility / osmotic (incl. Margules VLLECALC) | phase-eq driver on `activity_models` + `pr_mixture` | to build | 11 |
| Chemical equilibrium $K_a(T)$ (CHEMEQ); ionization / pH | **`chem_equilibrium.py`** (new) | to build | 13, 15 |
| Adiabatic flame / reaction temperature | application (chem_eq + energy balance) | to build | 14 |
| Electrochemistry — Nernst / cell potentials (fuel cells, batteries) | **`electrochem.py`** (new) | to build | 14 |
| Fitting utils: partial molar via R–K | **`fitting.py`** (new) | **done 2026-08-09** (`RedlichKister`, `tangent_intercepts`; reproduces printed Table 8.6-2 exactly) | 1, 8 |
| Fitting utils: Antoine ($P^{vap}$) fits | `fitting.py` | to add | 1, 7 |

Legend: done · partial · to build · next up, and blocking a notebook the chapter needs

## Target module layout

```
thermo/
  data.py              property + reaction databases (code/data/*.csv)
  charts.py            chart craft: weights, grids, label placement, book typography
  ph_chart.py          P-H charts from a cubic EOS (Figs. 3.3-2/3.3-3, 5.1-3, c06uf002)
  steam_chart.py       Mollier + T-S charts from Appendix A.III (Figs. 3.3-1a/b, c05uf001)
  cubic.py             CubicEOS base: roots, fugacity, spinodals, P^vap
  peng_robinson.py     pure-fluid PR + PRSV  (add P–H/P–S, VLE envelope)
  van_der_waals.py     pure-fluid vdW (Fig. 7.5-2 curve a)
                       SRK / RK: two more subclasses of CubicEOS
  fitting.py           Redlich-Kister correlation + partial molar properties (Sec. 8.6)
                       Antoine fit
  pr_mixture.py        vdW mixing rules, mixture fugacity, bubble/dew/flash
                       iso-H / iso-S flash
  activity_models.py   Margules 1&2, Redlich-Kister, van Laar, Wilson, NRTL,
                          Flory-Huggins, UNIQUAC, regular solution (2026-08-12)
  unifac.py            modified done (the book's model — D1); original retired; VLE driver to add
  electrolytes.py      Debye-Huckel + extensions, Table 9.10-1 (Sec. 9.10); ch15 reuses it
  wong_sandler.py      Wong-Sandler mixing rule + GexFromUNIFAC (Sec. 9.9); ch10 uses it
  chem_equilibrium.py  Ka(T) from ΔfG, ΔfH, Cp; ionization / pH
  electrochem.py       Nernst / cell potentials (fuel cells, batteries)
```

**Both names settled 2026-08-12.** `electrolytes.py`, not `debye_huckel.py`: ch15 layers
ionic strength, pH and Gibbs–Donnan on the same module, and the *extensions* of §9.10 are
not Debye–Hückel proper. `wong_sandler.py` keeps the name the book gives the rule.

**`electrolytes` does not inherit `ActivityModel`, and that is the precedent to carry.**
Sec. 9.5's models map mole fractions to one coefficient per species; this one maps ionic
strength to a single mean for the salt, because electroneutrality makes the individual ion
coefficients unmeasurable. **Where the physics changes the independent variable, the
interface changes too** — a shared interface would be a lie about what is computed. Expect
the same question at `chem_equilibrium` (extent of reaction) and `electrochem` (potential).

## Per-chapter build-out (ch8–15)

Titles below are the book's actual chapter titles (from the 5e table of contents,
`book_markup/ftoc.pdf`; carried into 6e unchanged). Section numbers pinpoint where each
capability is introduced.

- **ch8 — The Thermodynamics of Multicomponent Mixtures.** **Done 2026-08-09.** Partial
  molar properties (§8.6 experimental determination of partial molar volume & enthalpy);
  Redlich–Kister fitting of mixing data (the 4e `PRTLMOLR` worksheet) → `fitting.py`. Three
  notebooks in `code/ch8/`, six printed QR keys, `codes:`/`figures:`/`tables:` wired into
  `ch08.yaml`. The chapter's two data tables are digitized in `code/data/mixing_*.csv`;
  Tables 8.6-2 and 8.6-4 are **generated**, not stored.
- **ch9 — Estimation of the Gibbs Energy and Fugacity of a Component in a Mixture.**
  **The model chapter, and two thirds done as of 2026-08-12.** Built: `pr_mixture.py`
  (§9.4/§9.7 fugacity of a species in a mixture via an EOS + mixing rules, with the $k_{ij}$
  table), `activity_models.py` (§9.5 correlative γ models), UNIFAC (§9.6 predictive —
  modified/Dortmund only, D1), `electrolytes.py` (§9.10), `wong_sandler.py` (§9.9), and
  **all nine notebooks in `code/ch9/`**, executing clean, with eight recomputed figures and
  **all eleven QR keys** in the manuscript. **The chapter's code is done.**
  **What this chapter proved about the method.** Making two printed expressions for the
  same quantity compute against each other found **two equation errors in eleven pages that
  had survived five editions** (Eqs. 9.5-12b, 9.5-18), and digitizing Table 9.5-2 found three
  more defects plus a live bug in our own legacy-derived data. Run the same sweep in ch10–15
  as each module lands. See [`../../revision_notes/c09.md`](../../revision_notes/c09.md) §§12–13.
- **ch10 — Vapor-Liquid Equilibrium in Mixtures.** **Next.** Bubble/dew/flash drivers
  along the chapter's three tracks: ideal (§10.1), low-pressure γ-based (§10.2, on
  `activity_models`/UNIFAC), high-pressure φ–φ EOS method (§10.3, on `pr_mixture` **and
  `wong_sandler`**). P–xy diagrams. §10.3's four demonstrations are all acetone/water and
  all use the Wong–Sandler rule, including the two falsifiable claims T14 left untested: that
  the van der Waals rules produce a **false liquid–liquid split** at 298 K where Wong–Sandler
  does not, and that $k_{12}=0$ is not ideal-solution behavior. Both need the bubble-point
  driver ch10 builds.
- **ch11 — Other Types of Phase Equilibria in Fluid Mixtures.** Gas solubility (§11.1),
  LLE (§11.2), VLLE (§11.3), distribution coefficient (§11.4), osmotic equilibrium
  (§11.5) — phase-eq drivers on the same models (incl. Margules VLLE).
- **ch12 — Mixture Phase Equilibria Involving Solids.** Solid solubility / SLE (§12.1),
  freezing-point depression (§12.3), solid mixtures (§12.4), environmental phase
  behavior (§12.5).
- **ch13 — Chemical Equilibrium.** `chem_equilibrium.py` (Ka(T) from `react_property.csv`):
  single-phase (§13.1), heterogeneous (§13.2), multiple reactions (§13.3), combined
  chemical + phase (§13.4), ionization / acidity (§13.5–§13.6).
- **ch14 — The Balance Equations for Chemical Reactors, Availability, and
  Electrochemistry.** Adiabatic reaction / flame temperature (§14.3, the 4e
  `ADIABATIC FLAME` worksheet), availability in reacting systems (§14.5), and
  **electrochemistry — fuel cells & batteries (§14.6–§14.7), a new capability**
  (`electrochem.py`: Nernst / cell potentials).
- **ch15 — Some Additional Biochemical Applications of Thermodynamics.** Solubility vs pH
  (§15.1), ionic strength (§15.2), ligand binding (§15.3), denaturation (§15.5),
  ATP-ADP coupling (§15.6), Gibbs–Donnan / membrane potentials (§15.8) — layered on
  `chem_equilibrium` (ionization), `activity_models`, and **`electrolytes` (T13)**, whose
  ionic-strength and mean-ionic-activity machinery §15.2 and §15.8 use directly.

## A note on where the chart modules landed

This roadmap originally filed *"PR pure: P–H / P–S plots"* under `peng_robinson.py`. They
went to **dedicated modules** instead, for two reasons. `peng_robinson.py` is 179 lines of
equation-of-state mathematics and adding 400 lines of matplotlib would bury it — and the
**drawing craft is shared with the steam charts, which have no equation of state at all**,
so it could not live under PR in any case. `ChartFluid` (the per-kilogram wrapper on the
chart's datum) is arguably EOS work, but it is a *chart* convention — kJ/kg, and H = S = 0
at the triple-point liquid — so it sits with the chart it serves.

**The organizing principle is unchanged and is in fact better served**: the package holds
the method, the chapter notebooks are the application sites. What forced this was the
decision that a figure's notebook lives in the figure's own chapter, which made three
derived figures (5.1-3, `c05uf001`, `c06uf002`) into ch5 and ch6 applications of chart
machinery that had been trapped inside a ch3 notebook.

## Conventions to carry forward (from ch6)

- **Notebooks apply the toolkit; the package holds the method.** Each example/homework
  gets a self-contained version (method shown inline) and a `_thermo` twin
  (`from thermo import ...`), per [`../ch6/README.md`](../ch6/README.md).
- **`from_database` for parameters.** Constants, $\omega$, $C_p$, $\Delta_f G/H$ come
  from `code/data/*.csv` — one source of truth (numbers may differ slightly from
  hand-entered book-table values; regenerate answer keys from the assigned version).
- **SI units throughout** (T in K, P in Pa, molar quantities per mol).
- **Validate every module** against a known limit or book illustration before wiring
  it into notebooks (as done for PR pure and modified UNIFAC — see `README.md`).
- **Wire notebooks into the companion site** via `6e_companion_site/source/chapters/chNN.yaml`
  `codes:` blocks, then rerun `6e_companion_site/build.py`.

## Related

- [`README.md`](README.md) — current package API and validation notes.
- [`../ch6/README.md`](../ch6/README.md) — the notebook naming + twin convention.
- [`../../TEXTBOOK_PLAN.md`](../../TEXTBOOK_PLAN.md) — overall revision plan (Appendix B, backlog).
- [`../../revision_notes/bapp02.md`](../../revision_notes/bapp02.md) — **Appendix B** drafting spec: this positioning (general-purpose computing vs specialized process modeling) and the data provenance are the print-facing narrative of what this roadmap builds.
- [`../ch9/README.md`](../ch9/README.md) — the ch9 notebook map (the mixture-model chapter).
- [`../../revision_notes/c09.md`](../../revision_notes/c09.md) — where T13 and T14 are specified.
- `legacy-4e/Sandler_4e_app_B.pdf` — the source catalog this roadmap is derived from.
