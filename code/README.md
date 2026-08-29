# Scientific computing stack — code for *Chemical, Biochemical, and Engineering Thermodynamics*, 6e

Python and Jupyter notebooks that accompany the 6th edition (Sandler & Furst). They
modernize the legacy Visual Basic, MATLAB, and Mathcad companion programs of earlier
editions into an open, license-free stack. This repository is one of the book's
**companion-website deliverables**; each chapter's notebooks are also linked from that
chapter's page on the companion site.

- **Distribution:** primary copy in the Gitea repo `thermohub`
  (`https://lem.che.udel.edu/git/furst/thermohub`), with a GitHub mirror. The Wiley
  companion website hosts the notebooks for direct download.

## Requirements

Python 3.10+ with Jupyter and the scientific stack — `numpy`, `scipy`, `pandas`,
`matplotlib`, `jupyterlab`, `openpyxl`. These are declared in [`pyproject.toml`](pyproject.toml),
so you do not have to install them one at a time.

### Recommended: uv

[uv](https://docs.astral.sh/uv/) is a single tool that installs Python itself, creates the
virtual environment, and installs the libraries — replacing the older
`conda`/`venv` + `pip` routine. Install it once:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh          # macOS / Linux
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"   # Windows
```

Then, from this directory, one command does everything:

```bash
uv run jupyter lab
```

The first run creates `.venv/` and installs the stack (a few seconds); later runs just
start Jupyter. You never activate the environment by hand — `uv run <command>` uses it.
To add a library for your own work: `uv add <package>`.

<details>
<summary>Without uv (python.org + venv + pip)</summary>

```bash
python3 -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install numpy scipy pandas matplotlib jupyterlab openpyxl
jupyter lab
```

</details>

Nothing to install at all? Self-contained notebooks run in **Google Colab** — see the
chapter pages on the companion website.

## Design: two computational levels

The code sits at two deliberately different levels — a ladder from a one-off script to a
reusable tool:

- **Notebooks — elementary and procedural.** Functions and scripts, standard numerical and
  plotting libraries, read top to bottom, mirroring the book's equations. This is the level
  most examples and homework use, and where a student learning both thermodynamics *and*
  Python should start.
- **The `thermo` package — a step up, object-oriented.** It uses classes (`PengRobinson`,
  `UNIFAC`) to bundle a compound's parameters with the methods that act on them:
  `pr = PengRobinson.from_database("ethane")` returns an object that carries ethane's
  $T_c$, $P_c$, $\omega$, and $C_p$ and knows how to compute its molar volume, fugacity,
  and departure functions. It behaves like a small reusable library — the natural step from
  a one-off notebook to a tool you reuse across problems, and a concrete example of *when*
  object orientation earns its extra abstraction. It stands as a more advanced computational
  tool than the straight-up notebooks.

The `_thermo` twin notebooks in each chapter *apply* the package; their self-contained
counterparts build the same method procedurally — see [`thermo/README.md`](thermo/README.md)
for the package's design.

## Notation

Molar quantities use an **underbar** ($\underline{V}$, $\underline{H}$), specific
(per-mass) quantities a **caret** ($\hat{V}$), total quantities a **plain** symbol, and
an **overbar** is reserved for a partial molar property ($\bar{G}_i$, introduced in
Ch. 8). This matches the book. Per-chapter notation files for an AI coding assistant are
published on the companion website.

## Prose style in the code

Everything a reader sees — markdown cells, docstrings, comments, and the text a cell
prints — is plain prose with **no emoji**. Triage markers must not cross into a
notebook or a module: a notebook is read on the
companion site and on GitHub, and its stored output is part of what a student reads. Where
a warning needs weight, say so in words (`NOT`, `must`, bold) rather than with a glyph.

## Layout

### Chapter notebooks

**`ch1/` — getting started with Python (a skills ladder)**

Work through these **in order**. The first two are read rather than run; each of the last
four is the anchor for one printed problem (1.9–1.12) and ends with a *Your turn* cell.

| Notebook | What it teaches |
|---|---|
| `Math_review.ipynb` | the calculus and algebra the book leans on, with a note beside each result saying where it turns up. Markdown only — no code to break |
| `Getting_started_Python_Jupyter.ipynb` | the bridge from an introductory CS Python course to the scientific stack: how a notebook executes, why a list will not do arithmetic and a NumPy array will, the four steps behind every plot in this book, and the mistakes that cost the most time |
| `vdW_molar_volume_CO2.ipynb` | root finding — rearrange to $f(\underline{V})=0$ and hand it to `brentq`. Finds where CO₂ stops behaving ideally (≈2.4 bar at 300 K), checks it against a hand calculation, and shows what three roots do to a root finder |
| `Thermometer_water_mercury.ipynb` | tabulated data, plotting, and numerical differentiation on **unevenly spaced** points (water/mercury volumes, Perry's Handbook). Why mercury works and water cannot, and how many digits survive differencing a table |
| `Heat_capacity_fitting_CO2.ipynb` | curve-fit the Appendix A.II ideal-gas $C_P^*$ polynomial to NIST CO₂ data and compare with the book's coefficients. Residuals, and why extrapolation fails. Paired with ch3's `Heat_capacity_from_molecular_structure.ipynb`, which measures what this one asserts |
| `Heat_capacity_integration_N2.ipynb` | write your own trapezoidal rule, watch it converge against the exact polynomial integral, and measure that its error falls as $h^2$. Ends with $\Delta\underline{H}$ and $\Delta\underline{U}$ |

**`ch3/` — energy, and the ideal-gas heat capacity**
| Notebook | What it teaches |
|---|---|
| `Heat_capacity_range_of_validity.ipynb` | why a correlation's temperature range is part of its data: Appendix A.II's two O₂ rows against NIST-JANAF over 100–1800 K, and what each costs once you integrate it into $\underline{H}$ and $\underline{S}$. **Also derives** the cryogenic-range rows in `thermo.APPENDIX_A2_CP_CRYO` — run it to check those coefficients |
| `LJ_interaction_energy_figure.ipynb` | Figure 3.3-5: the Lennard-Jones interaction energy between two molecules, for argon and methane |
| `Steam_charts_from_appendix_A3.ipynb` | Figure 3.3-1: the Mollier and $T$–$S$ charts for steam, drawn from Appendix A.III (`data/steam_*.csv`) rather than from an equation of state, so chart and table agree by construction. Every line family is an argument — zoom in, redraw at your own scale, read values numerically. Ends by reading Illustration 3.4-1 off both charts |
| `Heat_capacity_from_molecular_structure.ipynb` | where $C_P^*$ comes from: equipartition, and an Einstein term per vibrational mode. Methane from four spectroscopic frequencies, matching NIST-JANAF to 0.03 J/(mol·K), with the $4R$ floor no nonlinear gas can go under. Then the derived curve becomes noisy "experimental" data and is fitted 400 times — what the coefficients do, what the curve does, and why the fitted *range* is part of the correlation. Needs no data file |
| `PH_charts_methane_and_nitrogen.ipynb` | Figures 3.3-2 and 3.3-3: pressure–enthalpy charts for methane and nitrogen, drawn entirely from Peng–Robinson — three constants, four heat-capacity coefficients, one cubic. Derives methane's cryogenic $C_P^*$ from its vibrational spectrum (the printed Appendix A.II row is 25% low at 100 K), checks the equation against the book's own Wagner constants before drawing, and ends by working Illustration 6.5-1 on the nitrogen chart |
| `refrigerant_comparison.ipynb` | can R-1234yf replace HFC-134a? Three constants per fluid and Peng–Robinson: why matching *volumetric* refrigeration capacity, not heat of vaporization, is the test a drop-in has to pass |

**`ch5/` — liquefaction and power cycles**
| Notebook | What it teaches |
|---|---|
| `Linde_liquefaction_CH4_figure.ipynb` | Figure 5.1-3: the liquefaction path of Illustration 5.1-1 |
| `Rankine_cycle_steam_figure.ipynb` | The Rankine cycle of Illustration 5.2-1, drawn on the steam chart |

**`ch6/` — thermodynamic properties of real substances (Peng–Robinson EOS)**
| Notebook | What it teaches |
|---|---|
| `PR_eos_reference.ipynb` | the generalized PR EOS and its compressibility/molar-volume roots (from-scratch reference; refactored into the `thermo` package) |
| `PR_isotherms_N2_example.ipynb` | PR isotherms for nitrogen |
| `PR_isotherms_O2_example.ipynb` | **Fig. 6.4-3** (Illustration 6.4-1): oxygen $P$–$\underline V$ diagram with the vapor–liquid saturation envelope |
| `PR_enthalpy_O2_example.ipynb` | **Fig. 6.4-4** (Illustration 6.4-1): oxygen $P$–$\underline H$ diagram (ideal-gas $C_P$ integral + PR enthalpy departure) |
| `PR_entropy_O2_example.ipynb` | **Fig. 6.4-5** (Illustration 6.4-1): oxygen $T$–$\underline S$ diagram (ideal-gas $C_P/T$ integral + PR entropy departure) |
| `PR_throttle_CH4_example.ipynb` | enthalpy/entropy departure functions applied to a throttling (Joule–Thomson) process (methane) |
| `PR_throttle_C2H6_homework.ipynb` | homework: outlet temperature of throttled ethane |
| `PR_heat_capacity_C2H6_homework.ipynb` | homework: pressure dependence of the heat capacity of ethane |

Each example/homework above has a `_thermo` twin (`PR_*_thermo.ipynb`) that does the same
calculation using the `thermo` package (`PengRobinson.from_database`) rather than an inline
EOS — the self-contained version shows the method, the twin shows the reuse.

**`ch7/` — equilibrium and stability in one-component systems**
| Notebook | What it teaches |
|---|---|
| `van_der_waals_EOS.ipynb` | van der Waals EOS solve |
| `Peng_Robinson_EOS_isotherms_CO2.ipynb` | PR isotherms for carbon dioxide |
| `CO2_phases.ipynb` | pressure–temperature phase diagram of CO₂ |
| `N2_phases.ipynb` | pressure–temperature phase diagram of nitrogen |
| `VLE from fugacity.ipynb` | pure-component vapor–liquid equilibrium from equal fugacities (PR EOS) |

`ch7/README.md` documents exporting a notebook to PDF via `nbconvert`.

**`ch8/` — the thermodynamics of multicomponent mixtures**
| Notebook | What it teaches |
|---|---|
| `RK_partial_molar_volume_water_methanol_example.ipynb` | Table 8.6-2: partial molar volumes from density data |
| `RK_partial_molar_enthalpy_water_methanol_example.ipynb` | Table 8.6-4: partial molar enthalpies from heat-of-mixing data |
| `Hrxn_temperature_N2O4_example.ipynb` | Illustration 8.5-2: the standard heat of reaction versus temperature |

`ch8/README.md` documents how the two partial-molar tables were rebuilt from the
printed density and heat-of-mixing data.

**`ch9/` — the Gibbs energy and fugacity of a component in a mixture**
| Notebook | What it teaches |
|---|---|
| `PR_species_fugacity_ethane_butane_example.ipynb` | Species fugacity in a mixture from the Peng-Robinson EOS (ethane / n-butane) |
| `Lewis_Randall_vs_PR_CO2_methane_example.ipynb` | Lewis-Randall rule against the Peng-Robinson EOS (CO2 / methane) |
| `activity_coefficient_correlation_benzene_TMP_example.ipynb` | Correlating activity coefficients (benzene / 2,2,4-trimethyl pentane) |
| `regular_solution_benzene_TMP_example.ipynb` | Regular solution theory (benzene / 2,2,4-trimethyl pentane) |
| `UNIFAC_benzene_TMP_example.ipynb` | Modified (Dortmund) UNIFAC (benzene / 2,2,4-trimethyl pentane) |
| `Margules_activity_coefficients_figure.ipynb` | The one-constant Margules activity coefficients |
| `Gex_composition_curves_figure.ipynb` | Excess Gibbs energy against composition |
| `Debye_Huckel_HCl_NaCl_example.ipynb` | Electrolyte solutions — Debye-Hückel and its extensions (HCl, NaCl) |
| `Wong_Sandler_mixing_rule_acetone_water_example.ipynb` | The Wong-Sandler mixing rule — an activity coefficient model inside an equation of state (acetone / water) |
| `hemoglobin_activity_coefficient_example.ipynb` | Activity coefficient of a protein on a Henry's law basis (hemoglobin in water) |

`ch9/README.md` covers the activity-coefficient models and which of them the `thermo`
package implements.

**`ch10/` — vapor–liquid equilibrium in mixtures**
| Notebook | What it teaches |
|---|---|
| `raoult_diagrams_pentane_heptane_example.ipynb` | Vapor-liquid equilibrium diagrams for an ideal mixture (n-pentane / n-heptane) |
| `hexane_triethylamine_diagrams_figure.ipynb` | Phase diagrams for a nearly ideal mixture (hexane / triethylamine) |
| `vle_correlation_benzene_TMP_example.ipynb` | Correlating vapor-liquid equilibrium data (benzene / 2,2,4-trimethyl pentane) |
| `excess_properties_benzene_TMP_figure.ipynb` | Excess enthalpy and entropy from the temperature dependence of $G^{ex}$ (benzene / 2,2,4-trimethyl pentane) |
| `azeotrope_to_vanlaar_benzene_cyclohexane_example.ipynb` | A phase diagram from one azeotropic point (benzene / cyclohexane) |
| `pressure_swing_methyl_acetate_methanol_example.ipynb` | Nonideal phase diagrams, an azeotrope, and how to separate it (methyl acetate / methanol) |
| `pr_vle_co2_isopentane_figure.ipynb` | One binary parameter, and what it is worth (carbon dioxide / isopentane) |
| `acetone_water_mixing_rules_figure.ipynb` | Four ways to model a very nonideal mixture (acetone / water) |
| `azeotrope_test_ethyl_acetate_benzene_homework.ipynb` | Does this system have an azeotrope? Two models that disagree (ethyl acetate / benzene) |
| `unifac_pentane_propionaldehyde_example.ipynb` | Three ways to get a phase diagram, and how little data each needs (n-pentane / propionaldehyde) |
| `flory_huggins_benzene_PIB_example.ipynb` | A solvent above a polymer, where mole fraction stops being useful (benzene / polyisobutylene) |

`ch10/README.md` lists the four mixing rules the chapter compares.

**`ch11/` — other phase equilibria in fluid mixtures**
| Notebook | What it teaches |
|---|---|
| `gas_solubility_co2_toluene_example.ipynb` | Solubility of a gas in a liquid from an equation of state (CO2 in toluene) |
| `lle_isobutane_furfural_example.ipynb` | Liquid-liquid equilibrium from an activity coefficient model (isobutane / furfural) |
| `polymer_compatibility_PS_PMMA_example.ipynb` | Are two molten polymers compatible? (polystyrene / PMMA) |
| `vlle_eos_co2_decane_example.ipynb` | Three-phase equilibrium from an equation of state (CO2 / n-decane) |
| `lle_eos_mixing_rules_co2_decane_example.ipynb` | Two mixing rules against one measured phase diagram (CO2 / n-decane) |
| `fluorocarbon_solubility_figure.ipynb` | Consolute temperatures of fluorocarbon + hydrocarbon mixtures (Figure 11.2-3) |
| `margules_azeotrope_vs_phase_split_homework.ipynb` | When does an azeotrope become a phase split? (Problem 11.2-17) |
| `vlle_isobutane_furfural_figure.ipynb` | When should you suspect a second liquid phase? (Figures 11.3-1 and 11.3-2) |
| `ternary_extraction_MIK_acetone_water_example.ipynb` | Liquid-liquid extraction on a triangular diagram (MIK / acetone / water) |
| `staged_extraction_MIK_acetone_water_example.ipynb` | Two-stage liquid-liquid extraction (MIK / acetone / water) |
| `ternary_mass_balance_example.ipynb` | Reading a triangular diagram (a mass balance on three components) |
| `distribution_coefficient_bromine_CCl4_water_example.ipynb` | An activity coefficient measured with a separatory funnel (bromine / CCl4 / water) |
| `octanol_water_partition_coefficient_example.ipynb` | The octanol-water partition coefficient: estimating it and using it |
| `partition_coefficient_correlation_pollutants_example.ipynb` | Testing a correlation against 160 printed numbers (Table 11.4-1) |
| `gibbs_energy_of_transfer_amino_acids_example.ipynb` | The hydrophobic effect, from six solubilities (amino acids, water to ethanol) |
| `osmometry_molecular_weight_PVC_albumin_example.ipynb` | A molecular weight from a column of liquid (PVC, and serum albumin) |
| `osmotic_virial_chymotrypsin_example.ipynb` | A molecular weight and a second virial coefficient from one osmometer (alpha-chymotrypsin) |
| `osmotic_pressure_blood_example.ipynb` | Why saline is 0.9 % (the osmotic pressure of blood) |
| `reverse_osmosis_seawater_example.ipynb` | Pushing water the other way (reverse osmosis of seawater) |

**`ch12/` — phase equilibria involving solids**
| Notebook | What it teaches |
|---|---|
| `solid_solubility_naphthalene_hexane_example.ipynb` | How much solid dissolves, and what the activity coefficient is worth (naphthalene / n-hexane) |
| `supercritical_solubility_naphthalene_co2_figure.ipynb` | Naphthalene in supercritical carbon dioxide, three ways (Figure 12.1-1) |
| `heat_of_fusion_from_solubility_insulin_example.ipynb` | A heat of fusion you cannot measure in a calorimeter (insulin hexamer) |
| `freezing_point_and_eutectic_figure.ipynb` | Freezing points, and the eutectic as a crossing (Figure 12.3-1, ethyl benzene / toluene) |
| `environmental_partitioning_benzo_a_pyrene_example.ipynb` | One chemical, one measurement, five compartments (benzo[a]pyrene) |
| `pcb_bioconcentration_example.ipynb` | Bioconcentration, with a field measurement to check it against (PCBs in St. Lawrence eels) |

**`ch13/` — chemical equilibrium**
| Notebook | What it teaches |
|---|---|
| `nitrogen_tetroxide_dissociation_example.ipynb` | Nitrogen tetroxide, three questions about one reaction (Figures 13.1-3 and 13.1-4) |
| `high_pressure_ammonia_equilibrium_example.ipynb` | Ammonia synthesis at high pressure (Figures 13.1-5 and 13.1-6) |
| `acetic_acid_ionization_constant_example.ipynb` | Reaction in the liquid phase: an esterification and an ionization (Figure 13.1-7) |
| `ideal_gas_pressure_and_equilibrium_ratios_example.ipynb` | Pressure, dilution, and the three equilibrium ratios (Figure 13.1-2) |
| `measured_equilibrium_constants_example.ipynb` | When the equilibrium constant is measured, not computed (Figure 13.1-1) |
| `heterogeneous_gas_solid_equilibrium_example.ipynb` | When a solid is one of the reactants: four heterogeneous equilibria (Figure 13.2-1) |
| `solubility_product_silver_chloride_example.ipynb` | The solubility product of silver chloride, and what an inert salt does to it (Figure 13.2-2) |
| `steam_carbon_multireaction_example.ipynb` | Coal and steam: three reactions at once (Figures 13.3-1 and 13.3-2) |
| `coupled_reactions_common_ion_and_atp_example.ipynb` | Two reactions at once, twice over: a common ion and a coupled hydrolysis |
| `combined_chemical_phase_ammonia_example.ipynb` | Reaction and condensation at once: ammonia over water |
| `titration_strong_and_weak_acids_example.ipynb` | Two titration curves, and why only one of them is interesting |
| `deprotonation_benzoyl_tyrosine_figure.ipynb` | When does an acid group carry its charge? (Figure 13.5-1) |
| `dibasic_acid_phthalic_figure.ipynb` | A molecule with two places to lose a proton (Figures 13.6-1 and 13.6-2) |
| `amino_acid_charge_glycine_example.ipynb` | The charge on an amino acid, and how it moves with temperature (Figures 13.6-5, 13.6-6 and 13.6-11) |
| `amino_acid_titration_glycine_example.ipynb` | An amino acid as its own buffer (Figures 13.6-8 and 13.6-9) |
| `protein_charge_lysozyme_figure.ipynb` | Thirty-two ionizable groups on one molecule (Figure 13.6-7) |

**`ch14/` — reactors, availability, and electrochemistry**
| Notebook | What it teaches |
|---|---|
| `adiabatic_reaction_temperature_ethylbenzene_example.ipynb` | The adiabatic reaction temperature, and why the phase of the feed decides it |
| `tank_batch_and_tubular_reactor_ethyl_acetate_example.ipynb` | One reaction, three reactors, and the same heat load |
| `maximum_useful_work_availability_example.ipynb` | The most work a fuel can do: methane, gasoline, and a spoonful of sugar |
| `electrochemical_cell_potentials_example.ipynb` | Cell potentials, equilibrium constants, and pH |

**`ch15/` — biochemical applications**
| Notebook | What it teaches |
|---|---|
| `protein_denaturation_examples.ipynb` | Proteins unfold when they are cold, too: six figures from one stability curve |
| `protein_solubility_salt_and_temperature_examples.ipynb` | Salting in, salting out, and a precipitation temperature |
| `solubility_vs_pH_examples.ipynb` | Four ionizable solids, one equation, and where each curve breaks |
| `ligand_binding_cooperative_vs_single_site_example.ipynb` | One binding site or four, and what cooperativity buys an animal |
| `gibbs_donnan_and_ultracentrifuge_examples.ipynb` | A protein that cannot cross a membrane rearranges everything that can |
| `bioreactor_yield_factors_examples.ipynb` | Balancing a reactor whose products have no molecular formula |
| `bioreactor_energy_and_second_law_examples.ipynb` | The second law decides how much a fermenter can make |

### `data/` — reference data

Pure-component and reaction property tables exported from the legacy companion database
(618 compounds and 99 reaction species; critical constants, acentric factor, ideal-gas
$C_P$ coefficients = Appendix A.II, formation properties, vapor-pressure constants). See
[`data/README.md`](data/README.md) for the schema and units.

### `thermo/` — the reusable library (object-oriented)

Twenty-four modules reading the `data/` tables: the cubic equations of state and their
departure functions, the activity-coefficient models and their fitting, UNIFAC, mixture
mixing rules, the vapor–liquid, liquid–liquid and solid–liquid equilibrium solvers, chemical
reaction equilibrium, electrolytes, osmotic equilibrium, and the chart-drawing layer. It is
the reusable core of the Python **general-purpose computing** layer (see Appendix B) — open
and inspectable, in contrast to a licensed process simulator. Both UNIFAC parameter sets are
complete and selectable: the original of Fredenslund, Jones and Prausnitz and the later
temperature-dependent (Dortmund) revision. API in [`thermo/README.md`](thermo/README.md).

## Generated output

- `ch*/pdf/` — figure exports, written by the notebooks when they run.

## Provenance and scope

These notebooks and the `thermo` package modernize the pedagogical numerical routines of
the 5e Visual Basic, MATLAB, and Mathcad programs (pure-fluid EOS, departure functions,
vapor pressure, phase diagrams, UNIFAC). In the 6e the Python stack is the **general-purpose
computing** layer — for both learning and everyday engineering analysis, and it now covers
the mixture VLE/flash, activity-coefficient, and chemical-equilibrium calculations across
ch8–15; **Aspen Plus** is
the complementary tool for specialized process modeling. See Appendix B of the text for the
design rationale and the property-data provenance.
