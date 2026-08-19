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
counterparts build the same method procedurally — see [`thermo/ROADMAP.md`](thermo/ROADMAP.md)
for the package's design and build-out plan.

## Notation

Molar quantities use an **underbar** ($\underline{V}$, $\underline{H}$), specific
(per-mass) quantities a **caret** ($\hat{V}$), total quantities a **plain** symbol, and
an **overbar** is reserved for a partial molar property ($\bar{G}_i$, introduced in
Ch. 8). This matches the book. Per-chapter notation files for an AI coding assistant are
generated alongside the companion site (`6e_companion_site/build/grounding/chNN-notation.md`).

## Prose style in the code

Everything a reader sees — markdown cells, docstrings, comments, and the text a cell
prints — is plain prose with **no emoji**. Markers like the ones used for triage in the
revision notes must not cross into a notebook or a module: a notebook is read on the
companion site and on GitHub, and its stored output is part of what a student reads. Where
a warning needs weight, say so in words (`NOT`, `must`, bold) rather than with a glyph.
This was reset across `code/` on 2026-08-19, when 596 markers were removed from 75 files;
`archive/` was left frozen on purpose.

## Layout

### Chapter notebooks

**`ch1/` — getting started with Python (a skills ladder)**

Work through these **in order**. The first two are read rather than run; each of the last
four is the anchor for one printed problem (1.9–1.12) and ends with a *Your turn* cell.
Rewritten 2026-08-05 — the previous course-era versions are in `ch1/archive/`.

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
| `Steam_charts_from_appendix_A3.ipynb` | Figure 3.3-1: the Mollier and $T$–$S$ charts for steam, drawn from Appendix A.III (`../data/steam_*.csv`) rather than from an equation of state, so chart and table agree by construction. Every line family is an argument — zoom in, redraw at your own scale, read values numerically. Ends by reading Illustration 3.4-1 off both charts |
| `Heat_capacity_from_molecular_structure.ipynb` | where $C_P^*$ comes from: equipartition, and an Einstein term per vibrational mode. Methane from four spectroscopic frequencies, matching NIST-JANAF to 0.03 J/(mol·K), with the $4R$ floor no nonlinear gas can go under. Then the derived curve becomes noisy "experimental" data and is fitted 400 times — what the coefficients do, what the curve does, and why the fitted *range* is part of the correlation. Needs no data file |
| `PH_charts_methane_and_nitrogen.ipynb` | Figures 3.3-2 and 3.3-3: pressure–enthalpy charts for methane and nitrogen, drawn entirely from Peng–Robinson — three constants, four heat-capacity coefficients, one cubic. Derives methane's cryogenic $C_P^*$ from its vibrational spectrum (the printed Appendix A.II row is 25% low at 100 K), checks the equation against the book's own Wagner constants before drawing, and ends by working Illustration 6.5-1 on the nitrogen chart |
| `refrigerant_comparison.ipynb` | can R-1234yf replace HFC-134a? Three constants per fluid and Peng–Robinson: why matching *volumetric* refrigeration capacity, not heat of vaporization, is the test a drop-in has to pass |

**`ch4/` — entropy**
| Notebook | What it teaches |
|---|---|
| `water_cp.ipynb` | plot the heat capacity of liquid water from triple point to critical point (NIST/IAPWS data) |

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

### `data/` — reference data

Pure-component and reaction property tables exported from the legacy companion database
(618 compounds and 99 reaction species; critical constants, acentric factor, ideal-gas
$C_P$ coefficients = Appendix A.II, formation properties, vapor-pressure constants). See
[`data/README.md`](data/README.md) for the schema and units.

### `thermo/` — the reusable library (object-oriented)

A small, class-based package (`PengRobinson`, `UNIFAC`) providing the pure-fluid
Peng–Robinson EOS and UNIFAC activity coefficients, reading the `data/` tables. It is the
reusable core of the Python **general-purpose computing** layer (see Appendix B) — open and
inspectable, in contrast to a licensed process simulator. PR and modified (Dortmund) UNIFAC
are complete and verified; original UNIFAC awaits its R/Q constants. The package is being
built out to cover the mixture, phase-equilibrium, and chemical-equilibrium calculations of
ch8–15 — see [`thermo/ROADMAP.md`](thermo/ROADMAP.md). API in [`thermo/README.md`](thermo/README.md).

### `honors/` — advanced / molecular thermodynamics

Optional material for advanced students and the emerging-topics track:
- `honors/MD/` — a simple molecular-dynamics simulation in Python (see also
  `https://github.com/emfurst/SimpleMD`).
- `honors/honors maxwell boltzmann/` — the Maxwell–Boltzmann speed distribution.

## Not part of the deliverable

Working and superseded material, excluded from the companion site:
- `ch*/archive/` — superseded notebooks; `ch*/pdf/` — generated PDF exports (author's reference, not distributed).
- `test code/` — experimental scratch (e.g. ternary diagrams).
- `CHEG231_python_2025.ipynb` — a full-course working notebook.
- `PLANNING_2024.txt` — the original 2024 planning outline for this code package.

## Provenance and scope

These notebooks and the `thermo` package modernize the pedagogical numerical routines of
the 5e Visual Basic, MATLAB, and Mathcad programs (pure-fluid EOS, departure functions,
vapor pressure, phase diagrams, UNIFAC). In the 6e the Python stack is the **general-purpose
computing** layer — for both learning and everyday engineering analysis — and is being
extended to the mixture VLE/flash, activity-coefficient, and chemical-equilibrium
calculations across ch8–15 (see [`thermo/ROADMAP.md`](thermo/ROADMAP.md)); **Aspen Plus** is
the complementary tool for specialized process modeling. See Appendix B of the text for the
design rationale and the property-data provenance.
