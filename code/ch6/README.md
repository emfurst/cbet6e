# Chapter 6 — thermodynamic properties of real substances (Peng–Robinson EOS)

## Naming convention

Notebooks are named `PR_<topic>_<substance>_<role>.ipynb`:

- **topic** — `eos`, `isotherms`, `throttle`, `heat_capacity`, …
- **substance** — `CH4`, `C2H6`, `N2`, … (omitted for the general reference notebook)
- **role** — `reference` (the general method), `example` (worked), or `homework`

A `_thermo` suffix marks the **package version** of a notebook (see below).

## The notebooks

| Notebook | Role | Uses |
|---|---|---|
| `PR_eos_reference.ipynb` | reference | inline — builds the generalized PR EOS and its compressibility roots from scratch. This is the method the `thermo` package was refactored from. |
| `vdW_isotherms_example.ipynb` | example | inline — recreates **Figure 6.6-1**: isotherms of the van der Waals EOS in the *reduced* $P$–$\underline V$ plane, above and below the critical isotherm. Split out of `ch7/van_der_waals_EOS.ipynb` (2026-07-28) so the figure's QR code lands on a chapter-6 file; that notebook keeps Figures 7.3-1 to 7.3-4. The reduced-EOS routines `rvdw_P`/`rvdw_V` are shared **by copy** between the two — change one, change the other. |
| `PR_isotherms_N2_example.ipynb` | example | inline — nitrogen isotherms in the $P$–$\underline V$ plane |
| `PR_isotherms_O2_example.ipynb` | example | inline — recreates **Figure 6.4-3** (Illustration 6.4-1): oxygen isotherms on log–log $P$–$\underline V$ axes with the vapor–liquid saturation envelope and sub-critical tie lines (equal-fugacity), the same construction as the vdW Figure 7.3-4 notebook |
| `PR_enthalpy_O2_example.ipynb` | example | inline — recreates **Figure 6.4-4** (Illustration 6.4-1): oxygen isotherms in the $P$–$\underline H$ plane; enthalpy = ideal-gas $C_P$ integral + PR enthalpy departure (Eq. 6.4-29), ref. ideal gas 25 °C/1 bar |
| `PR_entropy_O2_example.ipynb` | example | inline — recreates **Figure 6.4-5** (Illustration 6.4-1): oxygen isobars in the $T$–$\underline S$ plane; entropy = ideal-gas $C_P/T$ integral $-R\ln(P/P_0)$ + PR entropy departure (Eq. 6.4-30) |
| `PR_throttle_CH4_example.ipynb` | example | inline — enthalpy/entropy departures for an isenthalpic (Joule–Thomson) throttle of methane |
| `PR_throttle_C2H6_homework.ipynb` | homework | inline — outlet temperature of throttled ethane |
| `PR_heat_capacity_C2H6_homework.ipynb` | homework | inline — pressure dependence of the heat capacity of ethane |

## Self-contained vs. `_thermo` twins

Each **example** and **homework** above has a `_thermo` twin that does the *same*
calculation but imports the shared package instead of re-deriving the EOS:

```python
import sys; sys.path.append("..")
from thermo import PengRobinson
pr = PengRobinson.from_database("methane")   # Tc, Pc, omega, Cp from code/data/pure_property.csv
```

| Self-contained (inline math) | Package twin |
|---|---|
| `PR_isotherms_N2_example.ipynb` | `PR_isotherms_N2_example_thermo.ipynb` |
| `PR_isotherms_O2_example.ipynb` | `PR_isotherms_O2_example_thermo.ipynb` |
| `PR_enthalpy_O2_example.ipynb` | `PR_enthalpy_O2_example_thermo.ipynb` |
| `PR_entropy_O2_example.ipynb` | `PR_entropy_O2_example_thermo.ipynb` |
| `PR_throttle_CH4_example.ipynb` | `PR_throttle_CH4_example_thermo.ipynb` |
| `PR_throttle_C2H6_homework.ipynb` | `PR_throttle_C2H6_homework_thermo.ipynb` |
| `PR_heat_capacity_C2H6_homework.ipynb` | `PR_heat_capacity_C2H6_homework_thermo.ipynb` |

The **self-contained** version shows the method (the cubic, the departure-function
formulas coded by hand); the **twin** shows the reuse — the book's "two ways to work."
The reference notebook has no twin: it *is* the from-scratch method.

> **Numbers differ slightly between a notebook and its twin.** The twins pull
> $T_c$, $P_c$, $\omega$, and $C_p$ coefficients from `code/data/pure_property.csv`
> via `from_database`, whereas the self-contained notebooks hard-code the SIS
> Table 6.6-1 values. In practice the difference is tiny (e.g. the methane throttle
> gives $P_2 = 41.46$ vs $41.45$ bar), but regenerate homework answer keys from the
> twin if that is the version you assign.

## Figure typography

Notebooks that generate **book art** set their type in Computer Modern, to match the
book, rather than matplotlib's default DejaVu Sans. The block lives next to the
matplotlib import:

```python
from shutil import which
mpl.rcParams["text.usetex"] = which("latex") is not None
mpl.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "pdf.fonttype": 42, "ps.fonttype": 42})
```

`usetex` is switched on only when a LaTeX binary is present, so the notebooks still run
in Colab, which has no TeX; the fallback keeps a serif face with CM math.

**Two consequences for label text, once usetex is on.** Every string is handed to LaTeX,
so:

- write degree signs as `$^\circ$C`, never a literal `°` — a bare `°` aborts the LaTeX run;
- write molar quantities as `$\underline{V}$`, not with the combining low line `V̲`;
- `style="italic"` is ignored — set italic symbols in math mode instead;
- bare `%`, `_`, `#`, `&` in any label will break the build.

Applied to the Fig. 6.4-3/-4/-5 notebooks and their `_thermo` twins, plus
`vdW_isotherms_example.ipynb` (Fig. 6.6-1). Verify with:

```python
import pymupdf; sorted({x[3] for p in pymupdf.open("pdf/Fig_6.4-3.pdf") for x in p.get_fonts()})
```

which should return `CM*` names only.

## Other folders

- `archive/` — earlier versions (was `old versions/`)
- `pdf/` — rendered PDF exports (stale filenames; regenerate as needed)

These notebooks are wired into the companion site via
`6e_companion_site/source/chapters/ch06.yaml` (`codes:` block). Each topic is **one
entry** that pairs its self-contained notebook (`file:`) with its `_thermo` twin
(`file_thermo:`); the site collapses the pair into a single card offering two labeled
download links (self-contained + thermo package). Update those paths if you rename
anything here, then rerun `6e_companion_site/build.py`.
