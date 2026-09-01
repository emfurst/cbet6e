# Chapter 6 — thermodynamic properties of real substances (Peng–Robinson EOS)

## Naming convention

Notebooks are named `PR_<topic>_<substance>_<role>.ipynb`:

- **topic** — `eos`, `isotherms`, `throttle`, `heat_capacity`, `discharge`, …
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
| `PR_properties_table_O2_example.ipynb` | example | **package** — generates **Table 6.4-4** (Illustration 6.4-1): $Z$, $\underline V$, $\underline H$, $\underline S$ for oxygen at 13 pressures × 11 temperatures, written to `output/Table_6.4-4.{txt,csv}`. The tabulated form of the same calculation the three chart notebooks plot. **Unpaired** — see below |
| `Helmholtz_fundamental_eos_O2_example.ipynb` | example | **package** — the chapter's only *fundamental* equation of state. Assembles $\underline A(T,\underline V)$ for oxygen from the book's own data (ideal-gas part from Appendix A.II $C_P^*$, residual part from one integration of the PR equation) and takes $P$, $Z$, $\underline S$, $\underline U$, $\underline H$, $\underline G$, $C_V$, $C_P$, $\phi$, $\alpha$, $\kappa_T$, $\mu$ from it **by finite differences**, quoting no property formula. Eqs. 6.4-2, 6.4-29, 6.4-30 and 6.2-35 are the *check*, as is all of Table 6.4-4 (worst 1.3e-3 J/mol vs. a printed 0.01). Backs **Illustration 6.2-4**. **Unpaired** — see below |
| `PR_discharge_N2_example.ipynb` | example | **package** — works **Illustration 6.7-1**: nitrogen withdrawn from an insulated 0.15 m³ cylinder, solved from the isentropic condition $\underline S(t{=}50)=\underline S(0)$ at the molar volume the mass balance fixes. Reproduces the printed 134.66 K / 40.56 bar, runs the text's hand iteration alongside the direct solve, recomputes four of the five rows of Table 6.5-1, and writes `output/Illustration_6.7-1.txt`. **Unpaired** — see below |
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
Four notebooks have no twin. The reference notebook *is* the from-scratch method.
`PR_properties_table_O2_example.ipynb` and `PR_discharge_N2_example.ipynb` are
**production** notebooks — they generate numbers the book prints (Table 6.4-4; and
Illustration 6.7-1's results, which also fill a row of Table 6.5-1 and of the Summary
table in Illustration 6.5-1) — so each needs one authoritative implementation rather
than two that differ in the last digit. `Helmholtz_fundamental_eos_O2_example.ipynb` has no
twin for a different reason: a self-contained/`_thermo` pair contrasts *inline math* against
*package reuse*, and this notebook's entire point is that the package's closed-form
`departure_H`/`departure_S` are the **check** on a calculation that derives them. Splitting
it would put the two halves of one argument in two files.

> **Numbers can differ slightly between a notebook and its twin.** The twins pull
> $T_c$, $P_c$, $\omega$, and $C_p$ coefficients from `code/data/pure_property.csv`
> via `from_database`, whereas the self-contained notebooks hard-code the SIS
> Table 6.6-1 values. In practice the difference is tiny (e.g. the methane throttle
> gives $P_2 = 41.46$ vs $41.45$ bar), but regenerate homework answer keys from the
> twin if that is the version you assign.
>
> **The oxygen pair is fixed (2026-08-04).** `pure_property.csv` carries
> $\omega = 0.025$ and $P_c = 50.4$ bar, against Table 6.6-1's $\omega = 0.021$
> ($\kappa = 0.4069$, the value Illustration 6.4-1 prints) and $P_c = 5.046$ MPa —
> a different parameter set, not a rounding, so `from_database("oxygen")` silently
> worked on a different basis from the illustration. `thermo.data` now carries
> **`TABLE_6_6_1`** beside `APPENDIX_A2_CP`, and the Figs. 6.4-3/-4/-5 twins build
>
> ```python
> pr = PengRobinson(**TABLE_6_6_1["oxygen"], cp=APPENDIX_A2_CP_CRYO["oxygen"])
> ```
>
> so each twin now matches its self-contained partner exactly. Any new notebook that
> reproduces printed numbers should do the same — as `PR_discharge_N2_example.ipynb`
> does with `TABLE_6_6_1["nitrogen"]`.

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

**The fallback needs matplotlib ≥ 3.11** (2026-08-09). Without TeX, math is rendered by
matplotlib's own `mathtext` engine, which implements a *subset* of LaTeX — and that
subset only gained `\underline` in 3.11.0. Colab still ships 3.10, where every
`$\underline{V}$` label raised

```
ParseFatalException: Unknown symbol: \underline
```

and killed the cell on the first draw, `tight_layout` included. Notebooks whose labels
use the molar underbar therefore open with a guard cell that upgrades matplotlib only
where it is too old, and does nothing on a machine that already has 3.11+:

```python
import sys
from importlib.metadata import version

if tuple(map(int, version("matplotlib").split(".")[:2])) < (3, 11):
    !pip install -q "matplotlib>=3.11"
    if "matplotlib" in sys.modules:
        print(...)   # the RESTART NOW banner
```

**Colab needs a session restart, and there is no way around it.** Colab imports
matplotlib as the session starts, so by the time the guard runs, `matplotlib` is
already in `sys.modules`. pip replaces the files on disk, but the kernel goes on
running the module it loaded at startup — so without a restart the labels fail exactly
as before, and the traceback is *confusing rather than obvious*: the line numbers come
from the old code objects in memory while the source text is read from the new file on
disk, so it points at lines whose content does not match the error. That mismatch is
the tell. The cell therefore prints a boxed **RESTART NOW: Runtime > Restart session,
then Run all** whenever it upgrades. One click, once per session.

The version is read through `importlib.metadata`, which inspects the installed metadata
*without* importing matplotlib — kept even though Colab has already imported it,
so the guard stays a true no-op everywhere else. And the whole thing self-retires: once
Colab's image moves to 3.11 the check stops firing, nothing installs, and no restart is
ever asked for again.

The alternative — rewriting all 37 label sites to emit `V̲` (combining low line U+0332)
under mathtext and `\underline{V}` under LaTeX — was considered and rejected 2026-08-10:
it needs no restart, but it puts a second notation spelling in the source permanently to
work around a gap that Colab will close on its own.

Also note `\underline` needs its braces. Real LaTeX accepts `\underline V`, but
mathtext declares a required group and rejects it — write `$\underline{V}$` always.

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

- `pdf/` — rendered PDF exports, written by the notebooks; regenerate as needed

Each topic appears on the companion site as a single card offering two labeled downloads:
the self-contained notebook and its `_thermo` twin.
