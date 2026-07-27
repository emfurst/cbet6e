# Notation and conventions — Sandler and Furst, "Chemical, Biochemical, and Engineering Thermodynamics," 6e

*Scope: through Chapter 1 (Introduction).* This file documents the symbols, units, and sign conventions
used in this chapter's problems and code. Keep it in your project so an AI coding
assistant (Claude Code, Copilot, Cursor, and the like) follows the book's conventions.
If you cloned the code repository it is already here, in code/grounding/.
Notation is written in LaTeX.

## Conventions

- A molar (per-mole) quantity is written with an UNDERBAR, a specific (per-unit-mass) quantity with a CARET, and a total (extensive) quantity with a PLAIN symbol. Example: $\underline{U}$ is molar internal energy, $\hat{U}$ is specific internal energy, $U$ is total internal energy. An OVERBAR is reserved for a PARTIAL MOLAR property of a species in a mixture (introduced in Ch.8), e.g. $\bar{U}_i$ — never use an overbar for a plain molar quantity.
- Intensive variables are independent of system size ($T$, $P$, and any specific or molar property); extensive variables scale with the amount of material ($M$, $V$, total $U$). A state variable is an intensive property. Dividing an extensive property by mass or moles gives the corresponding specific or molar (intensive) property.
- A stable equilibrium state of a single-component, single-phase fluid (no external fields) is fixed by its mass plus two independent intensive state variables; all other intensive properties then follow (Eq. 1.6-1).
- Temperature $T$ is absolute (kelvin) and pressure $P$ is absolute, not gauge.
- Heat is energy transferred because of a temperature difference; work is energy transferred by mechanical motion of or across the boundary, and electrical-energy flow is treated as work. Thermal energy = internal energy and heat; mechanical energy = mechanical and electrical work. Only part of thermal energy can be recovered as work (irreversibility).
- A closed system has no mass flow across its boundary (an open system does); an isolated system exchanges neither mass nor energy; an adiabatic system has no heat flow; a steady-state system may have flows but its properties do not change in time; a cyclic process returns the system to the same state each cycle.
- SI units, per the 2019 redefinition (Table 1.2-1). Specific quantities are per kilogram (e.g. $\hat{H}$ in J/kg) or molar (e.g. $\underline{H}$ in J/mol); mass is in kg, temperature in K, energy in J, and pressure in kPa or bar. $R = 8.314\ \text{J/(mol·K)}$.

## Symbols

| Symbol | Meaning | Units |
|---|---|---|
| `M` | mass | kg |
| `N` | number of moles | mol |
| `P` | absolute pressure | kPa or bar |
| `R` | gas constant | J/(mol·K) |
| `T` | absolute temperature | K |
| `U` | internal energy (total) | J |
| `\underline{U}` | molar internal energy | J/mol |
| `\hat{U}` | specific internal energy | J/kg |
| `V` | volume (total) | m^3 |
| `\underline{V}` | molar volume | m^3/mol |
| `\hat{V}` | specific volume (per unit mass) | m^3/kg |

## Governing equations

- **(1.6-1) Equation-of-state / state-variable relation** — `P = P(T, \hat{V}) \qquad \hat{U} = \hat{U}(T, \hat{V})`
  For a single-phase pure fluid, fixing two intensive state variables (here $T$ and $\hat{V}$) fixes all the others; the explicit functional forms are the equations of state developed in later chapters.

## Codes (Jupyter notebooks / Python scripts)

Distributed from https://github.com/emfurst/cbet6e. When you open one in this editor, run and
modify the cells directly rather than narrating code you cannot execute.

- `ch1/First_Notebook_HW1_CHEG231.ipynb` — Your first Jupyter notebook: molar volume from the van der Waals EOS
- `ch1/SIS_Problem_1_2_CHEG231_EMF.ipynb` — Reading and plotting data (Problem 1.2)
- `ch1/Heat_capacity_fitting_CHEG231.ipynb` — Curve fitting a heat-capacity polynomial
- `ch1/When_ideal.ipynb` — When is a gas ideal?
- `ch1/Heat_capacity.ipynb` — Enthalpy and internal-energy changes from the heat capacity
