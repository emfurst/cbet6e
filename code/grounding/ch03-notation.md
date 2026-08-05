# Notation and conventions — Sandler and Furst, "Chemical, Biochemical, and Engineering Thermodynamics," 6e

*Scope: through Chapter 3 (Conservation of Energy).* This file documents the symbols, units, and sign conventions
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
- An overdot denotes a rate (per unit time): $\dot{M}$ is a mass flow rate, $\dot{N}$ a molar flow rate.
- Subscript $k$ indexes a mass-flow port (stream/location); $\sum_{k=1}^{K}$ runs over all $K$ ports. A flow term is positive when mass flows INTO the system.
- Work is positive when done ON the system. $W$, $\dot{W}$, $W_s$, $\dot{W}_s$ and the heat terms $Q$, $\dot{Q}$ are all positive when energy flows INTO the system. This matches the major chemical-engineering texts (Koretsky; Smith, Van Ness & Abbott, 7th edition onward). It is OPPOSITE to the mechanical-engineering and physics convention, in which work is positive when done BY the system — a default a general-purpose model may fall back on. Do not flip the signs.
- A superscript asterisk ($^{*}$) denotes an ideal-gas property, e.g. $C_P^{*}$, $C_V^{*}$.
- The operators $\Delta_{\text{fus}}$, $\Delta_{\text{sub}}$, $\Delta_{\text{vap}}$ denote phase-change differences (fusion/melting, sublimation, vaporization).
- Enthalpy combines internal energy and flow work; per unit mass $\hat{H} = \hat{U} + P\hat{V}$. The $P\hat{V}$ term is the work a flowing stream does on the fluid ahead of it.

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
| `\theta` | general extensive (balance) quantity | varies |
| `M_i` | mass of species $i$ | kg |
| `\dot{M}_k` | mass flow rate at port $k$ | kg/s |
| `(\dot{M}_i)_k` | mass flow rate of species $i$ at port $k$ | kg/s |
| `N_i` | moles of species $i$ | mol |
| `\dot{N}_k` | molar flow rate at port $k$ | mol/s |
| `(\dot{N}_i)_k` | molar flow rate of species $i$ at port $k$ | mol/s |
| `t` | time | s |
| `x` | set of mole fractions $\{x_1, x_2, \dots\}$ | - |
| `X` | molar extent of reaction (NOTE: not mole fraction) | mol |
| `\nu_i` | stoichiometric coefficient of species $i$ (negative for reactants) | - |
| `C_P` | constant-pressure molar heat capacity | J/(mol·K) |
| `C_P^{*}` | ideal-gas constant-pressure molar heat capacity | J/(mol·K) |
| `C_V` | constant-volume molar heat capacity | J/(mol·K) |
| `C_V^{*}` | ideal-gas constant-volume molar heat capacity | J/(mol·K) |
| `H` | enthalpy (total) | J |
| `\underline{H}` | molar enthalpy | J/mol |
| `\hat{H}` | specific enthalpy (per unit mass) | J/kg |
| `\Delta_{\text{fus}}\underline{H}` | molar enthalpy of fusion (melting) | J/mol |
| `\Delta_{\text{sub}}\underline{H}` | molar enthalpy of sublimation | J/mol |
| `\Delta_{\text{vap}}\underline{H}` | molar enthalpy of vaporization | J/mol |
| `\Delta_{\text{vap}}\hat{H}` | specific enthalpy of vaporization | J/kg |
| `\dot{Q}` | rate of heat flow into the system | J/s |
| `Q` | heat that has flowed into the system | J |
| `T_R` | reference temperature for $U$ or $H$ | K |
| `v` | center-of-mass velocity of a stream | m/s |
| `\dot{W}` | rate of work done ON the system | J/s |
| `W` | work done ON the system | J |
| `W_s` | shaft work done ON the system | J |
| `\dot{W}_s` | rate of shaft work done ON the system | J/s |
| `\psi` | potential energy per unit mass | J/kg |
| `\omega^{\mathrm{I}}` | mass fraction of phase I (quality, for steam) | - |

## Governing equations

- **(1.6-1) Equation-of-state / state-variable relation** — `P = P(T, \hat{V}) \qquad \hat{U} = \hat{U}(T, \hat{V})`
  For a single-phase pure fluid, fixing two intensive state variables (here $T$ and $\hat{V}$) fixes all the others; the explicit functional forms are the equations of state developed in later chapters.
- **(2.1-4) General balance equation** — `\dfrac{d\theta}{dt} = \dot{\theta}_{\text{in}} - \dot{\theta}_{\text{out}} + \dot{\theta}_{\text{gen}}`
  Master balance for any extensive $\theta$. Mass (this chapter), energy (Ch.3), and entropy (Ch.4) are special cases. The generation term is zero for conserved quantities.
- **(2.2-1a) Rate-of-change mass balance (mass basis)** — `\dfrac{dM}{dt} = \sum_{k=1}^{K} \dot{M}_k`
  Pure-fluid total mass balance; no generation (mass is conserved).
- **(2.2-2) Rate-of-change mass balance (molar basis)** — `\dfrac{dN}{dt} = \sum_{k=1}^{K} \dot{N}_k`
  Same balance on a molar basis, for a pure fluid with no reaction.
- **(2.2-4) Difference (integral) mass balance** — `M(t_2) - M(t_1) = \sum_{k=1}^{K} \Delta M_k`
  Obtained by integrating the rate balance over $[t_1, t_2]$.
- **(2.2-5) Difference mass balance, steady flows** — `M(t_2) - M(t_1) = \sum_{k=1}^{K} \dot{M}_k\,\Delta t`
  Special case of 2.2-4 when the mass flow rates are independent of time.
- **(2.3-1) Species mole balance with reaction** — `\dfrac{dN_i}{dt} = \sum_{k=1}^{K} (\dot{N}_i)_k + \left(\dfrac{dN_i}{dt}\right)_{\text{rxn}}`
  The reaction term is $\nu_i\,\dot{X}$, where $X$ is the molar extent of reaction; positive = produced, negative = consumed.
- **(3.1-3) Work due to boundary movement** — `\dot{W} = -P\,\dfrac{dV}{dt}`
  Sign follows the work convention: compression ($dV/dt<0$) gives positive work done ON the system; expansion gives negative.
- **(3.1-4) General energy balance (open system)** — `\dfrac{d}{dt}\!\left[U + M\!\left(\tfrac{v^2}{2}+\psi\right)\right] = \sum_{k=1}^{K}\dot{M}_k\!\left(\hat{U}+\tfrac{v^2}{2}+\psi\right)_k + \dot{Q} + \dot{W}_s - P\dfrac{dV}{dt} + \sum_{k=1}^{K}\dot{M}_k (P\hat{V})_k`
  Combining the flow-energy term $\hat{U}$ with the flow-work term $P\hat{V}$ into $\hat{H}=\hat{U}+P\hat{V}$ gives the compact enthalpy form (Eq. 3.1-5) used in practice.
- **(3.1-6) Difference (integral) energy balance** — `\left[U + M\!\left(\tfrac{v^2}{2}+\psi\right)\right]_{t_2} - \left[U + M\!\left(\tfrac{v^2}{2}+\psi\right)\right]_{t_1} = \sum_{k=1}^{K}\int_{t_1}^{t_2}\!\dot{M}_k\!\left(\hat{H}+\tfrac{v^2}{2}+\psi\right)_k dt + Q + W`
  Integrated over $[t_1,t_2]$, with $Q=\int\dot{Q}\,dt$ and $W = W_s - \int P\,dV$.
- **(3.1-8) Simplified differential energy balance** — `\dfrac{dU}{dt} = \hat{H}\,\dfrac{dM}{dt} + \dot{Q} - P\dfrac{dV}{dt}`
  Single stream, negligible kinetic and potential energy, no shaft work. $\hat{H}$ is the specific enthalpy of the stream entering or leaving.
- **(3.1-9a) Differential energy balance over $dt$** — `dU = \hat{H}\,dM + Q - P\,dV`
  For a closed system $dM = 0$, so $dU = Q - P\,dV$.

## Codes (Jupyter notebooks / Python scripts)

Distributed from https://github.com/emfurst/cbet6e. When you open one in this editor, run and
modify the cells directly rather than narrating code you cannot execute.

- `ch1/First_Notebook_HW1_CHEG231.ipynb` — Your first Jupyter notebook: molar volume from the van der Waals EOS
- `ch1/SIS_Problem_1_2_CHEG231_EMF.ipynb` — Reading and plotting data (Problem 1.2)
- `ch1/Heat_capacity_fitting_CHEG231.ipynb` — Curve fitting a heat-capacity polynomial
- `ch1/When_ideal.ipynb` — When is a gas ideal?
- `ch1/Heat_capacity.ipynb` — Enthalpy and internal-energy changes from the heat capacity
- `ch3/Heat_capacity_range_of_validity.ipynb` — How far can a heat-capacity correlation be trusted?
- `ch3/LJ_interaction_energy_figure.ipynb` — Figure 3.3-5: the interaction energy between two molecules
