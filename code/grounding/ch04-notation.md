# Notation and conventions — Sandler and Furst, "Chemical, Biochemical, and Engineering Thermodynamics," 6e

*Scope: through Chapter 4 (Entropy: An Additional Balance Equation).* This file documents the symbols, units, and sign conventions
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
- The second law is stated in this book as $S_{\text{gen}} \ge 0$: the rate of entropy generation $\dot{S}_{\text{gen}} \ge 0$ in every real process, equal to zero only for a reversible process and at equilibrium. This entropy-generation form is the axiom the book uses — NOT "the entropy of the universe increases" and not a word-statement of Clausius or Kelvin–Planck. The Clausius and Kelvin–Planck statements are shown to FOLLOW from $S_{\text{gen}} \ge 0$ (Illustrations 4.1-1, 4.1-2), so do not treat them as the starting axiom.
- In the entropy balance the entropy carried by heat is $\dot{Q}/T$, where $T$ is the absolute temperature of the system AT the boundary where the heat crosses (not the temperature of a distant reservoir). For several heat-flow ports, replace $\dot{Q}/T$ by $\sum_j \dot{Q}_j/T_j$. The sign of $\dot{Q}$ follows Chapter 3: positive when heat flows INTO the system.
- A process is reversible when $\dot{S}_{\text{gen}} = 0$ (infinitesimal internal gradients; no viscous dissipation and no finite-temperature-difference heat flow). Reversible operation gives the MAXIMUM work obtainable from — or the MINIMUM work required for — a given change of state; it is the bound against which real devices are measured. "Reversible" is not a synonym for "adiabatic" ($\dot{Q}=0$) or "isentropic": an adiabatic reversible process is isentropic, but an adiabatic irreversible one is not.
- Helmholtz energy $A \equiv U - TS$ (Eq. 4.2-6) and Gibbs energy $G \equiv H - TS = U + PV - TS$ (Eq. 4.2-8); both are state functions. Use these definitions — $A$ is $U-TS$, not $U+PV$.
- This book writes molar (per-mole) quantities with an UNDERBAR ($\underline{S}$, $\underline{G}$, $\underline{A}$). An OVERBAR is NOT a molar quantity here — it is reserved for the partial molar property of a species in a mixture (Ch.8). If another source writes a plain per-mole quantity with an overbar, map it to this book's underbar; it means the same thing (per mole), just a different accent.
- The symbol $T_R$ here denotes the absolute temperature of a radiation-EMITTING body (Eq. 4.1-10). In Chapter 3, $T_R$ meant a reference temperature for $U$ or $H$. Same symbol, different meaning — disambiguate by context. Radiative heat transfer is rarely used in this book.

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
| `S` | entropy (total) | J/K |
| `\underline{S}` | molar entropy | J/(mol·K) |
| `\hat{S}` | specific entropy (per unit mass) | J/(kg·K) |
| `\dot{S}_{\text{gen}}` | rate of entropy generation within the system ($\ge 0$) | J/(K·s) |
| `S_{\text{gen}}` | total entropy generated over the interval ($\ge 0$) | J/K |
| `A` | Helmholtz energy (total), $A \equiv U - TS$ | J |
| `\underline{A}` | molar Helmholtz energy | J/mol |
| `G` | Gibbs energy (total), $G \equiv H - TS = U + PV - TS$ | J |
| `\underline{G}` | molar Gibbs energy | J/mol |
| `\hat{G}` | specific Gibbs energy (per unit mass) | J/kg |
| `\mathcal{A},\ \mathcal{B}` | availability (maximum useful shaft work for a change of state) | J |
| `W_s^{\text{rev}}` | shaft work in a reversible process | J |
| `W^{\text{rev}}` | work in a reversible process | J |
| `\gamma` | heat-capacity ratio $C_P/C_V$ | - |
| `\dot{Q}_R` | radiant heat flux (Stefan–Boltzmann radiation) | J/s |
| `T_R` | temperature of the body emitting radiation (reuses the Ch.3 symbol $T_R$) | K |

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
- **(4.1-5a) Entropy balance, open system (the second law)** — `\dfrac{dS}{dt} = \sum_{k=1}^{K}\dot{M}_k\hat{S}_k + \dfrac{\dot{Q}}{T} + \dot{S}_{\text{gen}}`
  The third balance equation, specializing the general balance (Eq. 2.1-4) with a non-zero generation term. Terms: convected entropy of the mass flows, entropy flow $\dot{Q}/T$ accompanying heat, and internal generation $\dot{S}_{\text{gen}} \ge 0$. On a molar basis replace $\dot{M}_k\hat{S}_k$ by $\dot{N}_k\underline{S}_k$.
- **(4.1-5b) Entropy balance, closed system** — `\dfrac{dS}{dt} = \dfrac{\dot{Q}}{T} + \dot{S}_{\text{gen}}`
  Set all $\dot{M}_k = 0$ in Eq. 4.1-5a. For a closed system entropy changes only through heat flow and internal generation.
- **(4.1-5c) Second-law axiom and equilibrium condition** — `S_{\text{gen}} \ge 0, \qquad \dfrac{dS}{dt} = 0 \ \text{at equilibrium}`
  The defining statement of the second law used throughout the book. $\dot{S}_{\text{gen}} = 0$ only for a reversible process; entropy is a maximum at equilibrium in an isolated, constant-volume system.
- **(4.1-9) Difference (integral) entropy balance** — `S_2 - S_1 = \sum_{k}\int_{t_1}^{t_2}\dot{M}_k\hat{S}_k\,dt + \int_{t_1}^{t_2}\dfrac{\dot{Q}}{T}\,dt + S_{\text{gen}}`
  Integrated over $[t_1,t_2]$, with $S_{\text{gen}} = \int_{t_1}^{t_2}\dot{S}_{\text{gen}}\,dt \ge 0$. If each stream's entropy is constant the flow term becomes $\sum_k \underline{S}_k\,\Delta N_k$; if $T$ is constant at the heat-flow port the heat term becomes $Q/T$.
- **(4.2-6) Definition of the Helmholtz energy** — `A \equiv U - TS`
  A state function. The reversible work at constant $N$, $V$, $T$ equals $\Delta A$ (Eq. 4.2-7); a real (irreversible) process needs $\Delta A + TS_{\text{gen}}$.
- **(4.2-8) Definition of the Gibbs energy** — `G \equiv U + PV - TS = H - TS`
  A state function. The reversible shaft work at constant $N$, $P$, $T$ equals $\Delta G$. $TS_{\text{gen}}$ is the mechanical energy dissipated to thermal energy by irreversibilities.
- **(4.2-13b) Fundamental relation (closed system)** — `dU = T\,dS - P\,dV`
  The differential internal-energy change for a closed system, valid for any process (reversible or not) since it relates state functions only. The open-system form is $dU = T\,dS - P\,dV + \hat{G}\,dM$ (Eq. 4.2-13a). Basis for computing entropy changes of real fluids in Ch.6.
- **(4.3-4) Maximum work from a heat engine (reversible)** — `-W = Q_1\left(\dfrac{T_1 - T_2}{T_1}\right)`
  Work done BY a cyclic/steady engine drawing heat $Q_1$ at $T_1$ and rejecting at $T_2$, in the reversible limit $S_{\text{gen}}=0$. Real engines produce less: $-W = Q_1(T_1-T_2)/T_1 - T_2 S_{\text{gen}}$ (Eq. 4.3-3).
- **(4.3-5) Carnot (maximum) engine efficiency** — `\eta = \dfrac{-W}{Q_1} = \dfrac{T_1 - T_2}{T_1}`
  Depends only on the absolute temperature levels, not the working fluid or design. An upper bound; real engines reach roughly half this. $T_2 = 0$ would be required to convert all heat to work — impossible.
- **(4.4-3) Ideal-gas entropy change with $T$ and $P$ (constant $C_P^{*}$)** — `\underline{S}(T_2,P_2) - \underline{S}(T_1,P_1) = C_P^{*}\ln\!\left(\dfrac{T_2}{T_1}\right) - R\ln\!\left(\dfrac{P_2}{P_1}\right)`
  The most-used entropy-change formula. Note the MINUS sign on the pressure term (entropy falls as pressure rises at fixed $T$). Uses the ideal-gas heat capacity $C_P^{*}$ from Ch.3; for temperature-dependent $C_P^{*}$ the term becomes $\int C_P^{*}\,dT/T$.
- **(4.4-2) Ideal-gas entropy change with $T$ and $V$ (constant $C_V^{*}$)** — `\underline{S}(T_2,\underline{V}_2) - \underline{S}(T_1,\underline{V}_1) = C_V^{*}\ln\!\left(\dfrac{T_2}{T_1}\right) + R\ln\!\left(\dfrac{\underline{V}_2}{\underline{V}_1}\right)`
  Companion of Eq. 4.4-3 in the $(T,\underline{V})$ variables; here the volume term is POSITIVE. Related by the ideal-gas law and $C_P^{*} = C_V^{*} + R$.
- **(4.4-6) Entropy change of a solid or liquid** — `\underline{S}(T_2) - \underline{S}(T_1) = \int_{T_1}^{T_2} C_P\,\dfrac{dT}{T}`
  For a solid or liquid the molar volume is nearly independent of $T$ and $P$ (so the $P\,d\underline{V}$ term drops) and $C_V \approx C_P$. For constant $C_P$ this is $C_P\ln(T_2/T_1)$.

## Codes (Jupyter notebooks / Python scripts)

Distributed from https://github.com/emfurst/cbet6e. When you open one in this editor, run and
modify the cells directly rather than narrating code you cannot execute.

- `ch1/First_Notebook_HW1_CHEG231.ipynb` — Your first Jupyter notebook: molar volume from the van der Waals EOS
- `ch1/SIS_Problem_1_2_CHEG231_EMF.ipynb` — Reading and plotting data (Problem 1.2)
- `ch1/Heat_capacity_fitting_CHEG231.ipynb` — Curve fitting a heat-capacity polynomial
- `ch1/When_ideal.ipynb` — When is a gas ideal?
- `ch1/Heat_capacity.ipynb` — Enthalpy and internal-energy changes from the heat capacity
- `ch3/LJ_interaction_energy_figure.ipynb` — Figure 3.3-5: the interaction energy between two molecules
