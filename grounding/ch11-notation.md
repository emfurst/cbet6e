# Notation and conventions — Sandler and Furst, "Chemical, Biochemical, and Engineering Thermodynamics," 6e

*Scope: through Chapter 11 (Other Types of Phase Equilibria in Fluid Mixtures).* This file documents the symbols, units, and sign conventions
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
- In cycle analysis, integer subscripts 1, 2, 3, 4, … label the successive STATES around the cycle, and a DOUBLE subscript labels the flow on the PATH between two states: $\dot{W}_{23}$ is the work on the path from state 2 to state 3, and $\dot{Q}_{12}$ is the heat on the path from state 1 to state 2. Do not read $\dot{W}_{23}$ as a single indexed stream $k$ (Ch.2) — it is a path quantity between two numbered states.
- Work and heat also carry a device label: $\dot{W}_T$ (turbine), $\dot{W}_P$ (pump or compressor), $\dot{Q}_B$ (boiler/evaporator), $\dot{Q}_C$ (condenser). The Chapter 3 sign convention still holds — flows INTO the system and work ON the system are positive — so a work-PRODUCING turbine has $\dot{W}_T < 0$ and a condenser REJECTING heat has $\dot{Q}_C < 0$.
- Two DIFFERENT figures of merit, easily confused. Thermal efficiency $\eta$ (power cycles) = net work out / heat in, always $< 1$. Coefficient of performance C.O.P. (refrigerators, heat pumps) = heat moved / net work in, routinely $> 1$. A C.O.P. above 1 is normal and is NOT a second-law violation — moving heat with work is not creating energy. Do not report a C.O.P. as an "efficiency" or cap it at 1.
- A subscript $S$ marks the ISENTROPIC (reversible-adiabatic) reference state: $T_{2,S}$ is the outlet temperature of a reversible-adiabatic path to $P_2$. Isentropic efficiency rates a real device against this ideal. Symbol collision to watch: $W_S$ (isentropic work, capital-$S$ subscript) is distinct from $W_s$ (shaft work, Ch.3, lowercase-$s$) and from $S$ (entropy).
- This chapter uses TWO different ideal-gas markers, and they are not interchangeable. The asterisk ($^{*}$, carried from Ch.3) is reserved for the ideal-gas HEAT CAPACITY only: $C_P^{*}(T)$, $C_V^{*}(T)$. The superscript "IG" marks an ideal-gas STATE property evaluated at the actual $T$ and $P$: $\underline{H}^{\mathrm{IG}}$, $\underline{S}^{\mathrm{IG}}$, $\underline{U}^{\mathrm{IG}}$, $\underline{V}^{\mathrm{IG}}$, $T^{\mathrm{IG}}$. Do not write $\underline{H}^{*}$ or $C_P^{\mathrm{IG}}$ — match the book: $^{*}$ for the heat capacity, $^{\mathrm{IG}}$ for the enthalpy/entropy/volume reference state.
- A departure (or residual) function is the REAL property minus the IDEAL-GAS property at the SAME temperature and pressure: $(\underline{H}-\underline{H}^{\mathrm{IG}})_{T,P}$ and $(\underline{S}-\underline{S}^{\mathrm{IG}})_{T,P}$ (Eqs. 6.4-22, 6.4-23). The book writes it out as $(\underline{H}-\underline{H}^{\mathrm{IG}})$ — it does NOT use a $\Delta \underline{H}^{\mathrm{dep}}$ or $\underline{H}^{R}$ symbol. Sign matters: real minus ideal. WATCH OUT — the generalized charts Fig. 6.6-4 and 6.6-5 plot the reverse, $(\underline{H}^{\mathrm{IG}}-\underline{H})/T_c$ and $(\underline{S}^{\mathrm{IG}}-\underline{S})$, so a sign flip is needed when reading them (book footnotes 15 and 16). A real-fluid change is then (ideal-gas change) + (departure at end state) - (departure at start state).
- Preferred equation-of-state (EOS) forms in this book, written on a molar basis (underbar). van der Waals: $P = RT/(\underline{V}-b) - a/\underline{V}^2$ (Eq. 6.2-38b). Redlich-Kwong: $P = RT/(\underline{V}-b) - a/[T^{1/2}\underline{V}(\underline{V}+b)]$ (Eq. 6.4-1). Soave-Redlich-Kwong (SRK): $P = RT/(\underline{V}-b) - a(T)/[\underline{V}(\underline{V}+b)]$ (Eq. 6.4-1b). Peng-Robinson (PR): $P = RT/(\underline{V}-b) - a(T)/[\underline{V}(\underline{V}+b)+b(\underline{V}-b)]$ (Eq. 6.4-2). PR is the DEFAULT working EOS of the text. $a$ is the attractive-energy parameter and $b$ the excluded-volume (co-volume) parameter; in vdW both are constants, in SRK/PR $a=a(T)$ depends on temperature through $\alpha(T)$.
- The symbol $\alpha$ has THREE distinct meanings in this chapter — disambiguate by context. (i) $\alpha$ = coefficient of thermal expansion $\dfrac{1}{\underline{V}}(\partial \underline{V}/\partial T)_P$ (Eq. 6.2-3), units K$^{-1}$. (ii) $\alpha(T)$ = the dimensionless temperature-dependent factor in the SRK/PR attractive term, with $\alpha(T_c)=1$ (Eqs. 6.7-3, 6.7-8). (iii) $\alpha$ = the coefficient of $Z^2$ in the cubic form $Z^3+\alpha Z^2+\beta Z+\gamma=0$ (Eq. 6.4-4), alongside $\beta,\gamma$. These are unrelated; never carry a value of one into another.
- Two different $\kappa$ (kappa) symbols. WITH a subscript, $\kappa_T$ is the isothermal compressibility $-\dfrac{1}{\underline{V}}(\partial \underline{V}/\partial P)_T$ (Eq. 6.2-4), units kPa$^{-1}$. WITHOUT a subscript, plain $\kappa$ is the dimensionless PR/SRK parameter fixing the temperature dependence of $\alpha(T)$, a function of the acentric factor $\omega$ (Eqs. 6.7-4, 6.7-9). Keep the subscript on the compressibility to tell them apart.
- The dimensionless EOS parameters $A$ and $B$ (used in the cubic-in-$Z$ form, Table 6.4-3) COLLIDE with earlier symbols. $A = aP/(RT)^2$ (van der Waals, Soave, PR) or $A = aP/(R^2T^{2.5})$ (Redlich-Kwong); $B = bP/RT$ — both dimensionless. Do NOT confuse the dimensionless $A$ with the Helmholtz energy $A=U-TS$ (Ch.4), nor the dimensionless $B$ with the second virial coefficient $B(T)$ (which has units m$^3$/mol). Judge by context: inside the cubic $Z^3+\alpha Z^2+\beta Z+\gamma=0$ they are the dimensionless parameters.
- Reduced properties are the actual property divided by its critical value: $T_r=T/T_c$, $P_r=P/P_c$, $\underline{V}_r=\underline{V}/\underline{V}_c$ (Eq. 6.6-5). A subscript $c$ marks a critical-point value ($T_c$, $P_c$, $\underline{V}_c$, $Z_c$). The critical and reduced molar volumes are written with the underbar ($\underline{V}_c$, $\underline{V}_r$), like every other molar quantity.
- The criterion for equilibrium DEPENDS ON THE CONSTRAINTS on the system (Table 7.1-1), and all three forms follow from $S_{\text{gen}} \ge 0$: for a closed ISOLATED system at constant $U$ and $V$, entropy is a MAXIMUM ($S=\text{max}$, $dS=0$, $d^2S<0$); for a closed system at constant $T$ and $V$, the Helmholtz energy is a MINIMUM ($A=\text{min}$, $dA=0$, $d^2A>0$); for a closed system at constant $T$ and $P$, the Gibbs energy is a MINIMUM ($G=\text{min}$, $dG=0$, $d^2G>0$). Do NOT state a blanket "entropy of the universe increases" or "Gibbs energy is always minimized" -- match the criterion to the constraints. The constant-$T$,$P$ Gibbs-energy form is the one used most in this book.
- At equilibrium between two phases of a PURE component the temperatures, pressures, and MOLAR GIBBS ENERGIES are equal: $T^{\mathrm{I}}=T^{\mathrm{II}}$, $P^{\mathrm{I}}=P^{\mathrm{II}}$, $\underline{G}^{\mathrm{I}}=\underline{G}^{\mathrm{II}}$ (Eqs. 7.2-15). For a ONE-component system the phase-equilibrium condition is written with the molar Gibbs energy $\underline{G}$, NOT the chemical potential $\mu$. The symbol $\mu$ (chemical potential / partial molar Gibbs energy) is introduced later for MIXTURES (Ch.8 onward); for a pure substance $\mu = \underline{G}$, but this book uses $\underline{G}$ here.
- Because the fugacity is defined from the Gibbs energy (Eq. 7.4-6), the phase-equilibrium condition $\underline{G}^{\mathrm{I}}=\underline{G}^{\mathrm{II}}$ is equivalent to EQUALITY OF FUGACITIES $f^{\mathrm{I}}(T,P)=f^{\mathrm{II}}(T,P)$ (Eq. 7.4-7a), or equivalently of fugacity coefficients $\phi^{\mathrm{I}}=\phi^{\mathrm{II}}$ (Eq. 7.4-7b), at the same $T$ and $P$. For vapor-liquid equilibrium this is $f^{L}=f^{V}$. This is the working criterion used for phase-equilibrium calculations throughout the book.
- A pure single phase is intrinsically STABLE only if it satisfies the thermal stability criterion $C_V>0$ (Eq. 7.2-12) AND the mechanical stability criterion $\left(\partial P/\partial \underline{V}\right)_T<0$, equivalently $\kappa_T>0$ (Eq. 7.2-13). These are consequences of $d^2S<0$ at constant $U,V$. Note the mechanical criterion is a STRICT inequality $<0$ (pressure must fall as molar volume rises at fixed $T$); the point where $\left(\partial P/\partial \underline{V}\right)_T=0$ is a limit of stability, not a stable state.
- The fluid CRITICAL POINT is the inflection point on the critical isotherm in the $P$-$\underline{V}$ plane (point $C$, the peak of the vapor-liquid coexistence dome and the terminus of the vapor-pressure curve in the $P$-$T$ plane). Mathematically it is where both the first and second isothermal volume derivatives of pressure vanish: $\left(\partial P/\partial \underline{V}\right)_T=0$ and $\left(\partial^2 P/\partial \underline{V}^2\right)_T=0$ at $(T_c,P_c)$. At the critical point the vapor and liquid become indistinguishable, so $\underline{G}^L=\underline{G}^V$, $\underline{S}^L=\underline{S}^V$, and all other properties coincide.
- Phases are labeled by SUPERSCRIPTS: $L$ (liquid), $V$ (vapor), $S$ (solid), and generic Roman numerals $\mathrm{I}$, $\mathrm{II}$ for unspecified phases; e.g. $\underline{H}^L$, $\underline{H}^V$, $\underline{H}^S$ are molar enthalpies of the liquid, vapor, and solid. WATCH THE COLLISION: a superscript $S$ means "solid", while a plain $S$ is the entropy -- disambiguate by position. Coexistence (saturation) pressures carry superscripts too: $P^{\text{vap}}(T)$ is the liquid vapor pressure, $P^{\text{sub}}(T)$ the solid sublimation pressure, and $P^{\text{sat}}$ a GENERAL coexistence pressure (equations in $P^{\text{sat}}$ apply to both). Phase-change differences are written $\Delta\theta=\theta^{\mathrm{I}}-\theta^{\mathrm{II}}$ (e.g. $\Delta_{\text{vap}}\underline{H}=\underline{H}^V-\underline{H}^L$).
- THE key notation of this chapter, and the FIRST use of the overbar in the book. A molar (per-mole) quantity is an UNDERBAR everywhere (Ch.1), and that does not change. New here: an overbar WITH A SPECIES SUBSCRIPT $i$, e.g. $\bar{H}_i$, $\bar{V}_i$, $\bar{G}_i$, denotes the PARTIAL MOLAR property of species $i$ in a mixture, defined by $\bar{\theta}_i = (\partial (N\theta)/\partial N_i)_{T,P,N_{j\neq i}}$ (Eq. 8.1-12). A partial molar property is a true mixture property (it must be measured for each mixture) and is generally NOT equal to the pure-component molar property of species $i$: $\bar{\theta}_i \neq \underline{\theta}_i$. Because both accents now appear WITH a species subscript, tell them apart by the ACCENT: underbar $\underline{V}_i$ = pure-component molar property of species $i$; overbar $\bar{V}_i$ = partial molar property of species $i$. (An overbar is never a plain molar quantity — that is always the underbar.)
- The chemical potential of species $i$ is DEFINED as its partial molar Gibbs energy: $\mu_i \equiv \bar{G}_i = (\partial G/\partial N_i)_{T,P,N_{j\neq i}}$ (Sec. 8.2). The symbols $\mu_i$ and $\bar{G}_i$ are interchangeable in this book. Numerically $\bar{G}_i$ also equals $(\partial U/\partial N_i)_{S,V,N_{j\neq i}} = (\partial A/\partial N_i)_{T,V,N_{j\neq i}} = (\partial H/\partial N_i)_{S,P,N_{j\neq i}}$, but ONLY the constant-$T,P$ derivative is a partial molar quantity; the others are equal in value yet are not "partial molar" properties (a partial molar derivative holds $T$ and $P$ fixed). SYMBOL COLLISION: in Ch.6 the bare symbol $\mu$ denotes the Joule-Thomson coefficient $(\partial T/\partial P)_H$ — a completely different quantity. The chemical potential ALWAYS carries a species subscript, $\mu_i$; a bare $\mu$ is the Joule-Thomson coefficient.
- Subscripts label SPECIES; superscripts label PHASES. A subscript $i$ (or $1,2,\dots$) is the chemical species; a superscript Roman numeral $\mathrm{I},\mathrm{II},\mathrm{III}$ (or $\mathrm{L},\mathrm{V}$ for liquid/vapor) is the phase. Thus $\bar{G}_i^{\mathrm{II}}$ is the partial molar Gibbs energy of species $i$ in phase II, and $N_i^{\mathrm{I}}$ is the moles of species $i$ in phase I. Mole fraction $x_i = N_i/N$; by convention $x_i$ is used for a liquid (or generic) phase and $y_i$ for the vapor phase in the later vapor-liquid-equilibrium chapters. Do not confuse the species subscript with the phase superscript.
- Calligraphic capitals are dimensionless COUNTS, not thermodynamic properties: $\mathcal{C}$ = number of components, $\mathcal{M}$ = number of independent chemical reactions, $\mathcal{P}$ = number of phases, $\mathcal{F}$ = number of degrees of freedom. They enter the Gibbs phase rule $\mathcal{F} = \mathcal{C} - \mathcal{M} - \mathcal{P} + 2$ (Eq. 8.9-6). The one-component phase rule of Ch.7 uses the same calligraphic $\mathcal{F}$ (there reducing to $\mathcal{F} = 3 - \mathcal{P}$).
- Reaction bookkeeping for multiple reactions. $\nu_{ij}$ is the stoichiometric coefficient of species $i$ in the $j$th independent reaction — the Ch.2 single-reaction $\nu_i$ generalized with a second subscript. Sign convention (from Ch.2): $\nu > 0$ for products, $\nu < 0$ for reactants, $\nu = 0$ for inerts. $X_j$ is the molar extent of the $j$th independent reaction, so the moles of species $i$ are $N_i = N_{i,0} + \sum_j \nu_{ij} X_j$ (Eq. 8.3-5). Track only a set of $\mathcal{M}$ INDEPENDENT reactions (the smallest set whose linear combinations give all others).
- Standard state and formation properties. A superscript $^{\circ}$ denotes the standard state: for a gas, the ideal gas at 1 bar and the temperature of interest (the $P=0$ state cannot be used); for solids and solutes, the conventions of Appendix A.IV. $\Delta_{\text{f}}\underline{H}_i^{\circ}$ and $\Delta_{\text{f}}\underline{G}_i^{\circ}$ are the standard-state molar enthalpy and Gibbs energy of FORMATION of species $i$ from its elements. The enthalpy (and Gibbs energy) zero is set once per element in its stable state at 25 C, 1 bar; every compound value is then referenced to that. Standard heats/Gibbs energies of reaction follow as $\Delta_{\text{rxn}}(\cdot)^{\circ} = \sum_i \nu_i\,\Delta_{\text{f}}(\cdot)_i^{\circ}$.
- Property changes on mixing. $\Delta_{\text{mix}}\theta \equiv \theta_{\text{mixture}} - \sum_i x_i\,\theta_i^{\text{pure}}$ (per mole) is the change in a molar property when the pure components are mixed at fixed $T$ and $P$; common cases are $\Delta_{\text{mix}}V$, $\Delta_{\text{mix}}H$ (the heat of mixing), and $\Delta_{\text{mix}}G$. An IDEAL mixture is one whose partial molar volume, internal energy, and enthalpy equal the pure-component molar values, so that $\Delta_{\text{mix}}V = \Delta_{\text{mix}}H = 0$; real mixtures generally deviate (e.g. sulfuric acid + water mixes strongly exothermically). The ideal SOLUTION used for the Gibbs energy and fugacity is developed in Ch.9.
- The mixture fugacity and mixture fugacity coefficient. $\bar{f}_i$ is the fugacity of species $i$ IN A MIXTURE and $\bar{\phi}_i = \bar{f}_i/(x_i P)$ its fugacity coefficient. BOTH CARRY AN OVERBAR — the same overbar Ch.8 introduced for partial molar properties, because both are mixture properties of species $i$. Distinguish them from the Ch.6-7 pure-fluid quantities, which carry NO accent: $f$ is the fugacity of a pure fluid and $\phi = f/P$ its fugacity coefficient. So $\phi$ (no accent) is a pure-fluid property and $\bar{\phi}_i$ (overbar, species subscript) is a mixture property. A caret is NEVER used for either; if you have seen $\hat{\phi}_i$ attributed to this book, it is wrong.
- THE most dangerous notation in the book: $\phi$ carries five distinct meanings across Ch.6-9, several separated only by an accent, a prime, or a case change. (1) $\phi$, no accent — fugacity coefficient of a PURE fluid, $f/P$ (Ch.6-7). (2) $\bar{\phi}_i$, overbar — fugacity coefficient of species $i$ IN A MIXTURE (Sec. 9.4). (3) $\phi_i$, bare with a species subscript — the VOLUME (segment) FRACTION of species $i$, $x_ir_i/\sum_j x_jr_j$, used in Flory-Huggins and UNIQUAC (Sec. 9.5). (4) $\phi_i^{\prime}$, primed — the volume fraction weighted by $r_i^{3/4}$, used ONLY in the modified-UNIFAC combinatorial term (Sec. 9.6). (5) $\Phi_i$, capital — the volume fraction defined from MOLAR VOLUMES, $x_i\underline{V}_i/\underline{V}_{\text{mix}}$, used only in regular solution theory (Sec. 9.6). Meanings (2) and (3) are declared nine lines apart in the same notation box and differ by one overbar. When reading or generating any expression in Secs. 9.5-9.6, establish which $\phi$ is meant before doing anything else.
- Subscripts carry the species-versus-group distinction, and it is load-bearing. Subscripts $i$ and $j$ denote SPECIES; $k$, $m$, and $n$ denote FUNCTIONAL GROUPS. Paired with the case convention (lowercase for species, uppercase for groups) this separates every group-contribution quantity from its species analog: $\theta_i$ / $\Theta_m$, $\gamma_i$ / $\Gamma_k$, $r_i,q_i$ / $R_k,Q_k$, $x_i$ / $X_m$, $\tau_{ij}$ / $\Psi_{mn}$. IT ALSO SEPARATES TWO UNRELATED USES OF $a$: $a_{ij}$ is an equation-of-state parameter (Sec. 9.4, Sec. 9.9) or a Wohl expansion parameter (Sec. 9.5, App. A9.2), while $a_{mn}$ is a UNIFAC group-interaction parameter (Sec. 9.6). They are different quantities with different units and different symmetry — see convention interaction_parameter_asymmetry.
- UNIFAC group-interaction parameters are NOT symmetric, and the chapter contains a symmetric look-alike that makes this easy to get wrong. In modified (Dortmund) UNIFAC, $\Psi_{mn} = \exp[-(a_{mn}+b_{mn}T+c_{mn}T^2)/T]$, and in general $a_{mn} \neq a_{nm}$, hence $\Psi_{mn} \neq \Psi_{nm}$. Every ordered pair of main groups carries its own three constants ($a$ in K, $b$ dimensionless, $c$ in K$^{-1}$), which is why the subscript ORDER in the residual term must be preserved. CONTRAST with Wohl: Appendix A9.2 states $a_{ij} = a_{ji}$ for the Wohl expansion parameter, and that is CORRECT — Wohl parameters genuinely are symmetric. So the chapter carries one symmetric $a$ and one asymmetric $a$, told apart only by whether the subscripts are species or group indices. Do not generalize the symmetry of one to the other.
- An activity coefficient is meaningless without its standard state. $\gamma_i$ is defined by $\bar{f}_i = x_i\gamma_i f_i$, where $f_i$ is the fugacity of species $i$ in the CHOSEN STANDARD STATE, and the numerical value of $\gamma_i$ depends entirely on that choice. The default in Secs. 9.5-9.6 is the pure-component (Lewis-Randall) standard state, for which $\gamma_i \to 1$ as $x_i \to 1$. Henry’s law uses a different one: $\gamma_i^{*}$ is normalized so $\gamma_i^{*} \to 1$ as $x_i \to 0$ (mole-fraction basis) and $\gamma_i^{\square}$ is its molality-based analog. A reported $\gamma_i$ without a stated standard state is not a number anyone can use — always ask which one.
- An EXCESS property is the departure from the IDEAL MIXTURE, not from the ideal gas: $\underline{\theta}^{\text{ex}} = \underline{\theta}_{\text{mix}} - \underline{\theta}_{\text{mix}}^{\text{IM}}$, where the ideal-mixture value is the mole-fraction-weighted sum of pure-component molar properties (plus the ideal entropy of mixing for $\underline{S}$ and $\underline{G}$). The chapter is built on $\underline{G}^{\text{ex}}$, which is related to the activity coefficients by $\underline{G}^{\text{ex}}/RT = \sum_i x_i \ln\gamma_i$ (Eq. 9.3-12). THAT IDENTITY IS THE CHAPTER’S BEST CHECK: any model that supplies both a $\underline{G}^{\text{ex}}$ and its activity coefficients must satisfy it exactly. It is what found the two errata corrected in this file. $\underline{A}^{\text{ex}}$ appears in Sec. 9.9 for the same reason at constant volume.
- Two Henry’s law constants, and one of them collides with enthalpy. $H_i(T,P)$ is the mole-fraction-based Henry’s law constant of species $i$, defined by $\bar{f}_i \to H_i x_i$ as $x_i \to 0$; $\mathcal{H}_i$ (script) is the MOLALITY-based analog, with units of pressure per molality. $H_i$ is typeset identically to a partial molar enthalpy in careless rendering and identically to the pure-component enthalpy $H$ of Ch.3-6; only context and the script $\mathcal{H}$ separate them. When generating text or code for Sec. 9.7, never abbreviate the Henry’s law constant to a bare $H$.
- Sec. 9.10 notation, reused in Ch.15. $\gamma_{\pm}$ is the MEAN IONIC activity coefficient — a single measurable quantity standing for the geometric mean of the individual ion activity coefficients, which are not separately measurable. $M_i$ is molality (moles of solute per kilogram of solvent), $M_{\pm}$ the mean ionic molality, and $I = \tfrac{1}{2}\sum_i M_i z_i^2$ the IONIC STRENGTH, where $z_+$ and $z_-$ are the cation and anion charges and $\nu_+$, $\nu_-$ the numbers of each produced by ionization. $M$ for molality collides with $M$ for molecular weight in Ch.1-3; in Secs. 9.10 and 15.x, $M$ with a species subscript is always molality.

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
| `\eta` | thermal efficiency of a cycle (also isentropic efficiency) | - |
| `\mathrm{C.O.P.}` | coefficient of performance of a refrigeration / heat-pump cycle | - |
| `\dot{W}_T` | rate of turbine work ($<0$ when producing work) | J/s |
| `\dot{W}_P` | rate of pump / compressor work ($>0$ when consuming work) | J/s |
| `\dot{Q}_B` | rate of heat into the boiler / evaporator | J/s |
| `\dot{Q}_C` | rate of heat from the condenser ($<0$) | J/s |
| `W_S` | isentropic (reversible-adiabatic) work — capital-$S$, NOT shaft work $W_s$ | J/kg |
| `K_T` | compression (pressure) ratio, $P_{\text{high}}/P_{\text{low}}$ | - |
| `T_b` | normal boiling temperature (at atmospheric pressure) | K |
| `\Delta_{\text{vap}}\underline{U}` | molar internal energy change on vaporization | J/mol |
| `\Delta_{\text{vap}}\hat{U}` | specific internal energy change on vaporization | J/kg |
| `a` | EOS attractive-energy parameter $a$ (or $a(T)$ for SRK/PR) | Pa·m⁶/mol² (vdW); depends on EOS |
| `b` | EOS excluded-volume (co-volume) parameter $b$ | m³/mol |
| `A` | dimensionless EOS parameter $A=aP/(RT)^2$ (collides with Helmholtz $A$) | - |
| `B` | dimensionless EOS parameter $B=bP/RT$ (collides with virial $B(T)$) | - |
| `B(T)` | second virial coefficient | m³/mol |
| `C(T)` | third virial coefficient | (m³/mol)² |
| `Z` | compressibility factor $Z=P\underline{V}/RT$ ($=1$ for an ideal gas) | - |
| `Z_c` | compressibility factor at the critical point | - |
| `P_c` | critical pressure | kPa (tables use MPa) |
| `T_c` | critical temperature | K |
| `\underline{V}_c` | critical molar volume | m³/mol |
| `T_r` | reduced temperature $T_r=T/T_c$ | - |
| `P_r` | reduced pressure $P_r=P/P_c$ | - |
| `\underline{V}_r` | reduced molar volume $\underline{V}_r=\underline{V}/\underline{V}_c$ | - |
| `\alpha` | coefficient of thermal expansion $\frac{1}{\underline{V}}(\partial \underline{V}/\partial T)_P$ | K⁻¹ |
| `\alpha(T)` | temperature-dependent factor in the SRK/PR attractive term ($\alpha(T_c)=1$) | - |
| `\kappa_T` | isothermal compressibility $-\frac{1}{\underline{V}}(\partial \underline{V}/\partial P)_T$ | kPa⁻¹ |
| `\kappa` | PR/SRK parameter setting the $T$-dependence of $\alpha(T)$ (function of $\omega$) | - |
| `\mu` | Joule-Thomson coefficient $(\partial T/\partial P)_H$ | K/kPa |
| `\omega` | acentric factor (Pitzer) | - |
| `(\underline{H}-\underline{H}^{\mathrm{IG}})_{T,P}` | molar enthalpy departure (real minus ideal gas at same $T,P$) | J/mol |
| `(\underline{S}-\underline{S}^{\mathrm{IG}})_{T,P}` | molar entropy departure (real minus ideal gas at same $T,P$) | J/(mol·K) |
| `f` | fugacity (defined in Sec. 7.4; the effective pressure that makes real-fluid relations take the ideal form) | kPa |
| `\phi` | fugacity coefficient, $\phi = f/P$ | - |
| `f_{\text{sat}}(T)` | fugacity of a phase at its phase-change (saturation) pressure at temperature $T$ | kPa |
| `\mathcal{F}` | number of degrees of freedom (Gibbs phase rule) | - |
| `P^{\text{vap}}(T)` | vapor pressure (liquid-vapor coexistence pressure) at $T$ | kPa |
| `P^{\text{sub}}(T)` | sublimation pressure (solid-vapor coexistence pressure) at $T$ | kPa |
| `P^{\text{sat}}` | general phase-coexistence (saturation) pressure | kPa |
| `P_t` | triple-point pressure | kPa |
| `T_t` | triple-point temperature | K |
| `\omega^{V}` | fraction of a two-phase mixture that is vapor (molar basis; $\omega^{L}=1-\omega^{V}$). COLLISION: $\omega$ alone is the Ch.6 acentric factor, and $\omega^{\mathrm{I}}$ in Ch.3 is the steam quality | - |
| `\bar{\theta}_i` | partial molar property of species $i$ (generic; e.g. $\bar{V}_i,\bar{H}_i,\bar{G}_i$) $= (\partial(N\theta)/\partial N_i)_{T,P,N_{j\neq i}}$ | varies (per mol) |
| `\bar{V}_i` | partial molar volume of species $i$ | m^3/mol |
| `\bar{H}_i` | partial molar enthalpy of species $i$ (NOT the pure-component molar enthalpy) | J/mol |
| `\bar{S}_i` | partial molar entropy of species $i$ | J/(mol·K) |
| `\bar{G}_i` | partial molar Gibbs energy of species $i$ (= chemical potential $\mu_i$) | J/mol |
| `\underline{V}_i` | pure-component molar volume of species $i$ (underbar); distinct from the partial molar $\bar{V}_i$ (overbar) | m^3/mol |
| `\mu_i` | chemical potential of species $i$ ($\equiv \bar{G}_i$) | J/mol |
| `x_i` | mole fraction of species $i$ ($= N_i/N$; liquid or generic phase) | - |
| `y_i` | mole fraction of species $i$ in the vapor phase (VLE convention) | - |
| `w_i` | mass fraction of species $i$ in a phase | - |
| `w^{\mathrm{I}}` | fraction of total system mass in phase I | - |
| `\mathcal{C}` | number of components | - |
| `\mathcal{M}` | number of independent chemical reactions | - |
| `\mathcal{P}` | number of phases | - |
| `\nu_{ij}` | stoichiometric coefficient of species $i$ in the $j$th independent reaction (generalizes Ch.2 $\nu_i$) | - |
| `X_j` | molar extent of the $j$th independent reaction | mol |
| `\Delta_{\text{mix}}\theta` | change in a molar property on mixing (e.g. $\Delta_{\text{mix}}V$, $\Delta_{\text{mix}}H$, $\Delta_{\text{mix}}G$) | varies |
| `\Delta_{\text{rxn}}H^{\circ}` | standard-state heat (enthalpy change) of reaction; $\Delta_{\text{rxn}}G^{\circ}$ is the analog for Gibbs energy | J/mol |
| `\Delta_{\text{f}}\underline{H}_i^{\circ}` | standard-state molar enthalpy of formation of species $i$ ($\Delta_{\text{f}}\underline{G}_i^{\circ}$ = Gibbs analog) | J/mol |
| `\bar{f}_i` | fugacity of species $i$ in a mixture (OVERBAR, never a caret) | kPa |
| `\bar{\phi}_i` | fugacity coefficient of species $i$ in a mixture, $= \bar{f}_i/(x_i P)$ (OVERBAR; distinct from the pure-fluid $\phi$ of Ch.6-7 and from the volume fraction $\phi_i$) | - |
| `\gamma_i` | activity coefficient of species $i$; meaningless without its standard state | - |
| `\underline{G}^{\text{ex}}` | molar excess Gibbs energy — departure from the IDEAL MIXTURE, not the ideal gas | J/mol |
| `\underline{A}^{\text{ex}}` | molar excess Helmholtz energy (Sec. 9.9, Wong-Sandler) | J/mol |
| `\phi_i` | VOLUME (segment) fraction of species $i$, $x_ir_i/\sum_j x_jr_j$ — NOT a fugacity coefficient (Secs. 9.5, 9.6) | - |
| `\phi_i^{\prime}` | volume fraction weighted by $r_i^{3/4}$, used only in the modified-UNIFAC combinatorial term (Eq. 9.6-12a) | - |
| `\Phi_i` | volume fraction from molar volumes, $x_i\underline{V}_i/\underline{V}_{\text{mix}}$ — regular solution theory ONLY (kept as $\Phi_i$ deliberately, for consistency with earlier editions) | - |
| `\theta_i` | SURFACE AREA fraction of species $i$, $x_iq_i/\sum_j x_jq_j$ (Secs. 9.5, 9.6); an undecorated $\theta$ elsewhere in the book is a generic property | - |
| `r_i,\;q_i` | volume and surface area parameters of SPECIES $i$, summed from its functional groups | - |
| `R_k,\;Q_k` | volume and surface area parameters of FUNCTIONAL GROUP $k$ (Table 9.5-2, modified/Dortmund set) | - |
| `\Theta_m` | surface area fraction of functional group $m$ in a mixture | - |
| `X_m` | mole fraction of functional group $m$ in a mixture | - |
| `\Gamma_k` | activity coefficient of functional group $k$ in a mixture; $\Gamma_k^{(i)}$ in a pure fluid of species $i$ | - |
| `\nu_k^{(i)}` | number of functional groups of type $k$ in a molecule of species $i$ | - |
| `\tau_{ij}` | UNIQUAC interaction parameter for the species pair $i$-$j$; $\tau_{ij} \neq \tau_{ji}$ | - |
| `\Psi_{mn}` | UNIFAC interaction parameter for group $m$ with group $n$; $\Psi_{mn} \neq \Psi_{nm}$ | - |
| `a_{mn},\;b_{mn},\;c_{mn}` | modified-UNIFAC group-interaction parameters determining $\Psi_{mn}$; ORDERED pair, $a_{mn}\neq a_{nm}$. Not printed in the book — they ship with the code | K, -, K^{-1} |
| `\Lambda_{ij}` | Wilson parameters; $\Lambda_{ii}=1$ and $\Lambda_{ij}\neq\Lambda_{ji}$ | - |
| `\alpha` | NRTL nonrandomness parameter (NOT the van Laar $\alpha$ of Eq. 9.5-9, which is a different quantity in the same chapter) | - |
| `\delta_i` | solubility parameter of species $i$ (Table 9.6-1, regular solution theory). tabulated in calorie-based units; the code converts on read | (J/m^3)^{1/2} |
| `k_{ij}` | binary interaction parameter in the van der Waals one-fluid mixing rule (Eq. 9.4-9); fitted, and often but not always symmetric | - |
| `z` | lattice coordination number in UNIQUAC/UNIFAC, taken as 10 (not a compressibility factor $Z$, and not an ionic charge $z_+$) | - |
| `l_i` | UNIQUAC combinatorial abbreviation, $l_i=\tfrac{z}{2}(r_i-q_i)-(r_i-1)$ | - |
| `\chi` | Flory interaction parameter (the "chi parameter") in the Flory-Huggins model | - |
| `H_i(T,P)` | Henry’s law constant of species $i$, mole-fraction basis (collides visually with enthalpy $H$) | kPa |
| `\mathcal{H}_i` | molality-based Henry’s law constant of species $i$ | kPa/M |
| `\gamma_i^{*}` | Henry’s law activity coefficient of species $i$, mole-fraction basis ($\to 1$ as $x_i \to 0$) | - |
| `\gamma_i^{\square}` | Henry’s law activity coefficient of species $i$, molality basis | - |
| `\gamma_{\pm}` | mean ionic activity coefficient in an electrolyte solution (Sec. 9.10; reused in Ch.15) | - |
| `I` | ionic strength, $I=\tfrac{1}{2}\sum_i M_i z_i^2$ | M |
| `M_i` | molality of species $i$ ($M$ is molecular weight in Ch.1-3; with a species subscript here it is molality) | mol/kg solvent |

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
- **(5.1 (flash)) Joule–Thomson flash / vapor fraction (liquefaction)** — `\hat{H}_{\text{in}} = (1-\omega)\,\hat{H}^{L} + \omega\,\hat{H}^{V}`
  Throttling is isenthalpic ($\hat{H}_{\text{in}} = \hat{H}_{\text{out}}$), so the feed enthalpy sets the vapor mass fraction $\omega$ leaving the flash drum; $(1-\omega)$ is the liquefied fraction. This lever-rule split is the basis of the liquefaction yield. ($\omega$ here is the flash vapor fraction; cf. the phase mass fraction $\omega^{\mathrm{I}}$ of Ch.3.)
- **(5.2 (Rankine eta)) Thermal efficiency of a power cycle** — `\eta = \dfrac{-(\dot{W}_T + \dot{W}_P)}{\dot{Q}_B}`
  Net work out over heat in. With the Ch.3 sign convention the turbine term $\dot{W}_T<0$ dominates the pump term $\dot{W}_P>0$, so the numerator is positive. Equals the Carnot efficiency only for a fully reversible cycle between two reservoirs.
- **(5.2 (COP)) Coefficient of performance (refrigeration / heat pump)** — `\mathrm{C.O.P.} = \dfrac{-\dot{Q}_B}{\dot{W}_P + \dot{W}_T}`
  Heat removed from the low-temperature region per unit net work supplied. Routinely greater than 1. This is the refrigeration figure of merit — reported instead of efficiency for refrigerators, air conditioners, and heat pumps.
- **(5.2 (pump work)) Reversible pump work (incompressible liquid)** — `\dot{W}_P = \dot{M}\,\hat{V}_1\,(P_2 - P_1)`
  Because a liquid's specific volume and internal energy are nearly independent of pressure. Pump work is small compared with turbine work — a defining advantage of the Rankine cycle over a gas cycle.
- **(5.3-1) Isentropic (ideal) work of a turbine or compressor** — `W_S = \hat{H}(T_{2,S}, P_2) - \hat{H}(T_1, P_1)`
  Work along a reversible-adiabatic (constant-entropy) path to the outlet pressure $P_2$; $T_{2,S}$ is the isentropic outlet temperature. The ideal against which a real device is rated.
- **(5.3-3) Isentropic efficiency (turbine)** — `\eta = \dfrac{\hat{H}(T_2, P_2) - \hat{H}(T_1, P_1)}{\hat{H}(T_{2,S}, P_2) - \hat{H}(T_1, P_1)}`
  Actual work / isentropic work for a work-PRODUCING device (turbine). For a work-CONSUMING device (compressor, pump) the ratio inverts to $W_S/W_A$ so that $\eta<1$. Distinct from the cycle thermal efficiency above.
- **(5.3-4) Actual work from isentropic efficiency** — `W_A = \eta\,W_S`
  For a turbine. For a compressor, $W_A = W_S/\eta$. Use with 5.3-1 and the outlet-pressure state to get the real work and exit temperature.
- **(6.2-3) Coefficient of thermal expansion (definition)** — `\alpha = \dfrac{1}{\underline{V}}\left(\dfrac{\partial \underline{V}}{\partial T}\right)_P`
  A measurable volumetric property of the fluid. With $\kappa_T$ it lets one reduce most thermodynamic derivatives to $C_P$, $\alpha$, and $\kappa_T$. Same $\alpha$ symbol is also used for the SRK/PR $\alpha(T)$ factor and the cubic coefficient — see the conventions.
- **(6.2-4) Isothermal compressibility (definition)** — `\kappa_T = -\dfrac{1}{\underline{V}}\left(\dfrac{\partial \underline{V}}{\partial P}\right)_T`
  Note the leading MINUS sign (volume falls as pressure rises, so $\kappa_T>0$). Keep the $T$ subscript to distinguish it from the unsubscripted PR/SRK parameter $\kappa$.
- **(6.2-27) Joule-Thomson coefficient** — `\mu = \left(\dfrac{\partial T}{\partial P}\right)_H = -\dfrac{\underline{V} - T\left(\dfrac{\partial \underline{V}}{\partial T}\right)_P}{C_P}`
  Temperature change per unit pressure change in a constant-enthalpy (throttling) expansion. Derived from $d\underline{H}=C_P\,dT+[\underline{V}-T(\partial\underline{V}/\partial T)_P]\,dP$ (Eq. 6.2-22) by setting $d\underline{H}=0$. For an ideal gas $\mu=0$. In terms of $\alpha$: $\mu = -\dfrac{\underline{V}}{C_P}(1-T\alpha)$.
- **(6.2-35) Relation between $C_P$ and $C_V$** — `C_P - C_V = -T\left(\dfrac{\partial P}{\partial \underline{V}}\right)_T\left(\dfrac{\partial \underline{V}}{\partial T}\right)_P^{2} = \dfrac{T\underline{V}\alpha^2}{\kappa_T}`
  Always $\ge 0$. Reduces to $C_P-C_V=R$ for an ideal gas. Ties together the four common derivatives $C_P,C_V,\alpha,\kappa_T$, of which only three are independent.
- **(6.2-38b) van der Waals equation of state** — `P = \dfrac{RT}{\underline{V}-b} - \dfrac{a}{\underline{V}^{2}}`
  The prototype cubic EOS ($a,b$ constants). Historically first to predict vapor-liquid transition; used in the text mainly as a teaching prototype, not for design. Parameters from critical data via Eq. 6.6-4a.
- **(6.4-2) Peng-Robinson equation of state** — `P = \dfrac{RT}{\underline{V}-b} - \dfrac{a(T)}{\underline{V}(\underline{V}+b)+b(\underline{V}-b)}`
  The book's default engineering EOS (Peng & Robinson, 1976). Compare Redlich-Kwong $P=RT/(\underline{V}-b)-a/[T^{1/2}\underline{V}(\underline{V}+b)]$ (Eq. 6.4-1) and Soave-Redlich-Kwong $P=RT/(\underline{V}-b)-a(T)/[\underline{V}(\underline{V}+b)]$ (Eq. 6.4-1b). Note the PR denominator $\underline{V}(\underline{V}+b)+b(\underline{V}-b)$ — a frequent source of transcription error.
- **(6.4-4) Cubic (compressibility) form of the EOS** — `Z^{3} + \alpha Z^{2} + \beta Z + \gamma = 0`
  Any of the cubic EOS can be recast in $Z=P\underline{V}/RT$ with dimensionless $A=aP/(RT)^2$ (or $aP/(R^2T^{2.5})$ for RK) and $B=bP/RT$; $\alpha,\beta,\gamma$ from Table 6.4-3 (e.g. PR: $\alpha=-1+B$, $\beta=A-3B^2-2B$, $\gamma=-AB+B^2+B^3$). Solve for $Z$, then $\underline{V}=ZRT/P$. Here $\alpha,\beta,\gamma$ are cubic coefficients, NOT the thermal-expansion $\alpha$.
- **(6.4-5) Virial equation of state** — `\dfrac{P\underline{V}}{RT} = 1 + \dfrac{B(T)}{\underline{V}} + \dfrac{C(T)}{\underline{V}^{2}} + \cdots`
  A power series in reciprocal molar volume (density) about the ideal-gas limit. $B(T)$ and $C(T)$ are the second and third virial coefficients. Truncated at $B(T)$ it is a low-density approximation only — do not use above about 10 bar for most fluids. Not applicable to the liquid phase.
- **(6.4-22) Enthalpy departure function** — `(\underline{H}-\underline{H}^{\mathrm{IG}})_{T,P} = \int_{P=0}^{P}\left[\underline{V} - T\left(\dfrac{\partial \underline{V}}{\partial T}\right)_P\right]dP`
  Real minus ideal-gas enthalpy at fixed $T,P$. Entropy departure (Eq. 6.4-23): $(\underline{S}-\underline{S}^{\mathrm{IG}})_{T,P} = -\int_{P=0}^{P}[(\partial\underline{V}/\partial T)_P - R/P]\,dP$. For $\underline{V},T$-explicit EOS (vdW/PR) the more convenient forms are $\underline{H}-\underline{H}^{\mathrm{IG}}=RT(Z-1)+\int_\infty^{\underline{V}}[T(\partial P/\partial T)_{\underline{V}}-P]\,d\underline{V}$ (Eq. 6.4-27) and $\underline{S}-\underline{S}^{\mathrm{IG}}=R\ln Z+\int_\infty^{\underline{V}}[(\partial P/\partial T)_{\underline{V}}-R/\underline{V}]\,d\underline{V}$ (Eq. 6.4-28). The full property change is $\underline{H}(T_2,P_2)-\underline{H}(T_1,P_1)=[\underline{H}^{\mathrm{IG}}(T_2,P_2)-\underline{H}^{\mathrm{IG}}(T_1,P_1)]+(\underline{H}-\underline{H}^{\mathrm{IG}})_{T_2,P_2}-(\underline{H}-\underline{H}^{\mathrm{IG}})_{T_1,P_1}$ (Eq. 6.4-19).
- **(6.6-1) Critical-point conditions** — `\left(\dfrac{\partial P}{\partial \underline{V}}\right)_{T_c} = 0 \quad\text{and}\quad \left(\dfrac{\partial^{2} P}{\partial \underline{V}^{2}}\right)_{T_c} = 0 \quad\text{at } P_c,\ \underline{V}_c`
  The critical isotherm has an inflection point with a horizontal tangent in the $P$-$\underline{V}$ plane. Applying these to an EOS fixes its parameters in terms of critical properties.
- **(6.6-4a) van der Waals parameters from critical properties** — `a = \dfrac{27 R^{2} T_c^{2}}{64 P_c}, \qquad b = \dfrac{R T_c}{8 P_c}`
  From Eq. 6.6-1 applied to van der Waals. Gives the universal (and inaccurate) prediction $Z_c=P_c\underline{V}_c/RT_c=3/8=0.375$ (Eq. 6.6-3c), whereas real fluids have $Z_c\approx0.23$-$0.31$. Because $Z_c$ is wrong, using $T_c$ and $P_c$ (not $\underline{V}_c$) to get $a,b$ is preferred.
- **(6.6-7) Two-parameter corresponding states** — `Z = \dfrac{P\underline{V}}{RT} = Z(T_r, P_r)`
  The principle of corresponding states: all fluids obeying it have the same $Z$ at the same reduced temperature and pressure (Fig. 6.6-3). Systematic failures (different $Z_c$) motivate the three-parameter extension $Z=Z(T_r,P_r,\omega)$ using the acentric factor.
- **(6.6-omega) Acentric factor (definition)** — `\omega = -1.0 - \log_{10}\!\left[\dfrac{P^{\mathrm{vap}}(T_r=0.7)}{P_c}\right]`
  Pitzer's third corresponding-states parameter, built from the reduced vapor pressure at $T_r=0.7$ (near the normal boiling point). Roughly zero for spherical molecules (Ar, $\approx-0.004$) and larger for elongated/polar ones. Tabulated in Table 6.6-1; feeds the generalized SRK/PR $\kappa$.
- **(6.7-1) Generalized Peng-Robinson parameters** — `\begin{aligned} a(T) &= 0.45724\,\dfrac{R^{2}T_c^{2}}{P_c}\,\alpha(T), \qquad b = 0.07780\,\dfrac{R T_c}{P_c} \\ \sqrt{\alpha(T)} &= 1 + \kappa\left(1 - \sqrt{T/T_c}\right), \qquad \kappa = 0.37464 + 1.54226\,\omega - 0.26992\,\omega^{2} \end{aligned}`
  PR as a three-parameter ($T_c,P_c,\omega$) generalized EOS (Eqs. 6.7-1 to 6.7-4); $\alpha(T_c)=1$. SRK uses the same structure with different constants (Eqs. 6.7-6 to 6.7-9): $a(T)=0.42748\,R^2T_c^2/P_c\,\alpha(T)$, $b=0.08664\,RT_c/P_c$, $\kappa=0.480+1.574\,\omega-0.176\,\omega^2$. The temperature derivative $da/dT=-0.45724\,(R^2T_c^2/P_c)\,\kappa\sqrt{\alpha/(TT_c)}$ is needed for the PR enthalpy/entropy departures (Eqs. 6.4-29, 6.4-30).
- **(7.1-5) Equilibrium criterion, closed isolated system (constant $U$, $V$)** — `S = \text{maximum}, \qquad dS = 0 \ \text{and}\ d^2S < 0 \ \text{at equilibrium (const } U,V)`
  Follows from $S_{\text{gen}}\ge 0$ (Ch.4). $dS=0$ locates a candidate state; $d^2S<0$ makes it a true (stable) maximum rather than a metastable or unstable one. The companion criteria under other constraints are $A=\text{min}$ at constant $T,V$ (Eq. 7.1-10) and $G=\text{min}$ at constant $T,P$ (Eq. 7.1-12).
- **(7.1-12) Equilibrium criterion, closed system at constant $T$ and $P$** — `G = \text{minimum}, \qquad dG = 0 \ \text{and}\ d^2G > 0 \ \text{at equilibrium (const } T,P)`
  The most-used equilibrium criterion in this book. Derived by eliminating $\dot{Q}$ between the energy and entropy balances to get $dG/dt = -T\dot{S}_{\text{gen}}\le 0$. For a nonuniform system it again yields $P$ and $\underline{G}$ equal across regions (Eq. 7.1-9c).
- **(7.2-12) Stability criteria for a pure phase** — `C_V > 0 \qquad \text{and} \qquad \left(\dfrac{\partial P}{\partial \underline{V}}\right)_T < 0 \ \Longleftrightarrow\ \kappa_T = -\dfrac{1}{\underline{V}}\left(\dfrac{\partial \underline{V}}{\partial P}\right)_T > 0`
  Thermal (Eq. 7.2-12) and mechanical (Eq. 7.2-13) stability. Consequences of $d^2S<0$. Any EOS region with $(\partial P/\partial \underline{V})_T>0$ (e.g. inside a van der Waals loop) is not physically realizable; it signals a two-phase region.
- **(7.2-15) Conditions for phase equilibrium (all constraints)** — `T^{\mathrm{I}} = T^{\mathrm{II}}, \qquad P^{\mathrm{I}} = P^{\mathrm{II}}, \qquad \underline{G}^{\mathrm{I}} = \underline{G}^{\mathrm{II}}`
  The central result of the chapter (same as Eq. 7.1-9c). Equal temperature and pressure are the "obvious" mechanical/thermal conditions; equality of the molar Gibbs energies is the nontrivial chemical condition that fixes which pressure the coexistence isobar is drawn at. Valid for vapor-liquid, solid-liquid, and solid-vapor equilibrium.
- **(7.3-1b) Lever rule (two-phase molar volume)** — `\omega^{V} = \dfrac{\underline{V} - \underline{V}^{L}}{\underline{V}^{V} - \underline{V}^{L}}`
  In the two-phase region $T$ and $P$ are fixed, so the phase molar volumes $\underline{V}^L,\underline{V}^V$ are fixed and only the vapor fraction $\omega^V$ varies from 0 to 1. Analogous lever rules hold for $\underline{H}$, $\underline{U}$, $\underline{G}$, $\underline{S}$, $\underline{A}$: $\underline{\theta}=\omega^V\underline{\theta}^V+(1-\omega^V)\underline{\theta}^L$.
- **(7.4-6) Definition of fugacity and fugacity coefficient** — `f = P\exp\!\left[\dfrac{\underline{G}(T,P)-\underline{G}^{\text{IG}}(T,P)}{RT}\right] = P\exp\!\left[\dfrac{1}{RT}\int_{0}^{P}\!\left(\underline{V}-\dfrac{RT}{P}\right)dP\right], \qquad \phi = \dfrac{f}{P}`
  DEFINING equation (Eqs. 7.4-6a and 7.4-6b). The fugacity has units of pressure and $f\to P$, $\phi\to 1$ as $P\to 0$ (ideal-gas limit). It is NOT simply a "corrected pressure" -- it is defined through the difference between the real and ideal-gas molar Gibbs energies. $\underline{G}^{\text{IG}}$ is the ideal-gas molar Gibbs energy at the same $T,P$.
- **(7.4-7) Fugacity form of the phase-equilibrium criterion** — `f^{\mathrm{I}}(T,P) = f^{\mathrm{II}}(T,P) \qquad \text{equivalently} \qquad \phi^{\mathrm{I}}(T,P) = \phi^{\mathrm{II}}(T,P)`
  Follows directly from $\underline{G}^{\mathrm{I}}=\underline{G}^{\mathrm{II}}$ (Eq. 7.2-15c) plus the definition of fugacity. For vapor-liquid equilibrium, $f^{L}=f^{V}$. This is the criterion actually used in EOS-based phase-equilibrium calculations.
- **(7.4-8) Fugacity coefficient from a volume-explicit equation of state** — `\ln\dfrac{f}{P} = \ln\phi = \dfrac{1}{RT}\int_{\underline{V}=\infty}^{\underline{V}}\!\left(\dfrac{RT}{\underline{V}} - P\right)d\underline{V} \; - \ln Z + (Z-1)`
  The working form for cubic/virial EOS, which are pressure-explicit in $\underline{V}$. Here $Z=P\underline{V}/RT$ (Ch.6). Use the vapor (large-$\underline{V}$) root for a gas, the liquid (small-$\underline{V}$) root for a liquid. Closed-form results follow for the virial (Eq. 7.4-12), van der Waals (Eq. 7.4-13), and Peng-Robinson (Eq. 7.4-14) equations.
- **(7.4-18) Fugacity of a pure liquid (Poynting correction)** — `f^{L}(T,P) = f_{\text{sat}}(T)\,\exp\!\left[\dfrac{1}{RT}\int_{P^{\text{vap}}(T)}^{P}\underline{V}^{L}\,dP\right], \qquad f_{\text{sat}}(T) = P^{\text{vap}}(T)\left(\dfrac{f}{P}\right)_{\text{sat},T}`
  The exponential is the Poynting pressure correction: it accounts for the system pressure exceeding the vapor pressure. Because $\underline{V}^{L}$ is small, it matters only at high $P$ (or cryogenic $T$). Common approximations: $f^{L}\approx P^{\text{vap}}(T)$ (Eq. 7.4-19, low pressure); $f^{L}\approx f_{\text{sat}}(T)$ (Eq. 7.4-20); or the incompressible-liquid form $f^{L}=P^{\text{vap}}(f/P)_{\text{sat}}\exp[\underline{V}^{L}(P-P^{\text{vap}})/RT]$ (Eq. 7.4-21).
- **(7.6-1) Gibbs phase rule for a one-component system** — `\mathcal{F} = 2\mathcal{P} - 3(\mathcal{P}-1) = 3 - \mathcal{P}`
  CAUTION: $\mathcal{P}$ here is the NUMBER OF PHASES, not pressure. Each of the $\mathcal{P}$ phases needs 2 state variables ($2\mathcal{P}$ total); equilibrium imposes $\mathcal{P}-1$ equalities each of $T$, $P$, and $\underline{G}$ ($3(\mathcal{P}-1)$ total). So single-phase $\Rightarrow \mathcal{F}=2$, two-phase $\Rightarrow \mathcal{F}=1$, triple point (three phases) $\Rightarrow \mathcal{F}=0$. The degrees of freedom must be INTENSIVE per-phase variables, not overall two-phase molar properties.
- **(7.7-4) Clapeyron equation** — `\left(\dfrac{\partial P^{\text{sat}}}{\partial T}\right)_{\text{coex}} = \dfrac{\Delta \underline{S}}{\Delta \underline{V}} = \dfrac{\Delta \underline{H}}{T\,\Delta \underline{V}}`
  Exact for ANY phase transition (solid-liquid, liquid-vapor, solid-vapor). $\Delta\theta=\theta^{\mathrm{I}}-\theta^{\mathrm{II}}$. Relates the coexistence-curve slope to the enthalpy and volume changes of the transition. Water is the classic exception: $\Delta_{\text{fus}}\underline{H}>0$ but $\Delta_{\text{fus}}\underline{V}<0$, so its ice-water line has negative slope.
- **(7.7-5a) Clausius-Clapeyron equation** — `\dfrac{dP^{\text{vap}}}{dT} = \dfrac{P^{\text{vap}}\,\Delta_{\text{vap}}\underline{H}}{RT^{2}} \qquad \text{or} \qquad \dfrac{d\ln P^{\text{vap}}}{dT} = \dfrac{\Delta_{\text{vap}}\underline{H}}{RT^{2}}`
  Clapeyron specialized to vapor-liquid (or vapor-solid) equilibrium using $\underline{V}^V\gg\underline{V}^L$ and an ideal vapor ($\Delta_{\text{vap}}\underline{V}\approx\underline{V}^V=RT/P$). Taking $\Delta_{\text{vap}}\underline{H}$ constant and integrating gives $\ln[P^{\text{vap}}(T_2)/P^{\text{vap}}(T_1)] = -(\Delta_{\text{vap}}\underline{H}/R)(1/T_2 - 1/T_1)$ (Eq. 7.7-6): $\ln P^{\text{vap}}$ is linear in $1/T$.
- **(8.1-12) Definition of a partial molar property** — `\bar{\theta}_i = \bar{\theta}_i(T,P,\underline{x}) = \left(\dfrac{\partial (N\theta)}{\partial N_i}\right)_{T,P,N_{j\neq i}}`
  The defining equation of the chapter. $\theta$ is any molar mixture property and $N\theta$ the corresponding total property; the derivative is taken at constant $T$, $P$, and all other mole numbers $N_{j\neq i}$. The result $\bar{\theta}_i$ is the partial molar property of species $i$ — a true mixture property, generally NOT the pure-component molar value $\underline{\theta}_i$. Derivatives in which $T$ and/or $P$ are NOT held fixed are not partial molar quantities.
- **(8.1-13) Mixture property from partial molar properties** — `\theta = \sum_{i=1}^{\mathcal{C}} x_i\,\bar{\theta}_i(T,P,\underline{x})`
  Any molar mixture property is the mole-fraction-weighted sum of the partial molar properties. Looks like the naive rule $\theta = \sum x_i \underline{\theta}_i$ (Eq. 8.1-1) but uses the partial molar $\bar{\theta}_i$, not the pure-component molar value; the two agree only for an ideal mixture. Applies term by term to $V,H,U,S,G,A$ (Table 8.1-1).
- **(8.1-14) Volume change on mixing (definition)** — `\Delta_{\text{mix}}V = V(T,P,\underline{x}) - \sum_{i=1}^{\mathcal{C}} N_i\,\underline{V}_i(T,P) = \sum_{i=1}^{\mathcal{C}} N_i\left[\bar{V}_i(T,P,\underline{x}) - \underline{V}_i(T,P)\right]`
  The isothermal volume change on forming the mixture from pure components. The heat of mixing is the identical construction for enthalpy, $\Delta_{\text{mix}}H = \sum_i N_i[\bar{H}_i - \underline{H}_i]$ (Eq. 8.1-15). Note the explicit contrast between the partial molar $\bar{V}_i$ (overbar) and the pure-component molar $\underline{V}_i$ (underbar).
- **(8.2-1) Total differential of $G$; chemical potential** — `dG = -S\,dT + V\,dP + \sum_{i=1}^{\mathcal{C}} \bar{G}_i\,dN_i, \qquad \mu_i \equiv \bar{G}_i = \left(\dfrac{\partial G}{\partial N_i}\right)_{T,P,N_{j\neq i}}`
  Extends the pure-fluid $dG = -S\,dT + V\,dP$ (Eq. 6.2-12) with a composition term. The partial molar Gibbs energy $\bar{G}_i$ is the chemical potential $\mu_i$. Parallel forms: $dU = T\,dS - P\,dV + \sum \bar{G}_i\,dN_i$ (8.2-4) and $dA = -S\,dT - P\,dV + \sum \bar{G}_i\,dN_i$ (8.2-5); and $\bar{G}_i = \bar{H}_i - T\bar{S}_i$ (from Eq. 8.1-16b).
- **(8.2-9b) Gibbs-Duhem equation (constant $T$, $P$)** — `\sum_{i=1}^{\mathcal{C}} x_i\,d\bar{\theta}_i \Big|_{T,P} = 0`
  A consistency relation among the partial molar properties: at fixed $T$ and $P$ their composition-weighted differentials sum to zero, so the $\bar{\theta}_i$ cannot all vary independently. The general form (variable $T,P$) is $-(\partial\theta/\partial T)_{P,N}\,dT - (\partial\theta/\partial P)_{T,N}\,dP + \sum_i x_i\,d\bar{\theta}_i = 0$ (Eq. 8.2-8b). Central to reducing the data needed to characterize a mixture and to testing data consistency.
- **(8.2-13b) Gibbs-Duhem equation for the Gibbs energy (constant $T$, $P$)** — `\sum_{i=1}^{\mathcal{C}} x_i\,d\bar{G}_i \Big|_{T,P} = 0`
  The $\theta = G$ specialization (equivalently $\sum_i N_i\,d\bar{G}_i = 0$, Eq. 8.2-13a; the mole-number form is 8.2-12a,b). For a binary mixture at constant $T,P$ it gives $x_1(\partial\bar{G}_1/\partial x_1) + x_2(\partial\bar{G}_2/\partial x_1) = 0$ (Eq. 8.2-20) — used repeatedly in the phase-equilibrium chapters.
- **(8.3-5) Mole numbers with multiple reactions** — `N_i = N_{i,0} + \sum_{j=1}^{\mathcal{M}} \nu_{ij}\,X_j`
  Species mole numbers in a closed reacting system, summed over the $\mathcal{M}$ independent reactions ($N_{i,0}$ = initial moles). Reduces to the single-reaction $N_i = N_{i,0} + \nu_i X$ (Eq. 8.3-2a) when $\mathcal{M}=1$. Printed p.400 mislabels the reappearance of this relation "(8.35)".
- **(8.4-3) Energy balance for a reacting mixture (molar basis)** — `\dfrac{dU}{dt} = \sum_{k=1}^{K} \left(\dot{N}\,\underline{H}\right)_k + \dot{Q} + \dot{W}_s - P\dfrac{dV}{dt}, \qquad \left(\dot{N}\,\underline{H}\right)_k = \sum_{i=1}^{\mathcal{C}} \left(\dot{N}_i\,\bar{H}_i\right)_k`
  Same form as the pure-fluid energy balance (Ch.3), but each flow stream is a mixture, so its enthalpy flow is built from the PARTIAL MOLAR enthalpies $\bar{H}_i$ of its species (Eq. 8.4-4). The heat of reaction is contained implicitly through the differing inlet/outlet mole numbers. Sign convention (Ch.3): $\dot{Q}$, $\dot{W}_s$ positive when done ON the system.
- **(8.5-3) Standard heat of reaction from heats of formation** — `\Delta_{\text{rxn}}H^{\circ}(T) = \sum_{i=1}^{\mathcal{C}} \nu_i\,\Delta_{\text{f}}\underline{H}_i^{\circ}(T)`
  The reference-state atomic-species enthalpies cancel by conservation of atoms, leaving a sum over formation enthalpies weighted by stoichiometric coefficients (products positive, reactants negative). The Gibbs-energy analog is $\Delta_{\text{rxn}}G^{\circ}(T) = \sum_i \nu_i\,\Delta_{\text{f}}\underline{G}_i^{\circ}(T)$ (Eq. 8.5-4). Temperature dependence via $\Delta_{\text{rxn}}H^{\circ}(T) = \Delta_{\text{rxn}}H^{\circ}(25^\circ\mathrm{C}) + \sum_i \nu_i\int C_{P,i}^{\circ}\,dT$ (Eq. 8.5-5).
- **(8.7-7) Criteria for phase equilibrium in a multicomponent system** — `T^{\mathrm{I}} = T^{\mathrm{II}}, \qquad P^{\mathrm{I}} = P^{\mathrm{II}}, \qquad \bar{G}_i^{\mathrm{I}} = \bar{G}_i^{\mathrm{II}} \;\;\left(\mu_i^{\mathrm{I}} = \mu_i^{\mathrm{II}}\right) \quad \text{for each species } i`
  The central result: thermal, mechanical, and species (diffusive) equilibrium. It is the partial molar Gibbs energy / chemical potential of each species — NOT its concentration — that must be equal across phases. Replacing $\bar{G}_i$ by the pure-component molar $G$ recovers the one-component criteria of Sec. 7.1. Generalizes to any number of phases: $\bar{G}_i^{\mathrm{I}} = \bar{G}_i^{\mathrm{II}} = \bar{G}_i^{\mathrm{III}} = \cdots$ (Eq. 8.7-9). In Ch.9 this is recast in the working form $\bar{f}_i^{\mathrm{I}} = \bar{f}_i^{\mathrm{II}}$.
- **(8.8-4) Criteria for chemical equilibrium** — `\sum_{i=1}^{\mathcal{C}} \nu_{ij}\,\bar{G}_i = 0 \qquad j = 1, 2, \dots, \mathcal{M}`
  At constant $T,P$ the Gibbs energy is a minimum; each independent reaction must itself be at equilibrium. Single-reaction case: $\sum_i \nu_i \bar{G}_i = 0$ (Eq. 8.8-1). For combined phase + chemical equilibrium, both Eq. 8.7-9 (phase) and this criterion (in each phase) must hold simultaneously (Eqs. 8.8-15, 8.8-18).
- **(8.9-6) Gibbs phase rule** — `\mathcal{F} = \mathcal{C} - \mathcal{M} - \mathcal{P} + 2`
  Degrees of freedom = number of intensive state variables one may independently fix. Counts unknowns $\mathcal{P}(\mathcal{C}+1)$ minus the interrelations $2(\mathcal{P}-1) + \mathcal{C}(\mathcal{P}-1) + \mathcal{M}$. For a nonreacting pure fluid ($\mathcal{C}=1,\mathcal{M}=0$) it reduces to $\mathcal{F} = 3 - \mathcal{P}$ (Sec. 7.6). Tells you whether an equilibrium problem is well posed; it does not solve it.
- **(9.3-12) Excess Gibbs energy and the activity coefficients** — `\dfrac{\underline{G}^{\text{ex}}}{RT} = \sum_i x_i \ln\gamma_i`
  THE CHAPTER’S SELF-CHECK. Any activity-coefficient model that supplies both an expression for $\underline{G}^{\text{ex}}$ and expressions for the $\gamma_i$ must satisfy this identity EXACTLY — it follows from the definition of a partial molar property (Ch.8) applied to $\underline{G}^{\text{ex}}$, and is independent of the model. Running it numerically on the book’s eight models found two printed equations that fail it (see 9.5-12b and 9.5-18 below). Use it whenever you are asked to produce or verify activity coefficients.
- **(9.4-8) van der Waals one-fluid mixing rules** — `a_{\text{mix}} = \sum_i\sum_j x_i x_j a_{ij}, \qquad b_{\text{mix}} = \sum_i x_i b_i`
  How a cubic equation of state is extended to a mixture: the mixture $a$ is a double sum over species pairs, the mixture $b$ a simple mole-fraction average. Here $a$ and $b$ are EQUATION-OF-STATE parameters — nothing to do with the UNIFAC $a_{mn}$ (see convention species_vs_group_indices).
- **(9.4-9) Combining rule and the binary interaction parameter** — `a_{ij} = \sqrt{a_i a_j}\,(1 - k_{ij})`
  The geometric mean, corrected by a fitted binary interaction parameter $k_{ij}$. $k_{ij}$ is small (order 0.1 or less) but the results are sensitive to it. A missing $k_{ij}$ is NOT safely assumed zero — say so rather than defaulting silently; the book tabulates values in Table 9.4-1.
- **(9.5-9) van Laar activity coefficients** — `\ln\gamma_1 = \dfrac{\alpha}{\left(1+\dfrac{\alpha x_1}{\beta x_2}\right)^{2}}, \qquad \ln\gamma_2 = \dfrac{\beta}{\left(1+\dfrac{\beta x_2}{\alpha x_1}\right)^{2}}`
  A two-parameter correlative model, derived from the Wohl expansion with $\alpha = 2q_1a_{12}$ and $\beta = 2q_2a_{12}$. This $\alpha$ is NOT the NRTL nonrandomness parameter $\alpha$, which appears four pages later in the same chapter.
- **(9.5-12a) Wilson activity coefficients (binary)** — `\ln\gamma_1 = -\ln(x_1+\Lambda_{12}x_2) + x_2\left[\dfrac{\Lambda_{12}}{x_1+\Lambda_{12}x_2} - \dfrac{\Lambda_{21}}{x_2+\Lambda_{21}x_1}\right]`
  The Wilson model, with $\Lambda_{ij}\neq\Lambda_{ji}$ and $\Lambda_{ii}=1$. It cannot predict liquid-liquid phase splitting, which is why Ch.11 does not use it.
- **(9.5-12b) Wilson infinite-dilution activity coefficients (CORRECTED)** — `\ln\gamma_1^{\infty} = -\ln\Lambda_{12} + 1 - \Lambda_{21}, \qquad \ln\gamma_2^{\infty} = -\ln\Lambda_{21} + 1 - \Lambda_{12}`
  CORRECTED ERRATUM — the 5e prints the SAME subscript twice in each expression ($\ldots+1-\Lambda_{12}$ in the first). The form above is the $x_1\to0$ and $x_2\to0$ limit of Eq. 9.5-12a. This is not cosmetic: for a Wilson fit to the chapter’s own benzene / 2,2,4-trimethyl pentane data the printed form returns $\gamma_1^{\infty}=0.61$ where the model itself gives $1.48$ — the SIGN of the deviation from ideality is wrong. If you have seen the printed form, discard it.
- **(9.5-14b) NRTL infinite-dilution activity coefficients** — `\ln\gamma_1^{\infty} = \tau_{21} + \tau_{12}e^{-\alpha\tau_{12}}`
  Correct as printed — verified as the limit of the full NRTL expression. Included here as the counterexample to 9.5-12b: same kind of equation, same chapter, no error.
- **(9.5-18) Flory-Huggins activity coefficient of the larger species (CORRECTED)** — `\ln\gamma_2 = \ln\dfrac{\phi_2}{x_2} - (m-1)\phi_1 + m\chi\phi_1^{2}`
  CORRECTED ERRATUM — the 5e prints $+(m-1)\phi_1$. The sign above follows from differentiating the book’s own $\underline{G}^{\text{ex}}$ (Eq. 9.5-17). Its partner $\ln\gamma_1 = \ln(\phi_1/x_1)+(1-1/m)\phi_2+\chi\phi_2^{2}$ IS correct, so only one sign moves. The discrepancy reaches $1.75$ in $\ln\gamma$ at $m=3$, $\chi=0.5$, $x_1=0.7$. Here $\phi_i$ is a VOLUME FRACTION and $m$ the size ratio $v_2/v_1$.
- **(9.5-23b) UNIQUAC residual activity coefficient** — `\ln\gamma_i^{\text{residual}} = q_i\left[1 - \ln\!\left(\sum_j \theta_j\tau_{ji}\right) - \sum_j \dfrac{\theta_j\tau_{ij}}{\sum_k \theta_k\tau_{kj}}\right]`
  SUBSCRIPT ORDER IS LOAD-BEARING: $\tau_{ji}$ inside the logarithm, $\tau_{ij}$ in the sum. $\tau$ is not symmetric, so swapping them changes every computed $\gamma_i$ — and the pure-component limit does NOT catch the error ($\theta_i=1$, $\tau_{ii}=1$ gives $\ln\gamma_i=0$ either way). It only shows up in a mixture. Matches Abrams and Prausnitz (1975).
- **(9.6-12a) Modified (Dortmund) UNIFAC combinatorial term** — `\phi_i^{\prime} = \dfrac{x_i r_i^{3/4}}{\sum_j x_j r_j^{3/4}}`
  THE $3/4$ EXPONENT IS THE DORTMUND MODIFICATION and is absent from original UNIFAC. Its presence, together with the Dortmund $R$/$Q$ values in Table 9.5-2, is what settles that this book teaches MODIFIED UNIFAC. Never mix Dortmund $R$/$Q$ with original-UNIFAC interaction parameters or vice versa.
- **(9.6-13b) Modified UNIFAC group-interaction parameter** — `\Psi_{mn} = \exp\left[-\dfrac{a_{mn}+b_{mn}T+c_{mn}T^{2}}{T}\right]`
  The temperature-dependent three-parameter form, equivalently $\exp[-(u_{mn}-u_{nn})/kT]$ with $(u_{mn}-u_{nn})/k = a_{mn}+b_{mn}T+c_{mn}T^2$. $a_{mn}\neq a_{nm}$, so $\Psi_{mn}\neq\Psi_{nm}$ — see convention interaction_parameter_asymmetry. Original UNIFAC uses the one-parameter $\exp(-a_{mn}/T)$; if you are given a single $a_{mn}$ per pair, you have original-UNIFAC parameters and they do not belong with Table 9.5-2’s $R$ and $Q$.

## Codes (Jupyter notebooks / Python scripts)

Distributed from https://github.com/emfurst/cbet6e. When you open one in this editor, run and
modify the cells directly rather than narrating code you cannot execute.

- `ch1/Math_review.ipynb` — Mathematics you will use in this book
- `ch1/Getting_started_Python_Jupyter.ipynb` — Getting started with Python and Jupyter
- `ch1/vdW_molar_volume_CO2.ipynb` — When is carbon dioxide an ideal gas?
- `ch1/Thermometer_water_mercury.ipynb` — A mercury thermometer, and why not a water one
- `ch1/Heat_capacity_fitting_CO2.ipynb` — Curve fitting a heat-capacity polynomial
- `ch1/Heat_capacity_integration_N2.ipynb` — Integrating a heat capacity: how many intervals is enough?
- `ch3/Heat_capacity_range_of_validity.ipynb` — How far can a heat-capacity correlation be trusted?
- `ch3/refrigerant_comparison.ipynb` — Can a new refrigerant replace an old one?
- `ch3/Steam_charts_from_appendix_A3.ipynb` — Figure 3.3-1: the steam charts, drawn from Appendix A.III
- `ch3/PH_charts_methane_and_nitrogen.ipynb` — Figures 3.3-2 and 3.3-3: methane and nitrogen from the Peng-Robinson equation
- `ch3/Heat_capacity_from_molecular_structure.ipynb` — Where a heat capacity comes from, and what fitting one costs
- `ch3/LJ_interaction_energy_figure.ipynb` — Figure 3.3-5: the interaction energy between two molecules
- `ch5/Linde_liquefaction_CH4_figure.ipynb` — Figure 5.1-3: the liquefaction path of Illustration 5.1-1
- `ch5/Rankine_cycle_steam_figure.ipynb` — Illustration 5.2-1: the Rankine cycle on the steam T-S chart
- `ch6/vdW_isotherms_example.ipynb` — Figure 6.6-1: pressure–volume behavior of the van der Waals EOS
- `ch6/PR_eos_reference.ipynb` — Generalized Peng-Robinson equation of state
- `ch6/PR_isotherms_N2_example.ipynb` — Plotting PR isotherms (nitrogen)
  - `ch6/PR_isotherms_N2_example_thermo.ipynb` — thermo-package version
- `ch6/PR_isotherms_O2_example.ipynb` — Figure 6.4-3: pressure–volume diagram for oxygen
  - `ch6/PR_isotherms_O2_example_thermo.ipynb` — thermo-package version
- `ch6/PR_properties_table_O2_example.ipynb` — Table 6.4-4: thermodynamic properties of oxygen
- `ch6/PR_enthalpy_O2_example.ipynb` — Figure 6.4-4: pressure–enthalpy diagram for oxygen
  - `ch6/PR_enthalpy_O2_example_thermo.ipynb` — thermo-package version
- `ch6/PR_entropy_O2_example.ipynb` — Figure 6.4-5: temperature–entropy diagram for oxygen
  - `ch6/PR_entropy_O2_example_thermo.ipynb` — thermo-package version
- `ch6/PR_discharge_N2_example.ipynb` — Illustration 6.7-1: discharging a nitrogen cylinder with the Peng-Robinson EOS
- `ch6/Helmholtz_fundamental_eos_O2_example.ipynb` — A fundamental equation of state: every property of oxygen from one function
- `ch6/PR_throttle_CH4_example.ipynb` — Departure functions and a throttling (Joule-Thomson) calculation for methane
  - `ch6/PR_throttle_CH4_example_thermo.ipynb` — thermo-package version
- `ch6/PR_throttle_C2H6_homework.ipynb` — Homework: outlet temperature of throttled ethane
  - `ch6/PR_throttle_C2H6_homework_thermo.ipynb` — thermo-package version
- `ch6/PR_heat_capacity_C2H6_homework.ipynb` — Homework: pressure dependence of the heat capacity of ethane
  - `ch6/PR_heat_capacity_C2H6_homework_thermo.ipynb` — thermo-package version
- `ch7/van_der_waals_EOS.ipynb` — van der Waals equation of state — Figures 7.3-1 to 7.3-4
- `ch7/mechanical_stability_CO2.ipynb` — Mechanical stability and the coexistence region — carbon dioxide
- `ch7/CO2_phases.ipynb` — Phase diagram of carbon dioxide — the quantitative Figure 7.3-6
- `ch7/N2_phases.ipynb` — Phase diagram of nitrogen — the same calculation, different choices
- `ch7/vapor_pressure_n_butane.ipynb` — Vapor pressure of n-butane — Figure 7.5-2
- `ch7/vapor_pressure_oxygen.ipynb` — Vapor pressure of oxygen — Figure 7.5-3 and Table 7.5-1
- `ch7/properties_chart_oxygen.ipynb` — Completing the oxygen properties chart — Table 7.5-2
- `ch7/clausius_clapeyron_isooctane.ipynb` — The Clausius-Clapeyron equation — Illustration 7.7-1 and Figure 7.7-1
- `ch7/vapor_pressure_water.ipynb` — Vapor pressure of water — PR vs. PRSV, Illustration 7.5-3
- `ch7/surface_effects_droplets.ipynb` — Surface effects on small drops — Figure 7.8-1
- `ch8/RK_partial_molar_volume_water_methanol_example.ipynb` — Table 8.6-2: partial molar volumes from density data
- `ch8/RK_partial_molar_enthalpy_water_methanol_example.ipynb` — Table 8.6-4: partial molar enthalpies from heat-of-mixing data
- `ch8/Hrxn_temperature_N2O4_example.ipynb` — Illustration 8.5-2: the standard heat of reaction versus temperature
- `ch9/PR_species_fugacity_ethane_butane_example.ipynb` — Species fugacity in a mixture from the Peng-Robinson EOS (ethane / n-butane)
- `ch9/Lewis_Randall_vs_PR_CO2_methane_example.ipynb` — Lewis-Randall rule against the Peng-Robinson EOS (CO2 / methane)
- `ch9/activity_coefficient_correlation_benzene_TMP_example.ipynb` — Correlating activity coefficients (benzene / 2,2,4-trimethyl pentane)
- `ch9/regular_solution_benzene_TMP_example.ipynb` — Regular solution theory (benzene / 2,2,4-trimethyl pentane)
- `ch9/UNIFAC_benzene_TMP_example.ipynb` — Modified (Dortmund) UNIFAC (benzene / 2,2,4-trimethyl pentane)
- `ch9/Margules_activity_coefficients_figure.ipynb` — The one-constant Margules activity coefficients
- `ch9/Gex_composition_curves_figure.ipynb` — Excess Gibbs energy against composition
- `ch9/Debye_Huckel_HCl_NaCl_example.ipynb` — Electrolyte solutions — Debye-Hückel and its extensions (HCl, NaCl)
- `ch9/Wong_Sandler_mixing_rule_acetone_water_example.ipynb` — The Wong-Sandler mixing rule — an activity coefficient model inside an equation of state (acetone / water)
- `ch9/hemoglobin_activity_coefficient_example.ipynb` — Activity coefficient of a protein on a Henry's law basis (hemoglobin in water)
- `ch10/raoult_diagrams_pentane_heptane_example.ipynb` — Vapor-liquid equilibrium diagrams for an ideal mixture (n-pentane / n-heptane)
- `ch10/hexane_triethylamine_diagrams_figure.ipynb` — Phase diagrams for a nearly ideal mixture (hexane / triethylamine)
- `ch10/vle_correlation_benzene_TMP_example.ipynb` — Correlating vapor-liquid equilibrium data (benzene / 2,2,4-trimethyl pentane)
- `ch10/excess_properties_benzene_TMP_figure.ipynb` — Excess enthalpy and entropy from the temperature dependence of $G^{ex}$ (benzene / 2,2,4-trimethyl pentane)
- `ch10/azeotrope_to_vanlaar_benzene_cyclohexane_example.ipynb` — A phase diagram from one azeotropic point (benzene / cyclohexane)
- `ch10/pressure_swing_methyl_acetate_methanol_example.ipynb` — Nonideal phase diagrams, an azeotrope, and how to separate it (methyl acetate / methanol)
- `ch10/pr_vle_co2_isopentane_figure.ipynb` — One binary parameter, and what it is worth (carbon dioxide / isopentane)
- `ch10/acetone_water_mixing_rules_figure.ipynb` — Four ways to model a very nonideal mixture (acetone / water)
- `ch10/azeotrope_test_ethyl_acetate_benzene_homework.ipynb` — Does this system have an azeotrope? Two models that disagree (ethyl acetate / benzene)
- `ch10/unifac_pentane_propionaldehyde_example.ipynb` — Three ways to get a phase diagram, and how little data each needs (n-pentane / propionaldehyde)
- `ch10/flory_huggins_benzene_PIB_example.ipynb` — A solvent above a polymer, where mole fraction stops being useful (benzene / polyisobutylene)
- `ch11/gas_solubility_co2_toluene_example.ipynb` — Solubility of a gas in a liquid from an equation of state (CO2 in toluene)
- `ch11/lle_isobutane_furfural_example.ipynb` — Liquid-liquid equilibrium from an activity coefficient model (isobutane / furfural)
- `ch11/polymer_compatibility_PS_PMMA_example.ipynb` — Are two molten polymers compatible? (polystyrene / PMMA)
- `ch11/vlle_eos_co2_decane_example.ipynb` — Three-phase equilibrium from an equation of state (CO2 / n-decane)
- `ch11/lle_eos_mixing_rules_co2_decane_example.ipynb` — Two mixing rules against one measured phase diagram (CO2 / n-decane)
- `ch11/fluorocarbon_solubility_figure.ipynb` — Consolute temperatures of fluorocarbon + hydrocarbon mixtures (Figure 11.2-3)
- `ch11/margules_azeotrope_vs_phase_split_homework.ipynb` — When does an azeotrope become a phase split? (Problem 11.2-17)
- `ch11/vlle_isobutane_furfural_figure.ipynb` — When should you suspect a second liquid phase? (Figures 11.3-1 and 11.3-2)
- `ch11/ternary_extraction_MIK_acetone_water_example.ipynb` — Liquid-liquid extraction on a triangular diagram (MIK / acetone / water)
- `ch11/staged_extraction_MIK_acetone_water_example.ipynb` — Two-stage liquid-liquid extraction (MIK / acetone / water)
- `ch11/ternary_mass_balance_example.ipynb` — Reading a triangular diagram (a mass balance on three components)
- `ch11/distribution_coefficient_bromine_CCl4_water_example.ipynb` — An activity coefficient measured with a separatory funnel (bromine / CCl4 / water)
- `ch11/octanol_water_partition_coefficient_example.ipynb` — The octanol-water partition coefficient: estimating it and using it
- `ch11/partition_coefficient_correlation_pollutants_example.ipynb` — Testing a correlation against 160 printed numbers (Table 11.4-1)
- `ch11/gibbs_energy_of_transfer_amino_acids_example.ipynb` — The hydrophobic effect, from six solubilities (amino acids, water to ethanol)
- `ch11/osmometry_molecular_weight_PVC_albumin_example.ipynb` — A molecular weight from a column of liquid (PVC, and serum albumin)
- `ch11/osmotic_virial_chymotrypsin_example.ipynb` — A molecular weight and a second virial coefficient from one osmometer (alpha-chymotrypsin)
- `ch11/osmotic_pressure_blood_example.ipynb` — Why saline is 0.9 % (the osmotic pressure of blood)
- `ch11/reverse_osmosis_seawater_example.ipynb` — Pushing water the other way (reverse osmosis of seawater)
