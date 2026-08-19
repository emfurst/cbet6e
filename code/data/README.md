# Property data

Reference data extracted from the legacy 5e companion databases (Microsoft Access
`.mdb`, read by the old Visual Basic programs) and exported to CSV so the numbers
survive independent of VB6, Access, MATLAB, and Mathcad. One-time export with
`mdb-tools`; the source `.mdb` files are the 5e companion archives.

## Source and credit

> The constants of pure fluids in `pure_property.csv` are **adapted from R. C. Reid,
> J. M. Prausnitz, and B. E. Poling, *The Properties of Gases and Liquids*, 4th ed.,
> McGraw-Hill, New York, 1986, Appendix A, with corrections and data from other
> sources.**

This is the same credit line carried in the printed Appendix B of earlier editions,
where it described the Visual Basic `PROPERTY` program that shipped this database; it
appears in the 6e Appendix B and is reproduced here so the attribution travels with the
data files themselves. The printed Appendix A tables — A.II (ideal-gas heat capacities)
and A.IV (enthalpies and Gibbs energies of formation) — are drawn from this database.

**License note.** The `thermo` package and the notebooks in this repository are released
under the repository's own license; **that license does not extend to the underlying
property data**, whose provenance is above. Keep this credit with the CSVs in any copy
or redistribution.

## `pure_property.csv` — 618 pure compounds

Pure-component constants. The `Cp*` columns are the **Appendix A.II ideal-gas
heat-capacity polynomial** coefficients.

| Column | Meaning | Units |
|---|---|---|
| `No` | index | — |
| `Formula`, `Name` | chemical formula, compound name | — |
| `Molwt` | molecular weight | g/mol |
| `Tfp`, `Tb` | normal freezing / boiling point | K |
| `Tc`, `Pc`, `Vc`, `Zc` | critical temperature / pressure / volume / compressibility | K, bar, cm³/mol, — |
| `Omega` | Pitzer acentric factor | — |
| `Dipm` | dipole moment | debye |
| `CpA, CpB, CpC, CpD` | ideal-gas heat capacity, $C_P^* = A + B\,T + C\,T^2 + D\,T^3$ | J/(mol·K) |
| `dHf`, `dGf` | ideal-gas standard enthalpy / Gibbs energy of formation at 298.2 K | J/mol |
| `Eq` | vapor-pressure equation selector (1/2/3, see below) | — |
| `VpA, VpB, VpC, VpD` | vapor-pressure constants | — |
| `Tmin`, `Tmax` | valid T-range for the vapor-pressure fit | K |
| `Lden`, `Tden` | liquid density at reference temperature | g/cm³, K |

Vapor-pressure correlations, selected by `Eq` (with $t = 1 - T/T_c$):
- `Eq=1` (Wagner): $\ln(P^{\text{vap}}/P_c) = (A t + B t^{1.5} + C t^3 + D t^6)/(1-t)$
- `Eq=2` (Riedel): $\ln P^{\text{vap}} = A - B/T + C\ln T + D\,P^{\text{vap}}/T^2$ (iterative)
- `Eq=3` (Antoine): $\ln P^{\text{vap}} = A - B/(T + C)$

## `react_property.csv` — 99 reaction species

Appendix A.IV (standard formation properties at 25 °C, 1 bar) joined to Appendix A.II
(ideal-gas heat capacities). For chemical-equilibrium and heat-of-reaction work.

| Column | Meaning | Units |
|---|---|---|
| `Name` | species, as a formula (`N2O4`, `H2O(g)`, `C6H6`) | — |
| `DG`, `DH` | Gibbs energy / enthalpy of formation | **kJ/mol** |
| `A, B, C, D, E` | heat-capacity coefficients, **as Appendix A.II prints them** | see below |
| `ID` | index | — |

> **The heat-capacity columns are stored scaled, and evaluating them raw is wrong by
> three orders of magnitude — silently.** Appendix A.II prints $b$, $c$ and $d$ multiplied
> by $10^{2}$, $10^{5}$ and $10^{9}$, and the CSV keeps that printed form. The formula is
>
> $$C_P^{*} = A + B\cdot10^{-2}\,T + C\cdot10^{-5}\,T^{2} + D\cdot10^{-9}\,T^{3} + E/T^{2}
> \qquad \mathrm{J/(mol\,K)}$$
>
> `E` is unscaled and is nonzero only for the fourteen solid species. **Use
> `thermo.data.reaction_cp(name)`**, which returns the coefficients already multiplied by
> `REACT_CP_SCALE`; `formation_enthalpy` and `formation_gibbs` likewise return **J/mol**,
> not the CSV's kJ/mol.
>
> Verified against the book's own arithmetic (2026-08-09): Illustration 8.5-2 prints the
> combined coefficients for $2\,\mathrm{NO_2} - \mathrm{N_2O_4}$ as $12.804$,
> $-7.239\times10^{-2}$, $4.301\times10^{-5}$, $1.5732\times10^{-8}$, and the scaled
> columns reproduce all four exactly — and then the printed heats of reaction at 200, 300,
> 400, 500 and 600 K to the last digit. Spot-checked against literature $C_P^{*}(298\ \mathrm{K})$
> for O₂, N₂, CO, CO₂, H₂O, CH₄, NH₃ and SO₂ as well.

## `mixing_*.csv` — property changes on mixing (Sec. 8.6)

The data behind the partial molar volume and enthalpy tables. Species 1 is water and
species 2 is methanol in both; **SI units**, where the book prints volumes multiplied by
$10^{6}$ and enthalpies in kJ/mol.

| file | is | columns |
|---|---|---|
| `mixing_water_methanol_volume.csv` | **Table 8.6-1** — density data at 298.15 K | `x1`, `rho_kg_m3`, `V_m3_mol`, `dmixV_m3_mol` |
| `mixing_water_methanol_enthalpy.csv` | **Table 8.6-3** — heat of mixing at 19.69 °C | `x1`, `Qplus_J_mol_MeOH`, `dmixH_J_mol` |

Load with `thermo.load_mixing_data("water-methanol-volume")` (or `…-enthalpy`).

**Each table carries an identity that checks the transcription**, and both pass at the
level of their printed rounding — which is the only audit available for a hand-typed table:

- volume: $\underline{V} = M_{\text{mix}}/\rho$ with $M_1 = 18.0153$, $M_2 = 32.042$ g/mol,
  and $\Delta_{\text{mix}}\underline{V} = \underline{V} - x_1\underline{V}_1 - x_2\underline{V}_2$
  from the two end rows (max deviation $9\times10^{-11}$ m³/mol);
- enthalpy: $\Delta_{\text{mix}}\underline{H} = (1-x_1)\,Q^{+}$, the relation printed in the
  table's own source note (max deviation 0.5 J/mol against data printed to 1 J/mol).

The `code/ch8/` notebooks run both checks in their opening cells. Note that neither of the
book's *derived* tables (8.6-2, 8.6-4) is stored here — they are **generated** by those
notebooks, per the production-notebook rule.

*Source:* Table 8.6-1 from the 5e Sec. 8.6; Table 8.6-3 from *International Critical
Tables*, Vol. 5, McGraw-Hill, New York, 1929, p. 159, as reprinted there.

## `pr_kij.csv` — Table 9.4-1, Peng-Robinson binary interaction parameters

**127 pairs over 20 species**, digitized from the printed Table 9.4-1 (5e p. 441):
`species_i, species_j, formula_i, formula_j, kij`. Stored as an upper-triangular **edge
list, not a matrix**, because the table is about 65% blank and *a blank is not a zero*.

**The blanks are the point.** Table 9.4-1's own footnote reads: *"Blanks indicate no
data are available from which the k12 could be evaluated. In such case use estimates from
mixtures of similar compounds."* `thermo.data.pr_kij_matrix` therefore returns the pairs
it could not find alongside the matrix, and `PRMixture.from_database(..., kij="table")`
**warns** rather than filling a zero in silently.

*Verified* against the book's own three worked values — $k$ = 0.010 for ethane/*n*-butane
(Illustration 9.4-3), 0.09 for methane/carbon dioxide (Illustration 9.4-4), 0.018 for
*n*-pentane/benzene (Illustration 9.4-5) — plus two structural checks: every entry lands
in the upper triangle, and the parser places exactly as many numbers as each printed row
contains. With these values `PRMixture` reproduces Illustration 9.4-3's fugacity table.

*Source:* Table 9.4-1, from H. Knapp, R. Döring, L. Oellrich, U. Plöcker and
J. M. Prausnitz, *Vapor-Liquid Equilibria for Mixtures of Low-Boiling Substances*,
DECHEMA Chemistry Data Series Vol. VI, Frankfurt/Main (1982), and other sources.

## UNIFAC group-contribution parameters

- **`unifac_subgroups.csv` — 92 subgroups over 46 main-group names**, digitized from the
  book's own **Table 9.5-2** (5e pp. 457–458):
  `subgroup_no, main_group_no, subgroup_name, main_group_name, R, Q, example`.
  These *R*/*Q* are the **modified (Dortmund)** values — the only set the book prints
  (`revision_notes/c09.md` D1).
- `unifac_interactions_modified.csv` — modified (Dortmund) main-group interactions:
  `main_i, main_j, a_ij, b_ij, c_ij` (temperature-dependent, $\Psi=\exp[-(a/T+b+cT)]$),
  1220 nonzero pairs over 56 main groups.
- `unifac_interactions_original.csv` — original UNIFAC main-group interactions:
  `main_i, main_j, a_ij` ($\Psi=\exp(-a/T)$), 838 nonzero pairs over 44 main groups.
  **In use**: `UNIFAC("original")` reads it. Superseded the 2026-08-12 "retained but
  unused" note (c09.md **D1**), which was written when the `original` branch raised.

**One asymmetry, and it is deliberate.** *R* and *Q* and the group inventory come from
the book; the **subgroup and main-group numbers do not**, because the chapter prints no
number column anywhere (which is why four references in Illustration 9.6-2 dangle). The
numbers come from the legacy MATLAB `UNIFAC_data.mat` extraction, and they have to,
because `unifac_interactions_modified.csv` is keyed on the main-group number — the two
files must agree or every $\Psi_{mn}$ is wrong.

**What the digitization found.** Full audit trail:
`code/ch9/validation/unifac_subgroups_table_9.5-2_validation.ipynb`.

| | |
|---|---|
| **the legacy file had duplicate subgroup numbers** | 37, 38 and 39 each appeared twice — the *original*-UNIFAC pyridine subgroups (C5H5N, C5H4N, C5H3N) and the Dortmund ones (AC2H2N, AC2HN, AC2N). `dict(zip(subgroup_no, R))` silently kept whichever row came last. `load_unifac_subgroups` now **refuses** a duplicate. |
| **six subgroups carried silicon names** | legacy main group 42 was `SiH2` with subgroups `SiH2`/`SiH`/`Si`, and 43 was `SiO` — but their *R* and *Q* are *exactly* the book's `cy-CH2`/`cy-CH`/`cy-C` and the three `cy-CH2 O` subgroups. In modified UNIFAC, main groups 42 and 43 **are** the cyclic groups. Recovering them by *R*/*Q* fingerprint is what let the cyclic groups keep the numbers the interaction table already uses. |
| **the book prints `CHO` twice** | the aldehyde and the ether subgroup. The ether is stored as `CH-O`. Filed as an erratum. |
| **two values gained digits** | `AC2N` *Q* = 0.3539 (legacy 0.3530), `cy-CON-CH3` *R* = 3.9819 (legacy 3.9810). The book wins. |
| **`CHCl3` is its own main group** | the book prints chloroform and CCl₃ under one `CCl3` heading, but the parameter set gives chloroform main group 45 and CCl₃ main group 23 — which is why their *R* differ where subgroups of a genuine main group always share *R*. This is the only place the printed heading and the parameter set's main groups disagree, and it accounts exactly for 46 names over 47 numbers. |

*Checked:* the group sums reproduce Illustration 9.5-2 exactly — benzene as 6 ACH gives
$r$ = 2.2578 and $q$ = 2.5926; 2,2,4-trimethyl pentane as 5 CH₃ + CH₂ + CH + C gives
$r$ = 5.0600 and $q$ = 6.3675. `thermo.data.unifac_groups(name)` derives group
assignments from the table's own *Example Assignments* column, so a notebook never
hand-transcribes a subgroup number.

> **That is true of the 6e, and it is worth knowing why it has to be said.** The **5e's**
> version of the same illustration printed 3.1878 / 2.4000 and 5.8463 / 5.0080 — **original**
> UNIFAC sums — while citing Table 9.5-2, which by then held the Dortmund set. The 6e
> recomputed it. Those 5e values are not a defect to fix; they are the **source** of
> `unifac_subgroups_original.csv` below.

## `unifac_subgroups_original.csv` — 85 subgroups over 44 main groups

The **original**-UNIFAC *R*/*Q* set, digitized from **Table 7.5-2 of the 2nd edition of this
book** (*Chemical and Engineering Thermodynamics*, pp. 333–334) — the parameter set the book
itself published when it taught original UNIFAC.
`research/references/ch10/SIS-CET-2e-table-7.5-2.pdf`, supplied by the author 2026-08-17.

Read from a **600-dpi render, not from `pdftotext`**: the OCR mangles the numerals
(`1..4457`, `).6908`, `.9031`). Rebuild and re-check with
`python3 tools/build_unifac_subgroups_original.py`, which carries the transcription and its
tests together.

**Why this file could not simply be extracted like the others.** The original *R*/*Q* did not
survive in the 5e's materials: the Visual Basic package (`UNNRQnew.asc`) and the MATLAB
`UNIFAC_data.mat` both carry the Dortmund values, and the Mathcad worksheets carry none. Only
the original *interaction* matrix survived, as `Unfa44` in the `.mat`, which is where
`unifac_interactions_original.csv` comes from. Until 2026-08-17 this file held **six** rows,
recovered from the 5e's Illustration 9.5-2 — which had worked its whole calculation on
original-UNIFAC values while citing Table 9.5-2. **Those six now serve as an independent check
on the 2e transcription, and all six agree to the last digit.**

### What arbitrates the transcription

Four checks, chosen because they fail differently. All are in the build script.

| check | what it catches |
|---|---|
| the **six recovered rows** reproduce | a wholly independent source, 5e print vs 2e print |
| **main groups = the interaction matrix's keys**, exactly 1–44 | a dropped or invented main group |
| **one hydrogen, one increment**: *R* falls 0.2274 per aliphatic H, 0.1661 per aromatic H | a transposed or misread digit |
| the **5e's printed group sums** — benzene *r* = 3.1878, *q* = 2.4000; 2,2,4-trimethyl pentane *r* = 5.8463, *q* = 5.0080 | any error in the rows those sums touch |

**The increment law fails on the chlorinated series, and that is the source's doing.** The 2e
prints CH₂Cl 1.4654, CHCl 1.2380, CCl 1.0060 — steps of 0.2274 then **0.2320** — and CH₂Cl₂
2.2564, CHCl₂ 2.0606, CCl₂ 1.8016 — **0.1958** then **0.2590**. Both were re-read off the
600-dpi render and are as printed. Recorded rather than quietly excluded, because the next
reader will notice them too.

**One printed example was moved.** The 2e prints *"Chloroform"* against **CCl₄**, one row
below where it belongs — chloroform is CHCl₃. Corrected in the row rather than in the consumer,
with the reason in that row's `source` column.

### The numbering trap — the two sets do not agree on what a subgroup number means

Subgroup numbers are **not in the table**; they are the standard UNIFAC numbering, which
`unifac_subgroups.csv` already carries. The two parameter sets share it **through 77**, and
then part company:

| number | original | modified (Dortmund) |
|---|---|---|
| 27 | `FCH2O` | `cy-CH2 OCH2` |
| 37–39 | `C5H5N`, `C5H4N`, `C5H3N` | `AC2H2N`, `AC2HN`, `AC2N` |
| 78–81 | `SiH3`, `SiH2`, `SiH`, `Si` | `cy-CH2`, `cy-CH`, `cy-C`, `OH(s)` |
| 82–85 | `SiH2O`, `SiHO`, `SiO`, `NMP` | `OH(t)`, two `cy-CH2 O` variants, `CNH2` |

**Both numbers exist in both tables, so a mismatched assignment does not raise — it returns a
plausible wrong answer.** Tetrahydrofuran is the trap in miniature: original UNIFAC calls it
1 `FCH2O` + 3 `CH2` = `{27: 1, 2: 3}`, Dortmund 1 `cy-CH2 OCH2` + 2 `cy-CH2` = `{27: 1, 78: 2}`.
Different numbers, different counts, same molecule. This is also why the legacy MATLAB file
carried *silicon* names on main groups 42 and 43 while holding Dortmund's cyclic values — the
same collision, seen from the other side (see the digitization audit above).

**So call `unifac_groups(name, kind)`**, which reads the matching table, instead of reusing a
group dict across kinds. `UNIFAC(kind=...)` pairs *R*/*Q* with $a_{mn}$ for the same reason.

### DECIDED 2026-08-17 — the original model is a supported option, not a one-off

**AUTHOR:** *"We made the decision to use new Dortmund UNIFAC versus old UNIFAC, but have built
the code now such that it can use both. That's a similar approach to packages like Aspen, and
we should just make it available. Our default will still be Dortmund UNIFAC, but the option to
try the other will be there."*

So `kind` is a first-class choice, exactly as Aspen Plus offers `UNIFAC` against `UNIF-DMD`,
with Dortmund the default everywhere in the book. With this table in place the option is real
across the whole group inventory — water, alcohols, ketones and the rest — and not just for the
hydrocarbons and aldehydes Figure 10.2-8 needed.

> **Never mix the two sets.** *R*, *Q* and $a_{mn}$ come from one regression together; the
> 6e's own footnote to Table 9.5-2 says so. `UNIFAC(kind=...)` selects both together, which is
> why neither table is meant to be loaded on its own.

## `steam_*.csv` — Appendix A.III, the thermodynamic properties of water and steam

A digitization of the book's **own Appendix A.III**, extracted from the page proofs by
`tools/parse_appendix_A3.py`. Figure 3.3-1(a) (the Mollier diagram) and 3.3-1(b) (the
$T$–$S$ diagram) are drawn from these files, which is why the charts and the appendix a
student reads them against cannot disagree.

| file | what it is | rows |
|---|---|---|
| `steam_saturation_T.csv` | saturation, entered by temperature, 0.01–374.14 °C | 71 |
| `steam_saturation_P.csv` | saturation, entered by pressure | 73 |
| `steam_superheat.csv` | superheated vapor, 36 isobars (0.01–60 MPa) × 25 temperature levels (50–1300 °C) | 545 |
| `steam_compressed_liquid.csv` | compressed liquid, 6 isobars (5–50 MPa) | 111 |
| `steam_solid_vapor.csv` | saturated solid–vapor (below the triple point) | 22 |

Columns: `T_C`, `P_MPa`, and then `V`/`U`/`H`/`S` — with the `l`/`v`/`d` suffixes
(`Hl`, `Hv`, `dH`) for the two saturation tables, which carry the liquid value, the
vapor value and their difference. In the superheat and compressed-liquid files a row
with `sat = True` is the table's own `Sat.` row: the saturated state at that pressure,
whose temperature is in `Tsat_C` rather than `T_C`.

**Verified, not merely parsed.** Every row satisfies $P\hat V = \hat H - \hat U$ to
printing precision, the saturation rows satisfy $\Delta\hat X = \hat X^{\rm V} -
\hat X^{\rm L}$, and ten values typed by hand off the printed page reproduce to every
printed digit.

**Five cells are corrected.** The 5e appendix carries five typeset errors, each
established two independent ways before being changed; the arithmetic is worked out in
`ch3/validation/Appendix_A3_corrections.ipynb`, which is the record for this table.
They are applied by `tools/parse_appendix_A3.py` (`DATA_FIXES`) as it writes these
files, so the CSVs remain exactly what the parser produces and a re-run cannot silently
revert them.

| table | row | column | printed | corrected |
|---|---|---|---|---|
| saturation, pressure | 1.4 MPa | $\hat H^{\rm V}$ | 2790.6 | **2790.0** |
| superheat | 3.5 MPa, 500 °C | $\hat U$ | 3103.0 | **3103.8** |
| superheat | 10 MPa, 1300 °C | $\hat U$ | 4460.5 | **4660.5** |
| compressed liquid | 50 MPa, 240 °C | $\hat V$ | 0.001 107 2 | **0.001 170 2** |
| compressed liquid | 50 MPa, 60 °C | $\hat S$ | 0.8502 | **0.8052** |

These same five corrections go into the printed 6e Appendix A.III. Until they do,
the CSVs and the printed page disagree in these five cells and nowhere else.

**EMF:** These five corrections made on 2026-08-09 11:41

**A.III is not IAPWS.** The appendix is an older generation of steam tables and
differs from current IAPWS-95/IF97 values in the last digit or two. That is the point —
consistency with the printed table — but nothing built from these files may be labeled
IAPWS.

> **Credit.** Appendix A.III is reproduced in the book *"From G. J. Van Wylen and R. E.
> Sontag, Fundamentals of Classical Thermodynamics, S.I. Version, 2nd ed., John Wiley &
> Sons, New York (1978). Used with permission."* Keep this credit with the CSVs in any
> copy or redistribution — the license note above applies to them as it does to
> `pure_property.csv`.
>
> **Permission to publish the digitized table itself is not yet confirmed.** Reusing
> the appendix in the book is a cleared grant; publishing it as a machine-readable CSV in
> a public repository is a different act. Wiley published both books, so this is likely
> straightforward — but confirm before the public release.

## Provenance
`pure_property.csv` / `react_property.csv` extracted from `pure_prop.mdb` (618 rows) and
`React.mdb` (99 rows) — the shared property database embedded in the 5e Visual Basic
Property / Peng-Robinson / ChemEq programs. `unifac_*.csv` extracted from the 5e MATLAB
`UNIFAC_data.mat`. `steam_*.csv` extracted from the 5e page proofs of Appendix A.III by
`tools/parse_appendix_A3.py`. `mixing_*.csv` transcribed from the 5e Tables 8.6-1 and
8.6-3 and verified against each table's own identity (2026-08-09). See Appendix B and
`revision_notes/bapp02.md` for the modernization decisions.
