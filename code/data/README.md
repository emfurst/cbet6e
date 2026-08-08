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

For chemical-equilibrium work. `Cp = A + B T + C T^2 + D T^3 + E/T^2`.

| Column | Meaning | Units |
|---|---|---|
| `Name` | species | — |
| `DG`, `DH` | Gibbs energy / enthalpy of formation | kJ/mol |
| `A, B, C, D, E` | 5-term heat-capacity polynomial coefficients | — |
| `ID` | index | — |

## UNIFAC group-contribution parameters

Extracted from the legacy MATLAB `UNIFAC_data.mat` for the Python UNIFAC substitute
(`code/thermo/`).

- `unifac_subgroups.csv` — 95 subgroups: `subgroup_no, main_group_no, subgroup_name,
  main_group_name, R, Q`. **These R/Q are the *modified* (Dortmund) values.**
- `unifac_interactions_modified.csv` — modified (Dortmund) UNIFAC main-group interactions:
  `main_i, main_j, a_ij, b_ij, c_ij` (temperature-dependent, $\Psi=\exp[-(a/T+b+cT)]$),
  1220 nonzero pairs over 56 main groups.
- `unifac_interactions_original.csv` — original UNIFAC main-group interactions:
  `main_i, main_j, a_ij` (single parameter, $\Psi=\exp(-a/T)$), 838 nonzero pairs over
  44 main groups (recovered from the `Unfa44` table).

**Gap for the *original* variant:** the `.mat` does not contain the original-UNIFAC
subgroup R/Q constants or the 44-main-group names (the Matlab GUI only ran modified
UNIFAC). Those must be sourced separately (book appendix or the published original-UNIFAC
tables) before the original variant is complete.

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
printing precision, and ten values typed by hand off the printed page reproduce to every
printed digit. Five rows do not — they are typos in the printed appendix, listed in
`legacy-5e/Errata_Sandler_5Ed.docx.md`. **The CSVs keep the printed values**: this is a
digitization of A.III, not a correction of it. The one place a correction is applied is
`CL_CORRECTIONS` in `ch3/Steam_charts_from_appendix_A3.ipynb`, for two compressed-liquid
cells from which a chart line is *derived*.

⚠️ **A.III is not IAPWS.** The appendix is an older generation of steam tables and
differs from current IAPWS-95/IF97 values in the last digit or two. That is the point —
consistency with the printed table — but nothing built from these files may be labeled
IAPWS.

> **Credit.** Appendix A.III is reproduced in the book *"From G. J. Van Wylen and R. E.
> Sontag, Fundamentals of Classical Thermodynamics, S.I. Version, 2nd ed., John Wiley &
> Sons, New York (1978). Used with permission."* Keep this credit with the CSVs in any
> copy or redistribution — the license note above applies to them as it does to
> `pure_property.csv`.
>
> ⬜ **Permission to publish the digitized table itself is not yet confirmed.** Reusing
> the appendix in the book is a cleared grant; publishing it as a machine-readable CSV in
> a public repository is a different act. Wiley published both books, so this is likely
> straightforward — but confirm before the public release.

## Provenance
`pure_property.csv` / `react_property.csv` extracted from `pure_prop.mdb` (618 rows) and
`React.mdb` (99 rows) — the shared property database embedded in the 5e Visual Basic
Property / Peng-Robinson / ChemEq programs. `unifac_*.csv` extracted from the 5e MATLAB
`UNIFAC_data.mat`. `steam_*.csv` extracted from the 5e page proofs of Appendix A.III by
`tools/parse_appendix_A3.py`. See Appendix B and `revision_notes/bapp02.md` for the
modernization decisions.
