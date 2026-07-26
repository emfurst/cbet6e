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

**Licence note.** The `thermo` package and the notebooks in this repository are released
under the repository's own licence; **that licence does not extend to the underlying
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

## Provenance
`pure_property.csv` / `react_property.csv` extracted from `pure_prop.mdb` (618 rows) and
`React.mdb` (99 rows) — the shared property database embedded in the 5e Visual Basic
Property / Peng-Robinson / ChemEq programs. `unifac_*.csv` extracted from the 5e MATLAB
`UNIFAC_data.mat`. See Appendix B and `revision_notes/bapp02.md` for the modernization
decisions.
