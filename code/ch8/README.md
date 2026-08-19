# Chapter 8 — the thermodynamics of multicomponent mixtures

Chapter 8 is the chapter of definitions and criteria — partial molar properties, the
Gibbs–Duhem equation, the equations of change for a mixture, the phase and chemical
equilibrium criteria, the phase rule. Most of it is derivation, and most of its
illustrations are read off a chart.

Two places in it are genuinely computational, and both are here.

## Naming convention

Notebooks are named `<method>_<topic>_<substance>_<role>.ipynb`, the same shape as
[`../ch6/README.md`](../ch6/README.md) but with the method named rather than the equation
of state (there is no equation of state in this chapter):

- **method** — `RK` (Redlich–Kister correlation), `Hrxn` (heat of reaction), …
- **topic** — `partial_molar_volume`, `partial_molar_enthalpy`, `temperature`, …
- **substance** — `water_methanol`, `N2O4`, …
- **role** — `example` (worked) or `homework`

## The notebooks

| Notebook | Role | Generates |
|---|---|---|
| `RK_partial_molar_volume_water_methanol_example.ipynb` | example, **production** | **Table 8.6-2** and **Figure 8.6-1** from Table 8.6-1 |
| `RK_partial_molar_enthalpy_water_methanol_example.ipynb` | example, **production** | **Table 8.6-4** and **Figure 8.6-3** from Table 8.6-3 |
| `Hrxn_temperature_N2O4_example.ipynb` | example | **Illustration 8.5-2** and its ΔrxnH°–T curve (`c08uf002`) |

All three are **unpaired** — no `_thermo` twin. The first two are *production* notebooks:
they generate numbers the book prints, so there is one authoritative implementation rather
than two that differ in the last digit (the ch6 rule, set by
`PR_properties_table_O2_example.ipynb`). The third is unpaired for a different reason —
its "method" is two integrals and a dictionary of stoichiometric coefficients, so a
self-contained version and a package version would be the same notebook.

## What they share

Both §8.6 notebooks are one idea applied twice. `thermo.RedlichKister` fits

$$\Delta_{\text{mix}}\theta = x_1 x_2 \sum_i a_i (x_1-x_2)^i \tag{Eq. 8.6-5a}$$

to a property change on mixing and differentiates it analytically, which is all
Eqs. 8.6-6a,b (volume) and 8.6-10a,b (any property) need. **The module never mentions
volume or enthalpy** — swapping one data set for the other is the whole difference between
the two notebooks, which is exactly the generalization Eq. 8.6-10 is making.

`thermo.tangent_intercepts` is the same answer by the book's graphical route: a point and a
slope give the two intercepts $A$ and $B$ of Fig. 8.6-1. The volume notebook runs both and
shows they agree.

## The thing worth teaching here

**A property change on mixing is easy to correlate; its derivative at the endpoints is
not.** Adding a term to the fit barely moves the curve and can move the infinite-dilution
partial molar properties by ten percent or more. Both notebooks show that with
`scan_order` rather than hiding it behind a default, and both say plainly that the two
infinite-dilution rows are the least certain numbers in their table.

The enthalpy notebook goes one step further and is the better of the two on this point:
its cross-validated error keeps improving to order 6, but orders 4, 6 and 7 all produce a
**non-monotone** partial molar excess near an endpoint — a wiggle in a region with no data.
Order 5 is chosen because it has the lowest cross-validated error *among the fits that
behave physically*. The best-fitting model and the best-behaved model are not the same
model, and no statistic in the scan can see the difference.

## Data

Both mixing data sets are in `code/data/`, in SI units:

| file | is | check it satisfies |
|---|---|---|
| `mixing_water_methanol_volume.csv` | Table 8.6-1 | $\underline{V} = M_{\text{mix}}/\rho$, and $\Delta_{\text{mix}}\underline V$ from the pure values |
| `mixing_water_methanol_enthalpy.csv` | Table 8.6-3 | $\Delta_{\text{mix}}\underline H = (1-x_1)Q^{+}$ |

Both notebooks run those checks in their opening cells, and both tables pass at the level
of their printed rounding.

**`react_property.csv`'s heat-capacity columns are stored as Appendix A.II *prints*
them** — $b$, $c$, $d$ multiplied by $10^{2}$, $10^{5}$, $10^{9}$. Use
`thermo.data.reaction_cp`, which returns them scaled. Reading the raw columns gives an
answer wrong by three orders of magnitude with no warning.

## Outputs

`output/` holds the generated tables (`.csv` for downstream use, `.txt` formatted for
reading); `pdf/` holds the figures. Neither is published — the allowlist in
`6e_companion_site/tools/publish_code.sh` copies notebooks and this README only.

## Validation

`validation/` compares the generated tables with the 5e as printed. It is author's
material and is excluded from the public repo automatically (anything in a subdirectory of
`code/chNN/` is). Its result, in one line: **Table 8.6-2 reproduces to the last digit and
the check finds a sign error in the printed table; Table 8.6-4 does not reproduce at any
order, because it was drawn rather than computed.**
