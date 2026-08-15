# `code/ch9/` — Estimation of the Gibbs energy and fugacity of a component in a mixture

Chapter 9 is where the mixture models live. Naming follows the ch6 convention,
`<method>_<topic>_<substance>_<role>`; see [`code/ch6/README.md`](../ch6/README.md).

## Notebooks

| notebook | backs | uses |
|---|---|---|
| [`PR_species_fugacity_ethane_butane_example.ipynb`](PR_species_fugacity_ethane_butane_example.ipynb) | Illustrations 9.4-2, **9.4-3** | `PRMixture`, Table 9.4-1 |
| [`Lewis_Randall_vs_PR_CO2_methane_example.ipynb`](Lewis_Randall_vs_PR_CO2_methane_example.ipynb) | **Illustration 9.4-4** | `PRMixture`, `PengRobinson` |
| [`activity_coefficient_correlation_benzene_TMP_example.ipynb`](activity_coefficient_correlation_benzene_TMP_example.ipynb) | **Illustration 9.5-1**, Figs. 9.5-4, 9.5-5 | `OneConstantMargules`, `VanLaar` |
| [`regular_solution_benzene_TMP_example.ipynb`](regular_solution_benzene_TMP_example.ipynb) | **Illustration 9.6-1**, Fig. 9.6-1 | `RegularSolution`, Table 9.6-1 |
| [`UNIFAC_benzene_TMP_example.ipynb`](UNIFAC_benzene_TMP_example.ipynb) | Illustrations **9.5-2**, **9.6-2**, Fig. 9.6-2 | `UNIFAC`, `UNIQUACModel`, Table 9.5-2 |
| [`Gex_composition_curves_figure.ipynb`](Gex_composition_curves_figure.ipynb) | Figs. **9.5-1**, **9.5-2** | all the correlative models |
| [`Margules_activity_coefficients_figure.ipynb`](Margules_activity_coefficients_figure.ipynb) | Fig. **9.5-3** | `OneConstantMargules` |
| [`Debye_Huckel_HCl_NaCl_example.ipynb`](Debye_Huckel_HCl_NaCl_example.ipynb) | Illustrations **9.10-1**, **9.10-2**, Fig. **9.10-2**, `c09uf003` | `DebyeHuckel`, Table 9.10-1 |
| [`hemoglobin_activity_coefficient_example.ipynb`](hemoglobin_activity_coefficient_example.ipynb) | **Illustration 9.7-1**, `c09uf002` | `LogPolynomialActivity` |
| [`Wong_Sandler_mixing_rule_acetone_water_example.ipynb`](Wong_Sandler_mixing_rule_acetone_water_example.ipynb) | **Sec. 9.9** (no illustration, no figure) | `WongSandler`, `GexFromUNIFAC` |

⭐ **Three of these are one system through three models.** The correlation, regular
solution and UNIFAC notebooks all use benzene / 2,2,4-trimethyl pentane at 55 °C — the
same seven data points, correlated, then predicted two different ways. Read them in that
order; the last one carries the comparison table. The models genuinely disagree, and the
ranking (fitted correlation > UNIFAC > regular solution) is the chapter's argument.

**Where the data come from.** Illustration 9.5-1 says only *"the points in Figs. 9.5-4 and
9.5-5."* The numbers are printed in **Illustration 10.2-4** — Weissman and Wood's
vapor–liquid equilibrium measurements — and that is the source used here. Digitizing the
figure art would be a worse source than the table the figure was drawn from.

## `validation/`

Not teaching material. [`unifac_subgroups_table_9.5-2_validation.ipynb`](validation/unifac_subgroups_table_9.5-2_validation.ipynb)
is the audit trail for digitizing **Table 9.5-2** into `code/data/unifac_subgroups.csv`:
the checks the table makes possible, and the five defects the comparison found. Excluded
from the public repository by the same rule as the other `validation/` folders.

⚠️ Its Colab bootstrap checks `../../thermo` as well as `../thermo`, because it sits one
level deeper than the chapter notebooks. Without that, the bootstrap clones a published
copy of the package over the local one and the notebook tests the wrong code.

⭐ **The electrolyte notebook is the odd one out, and says so.** Every other notebook here
computes activity coefficients from mole fractions. §9.10's independent variable is ionic
strength, and its answer is a single *mean* coefficient for the salt, because
electroneutrality makes the individual ion coefficients unmeasurable. It is also the only
ch9 notebook Chapter 15 reuses.

⭐ **The Wong–Sandler notebook is the other odd one out, for the opposite reason.** It
backs a *section*, not an illustration or a figure, because §9.9 has none of either —
it is prose from end to end, and it is the one place the book teaches a method its own
author invented. It is the only notebook in the book whose printed key is a **section
key**, `s9.9`. Having nothing printed to reproduce, it checks the derivation against
itself instead: the two boundary conditions the mixing rule is built from, verified to
machine precision, and Eq. 9.9-8 confirmed by taking the equation of state to
$10^{11}$ bar.

## All ten are built

✅ **2026-08-13.** The chapter's code is complete — **ten notebooks, twelve QR keys, nine
recomputed figures**, closed out by Illustration 9.7-1 (hemoglobin), which gave §9.7 its
first notebook and left §§9.1–9.3 as the only sections with no code behind them. The VLE demonstrations of the Wong–Sandler rule (§10.3, Figs. 10.3-9
to 10.3-13, all acetone/water) belong to chapter 10 and need its bubble-point drivers;
`WongSandler.ln_phi` is what those will call.

## Two figures in §9.10 that are not what they look like

- **Fig. 9.10-1 is not recomputed, and must not be** — the same rule as Fig. 9.5-2. Its
  three solid curves are Robinson and Stokes' measurements for NaCl, CaCl₂ and CuSO₄, and
  we hold neither the data nor the right to redraw them. It keeps its 5e art. Its
  *caption* is still load-bearing: "for NaCl $I = M$; for CaCl₂ $I = 3M$; and for CuSO₄
  $I = 4M$" is the check on Eq. 9.10-16 the notebook runs.
- **`c09uf003` — Illustration 9.10-2's plot — had no 5e art file at all.** `c09uf001` and
  `c09uf002` are in the 5e figure set and `c09uf003` is not, though the manuscript has the
  slot and the printed page (5e p. 489) has the figure. Second such gap in the chapter
  after Fig. 9.7-3, and the first supplied by recomputation rather than recovery.
  ⚠️ Its points come from the illustration's printed **table**, not from the printed art,
  which is misregistered against its own axes by about 0.03 in $\gamma_\pm$ and drops the
  $M = 6$ point entirely.
