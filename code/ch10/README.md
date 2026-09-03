# `code/ch10/` — Vapor-liquid equilibrium in mixtures

Chapter 10 is where the mixture models of Chapter 9 are put to work on phase diagrams.
Naming follows the ch6 convention, `<method>_<topic>_<substance>_<role>`; see
[`code/ch6/README.md`](../ch6/README.md).

## Notebooks

| notebook | backs | QR key | uses |
|---|---|---|---|
| [`raoult_diagrams_pentane_heptane_example.ipynb`](raoult_diagrams_pentane_heptane_example.ipynb) | **Illustration 10.1-1**, `c10uf001a`–`d` | `i10.1-1` | `GammaPhi`, `Ideal`, `ClausiusClapeyron` |
| [`hexane_triethylamine_diagrams_figure.ipynb`](hexane_triethylamine_diagrams_figure.ipynb) | **Figures 10.1-1 … 10.1-5** | `f10.1-1` … `f10.1-5` | `GammaPhi`, `Ideal`, `psat_from_database` |
| [`azeotrope_to_vanlaar_benzene_cyclohexane_example.ipynb`](azeotrope_to_vanlaar_benzene_cyclohexane_example.ipynb) | **Illustration 10.2-1**, `c10uf002`–`003` | `i10.2-1` | `VanLaar`, `RegularSolution`, `azeotrope` |
| [`pressure_swing_methyl_acetate_methanol_example.ipynb`](pressure_swing_methyl_acetate_methanol_example.ipynb) | **Illustrations 10.2-2 / 10.2-3**, `c10uf004`–`005`; **Figures 10.2-2 / 10.2-3** | `i10.2-2`, `i10.2-3`, `f10.2-2`, `f10.2-3` | `VanLaar`, `txy`, `azeotrope`, `total_reflux_steps` |
| [`vle_correlation_benzene_TMP_example.ipynb`](vle_correlation_benzene_TMP_example.ipynb) | **Illustration 10.2-4**, `c10uf006`–`007`; **Figure 10.2-7** | `i10.2-4`, `f10.2-7` | `RedlichKisterGex`, `VanLaar`, `fit_binary` |
| [`excess_properties_benzene_TMP_figure.ipynb`](excess_properties_benzene_TMP_figure.ipynb) | **Figure 10.2-6** | `f10.2-6` | `RedlichKisterGex` |
| [`azeotrope_test_ethyl_acetate_benzene_homework.ipynb`](azeotrope_test_ethyl_acetate_benzene_homework.ipynb) | **Problem 10.2-72** — no print art | `p10.2-72` | `UNIFAC`, `VanLaar`, `txy`, `azeotrope` |


## The two benzene / 2,2,4-trimethyl pentane notebooks are a pair

They share a data set and split it the way the chapter does. `vle_correlation` reduces
**one** isotherm at 55 °C to activity coefficients and a Redlich-Kister correlation;
`excess_properties` uses **all five** of Weissman and Wood's temperatures, because
$H^{ex}$ and $S^{ex}$ come from a temperature derivative and one isotherm cannot supply
one. Read together they carry the section's real argument: the excess Gibbs energy is
what a VLE measurement gives you, and everything else costs accuracy to get.

**Reading the source paper (*J. Chem. Phys.* **32**, 1153) overturned two things this
folder previously asserted**, and both are recorded in the notebooks so they are not
reinstated:

1. The Redlich-Kister constants Illustration 10.2-4 prints, (1389.0, 419.45, 109.83)
   J/mol, are **not** Weissman and Wood's. Theirs at 55 °C are (1354.7, 424.5, 101.0).
2. Figure 10.2-7's consistency test fails on the printed numbers and **passes on the same
   measurements** once the vapor-phase correction the illustration drops is restored. The
   data is consistent; the reduction is what fails.

## Section 10.2's worked example changed system (2026-08-14)

The 5e worked Illustrations 10.2-2 and 10.2-3, Figures 10.2-2 to 10.2-4 and a two-column
distillation design on **ethyl acetate / benzene**. That azeotrope is not real: it is an
artifact of Table 9.5-1's van Laar constants for the pair, which imply γ<sup>∞</sup> = 3.16
and 2.51 for two species boiling 3 °C apart, and standard azeotrope compilations do not list
it. UNIFAC, predicting from structure alone, gives γ within a few percent of 1 and no
azeotrope; the 5e's own bracketed Aspen note already concedes *"the van Laar model is not a
good one for this system."*

The worked material is now **methyl acetate / methanol**, whose azeotrope is real (computed
0.655 at 53.35 °C against a reported 0.654 at 53.5 °C, from Table 9.5-1's own constants with
no tuning), and the ethyl acetate / benzene notebook is **kept and repointed at Problem
10.2-72**, which asks the very question it answers. The model-limits lesson is stronger as a
problem the student works than as an illustration built on the wrong answer.

**Two things about the new system are not the 5e's, and both matter.**

1. **The pressure swing runs upward.** Lowering the pressure moves *this* azeotrope away
   from the feed (x_az = 0.655 at 1.013 bar, 0.713 at 0.5 bar, 0.831 at 0.1013 bar), so a
   low-pressure second column cannot work at any pressure. The second column is at **3 bar**
   — the smallest extrapolation of Table 9.5-1's 53.7–64.6 °C fit range that separates, and
   less of a stretch than the 5e's own 15 °C column was.
2. **Both products leave as bottoms**, and which product leaves which port is *computed*, not
   assumed: `total_reflux_steps` applies the equilibrium curve or its inverse rather than
   presuming the staircase's direction. The fifth edition's own figure has column 2's
   distillate and bottoms swapped.

**Figures 10.2-2 and 10.2-3 were redrawn 2026-08-14, and the first version was wrong** — the
author caught it. Each column is now stepped off **from its own feed in both directions**: down
toward the reboiler, which reaches the 95 % product in 3 stages (column 1) and 4 (column 2) and
is clipped to end exactly on the specification line; and up toward the condenser, which
**pinches** against that pressure's azeotrope and stops. What was wrong before:

- The staircases spanned different things — column 1 from its distillate, column 2 from its feed
  — so neither figure showed the section its caption described, and the feed and product lines
  did not line up with the steps. One product line had no steps reaching it at all.
- `total_reflux_steps` could only ever *descend*: it applied the inverse curve unconditionally,
  so the rectifying section could not be drawn and asking for it silently produced the stripping
  section. It now takes `direction="up"`/`"down"` and reports whether it pinched.
- **The ascending staircase was then drawn reflected in the diagonal.** Its corners must lie
  **on the equilibrium curve** — ascending that corner is $(x_n, f(x_n))$, so the step goes *up*
  to the curve and then across; descending it is $(f^{-1}(x_n), x_n)$, across and then down.
  Reusing one order for both put the ascending steps on the far side of $x = y$, touching the
  curve nowhere. **The stage counts were unaffected**, which is why it passed every numerical
  check and had to be caught by eye. There is now a doctest asserting the invariant.
- The last stage steps **past** the product specification rather than being clipped to it,
  matching Illustration 10.1-6's own words (*"no further stages are needed"*) and keeping every
  corner on the curve.
- **The bottoms line marks where the staircase lands, not the specification** (author,
  2026-08-14). An integer stage count cannot meet a spec exactly — column 1 runs 0.400 → 0.209
  → 0.075 → 0.022, so two stages miss 0.05 and three overshoot — and a line at the spec floats
  between two corners in any correct construction. It is drawn at the achieved composition with
  the spec named in the label: `bottoms 0.022 (spec 0.05 met)`. Production must not "correct"
  these to 0.05 / 0.95 to match the prose; the figure reports **97.8 % methanol in three stages
  and 96.2 % methyl acetate in four**, and ¶769's 95 % is the requirement they satisfy.
- **The vertical labels are measured, not hand-placed.** A rotated label is as long as its text,
  so a start height that suits one panel runs off the axis in another — which is what happened
  to the 3 bar bottoms label. Each one is now drawn, measured with the renderer, and lifted if
  it would cross the frame. Checked in both renderings (the 3.5 in print panel and the narrower
  two-panel display figure): all six labels sit inside the axes.
- **Column 1's "22 stages" was an artifact of starting at the azeotrope**, which is a fixed
  point of the recursion: the first eleven steps moved the composition by less than 0.001. From
  the feed, as specified, it is **3 stages**. That is also why the count once flipped from 21 to
  22 on a change in the fourth decimal of a vapor-pressure constant — the number was measuring
  the arithmetic, not the chemistry. **A stage count needs two stated compositions**, and the
  distillate of an azeotropic column is not one of them.

**Figures 10.1-1 and 10.1-2 carry real data, transcribed from the source paper.**
Humphrey and Van Winkle, *J. Chem. Eng. Data* **12**, 526 (1967), Table II
— nine points for triethylamine / *n*-hexane
at 60 °C. The values are read from the paper's own
table, never digitized from the printed figure.
**The paper's component 1 is triethylamine and the book's abscissa is hexane**, so
every composition is flipped once, explicitly, on read.

**The transcription was checked against the table's own identity**, not by re-reading
the scan — the page is a two-column layout interleaving six binaries, so a slipped column
would be silent. Solving γᵢ = yᵢP/(xᵢPᵢˢᵃᵗ) row by row gives nine independent estimates of
each pure vapor pressure: triethylamine 0.3904 ± 0.0032 bar, hexane 0.7559 ± 0.0055.
Scatter under 1.5 % is what proves the columns are aligned.

Those imply pure vapor pressures 1.6 % and 0.3 % from the **0.3843 / 0.7583** SIS
prints. The Raoult lines are drawn from the book's values (what §10.1's prose and
Fig. 10.1-3 use) and the points plotted as measured; the ~0.002 bar offset is why the
points do not sit exactly on the lines at the ends, and the notebook states it.

The figures keep their "[Based on data of J. L. Humphrey and M. Van Winkle…]" credit —
recomputing the lines does not change whose the points are.

## QR keys

Printed keys resolve through `cbethermo.org/<key>` to the notebook page.

| goes beside | key | art |
|---|---|---|
| chapter opener | `c10` | `c10qf000.eps` |
| Illustration 10.1-1 | `i10.1-1` | `c10qx0101.eps` |
| Figure 10.1-1 | `f10.1-1` | `c10qf001.eps` |
| Figure 10.1-2 | `f10.1-2` | `c10qf002.eps` |
| Figure 10.1-3 | `f10.1-3` | `c10qf003.eps` |
| Figure 10.1-4 | `f10.1-4` | `c10qf004.eps` |
| Figure 10.1-5 | `f10.1-5` | `c10qf005.eps` |
| Illustration 10.2-1 | `i10.2-1` | `c10qx0201.eps` |
| Illustration 10.2-2 | `i10.2-2` | `c10qx0202.eps` |
| Illustration 10.2-3 | `i10.2-3` | `c10qx0203.eps` |
| Illustration 10.2-4 | `i10.2-4` | `c10qx0204.eps` |
| Figure 10.2-2 | `f10.2-2` | `c10qf014.eps` |
| Figure 10.2-3 | `f10.2-3` | `c10qf015.eps` |
| Problem 10.2-72 | `p10.2-72` | `c10qp272.eps` |



## Print art

| file | staged as | figure |
|---|---|---|
| `pdf/Illustration_10.1-1a.pdf` | `…/c10uf001a.pdf` | x-y, 50 °C |
| `pdf/Illustration_10.1-1b.pdf` | `…/c10uf001b.pdf` | P-x-y, 50 °C |
| `pdf/Illustration_10.1-1c.pdf` | `…/c10uf001c.pdf` | x-y, 1.013 bar |
| `pdf/Illustration_10.1-1d.pdf` | `…/c10uf001d.pdf` | **T-x-y**, 1.013 bar |
| `pdf/Fig_10.1-1.pdf` … `Fig_10.1-5.pdf` | `…/c10f001.pdf` … `c10f005.pdf` | hexane / triethylamine |
| `pdf/Illustration_10.2-1_fig1.pdf` | `…/c10uf002.pdf` | γ_B, γ_C vs. x_B |
| `pdf/Illustration_10.2-1_fig2.pdf` | `…/c10uf003.pdf` | P-x-y, benzene / cyclohexane |
| `pdf/Illustration_10.2-2_fig1.pdf` | `…/c10uf004.pdf` | P-x-y + x-y, 60 °C |
| `pdf/Illustration_10.2-2_fig2.pdf` | `…/c10uf005.pdf` | T-x-y + x-y, 1.013 bar |
| `pdf/Fig_10.2-2.pdf` | `…/c10f014.pdf` | column 1 at 1.013 bar, total reflux |
| `pdf/Fig_10.2-3.pdf` | `…/c10f015.pdf` | column 2 at 3 bar, total reflux |
| `pdf/Illustration_10.2-4_fig1.pdf` | `…/c10uf006.pdf` | Gᴱ and ΔG_mix |
| `pdf/Illustration_10.2-4_fig2.pdf` | `…/c10uf007.pdf` | partial pressures |

**The `c10uf00N` names above are this repository's, and are unique per figure.** The fifth
edition's own art filenames restart inside each illustration, so several of them collide;
the names in the table are the ones to use.

**Illustration 10.2-4's two figures.** Read off the fifth-edition page: **Figure 1** is Gᴱ and ΔG_mix against x_B; **Figure 2**
is the species partial pressures against x_B with the Raoult's law lines dashed.

**Figure 2 departs from the 5e deliberately.** The 5e draws it on twin ordinates at a
2:1 ratio, benzene left and 2,2,4-trimethyl pentane right. Both are pressures in bar, so
one axis carries both without distortion and removes the hazard of reading a curve against
the wrong scale. Flag to production; the caption needs no change.

Four separate files, not two pairs — the 5e sets these as four captioned panels, so the
art has to match that granularity. All four are Computer Modern throughout (checked with
`get_fonts`).

**Figure d is drawn as a T-x-y and the fifth edition's caption calls it a P-x-y.** The
caption is wrong; the art is right.

**Region labels are placed, not hard-coded.** `pxy_chart`/`txy_chart` put `Liquid` and
`Vapor` where each region is actually widest, by scanning the abscissa for the largest gap
between the bounding curve and its edge of the frame. The first version offset them a fixed
6 % from the axis limits, which on figure b put `Liquid` far above the bubble line and
`Vapor` straight through the dew line — wrong in opposite directions in one figure, which
is what a geometry-blind rule costs. Override with `region_pos=((x, v), (x, v))` in display
units, rename with `region_labels=`, or drop them with `regions=False`.

## Checked against the fifth edition's Aspen output

This chapter's first notebook was checked against **the Aspen Plus output Sandler shipped
with the fifth edition** — the spreadsheet the illustration's bracketed note points at,
82 rows.

**The result is sharper than "they roughly agree."** Driving Raoult's law with the
book's own correlation leaves a 1.76 % gap in bubble pressure. Driving it with **Aspen's
own endpoint vapor pressures** collapses that to **0.02 %**, and the vapor composition
from 4.6 × 10⁻³ to 1 × 10⁻⁵ — so the mixture thermodynamics are identical and the entire
visible difference is the pure-component correlations, exactly as the illustration claims.
It also proves Aspen was run on an **ideal** liquid here: no choice of pure vapor pressure
could collapse the difference otherwise.


## The two modules this chapter added

Chapter 10 needed machinery Chapter 9 did not have, and it went into the package rather
than into the notebooks.

**[`thermo/vle.py`](../thermo/vle.py)** — the γ-φ method. Activity coefficients in the
liquid, ideal gas in the vapor. Vapor-pressure correlations (`Antoine`,
`ClausiusClapeyron`, `TabulatedPsat`), the `Ideal` activity model, `GammaPhi` with
bubble/dew/flash, and the curve generators `pxy`, `txy`, `azeotrope`.

**[`thermo/vle_chart.py`](../thermo/vle_chart.py)** — the drawing layer. `pxy_chart`,
`txy_chart`, `xy_chart`, `tie_lines`, `mark_azeotrope`, `mccabe_thiele`. Lazily loaded,
so `from thermo import GammaPhi` stays free of matplotlib.

**Why the split is where it is.** The curve *generators* are in `vle.py`, in the
physics tier — the opposite of where `ph_chart.py` keeps `dome` and `isotherm`. A
saturation dome exists in order to be drawn; a VLE curve does not. This chapter prints
these numbers as **tables** (Illustration 10.1-1's *x*, *y*, *P*; Illustration 10.1-5's
L = 1.0 → 0.0 sweep; Table 10.2-3, which Illustration 10.2-6 generates), and a notebook
that wants numbers and draws nothing should not have to import matplotlib to get them.

**One generator serves both halves of the chapter.** `GammaPhi`'s solvers carry the
same names and signatures as `PRMixture`'s, so `pxy(model, T)` never learns which it was
handed. Sections 10.1 and 10.2 (γ-φ) and Section 10.3 (φ-φ) draw their diagrams through
the same code — verified, not assumed. It is the same trick as `ActivityModel.gamma(x, T)`
in Chapter 9.

## Validated against the book before any notebook used it

`vle.py` was checked against five illustrations first, per the ROADMAP rule:

| check | result |
|---|---|
| Illust. 10.1-1, Raoult's law | exact closed form |
| Illust. 10.1-2, ternary bubble/dew **T** | 334.9 / 350.6 K vs. printed 334.6 / 350.5 |
| Illust. 10.1-3, ternary bubble/dew **P** | 1.414 / 0.878 bar vs. printed 1.413 / 0.877 |
| Illust. 10.2-1, van Laar from one azeotropic point | α = 0.125, β = 0.0919 — exact |
| Illust. 10.2-2(b), isobaric azeotrope | x = 0.5149, T = 343.92 K vs. printed 0.515, 343.92 |
| Illust. 10.1-7, McCabe-Thiele | **4 stages**, and minimum reflux **q = 0.237** — both the book's |

**Two of the checks failed, and the book is what was wrong.** Illustration 10.2-2(a)'s
printed azeotrope (0.581, 0.518 bar) should read **0.518 and 1.165 bar** — the 5e's own
art agrees with the computation, not with its text — and Illustration 10.1-4's
**K₅ = 2.7406** contradicts the same page's vapor pressure and should be **2.686**. The
first of the two no longer applies — that illustration is methyl acetate / methanol in the
sixth edition, so there is no printed number left to correct — but the validation stands: those checks are
against the 5e as printed, and they are what caught the pair in the first place.

## Conventions

- **SI throughout** — T in K, P in Pa. The vapor-pressure correlations take the book's
  constants **exactly as printed in bar** and convert internally.
- **Three correlation forms, not one.** The chapter uses ln P = A − B/(T+C),
  ln P = A − ΔH_vap/RT, and log₁₀ P = −A/T + B. A module offering only "Antoine" cannot
  run the chapter.
- **Do not filter NaNs before plotting.** `pxy`/`txy` return NaN where the model has no
  solution, and the gap is the result — Figure 10.3-13 prints exactly such a gap
  ("region of nonconvergence with k₁₂ = 0"), which `pr_mixture` reproduces.
- **Black and white means BLACK — no gray.** Every line, marker and label in
  `vle_chart.py` is pure black (author, 2026-08-13); gray survives a laser printer and
  dies in offset printing at small sizes. Tie lines, the x = y diagonal and the
  McCabe-Thiele operating lines were all gray in the first cut and are now black.
  `charts.py` keeps a `GRAY` for the property charts' *grid rules* — a different job —
  and `vle_chart.py` does not import it. **Verified on the staged PDFs**, not just the
  source: all nine contain only gray levels 0 and 1.

  What separates one line from another is weight, dash and geometry:

  | element | style |
  |---|---|
  | envelope (bubble + dew) | solid, `LW["sat"]` — heaviest on the chart |
  | tie lines | solid, 0.6 — short horizontals inside the envelope |
  | x = y diagonal | **dashed**, 0.8 — the only dashed line on an x-y chart |
  | operating lines | solid, 0.9 — straight, against a curved equilibrium |
  | stage steps | solid, 0.7 — a staircase, against everything else |
  | measured data | open markers, never filled |
- **The .ipynb on disk is canonical.** Never regenerate a notebook from a builder script
  once it exists — edit it cell by cell.
