# Chemical, Biochemical, and Engineering Thermodynamics, 6th Edition — code and data

Python and Jupyter materials for Sandler & Furst, *Chemical, Biochemical, and Engineering
Thermodynamics*, 6th Edition (ISBN 9781394435449), replacing the Visual Basic, Mathcad,
and MATLAB programs of earlier editions.

**Companion website:** https://emfurst.github.io/cbet6e/ — per-chapter notation primers
for AI assistants, study prompts, and links to every notebook.

```
code/
  ch1/ ch4/ ch6/ ch7/   per-chapter notebooks (worked examples and homework)
  thermo/               reusable package: Peng-Robinson (pure + mixtures), UNIFAC
  data/                 reference property data as plain CSV (see data/README.md)
```

## Two ways to work

Many topics appear **twice**: a *self-contained* notebook that builds the method from the
governing equations, and a `_thermo` twin that imports the package instead. The first
shows how the method works; the second shows how it is reused.

```python
import sys; sys.path.append("..")
from thermo import PengRobinson
pr = PengRobinson.from_database("methane")   # Tc, Pc, omega, Cp from data/pure_property.csv
```

Self-contained notebooks need only `numpy`, `scipy`, `matplotlib`, and `pandas`, and run
in Google Colab with nothing installed. The `_thermo` twins additionally need this
repository's `thermo/` package and `data/` directory — so clone the repository:

```bash
git clone https://github.com/emfurst/cbet6e.git && cd cbet6e/code
git pull        # later, to pick up corrections
```

If you use an AI coding assistant, keep the chapter's notation file from the companion
website in your project so the assistant follows the book's conventions.

## Data credit

The constants of pure fluids in `data/pure_property.csv` are adapted from R. C. Reid,
J. M. Prausnitz, and B. E. Poling, *The Properties of Gases and Liquids*, 4th ed.,
McGraw-Hill, New York, 1986, Appendix A, with corrections and data from other sources.
Please keep this credit with the data — see [`data/README.md`](code/data/README.md).

## Licence

[MIT](LICENSE) for the code in this repository. The licence does **not** extend to the
third-party-derived property data above, nor to the text and figures of the printed book.

## Errata

Found a mistake? Please open an [issue](https://github.com/emfurst/cbet6e/issues).
