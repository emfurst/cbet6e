Title:  Chemical, Biochemical, and Engineering Thermodynamics, 6th Edition
Author: Stanley I. Sandler and Eric M. Furst
ISBN:   9781394435449

CODE AND DATA FOR THE SIXTH EDITION
===================================

This repository holds the Python and Jupyter materials for the sixth edition,
replacing the Visual Basic, Mathcad, and MATLAB programs distributed with earlier
editions. Companion website material (notation primers, study prompts, chapter
pages) is generated separately and linked from the book.

CONTENTS
  ch1/, ch4/, ch6/, ch7/   per-chapter Jupyter notebooks (worked examples and homework)
  thermo/                  reusable Python package (Peng-Robinson EOS, UNIFAC, mixtures)
  data/                    reference property data as plain CSV, with schema in data/README.md
  README.md                orientation and per-chapter guide
  pyproject.toml           the Python environment (library requirements), used by uv

HOW THE NOTEBOOKS ARE ORGANIZED
  Many topics appear twice: a self-contained notebook that builds the method from the
  governing equations (so you can see how it works), and a "_thermo" twin that imports
  the thermo package (so you can see how it is reused). Self-contained notebooks need
  only numpy, scipy, matplotlib, and pandas. The twins additionally need this
  repository's thermo package and data/ directory.

RUNNING THE CODE
  Requires Python 3.10+ with numpy, scipy, matplotlib, pandas, and Jupyter. The
  recommended tool is uv (https://docs.astral.sh/uv/), which installs Python, creates
  the virtual environment, and installs the libraries in one step:
    git clone https://github.com/emfurst/cbet6e.git
    cd cbet6e/code
    uv run jupyter lab
  The first run builds the environment from pyproject.toml; later runs just start
  Jupyter. If you prefer the older route, python3 -m venv plus pip install works too
  (see README.md). Notebooks can also be opened in an editor with an AI coding
  assistant.
  Self-contained notebooks also run in Google Colab with nothing installed; see the
  chapter pages on the companion website for one-click links.

DATA PROVENANCE AND CREDIT
  The constants of pure fluids in data/pure_property.csv are adapted from R. C. Reid,
  J. M. Prausnitz, and B. E. Poling, The Properties of Gases and Liquids, 4th ed.,
  McGraw-Hill, New York, 1986, Appendix A, with corrections and data from other
  sources. Please keep this credit with the data. See data/README.md.

LICENCE
  See LICENSE. The licence covers the code in this repository; it does not purport to
  relicense the third-party-derived property data described above.
