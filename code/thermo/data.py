"""Data loaders for the `thermo` package — reads the CSV tables in `code/data/`."""
from pathlib import Path
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def load_pure_properties():
    """The 618-compound pure-component property table (see code/data/README.md)."""
    return pd.read_csv(DATA_DIR / "pure_property.csv")


def get_compound(key):
    """Look up a compound by Formula or Name (case-insensitive). Returns a Series.

    Tries exact Formula, then exact Name, then a Name substring match.
    """
    df = load_pure_properties()
    k = str(key).strip().lower()
    for mask in (df["Formula"].str.lower() == k,
                 df["Name"].str.lower() == k,
                 df["Name"].str.lower().str.contains(k, na=False, regex=False)):
        hit = df[mask]
        if len(hit):
            return hit.iloc[0]
    raise KeyError(f"compound {key!r} not found in pure_property.csv")


# --- the book's own ideal-gas heat capacities -------------------------------
#
# `pure_property.csv` carries the **Reid-Prausnitz-Poling** ideal-gas Cp
# coefficients. Those are not the coefficients tabulated in the book's own
# Appendix A.II, and below room temperature the two sets part company: for oxygen
# they differ by 152 J/mol in H and 1.23 J/(mol K) in S at 73 K, referred to
# 298.15 K. Illustrations 6.4-1, 7.5-1 and 7.5-2 are computed with **Appendix
# A.II**, so a notebook that means to reproduce the printed tables must use these.
#
# Author decision 2026-08-03 (see revision_notes/c07.md): ch6 and ch7 both use
# Appendix A.II, so the book's tables reproduce and the two chapters agree.
#
# Form: Cp* = a + b T + c T^2 + d T^3, J/(mol K), T in K.
# Source: 5e Appendix A.II, "Combustion Gases (Low Temperature Range)", printed page
# 976 (= 9780470504796/pdf/bapp01.pdf p.4). The table gives b x 10^2, c x 10^5 and
# d x 10^9, already scaled below. Appendix A.II states its validity as 273-1800 K;
# Illustration 6.4-1 uses it from 73 K, an extrapolation that is the book's, not ours.
#
# Only the compounds the notebooks and problems need are transcribed. Add rows from
# Appendix A.II as later chapters require them.
APPENDIX_A2_CP = {
    "oxygen":   (25.460,  1.519e-2, -0.715e-5,  1.311e-9),
    # Added 2026-08-04 for Problems 7.33 and 7.34, which ask for Illustrations 6.4-1,
    # 7.5-1 and 7.5-2 redone for nitrogen and for water. Without these, a student
    # following the illustrations gets the database's Reid-Prausnitz-Poling set and
    # silently works on a different basis from the book.
    #   nitrogen MATTERS: over 77-300 K the two sets differ by 149 J/mol in H and
    #     1.12 J/(mol K) in S (2.3% and 2.8%) -- the same size as the oxygen gap that
    #     prompted the 2026-08-03 decision.
    #   water does NOT: the two agree to 0.07%. Transcribed anyway so that Problem
    #     7.34 rests on the book's own appendix rather than on a coincidence.
    "nitrogen": (28.883, -0.157e-2,  0.808e-5, -2.871e-9),
    "water":    (32.218,  0.192e-2,  1.055e-5, -3.593e-9),
}


def load_unifac_subgroups():
    """subgroup_no, main_group_no, subgroup_name, main_group_name, R, Q (modified R/Q)."""
    return pd.read_csv(DATA_DIR / "unifac_subgroups.csv")


def load_unifac_interactions(kind="modified"):
    """Main-group interaction parameters. `modified`: a,b,c; `original`: a only."""
    if kind == "modified":
        return pd.read_csv(DATA_DIR / "unifac_interactions_modified.csv")
    if kind == "original":
        return pd.read_csv(DATA_DIR / "unifac_interactions_original.csv")
    raise ValueError(f"kind must be 'modified' or 'original', got {kind!r}")
