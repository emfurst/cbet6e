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

# ---------------------------------------------------------------------------
# CRYOGENIC RANGE, 100-700 K. New in the 6e (author decision 2026-08-04).
#
# WHY THIS EXISTS. Illustration 6.4-1 tabulates oxygen from -100 C (173.15 K), and
# ch7's Illustrations 7.5-1/7.5-2 and Problems 7.33/7.34 follow it. Appendix A.II's
# "Low Temperature Range" row is valid 273-1800 K, so all of them extrapolate it, and
# the extrapolation is not benign: at 173.15 K the A.II cubic gives Cp* = 27.883
# against a true 29.118 J/(mol K) -- 1.235 low, 4.2% -- which is worth 67 J/mol in
# H and 0.32 J/(mol K) in S. That dwarfs every other error in Table 6.4-4, including
# the 5e's own float error (0.6 J/mol) and the choice of gas constant (0.02 J/mol).
#
# A.II's row is NOT a bad fit; it is a fit to the wrong range. Refitting 273-1800 K
# reproduces its accuracy almost exactly (0.48 vs 0.40 J/(mol K) max error), so
# nothing is gained by re-deriving it. What was missing was a row for the range the
# illustrations use -- and Appendix A.II already splits combustion gases into
# "Low Temperature Range" and "High Temperature Range", so a third range is the
# appendix's own idiom, not a departure from it.
#
# PROVENANCE. These four numbers per gas are NOT transcribed from a published table --
# they were FITTED. The fit is reproduced, and checked against the values below, in
#     code/ch3/Heat_capacity_range_of_validity.ipynb
# which is the audit trail: it is also the ch3 exercise on why a correlation's range is
# part of its data. If you doubt a coefficient, run that notebook.
#
# Least squares, in Appendix A.II's own cubic form, over 100-700 K -- which for oxygen
# is exactly the span of JANAF's first Shomate segment, so the fit never straddles a
# discontinuity in the reference. Reference data:
#
#   Chase, M. W., Jr., NIST-JANAF Thermochemical Tables, Fourth Edition,
#   J. Phys. Chem. Ref. Data Monograph 9 (1998). Shomate coefficients retrieved from
#   the NIST Chemistry WebBook, webbook.nist.gov (O2 species C7782447, N2 species
#   C7727379); data last reviewed March 1977, parameter fit January 2009.
#     O2  100-700 K:  A=31.32234  B=-20.23531  C=57.86644  D=-36.50624  E=-0.007374
#     N2  100-500 K:  A=28.98641  B=1.853978   C=-9.647459 D=16.63537   E=0.000117
#     N2  500-2000 K: A=19.50583  B=19.88705   C=-8.598535 D=1.369784   E=0.527601
#   Cp* = A + B t + C t^2 + D t^3 + E/t^2 with t = T/1000, J/(mol K).
#
# Accuracy over the span the illustrations need, against that reference:
#
#            max |Cp* err|      max |H err|            max |S err|
#   oxygen   0.173 J/(mol K)    2.8 J/mol  (was 67.1)  0.011  (was 0.323)
#   nitrogen 0.056 J/(mol K)    2.0 J/mol  (was 43.7)  0.009  (was 0.263)
#
# ⚠️ THESE ROWS ARE NOT DROP-IN REPLACEMENTS -- they are worse than useless above 700 K.
# Fitted over 100-700 K, the oxygen cubic returns Cp* = -6.5 J/(mol K) at 1800 K, a
# negative heat capacity. Below 700 K it beats the printed row everywhere, including by
# 8x inside the 273-700 K overlap (0.046 vs 0.397). The rule: below 700 K use CRYO,
# above 700 K use APPENDIX_A2_CP, whose fit over 700-1800 K is excellent (0.063).
#
# Water is deliberately absent: it is not a cryogenic fluid (it freezes at 273 K), so
# nothing in the book needs it below A.II's range, and the two agree to 0.07% anyway.
#
# ✅ IN PRINT as of 2026-08-04: the author added `deliverables/manuscript/bapp01.docx`
# and put these two rows into the "Combustion Gases (Low Temperature Range)" block,
# each one above its 273-1800 K sibling and distinguished by the Temperature Range
# column. No new block heading -- A.II's range column does the work.
#
# THE VALUES BELOW ARE EXACTLY WHAT THE APPENDIX PRINTS, at its own 3-decimal
# precision, NOT the full-precision fit. That is deliberate and it is the whole point:
# a reader who types the printed row into these notebooks must reproduce the printed
# Table 6.4-4. Rounding the fit to the printed precision costs at most 0.0023 J/(mol K)
# in Cp*, which is nothing beside the 0.17 the row already carries.
#
# ✅ VERIFIED AGAINST THE APPENDIX 2026-08-05: all eight coefficients below match what
# bapp01.docx prints, to the last digit. (A 234 -> 293 transposition in the nitrogen `a`
# was caught and fixed; it would have been a uniform +0.059 J/(mol K) offset in Cp*,
# worth -7.4 J/mol in H at 173 K -- four times larger than the 2 J/mol this row exists
# to achieve. Re-check this equality if either side is ever edited: it is what lets a
# reader type the printed row in here and reproduce Table 6.4-4.)
APPENDIX_A2_CP_CRYO = {
    "oxygen":   (30.171, -1.293e-2,  4.236e-5, -2.5828e-8),
    "nitrogen": (29.234, -0.102e-2,  0.025e-5,  6.339e-9),
}


# --- the book's own critical constants (SIS Table 6.6-1) --------------------
#
# WHY THIS EXISTS. `pure_property.csv` is the Reid-Prausnitz-Poling table, and for
# oxygen it does NOT agree with the book:
#
#             Table 6.6-1        pure_property.csv
#   Tc        154.6 K            154.6 K        same
#   Pc        5.046 MPa          50.4 bar       0.12% low
#   omega     0.021              0.025          -> kappa 0.4069 vs 0.4130
#
# That is a different EOS, not a rounding, and omega = 0.021 is the value that makes
# Eq. 6.7-4 give the kappa = 0.4069 Illustration 6.4-1 prints. So
# `PengRobinson.from_database("oxygen")` silently works on a different basis from the
# illustrations, which is why the ch6 `_thermo` twins used to disagree with their
# self-contained partners.
#
# Use these when reproducing anything the book prints:
#
#     pr = PengRobinson(**TABLE_6_6_1["oxygen"], cp=APPENDIX_A2_CP_CRYO["oxygen"])
#
# Pc in Pa, Tc in K, omega dimensionless. Only the compounds the notebooks need are
# transcribed; add rows from Table 6.6-1 as later chapters require them.
TABLE_6_6_1 = {
    "oxygen":   dict(Tc=154.6, Pc=5.046e6, omega=0.021, name="oxygen"),
    "nitrogen": dict(Tc=126.2, Pc=3.394e6, omega=0.040, name="nitrogen"),
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
