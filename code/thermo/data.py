"""Data loaders for the `thermo` package — reads the CSV tables in `code/data/`."""
import re
from pathlib import Path
import numpy as np
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
    # Added 2026-08-07 for the recomputed Fig. 3.3-2 (methane P-H chart). NOTE the
    # range: A.II prints methane under "Paraffinic Hydrocarbons" as 273-1500 K, not
    # 273-1800 K. It is excellent there -- 35.68 J/(mol K) at 298 K against a true
    # 35.65 -- and catastrophic below it. The chart starts at 111 K. Use the CRYO row.
    "methane":  (19.875,  5.021e-2,  1.268e-5, -11.004e-9),
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
# THESE ROWS ARE NOT DROP-IN REPLACEMENTS -- they are worse than useless above 700 K.
# Fitted over 100-700 K, the oxygen cubic returns Cp* = -6.5 J/(mol K) at 1800 K, a
# negative heat capacity. Below 700 K it beats the printed row everywhere, including by
# 8x inside the 273-700 K overlap (0.046 vs 0.397). The rule: below 700 K use CRYO,
# above 700 K use APPENDIX_A2_CP, whose fit over 700-1800 K is excellent (0.063).
#
# Water is deliberately absent: it is not a cryogenic fluid (it freezes at 273 K), so
# nothing in the book needs it below A.II's range, and the two agree to 0.07% anyway.
#
# IN PRINT as of 2026-08-04: the author added `deliverables/manuscript/bapp01.docx`
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
# VERIFIED AGAINST THE APPENDIX 2026-08-05: all eight coefficients below match what
# bapp01.docx prints, to the last digit. (A 234 -> 293 transposition in the nitrogen `a`
# was caught and fixed; it would have been a uniform +0.059 J/(mol K) offset in Cp*,
# worth -7.4 J/mol in H at 173 K -- four times larger than the 2 J/mol this row exists
# to achieve. Re-check this equality if either side is ever edited: it is what lets a
# reader type the printed row in here and reproduce Table 6.4-4.)
APPENDIX_A2_CP_CRYO = {
    "oxygen":   (30.171, -1.293e-2,  4.236e-5, -2.5828e-8),
    "nitrogen": (29.234, -0.102e-2,  0.025e-5,  6.339e-9),
    # ---------------------------------------------------------------------------
    # METHANE, 100-500 K. Added 2026-08-07 for the recomputed Fig. 3.3-2, author's
    # decision of the same day. NOT YET IN PRINT -- Appendix A.II must gain this
    # row as it gained the two above, or the methane chart rests on constants the
    # book does not contain.
    #
    # WHY IT IS THE WORST CASE IN THE BOOK. A.II's printed methane row is quoted
    # 273-1500 K. Extrapolated to 100 K it returns Cp* = 25.01 J/(mol K) against a
    # true 33.26 -- 25% low, and BELOW THE PHYSICAL FLOOR: a nonlinear polyatomic
    # cannot go under 4R = 33.258 J/(mol K), which is what methane's Cp* tends to
    # once its vibrations freeze out. Over 100-400 K, the span the chart covers,
    # that costs -635 J/mol in H and -4.5 J/(mol K) in S -- 8% of methane's entire
    # heat of vaporization, and six times the oxygen gap that prompted these rows.
    #
    # ITS PROVENANCE IS DIFFERENT FROM THE TWO ABOVE, AND BETTER. Those were
    # fitted to NIST-JANAF Shomate coefficients. THAT ROUTE DOES NOT EXIST HERE:
    # the WebBook carries no sub-298 K Shomate segment for methane (checked
    # 2026-08-07 -- JANAF gives 298-1300 K and 1300-6000 K only). It is not needed,
    # because methane's ideal-gas Cp* is DERIVABLE:
    #
    #     Cp* = 4R + R * sum_i g_i u_i^2 e^u_i / (e^u_i - 1)^2,  u_i = h c nu_i / kT
    #
    # -- translation (3/2 R) + rotation (3/2 R, nonlinear) + R, plus the harmonic
    # contribution of CH4's nine normal modes: nu = 2917 (x1), 1534 (x2), 3019 (x3),
    # 1306 (x3) cm^-1. That expression reproduces the JANAF Shomate fit to
    # +/-0.03 J/(mol K) over 298-1300 K, which is the validation; extending it BELOW
    # 298 K is then vibrational freeze-out toward an exact limit, not the
    # extrapolation of a curve fit. The derivation, the JANAF comparison and the fit
    # that produced the coefficients below are in
    #     code/ch3/PH_charts_methane_and_nitrogen.ipynb   (section 1)
    # which asserts them against this row, so a changed coefficient fails there. (The
    # O2 and N2 rows are audited in code/ch3/Heat_capacity_range_of_validity.ipynb
    # instead, because they were fitted to JANAF rather than derived.)
    #
    # RANGE. 100-500 K, not the 100-700 K of the two rows above: over 100-700 the
    # cubic reaches max |dCp*| = 0.496 J/(mol K) and +6.6 J/mol, where 100-500 gives
    # 0.362 and +2.2 -- matching the oxygen row already in print. A.II's own
    # Temperature Range column distinguishes them, so a third range is the appendix's
    # idiom rather than a departure from it. Methane's Cp* rises steeply past 500 K,
    # which is what a cubic anchored at 100 K cannot follow.
    #
    # Accuracy over 100-400 K, against the derivation above:
    #   max |dCp*| 0.362 J/(mol K)   dH +2.2 J/mol (was -635)   dS +0.010 (was -4.52)
    "methane":  (37.097, -5.282e-2, 18.976e-5, -9.2669e-8),
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
    # Added 2026-08-07 for Fig. 3.3-2. Table 6.6-1 and pure_property.csv disagree
    # here too, in the same direction as they do for oxygen and by about as much:
    #   Table 6.6-1   190.6 K   4.600 MPa   omega 0.008
    #   CSV (RPP)     190.4 K   4.60  MPa   omega 0.011
    # omega 0.008 vs 0.011 moves kappa by 0.004, which is small -- but the rule
    # stands: the printed figures use the printed constants.
    "methane":  dict(Tc=190.6, Pc=4.600e6, omega=0.008, name="methane"),
}


# --- Table 9.4-1: Peng-Robinson binary interaction parameters ---------------
#
# 127 pairs over 20 species, digitized from the printed Table 9.4-1 (5e p. 441).
# Stored as an upper-triangular edge list, not a matrix: the table is 65% blank, and
# a blank is *not* zero.
#
# THE BLANKS ARE THE POINT. Table 9.4-1's own footnote reads: "Blanks indicate no
# data are available from which the k12 could be evaluated. In such case use estimates
# from mixtures of similar compounds." A k_ij silently defaulted to zero is a
# different mixture, not a missing decimal -- so `pr_kij_matrix` returns the pairs it
# could not find alongside the matrix, and `PRMixture.from_database` warns rather than
# filling in quietly.
#
# Verified against the book's own three worked values: k = 0.010 for ethane/n-butane
# (Illustration 9.4-3), 0.09 for methane/carbon dioxide (Illustration 9.4-4), and
# 0.018 for n-pentane/benzene (Illustration 9.4-5). All 127 entries land in the upper
# triangle, and the parse places exactly as many numbers as the printed row contains.
def load_pr_kij():
    """Table 9.4-1 as a tidy edge list: species_i, species_j, formula_i, formula_j, kij."""
    return pd.read_csv(DATA_DIR / "pr_kij.csv")


def pr_kij_matrix(keys, strict=False):
    """Symmetric k_ij matrix for `keys`, plus the pairs Table 9.4-1 does not give.

    `keys` are matched against both the species names and the formulas as the table
    prints them, case-insensitively -- so ["methane", "CO2"] and ["CH4", "carbon
    dioxide"] both work.

    Returns
    -------
    kij : (n, n) array, zero on the diagonal and for any pair the table lacks.
    missing : list of (key_i, key_j) the table has no value for.

    Set `strict=True` to raise instead of returning `missing`.
    """
    df = load_pr_kij()
    lookup = {}
    for _, row in df.iterrows():
        a, b = row.species_i.lower(), row.species_j.lower()
        fa, fb = row.formula_i.lower(), row.formula_j.lower()
        for u in (a, fa):
            for v in (b, fb):
                lookup[(u, v)] = float(row.kij)
                lookup[(v, u)] = float(row.kij)
    n = len(keys)
    kij = np.zeros((n, n))
    missing = []
    for i in range(n):
        for j in range(i + 1, n):
            ki, kj = str(keys[i]).strip().lower(), str(keys[j]).strip().lower()
            if (ki, kj) in lookup:
                kij[i, j] = kij[j, i] = lookup[(ki, kj)]
            else:
                missing.append((keys[i], keys[j]))
    if missing and strict:
        raise KeyError(f"Table 9.4-1 has no k_ij for {missing}")
    return kij, missing


# --- Table 9.5-2: UNIFAC group volume and surface area parameters ------------
#
# 92 subgroups over 46 main-group names, digitized from the printed Table 9.5-2
# (5e pp. 457-458). The book's own table is the source of truth for R, Q and the group
# inventory; the subgroup and main-group *numbers* are not printed anywhere in the
# chapter and come from the legacy `UNIFAC_data.mat` extraction, because
# `unifac_interactions_modified.csv` is keyed on the main-group number and the two
# files have to agree. That asymmetry is recorded in code/data/README.md.
#
# Verified: the group sums reproduce Illustration 9.5-2 exactly -- benzene as 6 ACH
# gives r = 2.2578 and q = 2.5926, and 2,2,4-trimethyl pentane as 5 CH3 + CH2 + CH + C
# gives r = 5.0600 and q = 6.3675.
#
# That is true of the 6e, which recomputed the illustration on this table. The 5e's
# version of it printed 3.1878 / 2.4000 and 5.8463 / 5.0080 -- ORIGINAL-UNIFAC sums --
# while citing Table 9.5-2. Those 5e values are the source of
# code/data/unifac_subgroups_original.csv; see
# ch10/validation/unifac_original_vs_dortmund_validation.ipynb.
def _check_subgroups(df):
    """Integrity rules the table must satisfy. A silent violation here changes every
    activity coefficient, so it is checked on every load rather than documented.

    The legacy file this replaced violated the first rule: subgroup numbers 37, 38 and
    39 each appeared twice -- once for the original-UNIFAC pyridine subgroups
    (C5H5N, C5H4N, C5H3N) and once for the Dortmund ones (AC2H2N, AC2HN, AC2N) -- so
    `dict(zip(subgroup_no, R))` silently kept whichever row came last.
    """
    if not df.subgroup_no.is_unique:
        dup = sorted(df.loc[df.subgroup_no.duplicated(keep=False), "subgroup_no"].unique())
        raise ValueError(f"unifac_subgroups.csv: duplicate subgroup_no {dup} -- "
                         f"a dict keyed on it would silently drop rows")
    bad = df.groupby("main_group_no").main_group_name.nunique()
    if (bad > 1).any():
        raise ValueError(f"unifac_subgroups.csv: main_group_no with more than one name: "
                         f"{list(bad[bad > 1].index)}")
    if (df.R <= 0).any():
        raise ValueError("unifac_subgroups.csv: R must be positive")
    if (df.Q < 0).any():
        raise ValueError("unifac_subgroups.csv: Q must be non-negative")
    return df


def load_unifac_subgroups(kind="modified"):
    """Subgroup volume and surface-area parameters.

    `modified` is Table 9.5-2 in full: subgroup_no, main_group_no, subgroup_name,
    main_group_name, R, Q, example. This is the book's model (revision_notes/c09.md D1).

    `original` is Table 7.5-2 of **this book's 2nd edition** in full: 85 subgroups over
    44 main groups, which is exactly the main-group set
    `unifac_interactions_original.csv` is keyed on. It replaced a six-row stub on
    2026-08-17; those six rows, recovered from the 5e's Illustration 9.5-2, now serve
    as an independent check on the transcription and agree to the last digit.

    The two tables share subgroup numbers 1-26 and 40-77 and **disagree on 27,
    37-39 and 78-85** -- 27 is `FCH2O` here and `cy-CH2 OCH2` there. Both numbers
    exist in both tables, so a group assignment carried across kinds does not raise.
    Use `unifac_groups(name, kind)`.

    **Never mix the two sets.** The R/Q and the a_mn must come from the same
    parameter set; the 6e's own footnote to Table 9.5-2 says so. `UNIFAC(kind=...)`
    pairs them, which is why this loader is not meant to be called directly.
    """
    if kind == "modified":
        return _check_subgroups(pd.read_csv(DATA_DIR / "unifac_subgroups.csv"))
    if kind == "original":
        return _check_subgroups(pd.read_csv(DATA_DIR / "unifac_subgroups_original.csv"))
    raise ValueError(f"kind must be 'modified' or 'original', got {kind!r}")


def unifac_groups(name, kind="modified"):
    """Look up the group assignment the parameter set's Example Assignments column gives.

    Returns {subgroup_no: count} ready for `UNIFAC.gamma`, so a notebook can write
    `unifac_groups("benzene")` instead of hand-transcribing `{9: 6}` and risking the
    wrong subgroup number. Only the species the table uses as examples are available;
    anything else has to be assigned by hand from the table.

    **`kind` must match the `UNIFAC(kind=...)` the result is handed to.** The two
    parameter sets share subgroup numbers 1-26 and 40-77, but **27, 37-39 and 78-85
    name different groups in each** -- 27 is the original set's `FCH2O` and Dortmund's
    `cy-CH2 OCH2`, 78-80 are `SiH3`/`SiH2`/`SiH` against `cy-CH2`/`cy-CH`/`cy-C`. Both
    numbers exist in both tables, so a mismatched assignment does not raise; it returns
    a plausible wrong answer. Tetrahydrofuran is the trap in miniature: original UNIFAC
    calls it 1 FCH2O + 3 CH2, Dortmund 1 cy-CH2 OCH2 + 2 cy-CH2 -- different numbers,
    different counts, same species.
    """
    df = load_unifac_subgroups(kind)
    key = str(name).strip().lower()
    num = dict(zip(df.subgroup_name.str.replace(" ", ""), df.subgroup_no))
    for ex in df.example.dropna():
        if ":" not in ex:
            if ex.strip().lower() == key:            # e.g. "Water", "Chloroform"
                row = df[df.example == ex].iloc[0]
                return {int(row.subgroup_no): 1}
            continue
        species, assignment = ex.split(":", 1)
        if species.strip().lower() != key:
            continue
        groups = {}
        for part in assignment.split(","):
            part = part.strip()
            # Table 9.5-2 omits the count when it is 1 ("Methylamine: CH3 NH2",
            # "Dimethylamine: CH3 NH, 1 CH3"), so a missing leading integer means one.
            m = re.match(r"(?:(\d+)\s*)?(.+)", part)
            if not m:
                raise ValueError(f"cannot parse {part!r} in Table 9.5-2 example {ex!r}")
            count = int(m.group(1)) if m.group(1) else 1
            sub = m.group(2).replace(" ", "")
            if sub not in num:
                raise KeyError(f"example {ex!r} names subgroup {sub!r}, which is not a "
                               f"row of the {kind!r} table")
            groups[int(num[sub])] = groups.get(int(num[sub]), 0) + count
        return groups
    raise KeyError(f"{name!r} is not an example species in the {kind!r} table")


def load_unifac_interactions(kind="modified"):
    """Main-group interaction parameters. `modified`: a,b,c; `original`: a only."""
    if kind == "modified":
        return pd.read_csv(DATA_DIR / "unifac_interactions_modified.csv")
    if kind == "original":
        return pd.read_csv(DATA_DIR / "unifac_interactions_original.csv")
    raise ValueError(f"kind must be 'modified' or 'original', got {kind!r}")


# --- reaction species: Appendix A.IV formation data + Appendix A.II Cp* ------
#
# WATCH THE SCALING. `react_property.csv` stores the heat-capacity coefficients in
# the form Appendix A.II *prints* them, which is not the form you evaluate:
#
#     Cp* = A + B*1e-2 T + C*1e-5 T^2 + D*1e-9 T^3 + E/T^2       J/(mol K)
#
# Evaluating the raw columns gives Cp values wrong by three orders of magnitude,
# silently. `REACT_CP_SCALE` and `reaction_cp` below are the only places that
# factor should appear.
#
# Verified against the book's own arithmetic: Illustration 8.5-2 prints the
# combined coefficients for 2 NO2 - N2O4 as 12.804, -7.239e-2, 4.301e-5,
# 1.5732e-8, and the scaled columns reproduce all four exactly -- and then the
# printed heats of reaction at 200, 300, 400, 500 and 600 K to the last digit.
# (E is unscaled; it is nonzero only for the fourteen solid species.)
REACT_CP_SCALE = np.array([1.0, 1e-2, 1e-5, 1e-9, 1.0])


def load_reaction_properties():
    """The 99-species reaction table: Name, DG, DH (kJ/mol), A..E, ID.

    Appendix A.IV (standard Gibbs energies and enthalpies of formation at 25 C,
    1 bar) joined to Appendix A.II (ideal-gas heat capacities). See
    `reaction_cp` before using columns A..E.
    """
    return pd.read_csv(DATA_DIR / "react_property.csv")


def get_reaction_species(name):
    """One row of `react_property.csv` by exact `Name` (e.g. 'N2O4', 'H2O(g)')."""
    df = load_reaction_properties()
    hit = df[df["Name"] == str(name)]
    if not len(hit):
        hit = df[df["Name"].str.lower() == str(name).strip().lower()]
    if not len(hit):
        raise KeyError(f"species {name!r} not found in react_property.csv "
                       f"(names are formulas, e.g. 'N2O4', 'C6H6', 'H2O(g)')")
    return hit.iloc[0]


def reaction_cp(name):
    """Scaled ideal-gas Cp* coefficients (a, b, c, d, e) for one species, SI.

    Returns the coefficients already multiplied by `REACT_CP_SCALE`, so

        a, b, c, d, e = reaction_cp("NO2")
        Cp = a + b*T + c*T**2 + d*T**3 + e/T**2        # J/(mol K)
    """
    row = get_reaction_species(name)
    return np.array([float(row[k]) for k in "ABCDE"]) * REACT_CP_SCALE


def formation_enthalpy(name):
    """Standard enthalpy of formation at 25 C, 1 bar, in **J/mol** (CSV is kJ/mol)."""
    return float(get_reaction_species(name)["DH"]) * 1e3


def formation_gibbs(name):
    """Standard Gibbs energy of formation at 25 C, 1 bar, in **J/mol**."""
    return float(get_reaction_species(name)["DG"]) * 1e3


# --- mixing data (SIS Sec. 8.6) ---------------------------------------------

MIXING_DATASETS = {
    # key: (file, what it is, the identity that checks it)
    "water-methanol-volume": "mixing_water_methanol_volume.csv",
    "water-methanol-enthalpy": "mixing_water_methanol_enthalpy.csv",
}


def load_mixing_data(key):
    """A property-change-on-mixing data set from Sec. 8.6, in SI units.

    `water-methanol-volume`   -- Table 8.6-1: x1, rho (kg/m^3), V and dmixV (m^3/mol)
    `water-methanol-enthalpy` -- Table 8.6-3: x1, Q+ and dmixH (J/mol)

    Species 1 is water and species 2 is methanol in both. Note the book prints
    the volumes multiplied by 1e6 and the enthalpies in kJ/mol; the CSVs are SI,
    like everything else in the package.
    """
    try:
        fname = MIXING_DATASETS[key]
    except KeyError:
        raise KeyError(f"unknown mixing data set {key!r}; "
                       f"choose from {sorted(MIXING_DATASETS)}") from None
    return pd.read_csv(DATA_DIR / fname)
