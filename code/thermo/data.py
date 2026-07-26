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
