"""thermo — Aspen-optional Python substitutes for the cubic EOS models and UNIFAC.

    from thermo import PengRobinson, VanDerWaals, UNIFAC

See code/thermo/README.md. Reads reference data from code/data/.

The chart modules — `charts`, `ph_chart`, `steam_chart` — are loaded **lazily**, so
the line above stays matplotlib-free and costs nothing extra in a notebook that only
wants the equation of state. Import them by name when you need them:

    from thermo.charts import use_book_style, label_at
    from thermo.ph_chart import ChartFluid, ph_chart

`thermo.charts` also resolves as an attribute (`thermo.charts.label_at`), which is
what makes the lazy hook below worth having rather than just documenting the
submodule path.
"""
import importlib

from .cubic import CubicEOS
from .peng_robinson import PengRobinson
from .van_der_waals import VanDerWaals
from .pr_mixture import PRMixture
from .unifac import UNIFAC
from .fitting import RedlichKister, tangent_intercepts
from .data import (get_compound, load_pure_properties, APPENDIX_A2_CP,
                   APPENDIX_A2_CP_CRYO, TABLE_6_6_1,
                   load_mixing_data, load_reaction_properties, get_reaction_species)

_LAZY = {"charts", "ph_chart", "steam_chart"}


def __getattr__(name):
    """PEP 562 lazy submodule import — keeps matplotlib off the core import path."""
    if name in _LAZY:
        module = importlib.import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + list(_LAZY))


__all__ = ["CubicEOS", "PengRobinson", "VanDerWaals", "PRMixture", "UNIFAC",
           "RedlichKister", "tangent_intercepts",
           "get_compound", "load_pure_properties", "APPENDIX_A2_CP",
           "APPENDIX_A2_CP_CRYO", "TABLE_6_6_1",
           "load_mixing_data", "load_reaction_properties", "get_reaction_species",
           "charts", "ph_chart", "steam_chart"]
