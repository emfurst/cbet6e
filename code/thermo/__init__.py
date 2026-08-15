"""thermo — the book's models in Python: the cubic equations of state, the mixture
models, the activity coefficient models, UNIFAC, and the property charts.

    from thermo import PengRobinson, VanDerWaals, UNIFAC, VanLaar, NRTL

See code/thermo/README.md. Reads reference data from code/data/.

⚠️ **Two different UNIQUACs.** `thermo.UNIFAC` is the group-contribution model of
Sec. 9.6, which takes group counts. `thermo.UNIQUACModel` is the molecular model of
Sec. 9.5, which takes r and q -- exported under that name so the two do not collide;
it is `thermo.activity_models.UNIQUAC` in its own module.

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
from .activity_models import (ActivityModel, OneConstantMargules,
                              TwoConstantMargules, RedlichKisterGex, VanLaar,
                              Wilson, NRTL, FloryHuggins, RegularSolution,
                              fit_binary, TABLE_9_5_1_VAN_LAAR, TABLE_9_6_1)
from .activity_models import UNIQUAC as UNIQUACModel
from .electrolytes import (Electrolyte, DebyeHuckel, ionic_strength,
                           water_parameters, TABLE_9_10_1, ELECTROLYTES)
from .wong_sandler import (WongSandler, GexFromUNIFAC, C_STAR_PR, C_STAR_VDW)
from .vle import (VaporPressure, Antoine, ClausiusClapeyron, TabulatedPsat,
                  Ideal, GammaPhi, pxy, txy, azeotrope)
from .fitting import (RedlichKister, tangent_intercepts, LogPolynomialActivity,
                      ILLUSTRATION_9_7_1_HEMOGLOBIN, SIS_9_7_1_CONSTANTS)
from .data import (get_compound, load_pure_properties, APPENDIX_A2_CP,
                   APPENDIX_A2_CP_CRYO, TABLE_6_6_1,
                   load_mixing_data, load_reaction_properties, get_reaction_species,
                   load_pr_kij, pr_kij_matrix)

_LAZY = {"charts", "ph_chart", "steam_chart", "vle_chart"}


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
           "ActivityModel", "OneConstantMargules", "TwoConstantMargules",
           "RedlichKisterGex", "VanLaar", "Wilson", "NRTL", "FloryHuggins",
           "UNIQUACModel", "RegularSolution", "fit_binary",
           "Electrolyte", "DebyeHuckel", "ionic_strength", "water_parameters",
           "WongSandler", "GexFromUNIFAC", "C_STAR_PR", "C_STAR_VDW",
           "VaporPressure", "Antoine", "ClausiusClapeyron", "TabulatedPsat",
           "Ideal", "GammaPhi", "pxy", "txy", "azeotrope",
           "TABLE_9_10_1", "ELECTROLYTES",
           "TABLE_9_5_1_VAN_LAAR", "TABLE_9_6_1",
           "RedlichKister", "tangent_intercepts", "LogPolynomialActivity",
           "ILLUSTRATION_9_7_1_HEMOGLOBIN", "SIS_9_7_1_CONSTANTS",
           "get_compound", "load_pure_properties", "APPENDIX_A2_CP",
           "APPENDIX_A2_CP_CRYO", "TABLE_6_6_1",
           "load_mixing_data", "load_reaction_properties", "get_reaction_species",
           "load_pr_kij", "pr_kij_matrix",
           "charts", "ph_chart", "steam_chart", "vle_chart"]
