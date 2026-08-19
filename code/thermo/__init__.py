"""thermo — the book's models in Python: the cubic equations of state, the mixture
models, the activity coefficient models, UNIFAC, and the property charts.

Sections 11.4 and 11.5 add two small modules with no solver in either of them:
`partition` (distribution coefficients, K_OW, the Gibbs energy of transfer) and
`osmotic` (Eq. 11.5-4 and the osmotic virial expansion). Chapter 12 adds `sle`, which
is one equation -- the equality of fugacities with a pure solid -- solved for the
composition in Sec. 12.1 and for the temperature in Sec. 12.3. Chapter 13 adds
`reaction`, which is where the equilibrium constant, the extent of reaction, the
Ellingham construction and Gibbs minimization live -- and which replaces the 5e's
CHEMEQ and its Appendix B.I/B.II equilibrium-constant programs.

    from thermo import PengRobinson, VanDerWaals, UNIFAC, VanLaar, NRTL

See code/thermo/README.md. Reads reference data from code/data/.

⚠️ **Two different UNIQUACs.** `thermo.UNIFAC` is the group-contribution model of
Sec. 9.6, which takes group counts. `thermo.UNIQUACModel` is the molecular model of
Sec. 9.5, which takes r and q -- exported under that name so the two do not collide;
it is `thermo.activity_models.UNIQUAC` in its own module.

The chart modules — `charts`, `ph_chart`, `steam_chart`, `vle_chart`, `ternary` —
are loaded **lazily**, so
the line above stays matplotlib-free and costs nothing extra in a notebook that only
wants the equation of state. Import them by name when you need them:

    from thermo.charts import use_book_style, label_at
    from thermo.ph_chart import ChartFluid, ph_chart
    from thermo.ternary import ternary_axes, tie_line, check_labels

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
from .lle import (binary_lle, binary_lle_envelope, consolute_temperature,
                  lle_flash, vlle_binary, immiscible_pressure,
                  eos_binary_lle, eos_vlle_binary)
from .partition import (gamma_ratio_from_Kx, x_from_concentration, Kx_from_Kc,
                        kow_from_gamma, gamma_from_kow, solute_split,
                        gamma_ratio_from_solubility, x_from_molarity,
                        gibbs_energy_of_transfer, air_water_partition,
                        compartment_partition, compartment_concentrations,
                        COMPARTMENTS, EQ_12_5_8_CONSTANT)
from .sle import (ln_x_gamma, ideal_solubility, solubility, activity_coefficient,
                  heat_of_sublimation, solubility_parameter, heat_of_fusion,
                  ideal_solubility_in_gas, poynting, enhancement_factor,
                  solubility_in_gas, freezing_point, freezing_point_depression,
                  eutectic, Iteration)
from .osmotic import (osmotic_pressure, osmotic_pressure_dilute,
                      pressure_from_height, molecular_weight, virial_fit,
                      VirialFit, multicomponent_dilute, electrolyte_mole_fraction)
from .fitting import (RedlichKister, tangent_intercepts, LogPolynomialActivity,
                      ILLUSTRATION_9_7_1_HEMOGLOBIN, SIS_9_7_1_CONSTANTS)
from .data import (get_compound, load_pure_properties, APPENDIX_A2_CP,
                   APPENDIX_A2_CP_CRYO, TABLE_6_6_1,
                   load_mixing_data, load_reaction_properties, get_reaction_species,
                   load_pr_kij, pr_kij_matrix, formation_enthalpy, formation_gibbs,
                   reaction_cp)
from .reaction import (Reaction, Extent, equilibrium_extent, multireaction_extents,
                       gibbs_minimization, gibbs_curve, ellingham, activity,
                       K_nu_from_eos, K_gamma_from_model, elements,
                       formation_gibbs_T)

_LAZY = {"charts", "ph_chart", "steam_chart", "vle_chart", "ternary"}


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
           "binary_lle", "binary_lle_envelope", "consolute_temperature",
           "lle_flash", "vlle_binary", "immiscible_pressure",
           "eos_binary_lle", "eos_vlle_binary",
           "gamma_ratio_from_Kx", "x_from_concentration", "Kx_from_Kc",
           "kow_from_gamma", "gamma_from_kow", "solute_split",
           "gamma_ratio_from_solubility", "x_from_molarity",
           "gibbs_energy_of_transfer", "air_water_partition",
           "compartment_partition", "compartment_concentrations",
           "COMPARTMENTS", "EQ_12_5_8_CONSTANT",
           "ln_x_gamma", "ideal_solubility", "solubility", "activity_coefficient",
           "heat_of_sublimation", "solubility_parameter", "heat_of_fusion",
           "ideal_solubility_in_gas", "poynting", "enhancement_factor",
           "solubility_in_gas", "freezing_point", "freezing_point_depression",
           "eutectic", "Iteration",
           "osmotic_pressure", "osmotic_pressure_dilute", "pressure_from_height",
           "molecular_weight", "virial_fit", "VirialFit",
           "multicomponent_dilute", "electrolyte_mole_fraction",
           "TABLE_9_10_1", "ELECTROLYTES",
           "TABLE_9_5_1_VAN_LAAR", "TABLE_9_6_1",
           "RedlichKister", "tangent_intercepts", "LogPolynomialActivity",
           "ILLUSTRATION_9_7_1_HEMOGLOBIN", "SIS_9_7_1_CONSTANTS",
           "get_compound", "load_pure_properties", "APPENDIX_A2_CP",
           "APPENDIX_A2_CP_CRYO", "TABLE_6_6_1",
           "load_mixing_data", "load_reaction_properties", "get_reaction_species",
           "load_pr_kij", "pr_kij_matrix", "formation_enthalpy", "formation_gibbs",
           "reaction_cp",
           "Reaction", "Extent", "equilibrium_extent", "multireaction_extents",
           "gibbs_minimization", "gibbs_curve", "ellingham", "activity",
           "K_nu_from_eos", "K_gamma_from_model", "elements", "formation_gibbs_T",
           "charts", "ph_chart", "steam_chart", "vle_chart", "ternary"]
