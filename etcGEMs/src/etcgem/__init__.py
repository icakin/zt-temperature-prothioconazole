"""etcgem: enzyme- and temperature-constrained GEM thermal-performance
sensitivity analysis.

Vary proteome allocation and per-enzyme kcat temperature responses (MMRT) to see
the range of organismal thermal performance curves a given genome can generate.
"""
from .mmrt import MMRTParams, kcat_curve, solve_dH_dS
from .enzyme_cost import (
    EnzymeEntry, EnzymeCostTable, EnzymeConstrainedModel, Perturbation,
)
from .tpc import TPC, TPCDescriptors, compute_tpc
from .sensitivity import run_sensitivity, SensitivityResult
from . import providers
from . import dltkcat
from .dltkcat import (
    fit_mmrt, fit_predictions, export_targets, apply_fits_to_provider,
    write_csv_table, build_dltkcat_input, parse_dltkcat_output, select_substrate,
)
from . import decomposition
from .decomposition import run_decomposition, decompose_grid, DecompositionResult
from . import control
from .control import run_control, ControlResult
from . import sectors
from .sectors import add_proteome_sectors
from . import thermal_sampling
from .thermal_sampling import sample_thermal, apply_thermal_sample

__all__ = [
    "MMRTParams", "kcat_curve", "solve_dH_dS",
    "EnzymeEntry", "EnzymeCostTable", "EnzymeConstrainedModel", "Perturbation",
    "TPC", "TPCDescriptors", "compute_tpc",
    "run_sensitivity", "SensitivityResult", "providers",
    "dltkcat", "fit_mmrt", "fit_predictions", "export_targets",
    "apply_fits_to_provider", "write_csv_table",
    "decomposition", "run_decomposition", "decompose_grid", "DecompositionResult",
    "control", "run_control", "ControlResult",
    "sectors", "add_proteome_sectors",
    "thermal_sampling", "sample_thermal", "apply_thermal_sample",
]
__version__ = "0.1.0"
