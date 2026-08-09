"""Two-state (native <-> unfolded) thermal model, after Li et al. 2021 and the
MRes (Madkaikar 2023, github.com/adimadkaikar/MResProject, code/etcpy/etc.py).

In this model the high-temperature collapse of an enzyme's usable activity is set
by *denaturation* — the fraction of protein in the native (folded) state, keyed on
the melting temperature Tm — rather than by the curvature of a fictional MMRT
heat-capacity term. The effective per-flux enzyme cost is

    cost_i(T) = base_cost_i / ( rel_kcat_i(T) * f_N_i(T) )

where ``rel_kcat_i(T)`` is the transition-state (MMRT) turnover factor anchored at
the enzyme optimum Topt (Eyring form; NOT peak-normalised), and ``f_N_i(T)`` in
[0, 1] is the native fraction from a two-state unfolding equilibrium that goes to 0
above Tm, so cost inflates sharply near Tm and produces the falling limb.

All functions are direct ports of the reference ``etcpy/etc.py`` equations and use
its SI convention: energies in J/mol, entropies/heat-capacities in J/mol/K,
R = 8.314 J/mol/K, temperatures in K. They are vectorised over enzymes (arrays of
per-enzyme parameters at a scalar temperature T) for the sweep hot path.
"""
from __future__ import annotations

import numpy as np

# Convergence temperatures and fit constants from etcpy/etc.py (Li et al. 2021).
R_J = 8.314           # J / mol / K
TH = 373.5            # K, enthalpy convergence temperature
TS = 385.0            # K, entropy convergence temperature
T0_K = 273.15         # K, kcat reference temperature (0 C) used by calculate_kcatT
_SLOPE = 299.58
_INTERCEPT = 20008.0


# ---------------------------------------------------------------------------
# unfolding thermodynamics: (dHTH, dSTS, dCpu) from Tm and T90 or protein length
# ---------------------------------------------------------------------------
def dH_dS_dCpu_from_Tm_T90(Tm: float, T90: float):
    """Solve (dHTH, dSTS, dCpu) from Tm and T90 (both K).

    Solves the linear system encoding dHTH = slope*dSTS + intercept,
    dGu(Tm) = 0 and dGu(T90) = -R*T90*ln 9 (port of etc.py)."""
    a = np.array([[1, -_SLOPE, 0],
                  [1, -Tm, Tm - TH - Tm * np.log(Tm / TS)],
                  [1, -T90, T90 - TH - T90 * np.log(T90 / TS)]])
    b = np.array([_INTERCEPT, 0.0, -R_J * T90 * np.log(9)])
    dHTH, dSTS, dCpu = np.linalg.solve(a, b)
    return float(dHTH), float(dSTS), float(dCpu)


def dH_dS_dCpu_from_Tm_length(Tm: float, N: float):
    """Solve (dHTH, dSTS, dCpu) from Tm (K) and protein length N (residues).

    Uses the Sawle & Ghosh (2011) length relations for dHTH, dSTS, then solves
    dGu(Tm) = 0 for dCpu (port of etc.py get_dH_dS_dCpu_from_TmLength)."""
    from scipy.optimize import fsolve
    dHTH = (4 * N + 143) * 1000.0
    dSTS = 13.27 * N + 448.0

    def func(dCp):
        return dHTH + dCp * (Tm - TH) - Tm * dSTS - Tm * dCp * np.log(Tm / TS)

    dCpu = float(fsolve(func, 10000.0)[0])
    return float(dHTH), float(dSTS), dCpu


def unfolding_params(Tm: float, T90: float = None, length: float = None):
    """(dHTH, dSTS, dCpu) for one enzyme with the etc.py fallback logic.

    Prefer the Tm/T90 solution when a valid T90 (> Tm) is available and yields a
    non-negative dCpu; otherwise use the protein-length relation."""
    if T90 is not None and np.isfinite(T90) and T90 > Tm:
        dHTH, dSTS, dCpu = dH_dS_dCpu_from_Tm_T90(Tm, T90)
        if dCpu >= 0:
            return dHTH, dSTS, dCpu
    if length is None or not np.isfinite(length):
        length = 300.0   # generic bacterial protein length fallback
    return dH_dS_dCpu_from_Tm_length(Tm, length)


# ---------------------------------------------------------------------------
# native fraction f_N(T) (vectorised)
# ---------------------------------------------------------------------------
def dGu(T, dHTH, dSTS, dCpu):
    """Gibbs free energy of unfolding at T (K). Arrays broadcast over enzymes."""
    return dHTH + dCpu * (T - TH) - T * dSTS - T * dCpu * np.log(T / TS)


def native_fraction(T, dHTH, dSTS, dCpu):
    """Native (folded) fraction f_N(T) = 1/(1 + exp(-dGu/RT)) in [0, 1]."""
    g = dGu(T, dHTH, dSTS, dCpu)
    return 1.0 / (1.0 + np.exp(-g / (R_J * T)))


# ---------------------------------------------------------------------------
# transition-state kcat(T) factor, anchored at Topt (vectorised)
# ---------------------------------------------------------------------------
def rel_kcat(T, dHTH, dSTS, dCpu, dCpt, Topt):
    """Relative turnover kcat(T) / kcat(Topt) from transition-state theory.

    Direct port of etc.py ``calculate_kcatT`` with kcatTopt = 1 (the factor is
    linear in kcatTopt, so this returns kcat(T)/kcat(Topt)). ``dCpt`` is the
    transition-state heat capacity (J/mol/K); (dHTH, dSTS, dCpu) parameterise the
    unfolding term used in the activation enthalpy. Arrays broadcast over enzymes;
    T is a scalar in K. NOT peak-normalised — kcat(Topt) is the anchor."""
    dGuTopt = dHTH + dCpu * (Topt - TH) - Topt * dSTS - Topt * dCpu * np.log(Topt / TS)
    dHt = (dHTH + dCpu * (Topt - TH) - dCpt * (Topt - T0_K) - R_J * Topt
           - (dHTH + dCpu * (Topt - TH)) / (1.0 + np.exp(-dGuTopt / (R_J * Topt))))
    # kcat0 = kcat(T0) implied by anchoring kcat(Topt) = 1
    lnk0 = (np.log(Topt / T0_K) - (dHt + dCpt * (Topt - T0_K)) / (R_J * Topt)
            + dHt / (R_J * T0_K) + dCpt * np.log(Topt / T0_K) / R_J)
    lnkT = (np.log(T / T0_K) - (dHt + dCpt * (T - T0_K)) / (R_J * T)
            + dHt / (R_J * T0_K) + dCpt * np.log(T / T0_K) / R_J)
    return np.exp(lnkT - lnk0)


# ---------------------------------------------------------------------------
# temperature-dependent non-growth-associated maintenance (NGAM/ATPM)
# ---------------------------------------------------------------------------
def ngam_T(T, scale: float = 1.0, steepness: float = 1.0):
    """Temperature-dependent maintenance ATP demand at T (K), with a basal floor
    at 25 C (port of etc.py getNGAMT). Returns mmol ATP / gDW / h.

    ``scale`` multiplies the ~8.5 baseline amplitude; ``steepness`` multiplies the
    activation term (how fast maintenance rises with temperature) -- the
    peak-rounding lever. Both default to the ported Li getNGAMT values (1.0)."""
    def _f(Tk):
        return scale * 8.5 * (1.0 - 0.62 * np.exp((-0.5 * steepness / (8.617e-5))
                                                  * (1.0 / (273.15 + 25.0) - 1.0 / Tk)))
    val = _f(T)
    floor = _f(273.15 + 25.0)
    return float(max(val, floor))
