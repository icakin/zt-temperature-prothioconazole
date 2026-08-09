"""Thermal performance curve (TPC) computation and descriptors.

A TPC here is organismal growth rate (biomass flux, 1/h) as a function of
temperature, obtained by re-solving the enzyme-constrained FBA at each
temperature after rescaling every kcat with its MMRT response. Descriptors
summarise the emergent curve so sensitivity sweeps can be reduced to a handful
of ecologically meaningful numbers.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, Optional, Sequence

import numpy as np

from .enzyme_cost import Perturbation

K_BOLTZMANN_EV = 8.617333262e-5   # eV / K
R_KJ = 8.314462618e-3             # kJ / mol / K


def compute_tpc(pm, temps_C: Sequence[float], pert: Optional[Perturbation] = None,
                min_growth: float = 1e-6) -> "TPC":
    """Growth rate vs temperature for one parameter point.

    Parameters
    ----------
    pm       : ProvidedModel
    temps_C  : temperatures in degrees Celsius
    pert     : Perturbation applied to the whole enzyme table (default = none)
    """
    pert = pert or Perturbation()
    ecm = pm.ec
    budget = pert.budget if pert.budget is not None else ecm.default_budget
    use_alloc = pert.uses_allocation()
    temps_C = np.asarray(temps_C, dtype=float)
    growth = np.zeros_like(temps_C)
    for i, Tc in enumerate(temps_C):
        Tk = Tc + 273.15
        ecm.set_temperature(Tk, pert)
        if getattr(ecm, "_alloc_from_data", None) is not None and ecm._sectors is not None:
            # temperature-dependent allocation from measured proteomics (opt-in):
            # the measured f_sector(T) drives the sector split at each temperature.
            fm, fmaint = ecm._alloc_from_data.model_alloc(float(Tc))
            ecm.set_allocation(fm, fmaint, kappa_scale=pert.kappa_scale,
                               sigma_sat=pert.sigma_sat)
        elif use_alloc:
            # proteome-sector allocation path (opt-in; needs sectors wired)
            f_maint = pert.f_maint
            if pert.maint_to_bio is not None and ecm._sectors is not None:
                f_maint = ecm._sectors["f_maint_nom"] * (1.0 - pert.maint_to_bio)
            ecm.set_allocation(pert.f_metab, f_maint, kappa_scale=pert.kappa_scale,
                               sigma_sat=pert.sigma_sat)
        else:
            ecm.set_budget(budget, pert.group_alloc)
            # translation-efficiency lever on the biosynthesis cap (sectors mode)
            if ecm._sectors is not None and pert.kappa_scale != 1.0:
                s = ecm._sectors
                s["bio_constraint"].ub = s["f_bio_nom"] * s["P_total"] * float(pert.kappa_scale)
        g = ecm.model.slim_optimize()
        if g is None or not np.isfinite(g) or g < min_growth:
            g = 0.0
        growth[i] = g
    return TPC(temps_C=temps_C, growth=growth)


@dataclass
class TPCDescriptors:
    Topt_C: float          # temperature of peak growth (C)
    rmax: float            # peak growth rate (1/h)
    CTmin_C: float         # low-T growth cessation (C)
    CTmax_C: float         # high-T growth cessation (C)
    niche_width_C: float   # CTmax - CTmin
    B80_C: float           # width where growth >= 0.8 * rmax
    Ea_eV: float           # activation energy of rising limb (eV)
    Ea_kJmol: float        # same, kJ/mol
    skewness: float        # curve asymmetry (negative = long cold tail)

    def as_dict(self) -> Dict[str, float]:
        return asdict(self)


class TPC:
    def __init__(self, temps_C: np.ndarray, growth: np.ndarray):
        self.temps_C = np.asarray(temps_C, float)
        self.growth = np.asarray(growth, float)

    # -- crossings ----------------------------------------------------------
    def _cross(self, level: float, side: str, i_peak: int) -> float:
        T, g = self.temps_C, self.growth
        if side == "low":
            idx = range(i_peak, 0, -1)
        else:
            idx = range(i_peak, len(T) - 1)
        for i in idx:
            j = i - 1 if side == "low" else i + 1
            if (g[i] - level) * (g[j] - level) <= 0 and g[i] != g[j]:
                # linear interpolation for the crossing temperature
                f = (level - g[i]) / (g[j] - g[i])
                return float(T[i] + f * (T[j] - T[i]))
        return float(T[0] if side == "low" else T[-1])

    def descriptors(self, crit_frac: float = 0.05) -> TPCDescriptors:
        T, g = self.temps_C, self.growth
        if g.max() <= 0:
            nan = float("nan")
            return TPCDescriptors(nan, 0.0, nan, nan, nan, nan, nan, nan, nan)
        i_peak = int(np.argmax(g))
        rmax = float(g[i_peak])
        Topt = float(T[i_peak])
        crit = crit_frac * rmax
        CTmin = self._cross(crit, "low", i_peak)
        CTmax = self._cross(crit, "high", i_peak)
        B80_lo = self._cross(0.8 * rmax, "low", i_peak)
        B80_hi = self._cross(0.8 * rmax, "high", i_peak)
        return TPCDescriptors(
            Topt_C=Topt,
            rmax=rmax,
            CTmin_C=CTmin,
            CTmax_C=CTmax,
            niche_width_C=CTmax - CTmin,
            B80_C=B80_hi - B80_lo,
            Ea_eV=self._activation_energy_eV(i_peak),
            Ea_kJmol=self._activation_energy_eV(i_peak) * 96.485,
            skewness=self._skewness(),
        )

    def _activation_energy_eV(self, i_peak: int) -> float:
        """Boltzmann-Arrhenius slope of the rising limb: ln r vs 1/(k T)."""
        T, g = self.temps_C, self.growth
        rmax = g[i_peak]
        mask = (np.arange(len(T)) <= i_peak) & (g > 0.1 * rmax) & (g < 0.95 * rmax) & (g > 0)
        if mask.sum() < 3:
            return float("nan")
        Tk = T[mask] + 273.15
        x = 1.0 / (K_BOLTZMANN_EV * Tk)
        y = np.log(g[mask])
        slope = np.polyfit(x, y, 1)[0]
        return float(-slope)   # Ea in eV

    def _skewness(self) -> float:
        """Mass-weighted skewness of the curve (asymmetry of the TPC)."""
        g = self.growth
        if g.sum() <= 0:
            return float("nan")
        T = self.temps_C
        w = g / g.sum()
        mu = np.sum(w * T)
        var = np.sum(w * (T - mu) ** 2)
        if var <= 0:
            return 0.0
        return float(np.sum(w * (T - mu) ** 3) / var ** 1.5)

    def to_dict(self):
        return {"temps_C": self.temps_C.tolist(), "growth": self.growth.tolist()}
