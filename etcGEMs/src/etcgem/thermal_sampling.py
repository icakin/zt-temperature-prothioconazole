"""Correlated / calibrated per-enzyme thermal-parameter sampling (M1.2).

Instead of perturbing the thermal envelope with two global knobs (dTopt shift,
topt_scale), draw each enzyme's (Topt_i, dCp_i) from a one-factor model that
shares an organism-wide thermal regime plus per-enzyme idiosyncrasy:

    Topt_i = mean_Topt_i + sd_Topt_i * ( sqrt(rho)*Z  + sqrt(1-rho)*eps_i )
    dCp_i  = mean_dCp_i  + sd_dCp_i  * ( sqrt(rho)*Z2 + sqrt(1-rho)*eps2_i )

with shared Z, Z2 ~ N(0,1) per ensemble member and eps_i, eps2_i ~ N(0,1) per
enzyme. rho = shared_fraction in [0,1]: rho=1 -> a coherent whole-proteome shift
(like a global dTopt); rho=0 -> independent per-enzyme optima (explores
mismatched-optima space). Captures that thermostability co-varies across a genome.

Modes:
  * "correlated" : mean = nominal per-enzyme Topt/dCp; sd from config
    (topt_sd_K, dcp_sd_frac).
  * "posterior"  : mean + sd per enzyme from DLTKcat fits.csv (Topt_C/dCp +
    Topt_sd/dCp_sd from dltkcat.fit_*); enzymes without a fit fall back to
    correlated defaults.

Applied by mutating the enzyme entries and calling refresh_params (the same
mutate-then-refresh pattern as dltkcat.apply_fits_to_provider); always restore
the nominal parameters after an evaluation.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

DCP_MAX = -1e-3   # dCp must stay negative for a single-peaked MMRT curve


def nominal_thermal(pm):
    ents = pm.ec.table.entries
    return (np.array([e.Topt for e in ents], float),
            np.array([e.dCp for e in ents], float))


def apply_thermal_sample(pm, sample: Dict[str, np.ndarray]):
    """Mutate every enzyme entry's Topt/dCp from a sample, then refresh."""
    ents = pm.ec.table.entries
    Topt, dCp = sample["Topt"], sample["dCp"]
    for i, e in enumerate(ents):
        e.Topt = float(Topt[i])
        e.dCp = float(dCp[i])
    pm.ec.refresh_params()


def restore_thermal(pm, nom_Topt, nom_dCp):
    ents = pm.ec.table.entries
    for i, e in enumerate(ents):
        e.Topt = float(nom_Topt[i])
        e.dCp = float(nom_dCp[i])
    pm.ec.refresh_params()


def _posterior_means_sds(pm, posterior_df: pd.DataFrame, key: str,
                         topt_sd_K: float, dcp_sd_frac: float):
    """Per-enzyme (mean_Topt, mean_dCp, sd_Topt, sd_dCp) from a fits table,
    falling back to nominal mean + correlated sd where a fit is missing."""
    ents = pm.ec.table.entries
    nomT, nomC = nominal_thermal(pm)
    lut = {row[key]: row for _, row in posterior_df.iterrows()} if posterior_df is not None else {}
    meanT, meanC = nomT.copy(), nomC.copy()
    sdT = np.full(len(ents), float(topt_sd_K))
    sdC = np.abs(nomC) * float(dcp_sd_frac)
    for i, e in enumerate(ents):
        k = e.rxn_id if key == "rxn_id" else e.enzyme_id
        r = lut.get(k)
        if r is None or not bool(r.get("ok", True)):
            continue
        if np.isfinite(r.get("Topt_C", np.nan)):
            meanT[i] = float(r["Topt_C"]) + 273.15
        if np.isfinite(r.get("dCp", np.nan)):
            meanC[i] = float(r["dCp"])
        if np.isfinite(r.get("Topt_sd", np.nan)):
            sdT[i] = float(r["Topt_sd"])
        if np.isfinite(r.get("dCp_sd", np.nan)):
            sdC[i] = float(r["dCp_sd"])
    return meanT, meanC, sdT, sdC


def sample_thermal(pm, n: int, mode: str = "correlated", rho: float = 0.7,
                   topt_sd_K: float = 4.0, dcp_sd_frac: float = 0.3,
                   posterior_df: Optional[pd.DataFrame] = None, key: str = "rxn_id",
                   seed: int = 1) -> List[Dict[str, np.ndarray]]:
    """Return a list of n per-enzyme thermal samples (dicts with Topt, dCp, Z, Z2)."""
    ents = pm.ec.table.entries
    E = len(ents)
    if mode == "posterior":
        meanT, meanC, sdT, sdC = _posterior_means_sds(pm, posterior_df, key,
                                                      topt_sd_K, dcp_sd_frac)
    else:  # correlated
        meanT, meanC = nominal_thermal(pm)
        sdT = np.full(E, float(topt_sd_K))
        sdC = np.abs(meanC) * float(dcp_sd_frac)
    rho = float(np.clip(rho, 0.0, 1.0))
    a, b = np.sqrt(rho), np.sqrt(1.0 - rho)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        Z = float(rng.standard_normal())
        Z2 = float(rng.standard_normal())
        eps = rng.standard_normal(E)
        eps2 = rng.standard_normal(E)
        Topt = meanT + sdT * (a * Z + b * eps)
        dCp = np.minimum(meanC + sdC * (a * Z2 + b * eps2), DCP_MAX)
        out.append({"Topt": Topt, "dCp": dCp, "Z": Z, "Z2": Z2})
    return out
