"""Sensitivity analysis: sweep proteome allocation and kcat(T) responses.

Draws a Latin-hypercube sample over the perturbation parameters, computes a TPC
for each point, and reduces the ensemble to descriptor distributions plus
Spearman sensitivity indices (a PRCC-style measure of how strongly each input
moves each TPC feature).

Perturbation parameters (any subset can be swept):
    dTopt        shift (K) applied to every enzyme optimum
    topt_scale   scales the spread of enzyme optima about T0 (heterogeneity)
    dCp_scale    scales heat capacity of activation (curvature / breadth)
    budget_scale multiplies the total proteome pool (global allocation)
    alloc_<grp>  multiplies a group's sub-budget (allocation between pathways)
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

from .enzyme_cost import Perturbation
from .tpc import TPC, compute_tpc


def _lhs(n: int, d: int, seed: int) -> np.ndarray:
    """Latin-hypercube sample in the unit cube (scipy if present, else numpy)."""
    try:
        from scipy.stats.qmc import LatinHypercube
        return LatinHypercube(d=d, seed=seed).random(n)
    except Exception:
        rng = np.random.default_rng(seed)
        cut = (np.arange(n)[:, None] + rng.random((n, d))) / n
        return np.take_along_axis(cut, rng.random((n, d)).argsort(0), 0)


def _make_perturbation(sample: Dict[str, float], default_budget: float,
                       group_names: Sequence[str]) -> Perturbation:
    p = Perturbation(
        dTopt=sample.get("dTopt", 0.0),
        topt_scale=sample.get("topt_scale", 1.0),
        dCp_scale=sample.get("dCp_scale", 1.0),
        dTm=sample.get("dTm", 0.0),
    )
    if "budget_scale" in sample:
        p.budget = default_budget * sample["budget_scale"]
    for g in group_names:
        key = f"alloc_{g}"
        if key in sample:
            p.group_alloc[g] = sample[key]
    # proteome-sector allocation (opt-in; only used when sectors are wired)
    if "f_metab" in sample:
        p.f_metab = sample["f_metab"]
    if "f_maint" in sample:
        p.f_maint = sample["f_maint"]
    if "maint_to_bio" in sample:
        p.maint_to_bio = sample["maint_to_bio"]
    return p


@dataclass
class SensitivityResult:
    temps_C: np.ndarray
    curves: np.ndarray                 # (n_samples, n_temps)
    samples: pd.DataFrame              # sampled inputs
    descriptors: pd.DataFrame          # one row per sample
    sensitivity: pd.DataFrame          # Spearman rho, inputs x descriptors
    nominal: TPC                       # unperturbed reference curve

    def save(self, out_dir: str):
        os.makedirs(out_dir, exist_ok=True)
        np.save(os.path.join(out_dir, "temps_C.npy"), self.temps_C)
        np.save(os.path.join(out_dir, "curves.npy"), self.curves)
        self.samples.to_csv(os.path.join(out_dir, "samples.csv"), index=False)
        self.descriptors.to_csv(os.path.join(out_dir, "descriptors.csv"), index=False)
        self.sensitivity.to_csv(os.path.join(out_dir, "sensitivity_spearman.csv"))
        pd.DataFrame({"temp_C": self.nominal.temps_C,
                      "growth": self.nominal.growth}).to_csv(
            os.path.join(out_dir, "nominal_tpc.csv"), index=False)


def run_sensitivity(pm, temps_C: Sequence[float],
                    param_ranges: Dict[str, Tuple[float, float]],
                    n_samples: int = 200, seed: int = 1,
                    group_names: Sequence[str] = (),
                    crit_frac: float = 0.05,
                    envelope_samples=None,
                    progress: bool = True) -> SensitivityResult:
    """LHS sweep over param_ranges.

    ``envelope_samples`` (optional): a list of n_samples per-enzyme thermal draws
    (thermal_sampling.sample_thermal). When given, each sample's per-enzyme
    Topt/dCp is applied (mutating entries) before its TPC instead of the
    dTopt/topt_scale/dCp_scale knobs; the shared latent Z is recorded as a
    sensitivity input. When None, behaviour is exactly as before.
    """
    temps_C = np.asarray(temps_C, float)
    default_budget = pm.ec.default_budget
    names = list(param_ranges.keys())
    U = _lhs(n_samples, len(names), seed)

    # scale unit cube to parameter ranges
    sampled = {}
    for j, name in enumerate(names):
        lo, hi = param_ranges[name]
        sampled[name] = lo + U[:, j] * (hi - lo)
    samples_df = pd.DataFrame(sampled)

    use_env = envelope_samples is not None
    if use_env:
        from .thermal_sampling import (nominal_thermal, apply_thermal_sample,
                                        restore_thermal)
        nomT, nomC = nominal_thermal(pm)
        samples_df = samples_df.copy()
        samples_df["thermal_Z"] = [envelope_samples[i]["Z"] for i in range(n_samples)]

    curves = np.zeros((n_samples, len(temps_C)))
    desc_rows: List[dict] = []
    for i in range(n_samples):
        row = {k: float(samples_df.iloc[i][k]) for k in names}
        pert = _make_perturbation(row, default_budget, group_names)
        if use_env:
            apply_thermal_sample(pm, envelope_samples[i])
        tpc = compute_tpc(pm, temps_C, pert)
        curves[i] = tpc.growth
        desc_rows.append(tpc.descriptors(crit_frac).as_dict())
        if progress and (i + 1) % max(1, n_samples // 10) == 0:
            print(f"[sensitivity] {i + 1}/{n_samples}")
    if use_env:
        restore_thermal(pm, nomT, nomC)
    desc_df = pd.DataFrame(desc_rows)

    nominal = compute_tpc(pm, temps_C, Perturbation())
    sens = _spearman_matrix(samples_df, desc_df)
    return SensitivityResult(temps_C, curves, samples_df, desc_df, sens, nominal)


def _spearman_matrix(inputs: pd.DataFrame, outputs: pd.DataFrame) -> pd.DataFrame:
    from scipy.stats import spearmanr
    out = pd.DataFrame(index=inputs.columns, columns=outputs.columns, dtype=float)
    for a in inputs.columns:
        for b in outputs.columns:
            y = outputs[b].values
            m = np.isfinite(y)
            if m.sum() < 3 or np.nanstd(y[m]) == 0:
                out.loc[a, b] = np.nan
            else:
                out.loc[a, b] = spearmanr(inputs[a].values[m], y[m]).correlation
    return out


# ---------------------------------------------------------------------------
# equal-perturbation (standardised elasticity) sensitivity
# ---------------------------------------------------------------------------
# The LHS/Spearman sweep above scores influence with a RANK correlation over
# HAND-SET, UNEQUAL per-input ranges: a dial swept over a wider range can rank
# higher for that reason alone, and rank correlation measures monotonic
# *consistency*, not *magnitude*. The elasticity analysis instead moves EVERY
# input by the SAME standardised step h and reports a dimensionless, directly
# comparable magnitude (a local elasticity) per descriptor, so the ranking
# reflects the model's structural leverage rather than our range choices.


@dataclass
class ElasticityResult:
    inputs: List[str]
    descriptors: List[str]
    elasticity: pd.DataFrame   # inputs x descriptors: E[D,p] (dimensionless, signed)
    raw_delta: pd.DataFrame    # inputs x descriptors: raw D(p+) - D(p-)
    nominal: Dict[str, float]  # descriptor -> nominal value
    reference_scales: Dict     # h, dTopt reference scale, sector nominal, input list

    def save(self, out_dir: str):
        os.makedirs(out_dir, exist_ok=True)
        self.elasticity.to_csv(os.path.join(out_dir, "elasticity_table.csv"))
        self.raw_delta.to_csv(os.path.join(out_dir, "elasticity_raw_delta.csv"))
        np.sign(self.raw_delta).to_csv(os.path.join(out_dir, "elasticity_signs.csv"))
        import json
        with open(os.path.join(out_dir, "reference_scales.json"), "w") as fh:
            json.dump(self.reference_scales, fh, indent=2)
        with open(os.path.join(out_dir, "elasticity_nominal.json"), "w") as fh:
            json.dump(self.nominal, fh, indent=2)


def _elasticity_pert(inp: str, signed_h: float, default_budget: float,
                     dtopt_ref: float, f_metab_nom, f_maint_nom, group_names,
                     dtm_ref: float = 0.0):
    """Build the Perturbation for one input at a standardised signed step.

    Multiplicative inputs (nominal 1) move to 1+signed_h. The additive dTopt and
    dTm shifts (nominal 0) move by signed_h * (their reference scale) -- the sd of
    the per-enzyme Topt / Tm distribution respectively (so it is 'h of a natural
    temperature scale', comparable to an h fractional move of a multiplicative
    dial). Sector fractions move to nominal*(1+signed_h); set_allocation
    renormalises the complementary biosynthesis sector (f_bio = 1 - f_metab -
    f_maint) so the budget still sums."""
    if inp == "dTopt":
        return Perturbation(dTopt=signed_h * dtopt_ref)
    if inp == "dTm":
        return Perturbation(dTm=signed_h * dtm_ref)
    if inp == "topt_scale":
        return Perturbation(topt_scale=1.0 + signed_h)
    if inp == "dCp_scale":
        return Perturbation(dCp_scale=1.0 + signed_h)
    if inp == "budget_scale":
        return Perturbation(budget=default_budget * (1.0 + signed_h))
    if inp == "f_metab":
        return Perturbation(f_metab=f_metab_nom * (1.0 + signed_h), f_maint=f_maint_nom)
    if inp == "f_maint":
        return Perturbation(f_metab=f_metab_nom, f_maint=f_maint_nom * (1.0 + signed_h))
    if inp.startswith("alloc_"):
        return Perturbation(group_alloc={inp[len("alloc_"):]: 1.0 + signed_h})
    raise ValueError(f"unknown elasticity input: {inp}")


def run_elasticity(pm, temps_C: Sequence[float], inputs: Sequence[str],
                   h: float = 0.10, crit_frac: float = 0.05,
                   group_names: Sequence[str] = (),
                   descriptors=("Topt_C", "rmax", "CTmax_C", "CTmin_C",
                                "niche_width_C", "Ea_eV")) -> ElasticityResult:
    """Central-finite-difference elasticities at the nominal point: every input is
    perturbed by the SAME standardised step +/- h and each TPC descriptor D gets

        E[D,p] = ( D(p+) - D(p-) ) / ( 2 * h * D_nominal )

    a dimensionless, comparable sensitivity magnitude. Local (around the nominal)."""
    temps_C = np.asarray(temps_C, float)
    ec = pm.ec
    default_budget = ec.default_budget
    dtopt_ref = float(np.std(np.asarray(ec._Topt, float))) if getattr(ec, "_Topt", None) is not None \
        and len(ec._Topt) else 5.0
    dtm_ref = float(np.std(np.asarray(ec._Tm, float))) if getattr(ec, "_Tm", None) is not None \
        and len(ec._Tm) else 5.0
    sec = getattr(ec, "_sectors", None)
    f_metab_nom = sec["f_metab_nom"] if sec else None
    f_maint_nom = sec["f_maint_nom"] if sec else None

    nom = compute_tpc(pm, temps_C, Perturbation()).descriptors(crit_frac).as_dict()
    descs = [d for d in descriptors if d in nom]
    E = pd.DataFrame(index=list(inputs), columns=descs, dtype=float)
    R = pd.DataFrame(index=list(inputs), columns=descs, dtype=float)
    for inp in inputs:
        dp = compute_tpc(pm, temps_C, _elasticity_pert(inp, +h, default_budget, dtopt_ref,
                         f_metab_nom, f_maint_nom, group_names, dtm_ref)).descriptors(crit_frac).as_dict()
        dm = compute_tpc(pm, temps_C, _elasticity_pert(inp, -h, default_budget, dtopt_ref,
                         f_metab_nom, f_maint_nom, group_names, dtm_ref)).descriptors(crit_frac).as_dict()
        for D in descs:
            delta = dp.get(D, np.nan) - dm.get(D, np.nan)
            R.loc[inp, D] = delta
            Dn = nom.get(D, np.nan)
            E.loc[inp, D] = delta / (2.0 * h * Dn) if (Dn is not None and np.isfinite(Dn)
                                                       and abs(Dn) > 1e-9) else np.nan
    refs = {"h": h, "dTopt_reference_scale_K": round(dtopt_ref, 3),
            "dTm_reference_scale_K": round(dtm_ref, 3),
            "sector_f_metab_nom": f_metab_nom, "sector_f_maint_nom": f_maint_nom,
            "inputs": list(inputs),
            "note": "additive dTopt/dTm perturbed by +/- h*(their reference scale = sd of "
                    "per-enzyme Topt/Tm); multiplicative inputs by 1 +/- h; sector fractions "
                    "by nominal*(1 +/- h) with f_bio renormalised."}
    return ElasticityResult(list(inputs), descs, E, R,
                            {k: round(float(nom[k]), 4) for k in descs}, refs)
