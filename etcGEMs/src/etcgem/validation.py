"""Validate the emergent etc-GEM against two trusted, strain-matched E. coli TPCs.

Nothing is fit to growth: the model predicts each curve a priori under its medium
(availability, not pinned uptake) and we compare on RAW ABSOLUTE growth rate (1/h).

Curves (primary sources, replacing the retired Smith-derived 26-curve compilation):
  * MINIMAL -- Noll / Katipoglu-Yazan et al. 2023 (Data in Brief 48:109037),
    K-12 NCM3722, defined glucose-minimal (modified MOPS + 0.5 g/L glucose + 6 trace
    amino acids), 27-45 C, mean +/- SD over n wells. The quantitative anchor.
  * RICH   -- Erdos et al. 2026 (unpublished), K-12 MG1655 wt, LB, ~16-45 C.
    Figure-digitized (+/- ~0.1 h/1); lower precision, the rich-medium comparator.

The headline is the minimal-vs-rich magnitude contrast: does the medium-matched
proteome-sector allocation reproduce the observed rich > minimal peak?
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .enzyme_cost import Perturbation
from .tpc import TPC, compute_tpc

SRC = os.path.join("strains", "{strain}", "thermal", "sources")
LB_MEDIA_CSV = os.path.join("strains", "{strain}", "media", "LB_media.csv")
BHI_MEDIA_CSV = os.path.join("strains", "{strain}", "media", "BHI_media.csv")

# 6 trace amino acids supplemented in the Noll minimal medium (exchange REV ids)
NOLL_TRACE_AA = ["met__L", "his__L", "arg__L", "pro__L", "thr__L", "trp__L"]


@dataclass
class CurveSpec:
    key: str
    rel_path: str
    medium: str                 # "glucose_minimal" | "LB" | "BHI"
    label: str
    color: str
    carbon: str = "glc__D"
    digitized: bool = False


# Primary (and only main) validation curve: the exact GEM strain (K-12 MG1655) on
# a rich medium, spanning the cold rising limb through the hot collapse.
CURVES = [
    CurveSpec("vanderlinden_bhi",
              "vanderlinden2012_intfoodmicro/vanderlinden2012_mg1655_bhi_tpc.csv",
              "BHI",
              "Van Derlinden 2012 — K-12 MG1655, BHI (rich)", "tab:blue", digitized=True),
]

# The validation is Van Derlinden only (exact strain). Noll (NCM3722, wrong strain)
# and Erdos (unpublished LB cross-check) are dropped from the pipeline.
SECONDARY = []


# ---------------------------------------------------------------------------
def load_curve(strain: str, spec: CurveSpec) -> Dict:
    df = pd.read_csv(os.path.join(SRC.format(strain=strain), spec.rel_path))
    df = df[np.isfinite(pd.to_numeric(df["rate_per_h"], errors="coerce"))].copy()
    df["rate_per_h"] = df["rate_per_h"].astype(float)
    temps = df["temp_C"].to_numpy(float)
    rates = df["rate_per_h"].to_numpy(float)
    sd = df["rate_sd_per_h"].to_numpy(float) if "rate_sd_per_h" in df.columns else None
    n = df["n_wells"].to_numpy(float) if "n_wells" in df.columns else None
    if sd is not None and not np.isfinite(sd).any():
        sd = None
    return {"spec": spec, "temps_C": temps, "rate": rates, "sd": sd, "n": n,
            "study": str(df["study"].iloc[0]), "strain": str(df["strain"].iloc[0]),
            "medium_detail": str(df["medium_detail"].iloc[0]),
            "obs_rmax": float(np.max(rates)),
            "obs_Topt_C": float(temps[int(np.argmax(rates))])}


def _set_medium(pm, spec: CurveSpec, strain: str, extra_aa: Optional[List[str]] = None):
    from .providers import set_medium
    if spec.medium == "BHI":
        set_medium(pm, "BHI", bhi_media_csv=BHI_MEDIA_CSV.format(strain=strain))
    elif spec.medium == "LB":
        set_medium(pm, "LB", lb_media_csv=LB_MEDIA_CSV.format(strain=strain))
    else:
        set_medium(pm, "glucose_minimal", spec.carbon, True)
        if extra_aa:
            model = pm.ec.model
            for base in extra_aa:
                rev = f"EX_{base}_e_REV"
                if rev in model.reactions:
                    model.reactions.get_by_id(rev).upper_bound = 1000.0
            model.solver.update()


def _fit_stats(obs, pred):
    obs, pred = np.asarray(obs, float), np.asarray(pred, float)
    ss_res = float(np.sum((obs - pred) ** 2))
    ss_tot = float(np.sum((obs - np.mean(obs)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rmse = float(np.sqrt(np.mean((obs - pred) ** 2)))
    return r2, rmse


def _predict(pm, spec, strain, data_temps, extra_aa=None, pert=None, set_med=True):
    """Return (pred_at_data_temps, dense_T, dense_growth, descriptors)."""
    if set_med:
        _set_medium(pm, spec, strain, extra_aa)
    pert = pert or Perturbation()
    pred = compute_tpc(pm, data_temps, pert).growth
    dense = np.linspace(min(5.0, float(min(data_temps)) - 3),
                        max(52.0, float(max(data_temps)) + 6), 101)
    dg = compute_tpc(pm, dense, pert).growth
    desc = TPC(dense, dg).descriptors(0.05)
    return pred, dense, dg, desc


# ---------------------------------------------------------------------------
def _build_emergent_rich(strain):
    """Provider configured EXACTLY as the P2 tuning will use it: rich BHI,
    reconciled single proteome pool (P1c), coupled growth law ON (P1b), T-dependent
    maintenance restored. The emergent prediction sets the measured rich sector
    fractions and identity magnitude/envelope knobs."""
    from .calibration_multi import _build_pm_rich
    pm = _build_pm_rich(strain)   # BHI, growth law ON, static rich alloc, reconciled
    pert = Perturbation(f_metab=0.280, f_maint=0.360)   # measured rich (LB) sectors; else identity
    return pm, pert


def _score_curve(pm, spec, strain, growth_law_rich=False):
    cv = load_curve(strain, spec)
    if growth_law_rich and spec.medium == "BHI":
        pmr, pert = _build_emergent_rich(strain)
        pred, dense, dg, desc = _predict(pmr, spec, strain, cv["temps_C"], pert=pert, set_med=False)
    else:
        pred, dense, dg, desc = _predict(pm, spec, strain, cv["temps_C"])
    r2, rmse = _fit_stats(cv["rate"], pred)
    obs_desc = TPC(cv["temps_C"], cv["rate"]).descriptors(0.05)
    res = {
        "study": cv["study"], "strain": cv["strain"], "medium": spec.medium,
        "digitized": spec.digitized, "n": int(len(cv["temps_C"])),
        "T_range_C": [float(cv["temps_C"].min()), float(cv["temps_C"].max())],
        "abs_R2": round(r2, 3), "RMSE_per_h": round(rmse, 3),
        "obs_rmax": round(cv["obs_rmax"], 3), "pred_rmax": round(float(desc.rmax), 3),
        "obs_Topt_C": round(cv["obs_Topt_C"], 1), "pred_Topt_C": round(float(desc.Topt_C), 1),
        "obs_CTmax_C": round(float(obs_desc.CTmax_C), 1),
        "pred_CTmax_C": round(float(desc.CTmax_C), 1),
        "obs_Ea_eV": round(float(obs_desc.Ea_eV), 3),
        "pred_Ea_eV": round(float(desc.Ea_eV), 3),
    }
    cv["pred"] = pred; cv["dense_T"] = dense; cv["dense_g"] = dg; cv["desc"] = desc
    res["_cv"] = cv
    return res


def run(strain: str, out_dir: str, include_secondary: bool = True,
        growth_law_rich: bool = True, **_ignore) -> Dict:
    """Emergent (nothing-fit) validation. ``growth_law_rich`` (default True) predicts
    the rich BHI curve on the model configured as it will be tuned: reconciled single
    proteome pool + coupled growth law ON + T-dependent maintenance (rounded peak)."""
    from .config import resolve, build_provider
    pm = build_provider(resolve(strain))
    try:
        pm.ec.model.solver.configuration.timeout = 10
    except Exception:
        pass
    os.makedirs(out_dir, exist_ok=True)

    results, rows = {}, []
    specs = list(CURVES) + (list(SECONDARY) if include_secondary else [])
    for spec in specs:
        res = _score_curve(pm, spec, strain, growth_law_rich=growth_law_rich)
        res["role"] = "primary" if spec in CURVES else "secondary_crosscheck"
        results[spec.key] = res
        rows.append({"curve": spec.key, "role": res["role"],
                     **{k: v for k, v in res.items() if k not in ("_cv", "role")}})

    _plot_curves(results, out_dir)

    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "validation_trusted_table.csv"), index=False)
    clean = {k: {kk: vv for kk, vv in v.items() if kk != "_cv"} if isinstance(v, dict) else v
             for k, v in results.items()}
    with open(os.path.join(out_dir, "validation_trusted_summary.json"), "w") as fh:
        json.dump(clean, fh, indent=2)
    return clean


# ---------------------------------------------------------------------------
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _plot_curves(results, out_dir):
    """Emergent model vs the primary Van Derlinden curve (and any secondary
    cross-check), raw absolute rate. Primary first, one panel per curve."""
    plt = _mpl()
    keys = [s.key for s in CURVES if s.key in results] + \
           [s.key for s in SECONDARY if s.key in results]
    n = len(keys)
    fig, axes = plt.subplots(1, n, figsize=(6.0 * n, 5))
    axes = np.atleast_1d(axes)
    for ax, key in zip(axes, keys):
        r = results[key]; cv = r["_cv"]; spec = cv["spec"]
        if cv["sd"] is not None:
            ax.errorbar(cv["temps_C"], cv["rate"], yerr=cv["sd"], fmt="o", color="k",
                        ms=5, capsize=3, label="data (mean ± SD)", zorder=3)
        else:
            ax.plot(cv["temps_C"], cv["rate"], "o", color="k", ms=5,
                    label="data (digitized)", zorder=3)
        ax.plot(cv["dense_T"], cv["dense_g"], color=spec.color, lw=2.2,
                label="emergent model")
        tag = "primary" if spec in CURVES else "secondary cross-check"
        ax.set_title(f"{spec.label}\n({tag})", fontsize=10)
        ax.set_xlabel("Temperature (°C)"); ax.set_ylabel("Growth rate (1/h)")
        ax.text(0.03, 0.97, f"absolute $R^2$={r['abs_R2']}\nRMSE={r['RMSE_per_h']} 1/h\n"
                f"obs $r_{{max}}$={r['obs_rmax']} vs pred {r['pred_rmax']}\n"
                f"obs $T_{{opt}}$={r['obs_Topt_C']} vs pred {r['pred_Topt_C']} °C",
                transform=ax.transAxes, va="top", fontsize=8.5,
                bbox=dict(boxstyle="round", fc="0.96", ec="0.8"))
        ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("Emergent model vs Van Derlinden (K-12 MG1655, BHI; raw absolute rate, nothing fit)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p = os.path.join(out_dir, "validation_trusted_curves.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    return p
