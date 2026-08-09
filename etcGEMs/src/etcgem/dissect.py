"""Dissect the model around an arbitrary baseline operating point.

Given a provider at a chosen operating point (here the rich BHI, growth-law-ON,
reconciled single pool) and a BASELINE ``Perturbation`` (here the P2 v3 tuned
posterior medians), this module computes, around that baseline:

* equal-perturbation **elasticities** E[D,p] for every descriptor D and lever p
  (central finite difference at a standardised step h), with an optional
  **posterior-propagated** IQR band from a sample of posterior draws;
* a two-group **variance decomposition** (envelope vs magnitude) via a crossed
  functional-ANOVA / Shapley grid, plus the finer per-sub-group achievable-IQR
  magnitude and median+IQR reachable-TPC bands (kinetic / stability / allocation
  / catalysis / maintenance).

Nothing here is specific to the tuned model: pass any baseline Perturbation and
operating-point provider. The heavy TPC evaluations are farmed out to a Gurobi
process pool (each worker builds its own provider -- cobra models are not
picklable). Reuses the exact v3 free-parameter mapping (calibration_multi specs)
and the two-group Shapley primitive (decomposition.decompose_grid)."""
from __future__ import annotations

import json
import os
from dataclasses import replace
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .enzyme_cost import Perturbation
from .tpc import TPC, compute_tpc

# --- the unified P2 v3 lever set and its groupings ---------------------------
LEVERS = ["dTopt", "topt_scale", "dCp_scale", "dTm", "tm_scale",
          "kcat_scale", "kappa_scale", "sigma", "f_metab", "f_maint",
          "ngam_scale", "ngam_steepness"]
GROUPS = {                                   # finer 5-way split
    "kinetic":     ["dTopt", "topt_scale", "dCp_scale"],
    "stability":   ["dTm", "tm_scale"],
    "allocation":  ["f_metab", "f_maint"],
    "catalysis":   ["kcat_scale", "kappa_scale", "sigma"],
    "maintenance": ["ngam_scale", "ngam_steepness"],
}
ENVELOPE = GROUPS["kinetic"] + GROUPS["stability"]                    # 5 levers
MAGNITUDE = GROUPS["allocation"] + GROUPS["catalysis"] + GROUPS["maintenance"]  # 7
ADDITIVE = {"dTopt", "dTm"}                  # stepped by h * reference scale (K)
LEVER_ATTR = {lv: ("sigma_sat" if lv == "sigma" else lv) for lv in LEVERS}
DESCRIPTORS = ["Topt_C", "rmax", "CTmax_C", "niche_width_C", "Ea_eV"]


# --- baseline Perturbation + posterior draws from a v3 calibration dir --------
def tuned_pert_from_v3(v3_dir: str):
    """(Perturbation, medians dict) from a calibration_vanderlinden_v3 summary."""
    d = json.load(open(os.path.join(v3_dir, "summary.json")))
    med = {k: v["posterior_median"] for k, v in d["posterior"].items()}
    kw = {LEVER_ATTR[lv]: float(med[lv]) for lv in LEVERS if lv in med}
    return Perturbation(**kw), med


def load_draws(v3_dir: str, n: int = 150, seed: int = 0) -> List[Perturbation]:
    """Sample n posterior draws from the v3 flat chain -> Perturbations (identity
    on the noise term; sigma_disc is not a model lever)."""
    from .calibration_multi import build_vdl_specs, to_pert
    flat = np.load(os.path.join(v3_dir, "chain_flat.npy"))
    specs = build_vdl_specs()
    rng = np.random.default_rng(seed)
    idx = rng.choice(flat.shape[0], size=min(n, flat.shape[0]), replace=False)
    return [to_pert(flat[i], specs) for i in idx]


# --- per-lever standardised step / range about a baseline --------------------
def _refs(pm) -> Dict[str, float]:
    ec = pm.ec
    sd_topt = float(np.std(np.asarray(ec._Topt, float))) if getattr(ec, "_Topt", None) is not None \
        and len(ec._Topt) else 5.0
    sd_tm = float(np.std(np.asarray(ec._Tm, float))) if getattr(ec, "_Tm", None) is not None \
        and len(ec._Tm) else 6.0
    return {"dTopt": sd_topt, "dTm": sd_tm}


def _base_val(base: Perturbation, lever: str) -> float:
    v = getattr(base, LEVER_ATTR[lever])
    if v is not None:
        return float(v)
    return 0.0 if lever in ADDITIVE else 1.0


def _step(base: Perturbation, lever: str, signed_h: float, refs: Dict[str, float]) -> Perturbation:
    """Baseline with ONE lever moved by a standardised signed step: additive
    levers by signed_h * reference scale (K); the rest multiplicatively by
    (1+signed_h). sigma is clipped to its physical (0,1] support."""
    v0 = _base_val(base, lever)
    if lever in ADDITIVE:
        newv = v0 + signed_h * refs[lever]
    else:
        newv = v0 * (1.0 + signed_h)
        if lever == "sigma":
            newv = float(np.clip(newv, 1e-3, 1.0))
    return replace(base, **{LEVER_ATTR[lever]: newv})


def _group_ranges(base: Perturbation, levers: Sequence[str], H: float,
                  refs: Dict[str, float]) -> Dict[str, tuple]:
    out = {}
    for lv in levers:
        v0 = _base_val(base, lv)
        if lv in ADDITIVE:
            out[lv] = (v0 - H * refs[lv], v0 + H * refs[lv])
        elif lv == "sigma":
            out[lv] = (max(1e-3, v0 * (1 - H)), min(1.0, v0 * (1 + H)))
        else:
            out[lv] = (v0 * (1 - H), v0 * (1 + H))
    return out


def _pert_from_sample(base: Perturbation, sample: Dict[str, float]) -> Perturbation:
    return replace(base, **{LEVER_ATTR[lv]: float(v) for lv, v in sample.items()})


# --- Gurobi process pool: each worker builds its own rich provider -----------
_PMR = None
_GRID = None


def _winit(strain: str, grid, timeout: float = 30.0):
    global _PMR, _GRID
    from .calibration_multi import _set_default_solver, _build_pm_rich
    _set_default_solver("gurobi")
    _PMR = _build_pm_rich(strain)
    try:
        _PMR.ec.model.solver.configuration.timeout = timeout
    except Exception:
        pass
    _GRID = np.asarray(grid, float)


def _wdesc(pert: Perturbation) -> Dict[str, float]:
    g = compute_tpc(_PMR, _GRID, pert).growth
    return TPC(_GRID, g).descriptors(0.05).as_dict()


def _wcurve(pert: Perturbation) -> np.ndarray:
    return compute_tpc(_PMR, _GRID, pert).growth


# --- elasticity (base + posterior-propagated bands) --------------------------
def _elasticities(descs_plus, descs_minus, base_desc, h) -> pd.DataFrame:
    E = pd.DataFrame(index=LEVERS, columns=DESCRIPTORS, dtype=float)
    for k, lv in enumerate(LEVERS):
        for D in DESCRIPTORS:
            dp, dm = descs_plus[k].get(D, np.nan), descs_minus[k].get(D, np.nan)
            Dn = base_desc.get(D, np.nan)
            E.loc[lv, D] = (dp - dm) / (2.0 * h * Dn) if (np.isfinite(Dn) and abs(Dn) > 1e-9) else np.nan
    return E


def run_elasticity(pool, base: Perturbation, refs, h: float,
                   draws: Optional[List[Perturbation]] = None):
    """Central-difference elasticities at ``base`` plus, if ``draws`` given, the
    posterior IQR of each elasticity (recomputed at every draw)."""
    perts = [base] + [_step(base, lv, +h, refs) for lv in LEVERS] \
                   + [_step(base, lv, -h, refs) for lv in LEVERS]
    res = pool.map(_wdesc, perts)
    base_desc = res[0]
    E = _elasticities(res[1:1 + len(LEVERS)], res[1 + len(LEVERS):], base_desc, h)

    bands = None
    if draws:
        # one big pooled batch: for every draw, base + 2*L stepped perts
        batch, spans = [], []
        for b in draws:
            start = len(batch)
            batch.append(b)
            batch += [_step(b, lv, +h, refs) for lv in LEVERS]
            batch += [_step(b, lv, -h, refs) for lv in LEVERS]
            spans.append(start)
        alld = pool.map(_wdesc, batch)
        stack = {D: {lv: [] for lv in LEVERS} for D in DESCRIPTORS}
        L = len(LEVERS)
        for s in spans:
            bd = alld[s]
            Ei = _elasticities(alld[s + 1:s + 1 + L], alld[s + 1 + L:s + 1 + 2 * L], bd, h)
            for D in DESCRIPTORS:
                for lv in LEVERS:
                    stack[D][lv].append(Ei.loc[lv, D])
        bands = {}
        for D in DESCRIPTORS:
            rows = {}
            for lv in LEVERS:
                arr = np.asarray(stack[D][lv], float)
                arr = arr[np.isfinite(arr)]
                if arr.size:
                    rows[lv] = [float(np.percentile(arr, 5)), float(np.percentile(arr, 50)),
                                float(np.percentile(arr, 95))]
                else:
                    rows[lv] = [np.nan, np.nan, np.nan]
            bands[D] = rows
    return E, base_desc, bands


# --- two-group (envelope vs magnitude) variance + per-group IQR bands --------
def run_decomposition(pool, base: Perturbation, refs, H: float, grid,
                      n_env=24, n_mag=24, n_band=41, seed=1):
    from .sensitivity import _lhs
    from .decomposition import decompose_grid
    grid = np.asarray(grid, float)

    def _lhs_samples(levers, n, sd):
        rg = _group_ranges(base, levers, H, refs)
        U = _lhs(n, len(levers), sd)
        return [{lv: rg[lv][0] + U[i, k] * (rg[lv][1] - rg[lv][0])
                 for k, lv in enumerate(levers)} for i in range(n)]

    env_s = _lhs_samples(ENVELOPE, n_env, seed)
    mag_s = _lhs_samples(MAGNITUDE, n_mag, seed + 1)
    # crossed grid: env x mag
    pairs = [_pert_from_sample(base, {**e, **m}) for e in env_s for m in mag_s]
    flat = pool.map(_wdesc, pairs)
    grids = {D: np.full((n_env, n_mag), np.nan) for D in DESCRIPTORS}
    for idx, (i, j) in enumerate([(i, j) for i in range(n_env) for j in range(n_mag)]):
        for D in DESCRIPTORS:
            grids[D][i, j] = flat[idx].get(D, np.nan)

    var_rows = {}
    for D in DESCRIPTORS:
        g = decompose_grid(grids[D])          # axis0=env (A), axis1=mag (E)
        var_rows[D] = {"phi_envelope": g["phi_A"], "phi_magnitude": g["phi_E"],
                       "S_envelope": g["S_A"], "S_magnitude": g["S_E"],
                       "S_interaction": g["S_AE"]}
    var_df = pd.DataFrame(var_rows).T

    # per-sub-group marginal ensembles -> achievable IQR + reachable-TPC bands
    base_curve = np.asarray(pool.map(_wcurve, [base])[0], float)
    group_curves, mag_rows = {}, {}
    for gi, (gname, levers) in enumerate(GROUPS.items()):
        samp = _lhs_samples(levers, n_band, seed + 10 + gi)
        perts = [_pert_from_sample(base, s) for s in samp]
        curves = np.vstack(pool.map(_wcurve, perts))
        group_curves[gname] = curves
        dd = pd.DataFrame([TPC(grid, c).descriptors(0.05).as_dict() for c in curves])
        for D in DESCRIPTORS:
            mag_rows.setdefault(D, {})[gname] = float(
                np.nanpercentile(dd[D], 75) - np.nanpercentile(dd[D], 25)) if D in dd else np.nan
    mag_df = pd.DataFrame(mag_rows).T   # descriptor x group IQR
    return var_df, mag_df, group_curves, base_curve


# --- plots -------------------------------------------------------------------
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


_GROUP_OF = {lv: g for g, lvs in GROUPS.items() for lv in lvs}
_GROUP_COLOR = {"kinetic": "tab:blue", "stability": "tab:cyan", "allocation": "tab:orange",
                "catalysis": "tab:red", "maintenance": "tab:green"}


def plot_elasticity_heatmap(E: pd.DataFrame, out_dir, fname="elasticity_heatmap.png"):
    plt = _mpl()
    M = E.loc[LEVERS, DESCRIPTORS].astype(float).values
    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    vmax = np.nanmax(np.abs(M)) or 1.0
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(DESCRIPTORS))); ax.set_xticklabels(DESCRIPTORS, rotation=30, ha="right")
    ax.set_yticks(range(len(LEVERS))); ax.set_yticklabels(LEVERS)
    for i in range(len(LEVERS)):
        for j in range(len(DESCRIPTORS)):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7,
                        color="k" if abs(M[i, j]) < 0.6 * vmax else "w")
    fig.colorbar(im, ax=ax, label="standardised elasticity E[D,p]")
    ax.set_title("Tuned-model elasticity (equal step h; rich BHI)")
    fig.tight_layout(); p = os.path.join(out_dir, fname); fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_elasticity_tornado(E, bands, out_dir, D="rmax", fname=None):
    plt = _mpl()
    fname = fname or f"elasticity_tornado_{D}.png"
    vals = E[D].astype(float)
    order = vals.abs().sort_values(ascending=True).index
    y = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    xerr = None
    if bands and D in bands:
        lo = np.array([vals[lv] - bands[D][lv][0] for lv in order])
        hi = np.array([bands[D][lv][2] - vals[lv] for lv in order])
        xerr = np.vstack([np.clip(lo, 0, None), np.clip(hi, 0, None)])
    ax.barh(y, [vals[lv] for lv in order], color=[_GROUP_COLOR[_GROUP_OF[lv]] for lv in order],
            xerr=xerr, error_kw=dict(ecolor="0.3", lw=1, capsize=2), alpha=0.9)
    ax.set_yticks(y); ax.set_yticklabels(order)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel(f"elasticity E[{D},p]  (bars = posterior median; whiskers = 90% CI over draws)")
    ax.set_title(f"What moves {D} on the tuned model (rich BHI)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=g) for g, c in _GROUP_COLOR.items()],
              frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout(); p = os.path.join(out_dir, fname); fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_variance_partition(var_df, out_dir, fname="decomp_variance.png"):
    plt = _mpl()
    ds = [d for d in DESCRIPTORS if d in var_df.index]
    env = var_df.loc[ds, "phi_envelope"].astype(float).values
    mag = var_df.loc[ds, "phi_magnitude"].astype(float).values
    x = np.arange(len(ds))
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.bar(x, env, label="envelope (kinetic+stability)", color="tab:blue")
    ax.bar(x, mag, bottom=env, label="magnitude (allocation+catalysis+maintenance)", color="tab:orange")
    ax.set_xticks(x); ax.set_xticklabels(ds, rotation=20, ha="right")
    ax.set_ylabel("Shapley variance share"); ax.set_ylim(0, 1)
    ax.set_title("Tuned-model decomposition: envelope vs magnitude (rich BHI)")
    ax.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.32), ncol=1)
    fig.tight_layout(); p = os.path.join(out_dir, fname); fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_iqr_bands(group_curves, base_curve, grid, out_dir, fname="decomp_iqr_bands.png"):
    plt = _mpl()
    grid = np.asarray(grid, float)
    groups = list(group_curves)
    fig, axes = plt.subplots(1, len(groups), figsize=(3.1 * len(groups), 3.6), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, gname in zip(axes, groups):
        C = group_curves[gname]
        lo, md, hi = np.percentile(C, [25, 50, 75], axis=0)
        ax.fill_between(grid, lo, hi, color=_GROUP_COLOR[gname], alpha=0.3)
        ax.plot(grid, md, color=_GROUP_COLOR[gname], lw=2)
        ax.plot(grid, base_curve, "k--", lw=1.2)
        ax.set_title(gname, fontsize=10); ax.set_xlabel("T (°C)")
    axes[0].set_ylabel("Growth rate (1/h)")
    fig.suptitle("Reachable TPC per lever group (median + IQR; dashed = tuned baseline)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(out_dir, fname); fig.savefig(p, dpi=150); plt.close(fig)
    return p


# --- orchestrator ------------------------------------------------------------
def run(strain: str, v3_dir: str, out_base: str, *, grid=None, h=0.10, H=0.20,
        n_draws=150, n_env=24, n_mag=24, n_band=41, n_proc=0, allow_glpk=False,
        seed=1) -> Dict:
    """Elasticity (posterior-propagated) + envelope-vs-magnitude decomposition on
    the TUNED model. Writes elasticity_tuned/ and decompose_tuned/ under out_base."""
    from multiprocessing import Pool
    from .calibration_multi import _set_default_solver

    grid = np.linspace(5.0, 52.0, 48) if grid is None else np.asarray(grid, float)
    solver = _set_default_solver("gurobi")
    if solver != "gurobi" and not allow_glpk:
        raise SystemExit("[solver] Gurobi NOT active - stopping (set ALLOW_GLPK to override).")
    print(f"[solver] {solver}")

    base, med = tuned_pert_from_v3(v3_dir)
    draws = load_draws(v3_dir, n_draws, seed) if n_draws else None
    n_proc = n_proc or max(1, min(10, (os.cpu_count() or 2) - 2))

    el_dir = os.path.join(out_base, "elasticity_tuned")
    de_dir = os.path.join(out_base, "decompose_tuned")
    os.makedirs(el_dir, exist_ok=True); os.makedirs(de_dir, exist_ok=True)

    pool = Pool(processes=n_proc, initializer=_winit, initargs=(strain, grid, 30.0))
    try:
        refs = _refs_from_worker(pool)
        # pre-flight: base descriptors must reproduce the v3 headline
        base_desc = pool.map(_wdesc, [base])[0]
        print(f"[dissect] tuned baseline (rich BHI): rmax={base_desc['rmax']:.3f} "
              f"Topt={base_desc['Topt_C']:.1f} CTmax={base_desc['CTmax_C']:.1f} Ea={base_desc['Ea_eV']:.3f}")

        print(f"[dissect] elasticity: base + {len(draws or [])} posterior draws, {n_proc} workers")
        E, base_desc, bands = run_elasticity(pool, base, refs, h, draws=draws)
        E.to_csv(os.path.join(el_dir, "elasticity_table.csv"))
        with open(os.path.join(el_dir, "elasticity_bands.json"), "w") as fh:
            json.dump({"h": h, "median_descriptors": base_desc, "bands": bands,
                       "note": "bands = [5%,50%,95%] of E[D,p] over posterior draws"}, fh, indent=2)
        plot_elasticity_heatmap(E, el_dir)
        for D in ("rmax", "CTmax_C", "Topt_C", "Ea_eV"):
            plot_elasticity_tornado(E, bands, el_dir, D=D)

        print(f"[dissect] decomposition: {n_env}x{n_mag} crossed grid + {len(GROUPS)} marginal bands")
        var_df, mag_df, group_curves, base_curve = run_decomposition(
            pool, base, refs, H, grid, n_env=n_env, n_mag=n_mag, n_band=n_band, seed=seed)
    finally:
        pool.close(); pool.join()

    var_df.to_csv(os.path.join(de_dir, "decomposition_variance.csv"), index_label="descriptor")
    mag_df.to_csv(os.path.join(de_dir, "decomposition_iqr_magnitude.csv"), index_label="descriptor")
    np.save(os.path.join(de_dir, "recast_temps_C.npy"), grid)
    for gname, C in group_curves.items():
        np.save(os.path.join(de_dir, f"recast_curves_{gname}.npy"), C)
    np.save(os.path.join(de_dir, "recast_base_curve.npy"), base_curve)
    plot_variance_partition(var_df, de_dir)
    plot_iqr_bands(group_curves, base_curve, grid, de_dir)
    summ = {"operating_point": "rich (BHI), growth law ON, reconciled single pool; TUNED (v3 medians)",
            "h": h, "H": H, "tuned_medians": med,
            "baseline_descriptors": {k: round(float(base_desc[k]), 4) for k in base_desc},
            "variance_shares": json.loads(var_df.to_json(orient="index")),
            "iqr_magnitude": json.loads(mag_df.to_json(orient="index"))}
    with open(os.path.join(de_dir, "decompose_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    return {"elasticity_dir": el_dir, "decompose_dir": de_dir,
            "baseline_descriptors": base_desc, "elasticity": E, "variance": var_df, "iqr": mag_df}


def _refs_from_worker(pool):
    """Fetch the per-enzyme Topt/Tm reference scales from a worker's provider."""
    return pool.apply(_worker_refs)


def _worker_refs():
    return _refs(_PMR)


# --- identifiability (per-enzyme control) on the tuned model -----------------
def run_control_tuned(strain: str, v3_dir: str, out_dir: str, *, screen_top_k=80,
                      grid_n=95, allow_glpk=False) -> str:
    """Per-enzyme thermal control + identifiability about the TUNED rich baseline.
    Serial (mutates the shared model per enzyme); Gurobi keeps each solve ~ms."""
    from .calibration_multi import _set_default_solver, _build_pm_rich
    from . import control
    solver = _set_default_solver("gurobi")
    if solver != "gurobi" and not allow_glpk:
        raise SystemExit("[solver] Gurobi NOT active - stopping (set ALLOW_GLPK to override).")
    print(f"[solver] {solver}")
    base, _ = tuned_pert_from_v3(v3_dir)
    pm = _build_pm_rich(strain)
    try:
        pm.ec.model.solver.configuration.timeout = 30
    except Exception:
        pass
    grid = np.linspace(5.0, 52.0, grid_n)   # fine enough for argmax-Topt descriptors
    res = control.run_control(pm, grid, screen_top_k=screen_top_k, base_pert=base, progress=True)
    os.makedirs(out_dir, exist_ok=True)
    res.save(out_dir, no_plots=False)
    return out_dir

