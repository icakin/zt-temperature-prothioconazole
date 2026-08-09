"""Per-enzyme thermal control coefficients + identifiability (H1.1 / H1.2).

Which individual enzymes' thermal parameters set the organismal TPC envelope
(Topt, CT_max), which enzymes' capacity sets the baseline rate (rmax/B0), and --
the identifiability flip side -- which parameters the growth TPC is *insensitive*
to and therefore cannot be inferred from growth data alone (they need proteome /
flux data). Strain-level diagnostic; no sweep experiment required.

Two control-coefficient families
--------------------------------
* Thermal control (envelope), central finite difference on one enzyme's Topt/dCp:
      CC[D, Topt_i] = (D(Topt_i+dT) - D(Topt_i-dT)) / (2 dT)
  for organismal descriptors D in {Topt_C, CT_max_C, niche_width_C, ...}.
* Rate / flux control (magnitude) at temperature T:
      FCC_i(T) = d ln mu / d ln kcat_i
  computed by perturbing base_cost (lowering base_cost == raising kcat). Only
  enzymes with nonzero usage at T can have nonzero FCC.

Two-stage design (tractable on ~2500 enzymes)
---------------------------------------------
Stage A -- cheap screen over ALL enzymes from the nominal solution: usage share
u_i(T) = cost_i(T)|v_i(T)| / budget and analytic cost sensitivity
s_i(T) = d ln cost_i/dT = -d ln(relative_kcat_i)/dT. Rank; take top-k.
Stage B -- targeted finite differences (reusing the mutate-entry + refresh_params
pattern from dltkcat.apply_fits_to_provider) on the top candidates only.

Everything reuses EnzymeConstrainedModel, compute_tpc, TPC.descriptors and mmrt
unchanged. First-order control-magnitude proxy for identifiability -- NOT a full
Fisher-information / profile-likelihood analysis.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .enzyme_cost import Perturbation
from .mmrt import relative_kcat_vec
from .tpc import compute_tpc

TEMP_DESCRIPTORS = ["Topt_C", "CT_max_C", "niche_width_C"]
# TPC.descriptors names CT_max as CTmax_C; map the config-facing name to it.
_DESC_ALIAS = {"CT_max_C": "CTmax_C", "rmax": "rmax", "Topt_C": "Topt_C",
               "niche_width_C": "niche_width_C", "Ea_eV": "Ea_eV"}


def _desc(pm, tpc_temps_C, crit_frac, base=None):
    d = compute_tpc(pm, tpc_temps_C, base or Perturbation()).descriptors(crit_frac).as_dict()
    d["CT_max_C"] = d.get("CTmax_C", np.nan)   # expose the proposal-facing name
    return d


def _apply(ec, Tk, base=None):
    """Set the operating point at temperature Tk about the baseline perturbation:
    sector allocation (incl. kappa/sigma) when sectors are wired, else the scalar
    pool. Baseline None -> the identity/nominal point (unchanged behaviour)."""
    base = base or Perturbation()
    ec.set_temperature(Tk, base)
    if getattr(ec, "_sectors", None) is not None and base.uses_allocation():
        ec.set_allocation(base.f_metab, base.f_maint,
                          kappa_scale=base.kappa_scale, sigma_sat=base.sigma_sat)
    else:
        ec.set_budget(ec.default_budget)


def _analysis_temps(nom, tpc_temps_C):
    T = np.asarray(tpc_temps_C, float)
    topt = nom.get("Topt_C", np.nan)
    ctmax = nom.get("CT_max_C", nom.get("CTmax_C", np.nan))
    if not np.isfinite(topt):
        topt = float(np.median(T))
    supra = min(ctmax - 2.0, topt + 8.0) if np.isfinite(ctmax) else topt + 8.0
    cand = [topt - 10.0, topt, supra]
    lo, hi = float(T.min()), float(T.max())
    return [float(np.clip(round(t, 2), lo, hi)) for t in cand]


# ---------------------------------------------------------------------------
# result container
# ---------------------------------------------------------------------------
@dataclass
class ControlResult:
    analysis_temps_C: List[float]
    usage: pd.DataFrame            # u_i(T) per used enzyme, per analysis T
    thermal: pd.DataFrame          # thermal CCs + screen + rank
    rate: pd.DataFrame             # FCC_i(T) per enzyme
    identifiability: pd.DataFrame  # per (enzyme, param): ident, flag, top descriptor
    summary: dict

    def save(self, out_dir, no_plots=False):
        os.makedirs(out_dir, exist_ok=True)
        self.usage.to_csv(os.path.join(out_dir, "usage_by_temperature.csv"), index=False)
        self.thermal.to_csv(os.path.join(out_dir, "thermal_control.csv"), index=False)
        self.rate.to_csv(os.path.join(out_dir, "rate_control.csv"), index=False)
        self.identifiability.to_csv(os.path.join(out_dir, "identifiability.csv"), index=False)
        with open(os.path.join(out_dir, "summary.json"), "w") as fh:
            json.dump(self.summary, fh, indent=2)
        if not no_plots:
            try:
                plot_all(self, out_dir)
            except Exception as e:
                print(f"[control] plotting skipped ({e})")


# ---------------------------------------------------------------------------
# Stage A: usage + analytic thermal sensitivity (no per-enzyme re-solve)
# ---------------------------------------------------------------------------
def _usage_at(ecm, Tk, h=0.25, base=None):
    """Return (usage array u_i, cost sensitivity s_i = d ln cost/dT) at Tk, about
    the baseline operating point (sector allocation when wired)."""
    base = base or Perturbation()
    _apply(ecm, Tk, base)
    sol = ecm.model.optimize()
    ents = ecm.table.entries
    if sol.status != "optimal":
        return np.zeros(len(ents)), np.zeros(len(ents))
    v = np.array([abs(sol.fluxes.get(e.rxn_id, 0.0)) for e in ents])
    c = ecm._costs(Tk, base)
    budget = float(ecm._pool.ub) if getattr(ecm, "_sectors", None) is not None else ecm.default_budget
    u = c * v / max(budget, 1e-12)
    r1 = relative_kcat_vec(Tk + h, ecm._T0, ecm._Topt, ecm._dCp)
    r2 = relative_kcat_vec(Tk - h, ecm._T0, ecm._Topt, ecm._dCp)
    s = -(np.log(r1) - np.log(r2)) / (2 * h)   # d ln cost/dT
    return u, s


# ---------------------------------------------------------------------------
# Stage B: finite-difference control coefficients (mutate entry + refresh)
# ---------------------------------------------------------------------------
def _thermal_cc(pm, entry, attr, step, fractional, tpc_temps_C, descriptors, crit_frac, base=None):
    ec = pm.ec
    orig = getattr(entry, attr)
    plus = orig * (1 + step) if fractional else orig + step
    minus = orig * (1 - step) if fractional else orig - step
    setattr(entry, attr, plus); ec.refresh_params(); dP = _desc(pm, tpc_temps_C, crit_frac, base)
    setattr(entry, attr, minus); ec.refresh_params(); dM = _desc(pm, tpc_temps_C, crit_frac, base)
    setattr(entry, attr, orig); ec.refresh_params()
    denom = plus - minus
    return {D: (dP.get(D, np.nan) - dM.get(D, np.nan)) / denom for D in descriptors}


def _fcc(pm, entry, Tk, f, base=None):
    ec = pm.ec
    bc = entry.base_cost
    entry.base_cost = bc / (1 + f); ec.refresh_params()
    _apply(ec, Tk, base)
    gp = ec.model.slim_optimize()
    entry.base_cost = bc / (1 - f); ec.refresh_params()
    _apply(ec, Tk, base)
    gm = ec.model.slim_optimize()
    entry.base_cost = bc; ec.refresh_params()
    if gp and gm and np.isfinite(gp) and np.isfinite(gm) and gp > 0 and gm > 0:
        return (np.log(gp) - np.log(gm)) / (2 * f)
    return np.nan


# ---------------------------------------------------------------------------
# runner
# ---------------------------------------------------------------------------
def run_control(pm, tpc_temps_C, perturb=None, screen_top_k=100,
                descriptors=("Topt_C", "CT_max_C", "niche_width_C", "rmax", "Ea_eV"),
                identifiable_threshold=0.05, analysis_temps_C=None,
                crit_frac=0.05, progress=True, base_pert=None) -> ControlResult:
    perturb = perturb or {}
    dT = float(perturb.get("Topt_K", 1.0))
    dCp_frac = float(perturb.get("dCp_frac", 0.1))
    kcat_frac = float(perturb.get("kcat_frac", 0.02))
    ec = pm.ec
    ents = ec.table.entries
    ids = [e.rxn_id for e in ents]
    enz = [e.enzyme_id for e in ents]
    grp = [e.group for e in ents]

    nom = _desc(pm, tpc_temps_C, crit_frac, base_pert)
    temps = analysis_temps_C or _analysis_temps(nom, tpc_temps_C)
    temps = [float(t) for t in temps]
    labels = ["sub", "opt", "supra"][:len(temps)] + \
             [f"T{i}" for i in range(3, len(temps))]
    if progress:
        print(f"[control] nominal Topt={nom.get('Topt_C'):.1f}C "
              f"CTmax={nom.get('CT_max_C'):.1f}C; analysis T={temps}")

    # -- Stage A: usage + sensitivity at each analysis temperature --
    U, S = {}, {}
    for lab, tc in zip(labels, temps):
        u, s = _usage_at(ec, tc + 273.15, base=base_pert)
        U[lab], S[lab] = u, s
    hot = labels[-1]
    rate_screen = np.max(np.vstack([U[l] for l in labels]), axis=0)
    thermal_screen = U[hot] * np.clip(S[hot], 0, None)

    usage_rows = []
    for i, e in enumerate(ents):
        if max(U[l][i] for l in labels) <= 1e-9:
            continue
        row = {"rxn_id": ids[i], "enzyme_id": enz[i], "group": grp[i]}
        for l, tc in zip(labels, temps):
            row[f"u_{l}_{tc:g}C"] = float(U[l][i])
        usage_rows.append(row)
    usage_df = pd.DataFrame(usage_rows)

    thermal_rank = np.argsort(-thermal_screen)
    thermal_cand = [i for i in thermal_rank if thermal_screen[i] > 0][:screen_top_k]
    rate_rank = np.argsort(-rate_screen)
    rate_cand = [i for i in rate_rank if rate_screen[i] > 1e-9][:screen_top_k]
    if progress:
        print(f"[control] Stage B: {len(thermal_cand)} thermal candidates, "
              f"{len(rate_cand)} rate candidates")

    # -- Stage B thermal: CC of temperature descriptors to Topt_i, dCp_i --
    temp_descs = [d for d in descriptors if d in TEMP_DESCRIPTORS]
    trows = []
    for n, i in enumerate(thermal_cand):
        e = ents[i]
        cc_topt = _thermal_cc(pm, e, "Topt", dT, False, tpc_temps_C, temp_descs, crit_frac, base_pert)
        cc_dcp = _thermal_cc(pm, e, "dCp", dCp_frac, True, tpc_temps_C, temp_descs, crit_frac, base_pert)
        row = {"rxn_id": ids[i], "enzyme_id": enz[i], "group": grp[i],
               "thermal_screen": float(thermal_screen[i])}
        for D in temp_descs:
            row[f"CC[{D},Topt_i]"] = cc_topt[D]
            row[f"CC[{D},dCp_i]"] = cc_dcp[D]
        trows.append(row)
        if progress and (n + 1) % max(1, len(thermal_cand) // 5) == 0:
            print(f"[control] thermal CC {n+1}/{len(thermal_cand)}")
    thermal_df = pd.DataFrame(trows)
    if len(thermal_df):
        thermal_df = thermal_df.sort_values(
            f"CC[{temp_descs[0]},Topt_i]", key=lambda s: s.abs(), ascending=False)
        thermal_df.insert(0, "rank", range(1, len(thermal_df) + 1))

    # -- Stage B rate: FCC_i(T) for rate candidates at each analysis T --
    rrows = []
    for n, i in enumerate(rate_cand):
        e = ents[i]
        row = {"rxn_id": ids[i], "enzyme_id": enz[i], "group": grp[i]}
        for l, tc in zip(labels, temps):
            row[f"FCC_{l}_{tc:g}C"] = _fcc(pm, e, tc + 273.15, kcat_frac, base_pert) \
                if U[l][i] > 1e-9 else np.nan
        rrows.append(row)
        if progress and (n + 1) % max(1, len(rate_cand) // 5) == 0:
            print(f"[control] FCC {n+1}/{len(rate_cand)}")
    rate_df = pd.DataFrame(rrows)

    # -- summation checks --
    opt_lab = labels[1] if len(labels) > 1 else labels[0]
    fcc_opt_col = [c for c in rate_df.columns if c.startswith(f"FCC_{opt_lab}_")]
    sum_fcc = float(rate_df[fcc_opt_col[0]].sum()) if fcc_opt_col and len(rate_df) else np.nan
    topt_cc_col = f"CC[Topt_C,Topt_i]"
    sum_thermal = float(thermal_df[topt_cc_col].sum()) \
        if len(thermal_df) and topt_cc_col in thermal_df else np.nan

    # -- identifiability map (PROTEOME-WIDE: one row per enzyme x parameter) --
    ident_df = _identifiability(ents, ids, enz, thermal_screen, rate_screen,
                                thermal_df, rate_df, temp_descs, opt_lab,
                                identifiable_threshold)

    # Report BOTH the proteome-wide identifiable fraction (over all 3 params per
    # enzyme -- expected small: most parameters cannot be inferred from growth
    # alone, the etc-GEM under-determination that motivates omics) and, separately,
    # the fraction among the top-K finite-difference control enzymes (the
    # `refined` rows -- expected high, since those are selected as high-control).
    n_total = int(len(ident_df))                       # == 3 * n_enzymes
    n_ident = int(ident_df["identifiable_from_growth"].sum()) if n_total else 0
    refined = ident_df[ident_df["refined"]] if n_total else ident_df
    n_topk = int(len(refined))
    n_ident_topk = int(refined["identifiable_from_growth"].sum()) if n_topk else 0
    summary = {
        "analysis_temps_C": temps,
        "nominal": {k: nom.get(k) for k in ("Topt_C", "CT_max_C", "rmax", "niche_width_C")},
        "summation_check": {
            "sum_FCC_at_opt": sum_fcc,
            "sum_CC_Topt_org_wrt_Topt_i": sum_thermal,
        },
        "top_thermal_determinants": (thermal_df.head(10)[["rxn_id", "enzyme_id", topt_cc_col]]
                                     .to_dict("records") if len(thermal_df) else []),
        "top_rate_limiting_at_opt": (rate_df.reindex(
            rate_df[fcc_opt_col[0]].abs().sort_values(ascending=False).index)
            .head(10)[["rxn_id", "enzyme_id", fcc_opt_col[0]]].to_dict("records")
            if fcc_opt_col and len(rate_df) else []),
        "identifiability": {
            "method": "first-order control-magnitude proxy (screen), refined by "
                      "finite difference on top-K; NOT Fisher information",
            "n_params_total": n_total,
            # proteome-wide (headline: small)
            "n_identifiable_from_growth": n_ident,
            "frac_identifiable": (n_ident / n_total) if n_total else float("nan"),
            "frac_requires_omics": (1 - n_ident / n_total) if n_total else float("nan"),
            # among top-K finite-difference control enzymes (refined rows: high)
            "n_refined_topk": n_topk,
            "n_identifiable_topk": n_ident_topk,
            "frac_identifiable_topk": (n_ident_topk / n_topk) if n_topk else float("nan"),
            "threshold": identifiable_threshold,
        },
    }
    return ControlResult(temps, usage_df, thermal_df, rate_df, ident_df, summary)


def _norm_abs(series):
    m = series.abs().max()
    return series.abs() / m if m and np.isfinite(m) and m > 0 else series.abs() * 0.0


def _fd_ident_lookup(thermal_df, rate_df, temp_descs, opt_lab):
    """Finite-difference identifiability, keyed by (rxn_id, parameter).

    Reproduces the previous (top-K only) metric: each descriptor's CC is
    normalised by its max |CC| across the analysed enzymes (so descriptors are
    comparable) and ident = max_D |CC_norm|; for kcat_i the FCC at the optimum.
    Returns {(rxn_id, param): (ident, top_descriptor)}.
    """
    out = {}
    if len(thermal_df):
        for pname in ("Topt_i", "dCp_i"):
            cols = [f"CC[{D},{pname}]" for D in temp_descs if f"CC[{D},{pname}]" in thermal_df]
            if not cols:
                continue
            norm = pd.DataFrame({c: _norm_abs(thermal_df[c]) for c in cols})
            for k in range(len(thermal_df)):
                vals = norm.iloc[k]
                ident = float(vals.max()) if len(vals) else 0.0
                top = temp_descs[int(np.argmax(vals.values))] if len(vals) else None
                out[(thermal_df.iloc[k]["rxn_id"], pname)] = (ident, top)
    fcc_col = [c for c in rate_df.columns if c.startswith(f"FCC_{opt_lab}_")]
    if fcc_col and len(rate_df):
        norm = _norm_abs(rate_df[fcc_col[0]])
        for k in range(len(rate_df)):
            ident = float(norm.iloc[k]) if np.isfinite(norm.iloc[k]) else 0.0
            out[(rate_df.iloc[k]["rxn_id"], "kcat_i")] = (ident, "rmax")
    return out


def _identifiability(ents, ids, enz, thermal_screen, rate_screen,
                     thermal_df, rate_df, temp_descs, opt_lab, threshold):
    """PROTEOME-WIDE per (enzyme, parameter) identifiability proxy + flag.

    A first-order control-magnitude proxy (NOT a Fisher-information / profile-
    likelihood analysis): for every enzyme and each of its thermal/rate
    parameters, the identifiability score is the cheap Stage-A screen quantity
    (usage share x analytic parameter sensitivity), max-normalised across the
    proteome so it is comparable to the finite-difference metric. Where a top-K
    finite-difference control coefficient exists it replaces the proxy (more
    accurate) and the row is marked ``refined = True``.

    Parameters
    ----------
    thermal_screen : array over ALL enzymes -- usage x thermal sensitivity
        (shared first-order proxy for the thermal-envelope params Topt_i, dCp_i).
    rate_screen : array over ALL enzymes -- max usage share (proxy for kcat_i).
    """
    ts = np.asarray(thermal_screen, float)
    rs = np.asarray(rate_screen, float)
    ts_max = ts.max() if ts.size and ts.max() > 0 else 1.0
    rs_max = rs.max() if rs.size and rs.max() > 0 else 1.0
    t_norm = ts / ts_max
    r_norm = rs / rs_max
    fd = _fd_ident_lookup(thermal_df, rate_df, temp_descs, opt_lab)

    rows = []
    # thermal-envelope parameters: shared usage x thermal-sensitivity proxy
    for pname in ("Topt_i", "dCp_i"):
        for i, e in enumerate(ents):
            ident, top, refined = float(t_norm[i]), "CT_max_C", False
            if (ids[i], pname) in fd:
                ident, top = fd[(ids[i], pname)]
                refined = True
            rows.append({"rxn_id": ids[i], "enzyme_id": enz[i], "parameter": pname,
                         "ident": ident, "top_descriptor": top,
                         "identifiable_from_growth": ident > threshold,
                         "refined": refined})
    # catalytic capacity: usage-share proxy, refined by FCC where available
    for i, e in enumerate(ents):
        ident, top, refined = float(r_norm[i]), "rmax", False
        if (ids[i], "kcat_i") in fd:
            ident, top = fd[(ids[i], "kcat_i")]
            refined = True
        rows.append({"rxn_id": ids[i], "enzyme_id": enz[i], "parameter": "kcat_i",
                     "ident": ident, "top_descriptor": top,
                     "identifiable_from_growth": ident > threshold,
                     "refined": refined})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def _mpl():
    from .plotting import _mpl as m
    return m()


def plot_thermal_control(res, out_dir, fname="thermal_control_bar.png", topk=15):
    plt = _mpl()
    df = res.thermal
    if not len(df):
        return None
    c_topt = "CC[Topt_C,Topt_i]"
    c_ctmax = "CC[CT_max_C,Topt_i]"
    df = df.reindex(df[c_topt].abs().sort_values(ascending=False).index).head(topk)
    x = np.arange(len(df)); w = 0.4
    fig, ax = plt.subplots(figsize=(1.1 * len(df) + 2, 5))
    ax.bar(x - w / 2, df[c_topt].abs(), w, label="|CC[Topt_org, Topt_i]|", color="tab:blue")
    if c_ctmax in df:
        ax.bar(x + w / 2, df[c_ctmax].abs(), w, label="|CC[CT_max, Topt_i]|", color="tab:red")
    ax.set_xticks(x)
    ax.set_xticklabels(df["rxn_id"], rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("|thermal control coefficient|")
    ax.set_title("Top enzymes controlling the TPC envelope")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(out_dir, fname); fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_bottleneck_vs_temperature(res, out_dir, fname="bottleneck_vs_temperature.png", topn=6):
    plt = _mpl()
    df = res.usage
    if not len(df):
        return None
    ucols = [c for c in df.columns if c.startswith("u_")]
    temps = res.analysis_temps_C
    top = df.reindex(df[ucols].max(axis=1).sort_values(ascending=False).index).head(topn)
    fig, ax = plt.subplots(figsize=(7, 5))
    for _, r in top.iterrows():
        ax.plot(temps, [r[c] for c in ucols], marker="o", label=r["rxn_id"])
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("proteome usage share u_i(T)")
    ax.set_title("Limiting-enzyme usage across temperature")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    p = os.path.join(out_dir, fname); fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_identifiability_hist(res, out_dir, fname="identifiability_hist.png"):
    plt = _mpl()
    df = res.identifiability
    if not len(df):
        return None
    thr = res.summary["identifiability"]["threshold"]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(df["ident"].replace([np.inf, -np.inf], np.nan).dropna(), bins=30,
            color="tab:blue", alpha=0.8)
    ax.axvline(thr, color="k", ls="--", label=f"threshold={thr}")
    ax.set_xlabel("identifiability score  ident_i = max_D |CC_norm|")
    ax.set_ylabel("count (enzyme × parameter)")
    ax.set_title("Identifiability from the growth TPC")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(out_dir, fname); fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_all(res, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    return [f for f in (plot_thermal_control(res, out_dir),
                        plot_bottleneck_vs_temperature(res, out_dir),
                        plot_identifiability_hist(res, out_dir)) if f]


# ---------------------------------------------------------------------------
# CLI entry (called by cli.cmd_control)
# ---------------------------------------------------------------------------
def run(cfg, out_dir, no_plots=False):
    from .config import build_provider
    ctrl = cfg.get("control", {}) or {}
    pm = build_provider(cfg)
    try:
        pm.ec.model.solver.configuration.timeout = int(cfg.get("solver_timeout", 10))
    except Exception:
        pass
    # Envelope descriptors (esp. argmax-based Topt_C) need a fine TPC grid to
    # register single-enzyme shifts; refine the strain grid by grid_refine.
    g = cfg["temperature_grid"]
    refine = int(ctrl.get("grid_refine", 4))
    n = (int(g["n"]) - 1) * refine + 1
    tpc_temps = np.linspace(g["start_C"], g["stop_C"], n)
    res = run_control(
        pm, tpc_temps,
        perturb=ctrl.get("perturb"),
        screen_top_k=int(ctrl.get("screen_top_k", 100)),
        descriptors=tuple(ctrl.get("descriptors",
                          ["Topt_C", "CT_max_C", "niche_width_C", "rmax", "Ea_eV"])),
        identifiable_threshold=float(ctrl.get("identifiable_threshold", 0.05)),
        analysis_temps_C=ctrl.get("temperatures_C"),
        crit_frac=cfg.get("crit_frac", 0.05))
    res.save(out_dir, no_plots=no_plots)
    sc = res.summary["summation_check"]
    idf = res.summary["identifiability"]
    print(f"[control] summation: sum_FCC(opt)={sc['sum_FCC_at_opt']:.3g} "
          f"sum_CC_Topt={sc['sum_CC_Topt_org_wrt_Topt_i']:.3g}")
    print(f"[control] identifiable-from-growth (proteome-wide) "
          f"{idf['n_identifiable_from_growth']}/{idf['n_params_total']} params "
          f"({idf['frac_identifiable']:.1%}); among top-K control "
          f"{idf['n_identifiable_topk']}/{idf['n_refined_topk']} "
          f"({idf['frac_identifiable_topk']:.1%})")
    print(f"[control] wrote {out_dir}")
    return out_dir
