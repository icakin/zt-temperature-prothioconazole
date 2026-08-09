"""Allocation vs envelope variance decomposition of the TPC (hypothesis H1.3).

A single genome can generate a *range* of thermal performance curves by (a)
reshaping the per-enzyme thermal envelope (genome-set params: dTopt, topt_scale,
dCp_scale) and (b) reallocating its proteome (allocation-set params:
budget_scale, per-group allocations). H1.3 says these separate: the envelope
sets the *shape/position* of the TPC, allocation sets its *magnitude*.

This module quantifies that split with a two-group functional-ANOVA / Shapley
variance decomposition on a **crossed** design, reusing the existing engine
(`Perturbation`, `compute_tpc`, `TPC.descriptors`, `sensitivity._lhs`) unchanged.

Design
------
Draw M allocation samples (envelope held at nominal) and N envelope samples
(allocation at nominal) by Latin hypercube, then evaluate every crossed pair
(allocation_i, envelope_j) -> one Perturbation -> one TPC -> its descriptors.
That gives an M x N matrix ``f_ij`` per descriptor.

Decomposition (per descriptor, over finite cells)
-------------------------------------------------
    mu   = mean(f_ij)
    a_i  = mean_j f_ij - mu                 # allocation main effect
    e_j  = mean_i f_ij - mu                 # envelope main effect
    g_ij = f_ij - mu - a_i - e_j            # interaction
    V_A  = mean_i a_i^2 ; V_E = mean_j e_j^2 ; V_AE = mean_ij g_ij^2
    V    = V_A + V_E + V_AE
    S_A, S_E, S_AE = V_A/V, V_E/V, V_AE/V           # grouped Sobol fractions (sum 1)
    T_A, T_E       = S_A + S_AE, S_E + S_AE          # total-effect indices
    phi_A, phi_E   = S_A + S_AE/2, S_E + S_AE/2      # Shapley (exact for 2 groups)

The fractions are relative to the chosen input distributions (uniform over the
configured ranges) -- a structural, in-silico decomposition of what the model
*can* generate, not a claim about real cells.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .enzyme_cost import Perturbation
from .sensitivity import _lhs, _make_perturbation
from .tpc import compute_tpc

# Descriptors decomposed (magnitude: rmax; temperature/shape: the rest).
DESCRIPTORS = ["Topt_C", "rmax", "CTmin_C", "CTmax_C", "niche_width_C",
               "B80_C", "Ea_eV", "skewness"]
KEY_DESCRIPTORS = ["Topt_C", "rmax", "CTmax_C", "niche_width_C", "Ea_eV"]


# ---------------------------------------------------------------------------
# core math
# ---------------------------------------------------------------------------
def decompose_grid(F: np.ndarray) -> Dict[str, float]:
    """Two-group functional-ANOVA / Shapley split of an M x N descriptor grid.

    Drops fully-NaN rows/columns, then works over the remaining finite cells.
    Fractions sum to 1 by construction; returns NaNs if there is no variance.
    """
    F = np.asarray(F, float)
    keep_r = ~np.all(np.isnan(F), axis=1)
    keep_c = ~np.all(np.isnan(F), axis=0)
    F = F[np.ix_(keep_r, keep_c)]
    n_valid = int(np.isfinite(F).sum())
    out = dict(V_A=np.nan, V_E=np.nan, V_AE=np.nan, V=np.nan,
               S_A=np.nan, S_E=np.nan, S_AE=np.nan, T_A=np.nan, T_E=np.nan,
               phi_A=np.nan, phi_E=np.nan, n_valid=n_valid,
               n_alloc=int(F.shape[0]) if F.ndim == 2 else 0,
               n_env=int(F.shape[1]) if F.ndim == 2 else 0)
    if F.size == 0 or n_valid == 0:
        return out
    mu = float(np.nanmean(F))
    a = np.nanmean(F, axis=1) - mu
    e = np.nanmean(F, axis=0) - mu
    g = F - mu - a[:, None] - e[None, :]
    V_A = float(np.nanmean(a ** 2))
    V_E = float(np.nanmean(e ** 2))
    V_AE = float(np.nanmean(g ** 2))
    V = V_A + V_E + V_AE
    out.update(V_A=V_A, V_E=V_E, V_AE=V_AE, V=V)
    if V <= 0:
        return out
    S_A, S_E, S_AE = V_A / V, V_E / V, V_AE / V
    out.update(S_A=S_A, S_E=S_E, S_AE=S_AE,
               T_A=S_A + S_AE, T_E=S_E + S_AE,
               phi_A=S_A + S_AE / 2.0, phi_E=S_E + S_AE / 2.0)
    return out


def _lhs_table(ranges: Dict[str, Sequence[float]], n: int, seed: int) -> pd.DataFrame:
    names = list(ranges)
    U = _lhs(n, len(names), seed)
    cols = {nm: ranges[nm][0] + U[:, k] * (ranges[nm][1] - ranges[nm][0])
            for k, nm in enumerate(names)}
    return pd.DataFrame(cols, columns=names)


def _group_names(allocation_ranges: Dict[str, Sequence[float]]) -> List[str]:
    """Per-group allocation multipliers are named alloc_<grp>."""
    return [n[len("alloc_"):] for n in allocation_ranges if n.startswith("alloc_")]


# ---------------------------------------------------------------------------
# result container
# ---------------------------------------------------------------------------
@dataclass
class DecompositionResult:
    temps_C: np.ndarray
    alloc_samples: pd.DataFrame
    env_samples: pd.DataFrame
    grids: Dict[str, np.ndarray]              # descriptor -> (M, N)
    table: pd.DataFrame                       # one row per descriptor
    alloc_only_curves: np.ndarray             # (M, n_temps), envelope nominal
    env_only_curves: np.ndarray               # (N, n_temps), allocation nominal
    nominal_curve: np.ndarray                 # (n_temps,)

    # -- summary --------------------------------------------------------------
    def summary(self) -> Dict[str, dict]:
        out = {}
        for d, row in self.table.iterrows():
            phi_a, phi_e = row["phi_A"], row["phi_E"]
            dom = ("allocation" if phi_a > phi_e else "envelope") \
                if np.isfinite(phi_a) and np.isfinite(phi_e) else "undetermined"
            out[d] = {"dominant_axis": dom,
                      "S_allocation": _f(row["S_A"]), "S_envelope": _f(row["S_E"]),
                      "S_interaction": _f(row["S_AE"]),
                      "phi_allocation": _f(phi_a), "phi_envelope": _f(phi_e),
                      "n_valid": int(row["n_valid"])}
        return out

    # -- persistence ----------------------------------------------------------
    def save(self, out_dir: str, no_plots: bool = False) -> List[str]:
        os.makedirs(out_dir, exist_ok=True)
        written = []
        p = os.path.join(out_dir, "decomposition_table.csv")
        self.table.to_csv(p, index_label="descriptor")
        written.append(p)

        npz = os.path.join(out_dir, "grids.npz")
        payload = {f"grid_{d}": g for d, g in self.grids.items()}
        payload["allocation_samples"] = self.alloc_samples.to_numpy()
        payload["allocation_names"] = np.array(list(self.alloc_samples.columns))
        payload["envelope_samples"] = self.env_samples.to_numpy()
        payload["envelope_names"] = np.array(list(self.env_samples.columns))
        np.savez(npz, **payload)
        written.append(npz)

        np.save(os.path.join(out_dir, "allocation_only_curves.npy"), self.alloc_only_curves)
        np.save(os.path.join(out_dir, "envelope_only_curves.npy"), self.env_only_curves)
        np.save(os.path.join(out_dir, "temps_C.npy"), self.temps_C)

        summ = os.path.join(out_dir, "summary.json")
        with open(summ, "w") as fh:
            json.dump(self.summary(), fh, indent=2)
        written.append(summ)

        if not no_plots:
            try:
                written += plot_all(self, out_dir)
            except Exception as e:
                print(f"[decompose] plotting skipped ({e})")
        return written


def _f(x):
    return None if x is None or not np.isfinite(x) else float(x)


# ---------------------------------------------------------------------------
# runner
# ---------------------------------------------------------------------------
def run_decomposition(pm, temps_C, allocation_ranges, envelope_ranges,
                      n_alloc: int, n_env: int, seed: int = 1,
                      crit_frac: float = 0.05,
                      descriptors: Sequence[str] = DESCRIPTORS,
                      envelope_samples=None, progress: bool = True) -> DecompositionResult:
    """Crossed allocation x envelope evaluation + per-descriptor decomposition.

    Reuses Perturbation / compute_tpc / TPC.descriptors unchanged; the descriptor
    table it feeds the ANOVA is generic, so it can later target respiration/CUE
    curves instead of growth without touching the math.

    ``envelope_samples`` (optional): a list of per-enzyme thermal samples
    (from thermal_sampling.sample_thermal). When given, the ENVELOPE axis is those
    per-enzyme draws (applied by mutating entries) instead of the dTopt/topt_scale/
    dCp_scale knobs; the crossed design and the ANOVA/Shapley math are unchanged.
    When None, behaviour is exactly the knobs path.
    """
    temps_C = np.asarray(temps_C, float)
    default_budget = pm.ec.default_budget
    group_names = _group_names(allocation_ranges)

    alloc = _lhs_table(allocation_ranges, n_alloc, seed)
    a_names = list(alloc.columns)

    if envelope_samples is not None:
        from .thermal_sampling import (nominal_thermal, apply_thermal_sample,
                                        restore_thermal)
        n_env = len(envelope_samples)
        env = pd.DataFrame({"Z": [s["Z"] for s in envelope_samples],
                            "Z2": [s["Z2"] for s in envelope_samples]})
        nomT, nomC = nominal_thermal(pm)
        grids = {d: np.full((n_alloc, n_env), np.nan) for d in descriptors}
        total = n_alloc * n_env
        env_only = np.full((n_env, len(temps_C)), np.nan)
        for j, s in enumerate(envelope_samples):          # apply each draw once
            apply_thermal_sample(pm, s)
            env_only[j] = compute_tpc(pm, temps_C, Perturbation()).growth
            for i in range(n_alloc):
                ai = {k: float(alloc.iloc[i][k]) for k in a_names}
                pert = _make_perturbation(ai, default_budget, group_names)
                desc = compute_tpc(pm, temps_C, pert).descriptors(crit_frac).as_dict()
                for d in descriptors:
                    grids[d][i, j] = desc.get(d, np.nan)
            if progress and (j + 1) % max(1, n_env // 8) == 0:
                print(f"[decompose] crossed grid {(j + 1) * n_alloc}/{total}")
        restore_thermal(pm, nomT, nomC)
        alloc_only = np.vstack([
            compute_tpc(pm, temps_C,
                        _make_perturbation({k: float(alloc.iloc[i][k]) for k in a_names},
                                           default_budget, group_names)).growth
            for i in range(n_alloc)])
    else:
        env = _lhs_table(envelope_ranges, n_env, seed + 1)
        e_names = list(env.columns)
        grids = {d: np.full((n_alloc, n_env), np.nan) for d in descriptors}
        total = n_alloc * n_env
        for i in range(n_alloc):
            ai = {k: float(alloc.iloc[i][k]) for k in a_names}
            for j in range(n_env):
                ej = {k: float(env.iloc[j][k]) for k in e_names}
                pert = _make_perturbation({**ai, **ej}, default_budget, group_names)
                desc = compute_tpc(pm, temps_C, pert).descriptors(crit_frac).as_dict()
                for d in descriptors:
                    grids[d][i, j] = desc.get(d, np.nan)
            if progress and (i + 1) % max(1, n_alloc // 8) == 0:
                print(f"[decompose] crossed grid {(i + 1) * n_env}/{total}")
        # marginal ensembles (envelope nominal / allocation nominal)
        alloc_only = np.vstack([
            compute_tpc(pm, temps_C,
                        _make_perturbation({k: float(alloc.iloc[i][k]) for k in a_names},
                                           default_budget, group_names)).growth
            for i in range(n_alloc)])
        env_only = np.vstack([
            compute_tpc(pm, temps_C,
                        _make_perturbation({k: float(env.iloc[j][k]) for k in e_names},
                                           default_budget, group_names)).growth
            for j in range(n_env)])
    nominal = compute_tpc(pm, temps_C, Perturbation()).growth

    table = pd.DataFrame({d: decompose_grid(grids[d]) for d in descriptors}).T
    table = table[["V_A", "V_E", "V_AE", "V", "S_A", "S_E", "S_AE",
                   "T_A", "T_E", "phi_A", "phi_E", "n_valid", "n_alloc", "n_env"]]
    return DecompositionResult(temps_C, alloc, env, grids, table,
                               alloc_only, env_only, nominal)


# ---------------------------------------------------------------------------
# plots (match plotting.py: Agg, dpi 150)
# ---------------------------------------------------------------------------
def _mpl():
    from .plotting import _mpl as m
    return m()


def plot_achievable_ranges(res, out_dir, fname="achievable_ranges.png"):
    plt = _mpl()
    T = res.temps_C
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, C, title in ((axes[0], res.alloc_only_curves, "Allocation-only (envelope nominal)"),
                         (axes[1], res.env_only_curves, "Envelope-only (allocation nominal)")):
        for row in C:
            ax.plot(T, row, color="0.7", lw=0.4, alpha=0.4)
        med = np.median(C, axis=0)
        q1, q3 = np.percentile(C, [25, 75], axis=0)
        ax.fill_between(T, q1, q3, color="tab:blue", alpha=0.25, label="IQR")
        ax.plot(T, med, color="tab:blue", lw=2, label="median")
        ax.plot(T, res.nominal_curve, "k--", lw=2, label="nominal")
        ax.set_xlabel("Temperature (°C)")
        ax.set_title(title)
        ax.legend(frameon=False)
    axes[0].set_ylabel("Growth rate (1/h)")
    fig.suptitle("H1.3 achievable TPC ranges: allocation vs envelope")
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def plot_variance_partition(res, out_dir, fname="variance_partition.png"):
    plt = _mpl()
    descs = [d for d in KEY_DESCRIPTORS if d in res.table.index]
    S_A = res.table.loc[descs, "S_A"].astype(float).values
    S_E = res.table.loc[descs, "S_E"].astype(float).values
    S_AE = res.table.loc[descs, "S_AE"].astype(float).values
    x = np.arange(len(descs))
    fig, ax = plt.subplots(figsize=(1.4 * len(descs) + 2, 5))
    ax.bar(x, S_A, label="allocation (S_A)", color="tab:orange")
    ax.bar(x, S_E, bottom=S_A, label="envelope (S_E)", color="tab:blue")
    ax.bar(x, S_AE, bottom=S_A + S_E, label="interaction (S_AE)", color="0.6")
    ax.set_xticks(x)
    ax.set_xticklabels(descs, rotation=30, ha="right")
    ax.set_ylabel("variance fraction")
    ax.set_ylim(0, 1)
    ax.set_title("Grouped Sobol variance partition per descriptor")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def plot_shapley_effects(res, out_dir, fname="shapley_effects.png"):
    plt = _mpl()
    descs = [d for d in res.table.index if np.isfinite(res.table.loc[d, "phi_A"])]
    phi_A = res.table.loc[descs, "phi_A"].astype(float).values
    phi_E = res.table.loc[descs, "phi_E"].astype(float).values
    x = np.arange(len(descs))
    w = 0.4
    fig, ax = plt.subplots(figsize=(1.4 * len(descs) + 2, 5))
    ax.bar(x - w / 2, phi_A, w, label="allocation (φ_A)", color="tab:orange")
    ax.bar(x + w / 2, phi_E, w, label="envelope (φ_E)", color="tab:blue")
    ax.set_xticks(x)
    ax.set_xticklabels(descs, rotation=30, ha="right")
    ax.set_ylabel("Shapley effect (φ_A + φ_E = 1)")
    ax.set_ylim(0, 1)
    ax.set_title("Shapley variance attribution per descriptor")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def plot_all(res, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    return [plot_achievable_ranges(res, out_dir),
            plot_variance_partition(res, out_dir),
            plot_shapley_effects(res, out_dir)]


# ===========================================================================
# RECAST decomposition: range-fair, nominal-centred, magnitude-aware, +Tm
# ===========================================================================
# The runner above uses hand-set, per-axis UNEQUAL ranges (envelope wide, sector
# allocation narrow) and, in the shipped experiment, sampled f_metab 0.20-0.40
# while the measured nominal is ~0.48 -- so it (a) inflated the envelope's share
# by giving it more room and (b) was not centred on the operating point. It also
# reports only variance SHARES, so a large share on a negligible movement (CTmax)
# masquerades as importance, and it never perturbs Tm.
#
# The recast below fixes all four: every parameter (allocation AND envelope) is
# moved by the SAME standardised relative half-width H about the MEASURED NOMINAL
# (multiplicative -> nominal*(1 +/- H); additive dTopt/dTm -> +/- H*reference_scale
# = H*sd of the per-enzyme Topt/Tm distribution; sector fractions -> nominal*(1 +/-
# H) with biosynthesis renormalised). The envelope axis gains dTm and is split into
# a KINETIC sub-group {dTopt, topt_scale, dCp_scale} and a STABILITY sub-group
# {dTm} via a fully-crossed 3-way (allocation x kinetic x stability) functional
# ANOVA -- which yields BOTH the 2-group allocation-vs-envelope Shapley (grouping
# kinetic+stability as envelope) AND the 3-group allocation/kinetic/stability
# Shapley in one design. Alongside the variance shares it reports a MAGNITUDE
# measure (achievable IQR / range of each descriptor when each axis alone is
# perturbed), so share and magnitude are never conflated.


def decompose_grid3(F: np.ndarray) -> Dict[str, float]:
    """Fully-crossed 3-way functional-ANOVA / Shapley split of an
    (n_A, n_K, n_S) descriptor grid over the allocation (A), kinetic-envelope (K)
    and stability (S) axes.

    Returns the seven variance components (V_A,V_K,V_S,V_AK,V_AS,V_KS,V_AKS), the
    2-group allocation-vs-envelope shares/Shapley (envelope = K u S, so its internal
    K-S interaction counts as an envelope main effect), and the 3-group
    allocation/kinetic/stability Shapley (all three sum to 1).
    """
    F = np.asarray(F, float)
    n_valid = int(np.isfinite(F).sum())
    keys = ["V_A", "V_K", "V_S", "V_AK", "V_AS", "V_KS", "V_AKS", "V",
            # 2-group allocation vs envelope(=K+S)
            "S_A", "S_E", "S_AE", "phi_A", "phi_E",
            # 3-group allocation / kinetic / stability
            "S_kin", "S_stab", "phi_kin", "phi_stab",
            "n_valid"]
    out = {k: np.nan for k in keys}
    out["n_valid"] = n_valid
    if F.ndim != 3 or F.size == 0 or n_valid == 0:
        return out
    mu = float(np.nanmean(F))
    A = np.nanmean(F, axis=(1, 2)) - mu           # (n_A,)
    K = np.nanmean(F, axis=(0, 2)) - mu           # (n_K,)
    S = np.nanmean(F, axis=(0, 1)) - mu           # (n_S,)
    AK = np.nanmean(F, axis=2) - mu - A[:, None] - K[None, :]       # (n_A,n_K)
    AS = np.nanmean(F, axis=1) - mu - A[:, None] - S[None, :]       # (n_A,n_S)
    KS = np.nanmean(F, axis=0) - mu - K[:, None] - S[None, :]       # (n_K,n_S)
    AKS = (F - mu - A[:, None, None] - K[None, :, None] - S[None, None, :]
           - AK[:, :, None] - AS[:, None, :] - KS[None, :, :])
    V_A = float(np.nanmean(A ** 2)); V_K = float(np.nanmean(K ** 2))
    V_S = float(np.nanmean(S ** 2)); V_AK = float(np.nanmean(AK ** 2))
    V_AS = float(np.nanmean(AS ** 2)); V_KS = float(np.nanmean(KS ** 2))
    V_AKS = float(np.nanmean(AKS ** 2))
    V = V_A + V_K + V_S + V_AK + V_AS + V_KS + V_AKS
    out.update(V_A=V_A, V_K=V_K, V_S=V_S, V_AK=V_AK, V_AS=V_AS,
               V_KS=V_KS, V_AKS=V_AKS, V=V)
    if V <= 0:
        return out
    # 2-group: envelope = kinetic u stability
    V_Egrp = V_K + V_S + V_KS
    V_AEgrp = V_AK + V_AS + V_AKS
    out["S_A"] = V_A / V
    out["S_E"] = V_Egrp / V
    out["S_AE"] = V_AEgrp / V
    out["phi_A"] = (V_A + V_AEgrp / 2.0) / V
    out["phi_E"] = (V_Egrp + V_AEgrp / 2.0) / V
    # 3-group Shapley (allocation / kinetic / stability), each sums with the rest to 1
    out["S_kin"] = V_K / V
    out["S_stab"] = V_S / V
    out["phi_A3"] = (V_A + (V_AK + V_AS) / 2.0 + V_AKS / 3.0) / V
    out["phi_kin"] = (V_K + (V_AK + V_KS) / 2.0 + V_AKS / 3.0) / V
    out["phi_stab"] = (V_S + (V_AS + V_KS) / 2.0 + V_AKS / 3.0) / V
    return out


def _recast_ranges(pm, H: float):
    """Nominal-centred +/- H ranges for the three axes. Returns (alloc, kin, stab)
    dicts of name -> (lo, hi) plus the reference scales used."""
    ec = pm.ec
    sd_topt = float(np.std(np.asarray(ec._Topt, float))) if getattr(ec, "_Topt", None) is not None \
        and len(ec._Topt) else 5.0
    sd_tm = float(np.std(np.asarray(ec._Tm, float))) if getattr(ec, "_Tm", None) is not None \
        and len(ec._Tm) else 5.0
    sec = getattr(ec, "_sectors", None)
    if not sec:
        raise SystemExit("recast decomposition needs proteome sectors active "
                         "(the allocation axis is f_metab/f_maint)")
    fm, fma = float(sec["f_metab_nom"]), float(sec["f_maint_nom"])
    alloc = {"f_metab": (fm * (1 - H), fm * (1 + H)),
             "f_maint": (fma * (1 - H), fma * (1 + H))}
    kin = {"dTopt": (-H * sd_topt, H * sd_topt),
           "topt_scale": (1 - H, 1 + H),
           "dCp_scale": (1 - H, 1 + H)}
    stab = {"dTm": (-H * sd_tm, H * sd_tm)}
    refs = {"H": H, "sd_Topt_K": round(sd_topt, 3), "sd_Tm_K": round(sd_tm, 3),
            "f_metab_nom": fm, "f_maint_nom": fma}
    return alloc, kin, stab, refs


@dataclass
class RecastResult:
    temps_C: np.ndarray
    H: float
    refs: Dict
    grids: Dict[str, np.ndarray]              # descriptor -> (n_A, n_K, n_S)
    table: pd.DataFrame                       # variance decomposition per descriptor
    magnitude: pd.DataFrame                   # descriptor x axis -> IQR & range
    nominal_desc: Dict[str, float]
    curves: Dict[str, np.ndarray]             # axis -> (n_band, n_temps)
    desc_marg: Dict[str, pd.DataFrame]        # axis -> per-sample descriptors
    nominal_curve: np.ndarray

    def summary(self) -> Dict[str, dict]:
        out = {}
        for d, row in self.table.iterrows():
            phi_a, phi_e = row["phi_A"], row["phi_E"]
            dom = ("allocation" if phi_a > phi_e else "envelope") \
                if np.isfinite(phi_a) and np.isfinite(phi_e) else "undetermined"
            mag = self.magnitude.loc[d] if d in self.magnitude.index else {}
            out[d] = {"dominant_axis": dom,
                      "phi_allocation": _f(phi_a), "phi_envelope": _f(phi_e),
                      "S_interaction": _f(row["S_AE"]),
                      "phi_kinetic": _f(row.get("phi_kin")),
                      "phi_stability": _f(row.get("phi_stab")),
                      "nominal": _f(self.nominal_desc.get(d)),
                      "IQR_allocation": _f(mag.get("IQR_allocation")),
                      "IQR_kinetic": _f(mag.get("IQR_kinetic")),
                      "IQR_stability": _f(mag.get("IQR_stability")),
                      "IQR_envelope": _f(mag.get("IQR_envelope"))}
        return out

    def save(self, out_dir: str, no_plots: bool = False) -> List[str]:
        os.makedirs(out_dir, exist_ok=True)
        written = []
        p = os.path.join(out_dir, "decomposition_recast_table.csv")
        self.table.to_csv(p, index_label="descriptor"); written.append(p)
        pm_ = os.path.join(out_dir, "decomposition_recast_magnitude.csv")
        self.magnitude.to_csv(pm_, index_label="descriptor"); written.append(pm_)
        np.save(os.path.join(out_dir, "recast_temps_C.npy"), self.temps_C)
        for ax, C in self.curves.items():
            np.save(os.path.join(out_dir, f"recast_curves_{ax}.npy"), C)
        np.save(os.path.join(out_dir, "recast_nominal_curve.npy"), self.nominal_curve)
        summ = os.path.join(out_dir, "recast_summary.json")
        with open(summ, "w") as fh:
            json.dump({"H": self.H, "reference_scales": self.refs,
                       "descriptors": self.summary()}, fh, indent=2)
        written.append(summ)
        if not no_plots:
            try:
                written.append(plot_iqr_bands(self, out_dir))
                written.append(plot_recast_variance(self, out_dir))
            except Exception as e:
                print(f"[recast] plotting skipped ({e})")
        return written


def _descs_from_curves(curves: np.ndarray, temps_C, crit_frac):
    """Descriptor DataFrame for a (n, n_temps) growth-curve ensemble."""
    from .tpc import TPC
    rows = []
    for row in curves:
        rows.append(TPC(np.asarray(temps_C, float), np.asarray(row, float)).descriptors(crit_frac).as_dict())
    return pd.DataFrame(rows)


def run_decomposition_recast(pm, temps_C, H: float = 0.2,
                             n_alloc: int = 10, n_kin: int = 8, n_stab: int = 6,
                             n_band: int = 41, seed: int = 1, crit_frac: float = 0.05,
                             descriptors: Sequence[str] = DESCRIPTORS,
                             progress: bool = True) -> RecastResult:
    """Range-fair, nominal-centred, magnitude-aware allocation vs (kinetic,
    stability) envelope decomposition. One fully-crossed 3-way ANOVA grid gives the
    variance shares; four marginal ensembles (each axis perturbed alone over n_band
    LHS draws) give the achievable-IQR magnitude and the IQR-band TPC curves."""
    temps_C = np.asarray(temps_C, float)
    default_budget = pm.ec.default_budget
    alloc_r, kin_r, stab_r, refs = _recast_ranges(pm, H)

    A_tab = _lhs_table(alloc_r, n_alloc, seed)
    K_tab = _lhs_table(kin_r, n_kin, seed + 1)
    S_tab = _lhs_table(stab_r, n_stab, seed + 2)
    a_names, k_names, s_names = list(A_tab), list(K_tab), list(S_tab)

    grids = {d: np.full((n_alloc, n_kin, n_stab), np.nan) for d in descriptors}
    total = n_alloc * n_kin * n_stab
    done = 0
    for i in range(n_alloc):
        ai = {k: float(A_tab.iloc[i][k]) for k in a_names}
        for j in range(n_kin):
            kj = {k: float(K_tab.iloc[j][k]) for k in k_names}
            for l in range(n_stab):
                sl = {k: float(S_tab.iloc[l][k]) for k in s_names}
                pert = _make_perturbation({**ai, **kj, **sl}, default_budget, ())
                desc = compute_tpc(pm, temps_C, pert).descriptors(crit_frac).as_dict()
                for d in descriptors:
                    grids[d][i, j, l] = desc.get(d, np.nan)
                done += 1
        if progress and (i + 1) % max(1, n_alloc // 5) == 0:
            print(f"[recast] crossed 3-way grid {done}/{total}")

    table = pd.DataFrame({d: decompose_grid3(grids[d]) for d in descriptors}).T

    # marginal ensembles: perturb ONE axis (its multi-dim LHS), others nominal
    band = {"allocation": alloc_r, "kinetic": kin_r,
            "stability": stab_r, "envelope": {**kin_r, **stab_r}}
    curves, desc_marg = {}, {}
    for ax, ranges in band.items():
        tab = _lhs_table(ranges, n_band, seed + 7)
        nm = list(tab)
        C = np.vstack([
            compute_tpc(pm, temps_C,
                        _make_perturbation({k: float(tab.iloc[r][k]) for k in nm},
                                           default_budget, ())).growth
            for r in range(n_band)])
        curves[ax] = C
        desc_marg[ax] = _descs_from_curves(C, temps_C, crit_frac)
        if progress:
            print(f"[recast] marginal band '{ax}' ({n_band} draws)")

    nominal_tpc = compute_tpc(pm, temps_C, Perturbation())
    nominal_desc = nominal_tpc.descriptors(crit_frac).as_dict()

    # magnitude: achievable IQR & full range per descriptor per axis
    mag_rows = {}
    for d in descriptors:
        r = {}
        for ax in band:
            v = desc_marg[ax][d].replace([np.inf, -np.inf], np.nan).dropna().values \
                if d in desc_marg[ax].columns else np.array([])
            if v.size:
                q1, q3 = np.percentile(v, [25, 75])
                r[f"IQR_{ax}"] = float(q3 - q1)
                r[f"range_{ax}"] = float(np.max(v) - np.min(v))
            else:
                r[f"IQR_{ax}"] = np.nan; r[f"range_{ax}"] = np.nan
        mag_rows[d] = r
    magnitude = pd.DataFrame(mag_rows).T

    return RecastResult(temps_C, H, refs, grids, table, magnitude, nominal_desc,
                        curves, desc_marg, nominal_tpc.growth)


# ---------------------------------------------------------------------------
# recast plots (PART D centrepiece + variance)
# ---------------------------------------------------------------------------
def plot_iqr_bands(res: RecastResult, out_dir, fname="recast_iqr_bands.png"):
    """The centrepiece: per axis, the family of ABSOLUTE TPCs (growth 1/h) from
    perturbing that axis alone over +/-H, drawn as median + shaded IQR band, with
    the nominal overlaid. Allocation -> height; kinetic envelope -> cold side / Ea;
    Tm/stability -> hot side / CTmax."""
    plt = _mpl()
    T = res.temps_C
    panels = [("allocation", "Allocation (f_metab, f_maint)", "tab:orange"),
              ("kinetic", "Kinetic envelope (Topt, spread, curvature)", "tab:green"),
              ("stability", "Tm / stability (unfolding)", "tab:red"),
              ("envelope", "Full envelope (kinetic + stability)", "tab:blue")]
    panels = [p for p in panels if p[0] in res.curves]
    fig, axes = plt.subplots(1, len(panels), figsize=(4.4 * len(panels), 4.4),
                             sharey=True)
    if len(panels) == 1:
        axes = [axes]
    ymax = 0.0
    for ax, (key, title, col) in zip(axes, panels):
        C = res.curves[key]
        med = np.median(C, axis=0)
        q1, q3 = np.percentile(C, [25, 75], axis=0)
        ax.fill_between(T, q1, q3, color=col, alpha=0.30, label="IQR (±H)")
        ax.plot(T, med, color=col, lw=2, label="median")
        ax.plot(T, res.nominal_curve, "k--", lw=1.6, label="nominal")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Temperature (°C)")
        ax.legend(frameon=False, fontsize=8)
        ymax = max(ymax, float(np.nanmax(q3)))
    axes[0].set_ylabel("Growth rate (1/h)")
    for ax in axes:
        ax.set_ylim(0, ymax * 1.05)
    fig.suptitle(f"Achievable TPC bands per axis (nominal-centred, H={res.H}): "
                 "allocation moves height, kinetics the cold side, Tm the hot limit",
                 fontsize=11)
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_recast_variance(res: RecastResult, out_dir, fname="recast_variance.png"):
    """Two-panel: (left) 3-group Shapley (allocation/kinetic/stability) stacked per
    descriptor; (right) the achievable IQR magnitude per axis per descriptor -- so a
    large share on a negligible movement cannot masquerade as importance."""
    plt = _mpl()
    descs = [d for d in KEY_DESCRIPTORS if d in res.table.index]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.4))
    x = np.arange(len(descs))
    colors = ["tab:orange", "tab:green", "tab:red"]
    labels = ["allocation", "kinetic envelope", "Tm / stability"]
    pa = res.table.loc[descs, "phi_A3"].astype(float).values
    pk = res.table.loc[descs, "phi_kin"].astype(float).values
    ps = res.table.loc[descs, "phi_stab"].astype(float).values
    handles = [axL.bar(x, pa, color=colors[0])[0],
               axL.bar(x, pk, bottom=pa, color=colors[1])[0],
               axL.bar(x, ps, bottom=pa + pk, color=colors[2])[0]]
    axL.set_xticks(x); axL.set_xticklabels(descs, rotation=30, ha="right")
    axL.set_ylabel("Shapley variance share (sums to 1)")
    axL.set_ylim(0, 1); axL.set_title("3-group variance share (range-fair)")
    w = 0.27
    for off, ax_key, col in ((-w, "IQR_allocation", colors[0]),
                             (0.0, "IQR_kinetic", colors[1]),
                             (w, "IQR_stability", colors[2])):
        vals = res.magnitude.loc[descs, ax_key].astype(float).values \
            if ax_key in res.magnitude.columns else np.zeros(len(descs))
        axR.bar(x + off, vals, w, color=col)
    axR.set_xticks(x); axR.set_xticklabels(descs, rotation=30, ha="right")
    axR.set_ylabel("achievable IQR of descriptor (native units)")
    axR.set_title("Magnitude: IQR under each axis alone")
    # one shared legend, anchored BELOW both panels so it never overlaps the bars
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Recast decomposition: variance SHARE (left) vs MAGNITUDE (right)")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    return p


def run_recast(cfg: Dict, out_dir: str, no_plots: bool = False):
    """CLI entry for the recast decomposition (config kind: decompose with a
    `recast:` block, or defaults)."""
    from .config import build_provider, temperature_grid
    d = cfg.get("decomposition", {}) or {}
    rc = cfg.get("recast", d.get("recast", {})) or {}
    pm = build_provider(cfg)
    try:
        pm.ec.model.solver.configuration.timeout = int(cfg.get("solver_timeout", 10))
    except Exception:
        pass
    temps = temperature_grid(cfg)
    H = float(rc.get("H", 0.2))
    print(f"[recast] model={pm.name} grid={temps[0]:.0f}-{temps[-1]:.0f}°C H={H} "
          f"nA={rc.get('n_allocation', 10)} nK={rc.get('n_kinetic', 8)} "
          f"nS={rc.get('n_stability', 6)} band={rc.get('n_band', 41)}")
    res = run_decomposition_recast(
        pm, temps, H=H,
        n_alloc=int(rc.get("n_allocation", 10)), n_kin=int(rc.get("n_kinetic", 8)),
        n_stab=int(rc.get("n_stability", 6)), n_band=int(rc.get("n_band", 41)),
        seed=int(rc.get("seed", 1)), crit_frac=cfg.get("crit_frac", 0.05))
    res.save(out_dir, no_plots=no_plots)
    print(f"[recast] reference scales: {res.refs}")
    def _s(x):
        return "nan" if x is None or not np.isfinite(x) else f"{x:.2f}"
    for desc, info in res.summary().items():
        if desc in KEY_DESCRIPTORS:
            print(f"[recast] {desc:14s} dom={info['dominant_axis']:11s} "
                  f"phi_A={_s(info['phi_allocation'])} phi_E={_s(info['phi_envelope'])} "
                  f"(kin={_s(info['phi_kinetic'])} stab={_s(info['phi_stability'])}) "
                  f"| IQR alloc={_s(info['IQR_allocation'])} kin={_s(info['IQR_kinetic'])} "
                  f"stab={_s(info['IQR_stability'])}")
    print(f"[recast] wrote {out_dir}")
    return out_dir


# ---------------------------------------------------------------------------
# envelope-sampling helper (shared with sensitivity)
# ---------------------------------------------------------------------------
def build_envelope_samples(pm, cfg, n, seed):
    """Build per-enzyme thermal samples from an `envelope_sampling` config block,
    or return None (knobs mode / no block -> unchanged behaviour)."""
    es = cfg.get("envelope_sampling")
    if not es or es.get("mode", "knobs") == "knobs":
        return None
    from . import thermal_sampling as ts
    mode = es["mode"]
    posterior_df = None
    key = "rxn_id"
    if mode == "posterior":
        import os
        import pandas as _pd
        strain = cfg.get("_strain")
        src = es.get("posterior_from", "dltkcat")
        path = os.path.join("strains", str(strain), src, "fits.csv")
        if not os.path.exists(path):
            raise SystemExit(f"posterior envelope sampling needs {path} "
                             "(run `etcgem dltkcat parse` first)")
        posterior_df = _pd.read_csv(path)
        key = "rxn_id" if "rxn_id" in posterior_df.columns else "enzyme_id"
    return ts.sample_thermal(
        pm, n, mode=mode, rho=float(es.get("shared_fraction", 0.7)),
        topt_sd_K=float(es.get("topt_sd_K", 4.0)),
        dcp_sd_frac=float(es.get("dcp_sd_frac", 0.3)),
        posterior_df=posterior_df, key=key, seed=seed)


# keep a private alias used above
_build_envelope_samples = build_envelope_samples


# ---------------------------------------------------------------------------
# CLI entry (called by cli.cmd_decompose)
# ---------------------------------------------------------------------------
def run(cfg: Dict, out_dir: str, no_plots: bool = False):
    """Build the provider from a resolved config and run the decomposition."""
    from .config import build_provider, temperature_grid
    if "decomposition" not in cfg:
        raise SystemExit("config has no `decomposition:` block "
                         "(use an experiment of kind: decompose)")
    d = cfg["decomposition"]
    pm = build_provider(cfg)
    try:
        pm.ec.model.solver.configuration.timeout = int(cfg.get("solver_timeout", 10))
    except Exception:
        pass
    temps = temperature_grid(cfg)
    n_alloc, n_env = int(d.get("n_allocation", 24)), int(d.get("n_envelope", 24))
    # Optional per-enzyme envelope sampling (M1.2) replaces the 3-knob envelope axis.
    env_samples = _build_envelope_samples(pm, cfg, n_env, int(d.get("seed", 1)))
    mode = (cfg.get("envelope_sampling") or {}).get("mode", "knobs")
    print(f"[decompose] model={pm.name} grid={temps[0]:.0f}-{temps[-1]:.0f}°C "
          f"M={n_alloc} N={n_env} envelope={mode}")
    res = run_decomposition(
        pm, temps, d["allocation_params"], d["envelope_params"],
        n_alloc=n_alloc, n_env=n_env,
        seed=int(d.get("seed", 1)), crit_frac=cfg.get("crit_frac", 0.05),
        envelope_samples=env_samples)
    res.save(out_dir, no_plots=no_plots)
    for desc, info in res.summary().items():
        if desc in KEY_DESCRIPTORS:
            print(f"[decompose] {desc:14s} dominant={info['dominant_axis']:11s} "
                  f"phi_A={info['phi_allocation']} phi_E={info['phi_envelope']}")
    print(f"[decompose] wrote {out_dir}")
    return out_dir
