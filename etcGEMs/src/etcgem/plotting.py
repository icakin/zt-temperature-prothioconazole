"""Plots for TPC ensembles and sensitivity results. matplotlib only."""
from __future__ import annotations

import os

import numpy as np


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_ensemble(result, out_dir: str, fname: str = "tpc_ensemble.png"):
    """Spaghetti of all sampled TPCs, median +/- IQR band, nominal curve."""
    plt = _mpl()
    T, C = result.temps_C, result.curves
    fig, ax = plt.subplots(figsize=(7, 5))
    for row in C:
        ax.plot(T, row, color="0.7", lw=0.4, alpha=0.4)
    med = np.median(C, axis=0)
    q1, q3 = np.percentile(C, [25, 75], axis=0)
    ax.fill_between(T, q1, q3, color="tab:blue", alpha=0.25, label="IQR")
    ax.plot(T, med, color="tab:blue", lw=2, label="median")
    ax.plot(result.nominal.temps_C, result.nominal.growth, "k--", lw=2, label="nominal")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Growth rate (1/h)")
    ax.set_title("TPC ensemble across proteome allocation & kcat(T)")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def plot_descriptor_hist(result, out_dir: str, fname: str = "descriptor_distributions.png"):
    plt = _mpl()
    cols = ["Topt_C", "rmax", "CTmax_C", "niche_width_C", "B80_C", "Ea_eV"]
    cols = [c for c in cols if c in result.descriptors.columns]
    n = len(cols)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, c in zip(axes, cols):
        vals = result.descriptors[c].replace([np.inf, -np.inf], np.nan).dropna()
        ax.hist(vals, bins=25, color="tab:blue", alpha=0.8)
        ax.set_title(c)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("TPC descriptor distributions")
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def plot_sensitivity(result, out_dir: str, fname: str = "sensitivity_heatmap.png"):
    plt = _mpl()
    S = result.sensitivity.astype(float)
    fig, ax = plt.subplots(figsize=(1.2 * S.shape[1] + 2, 0.7 * S.shape[0] + 2))
    im = ax.imshow(S.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(S.shape[1]))
    ax.set_xticklabels(S.columns, rotation=45, ha="right")
    ax.set_yticks(range(S.shape[0]))
    ax.set_yticklabels(S.index)
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            v = S.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if abs(v) > 0.5 else "black", fontsize=8)
    ax.set_title("Spearman sensitivity (input × TPC descriptor)")
    fig.colorbar(im, ax=ax, label="rho")
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def plot_descriptor_intervals(result, out_dir: str,
                              fname: str = "descriptor_intervals.png"):
    """Headline calibrated-uncertainty view: per-descriptor median with the
    2.5-97.5 percentile interval. One panel per descriptor (independent scales)
    drawn as a horizontal point-interval. Works for any sweep."""
    plt = _mpl()
    cols = ["Topt_C", "rmax", "CTmax_C", "niche_width_C", "Ea_eV"]
    cols = [c for c in cols if c in result.descriptors.columns]
    if not cols:
        return None
    n = len(cols)
    fig, axes = plt.subplots(n, 1, figsize=(6, 0.9 * n + 1.2))
    axes = np.atleast_1d(axes).ravel()
    for ax, c in zip(axes, cols):
        vals = result.descriptors[c].replace([np.inf, -np.inf], np.nan).dropna()
        lo, med, hi = np.percentile(vals, [2.5, 50, 97.5])
        ax.plot([lo, hi], [0, 0], color="tab:blue", lw=3, alpha=0.4,
                solid_capstyle="round")
        ax.plot([med], [0], "o", color="tab:blue", ms=8, zorder=3)
        ax.annotate(f"{med:.3g}  [{lo:.3g}, {hi:.3g}]",
                    xy=(med, 0), xytext=(0, 8), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8)
        ax.set_yticks([])
        ax.set_ylabel(c, rotation=0, ha="right", va="center", fontsize=10)
        ax.margins(x=0.15)
        ax.margins(y=0.6)
    fig.suptitle("Calibrated descriptor intervals (median, 2.5–97.5%)")
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def plot_sector_tradeoff(result, out_dir: str, fname: str = "sector_tradeoff.png"):
    """Interior growth optimum from metabolic<->biosynthesis co-limitation.
    Only meaningful when the proteome-sector split was swept: returns None (and
    writes nothing) unless ``f_metab`` is a column of result.samples."""
    if "f_metab" not in result.samples.columns:
        return None
    plt = _mpl()
    fm = result.samples["f_metab"].values
    fmaint = (result.samples["f_maint"].values
              if "f_maint" in result.samples.columns
              else np.zeros_like(fm))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, c in zip(axes, ["rmax", "CTmax_C"]):
        if c not in result.descriptors.columns:
            ax.axis("off")
            continue
        y = result.descriptors[c].replace([np.inf, -np.inf], np.nan).values
        sc = ax.scatter(fm, y, c=fmaint, cmap="viridis", s=18, alpha=0.85)
        ax.set_xlabel("f_metab (metabolic sector)")
        ax.set_ylabel(c)
        ax.set_title(f"{c} vs metabolic allocation")
        fig.colorbar(sc, ax=ax, label="f_maint")
    fig.suptitle("Proteome-sector trade-off (interior optimum)")
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def plot_elasticity_heatmap(result, out_dir: str, fname: str = "elasticity_heatmap.png"):
    """Heatmap of standardised elasticities E[D,p] (inputs × descriptors), signed."""
    plt = _mpl()
    E = result.elasticity.astype(float)
    vmax = float(np.nanmax(np.abs(E.values))) or 1.0
    fig, ax = plt.subplots(figsize=(1.2 * E.shape[1] + 3, 0.7 * E.shape[0] + 2))
    im = ax.imshow(E.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(E.shape[1])); ax.set_xticklabels(E.columns, rotation=45, ha="right")
    ax.set_yticks(range(E.shape[0])); ax.set_yticklabels(E.index)
    for i in range(E.shape[0]):
        for j in range(E.shape[1]):
            v = E.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if abs(v) > 0.5 * vmax else "black", fontsize=8)
    ax.set_title(f"Standardised elasticity E[D,p]  (equal step h={result.reference_scales['h']})")
    fig.colorbar(im, ax=ax, label="elasticity")
    fig.tight_layout()
    p = os.path.join(out_dir, fname); fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_elasticity_tornado(result, out_dir: str,
                            descriptors=("rmax", "Ea_eV", "CTmax_C", "Topt_C", "niche_width_C"),
                            fname: str = "elasticity_tornado.png"):
    """Per-descriptor tornado bars: inputs ranked by |elasticity|, signed colour."""
    plt = _mpl()
    descs = [d for d in descriptors if d in result.elasticity.columns]
    n = len(descs)
    fig, axes = plt.subplots(1, n, figsize=(3.1 * n, 4.4))
    axes = np.atleast_1d(axes).ravel()
    for ax, D in zip(axes, descs):
        s = result.elasticity[D].astype(float).dropna()
        s = s.reindex(s.abs().sort_values().index)
        colors = ["tab:red" if v > 0 else "tab:blue" for v in s.values]
        ax.barh(range(len(s)), s.values, color=colors)
        ax.set_yticks(range(len(s))); ax.set_yticklabels(s.index, fontsize=8)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set_title(D, fontsize=10); ax.set_xlabel("elasticity")
    fig.suptitle("Which inputs drive each TPC descriptor (equal-perturbation elasticity; red +, blue −)")
    fig.tight_layout()
    p = os.path.join(out_dir, fname); fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_all(result, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    figs = [
        plot_ensemble(result, out_dir),
        plot_descriptor_hist(result, out_dir),
        plot_sensitivity(result, out_dir),
        plot_descriptor_intervals(result, out_dir),
    ]
    samples = getattr(result, "samples", None)
    if samples is not None and "f_metab" in samples.columns:
        figs.append(plot_sector_tradeoff(result, out_dir))
    return figs


# ===========================================================================
# Model anatomy: reference TPC, enzyme-parameter densities, example kcat(T)
# ===========================================================================
def plot_reference_tpc(pm, temps_C, out_dir: str, crit_frac: float = 0.05,
                       fname: str = "reference_tpc.png", op_label=None):
    """The model's growth TPC at the reference operating point, on RAW absolute
    rate (1/h), with the descriptors marked. ``op_label`` overrides the
    operating-point annotation (e.g. for a rich BHI reference); when None it falls
    back to the static glucose-minimal sector nominal."""
    plt = _mpl()
    from .tpc import compute_tpc
    from .enzyme_cost import Perturbation
    T = np.asarray(temps_C, float)
    tpc = compute_tpc(pm, T, Perturbation())
    g = tpc.growth
    d = tpc.descriptors(crit_frac)
    if op_label is None:
        sec = getattr(pm.ec, "_sectors", None) or {}
        fm = float(sec.get("f_metab_nom", float("nan")))
        fmn = float(sec.get("f_maint_nom", float("nan")))
        fb = 1.0 - fm - fmn if np.isfinite(fm) and np.isfinite(fmn) else float("nan")
        op_label = (f"Reference operating point: glucose-minimal\n"
                    f"$f_\\mathrm{{metab}}$={fm:.3f}, $f_\\mathrm{{bio}}$={fb:.3f}, "
                    f"$f_\\mathrm{{maint}}$={fmn:.3f}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(T, g, color="tab:blue", lw=2.2, zorder=3)
    # r_max / T_opt
    ax.axhline(d.rmax, color="0.7", ls=":", lw=1)
    ax.plot([d.Topt_C], [d.rmax], "o", color="tab:blue", ms=8, zorder=4)
    ax.annotate(f"$T_\\mathrm{{opt}}$ = {d.Topt_C:.0f} °C\n$r_\\mathrm{{max}}$ = {d.rmax:.2f} h⁻¹",
                xy=(d.Topt_C, d.rmax), xytext=(d.Topt_C - 15, d.rmax * 0.86),
                fontsize=9, ha="left")
    # CTmin / CTmax at crit_frac of rmax
    lvl = crit_frac * d.rmax
    for xc, lab, dx in ((d.CTmin_C, "$CT_\\mathrm{min}$", -1), (d.CTmax_C, "$CT_\\mathrm{max}$", 1)):
        if np.isfinite(xc):
            ax.axvline(xc, color="0.6", ls="--", lw=1)
            ax.annotate(f"{lab}\n{xc:.1f} °C", xy=(xc, lvl),
                        xytext=(xc + dx * 3.0, d.rmax * 0.12),
                        fontsize=9, ha="right" if dx < 0 else "left", color="0.3")
    # Ea marker on the rising limb
    i_peak = int(np.argmax(g))
    mask = (np.arange(len(T)) <= i_peak) & (g > 0.1 * d.rmax) & (g < 0.95 * d.rmax)
    if mask.sum() >= 3:
        Tm_mid = float(np.median(T[mask]))
        gm_mid = float(np.interp(Tm_mid, T, g))
        ax.annotate(f"rising limb\n$E_a$ = {d.Ea_eV:.2f} eV",
                    xy=(Tm_mid, gm_mid), xytext=(Tm_mid - 2, d.rmax * 0.55),
                    fontsize=9, ha="right",
                    arrowprops=dict(arrowstyle="->", color="0.5", lw=1))
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Growth rate (1/h)")
    ax.set_ylim(0, d.rmax * 1.18)
    ax.text(0.02, 0.97, op_label, transform=ax.transAxes, va="top", ha="left", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="0.96", ec="0.8"))
    ax.set_title("Reference growth thermal performance curve (complete model)")
    fig.tight_layout()
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_enzyme_param_densities(ec, out_dir: str,
                                fname: str = "enzyme_param_densities.png"):
    """Distributions of the grounded per-enzyme thermal parameters across all
    enzymes in the loaded model: optimum Topt (°C), melting temperature Tm (°C),
    and MMRT curvature dCp (kJ/mol/K). The Topt/Tm SDs are the reference scales the
    elasticity/decomposition use to set their standardised additive steps."""
    plt = _mpl()
    topt = np.asarray(ec._Topt, float) - 273.15
    tm = np.asarray(ec._Tm, float) - 273.15
    dcp = np.asarray(ec._uCpt, float) / 1000.0     # J/mol/K -> kJ/mol/K
    n = len(topt)
    panels = [(topt, "Enzyme optimum $T_\\mathrm{opt}$ (°C)", "tab:green"),
              (tm, "Melting temperature $T_m$ (°C)", "tab:red"),
              (dcp, "MMRT curvature $\\Delta C_p$ (kJ/mol/K)", "tab:purple")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, (v, title, col) in zip(axes, panels):
        vv = v[np.isfinite(v)]
        ax.hist(vv, bins=40, color=col, alpha=0.75)
        mu, sd = float(np.mean(vv)), float(np.std(vv))
        ax.axvline(mu, color="k", ls="--", lw=1.2)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("enzymes")
        ax.text(0.97, 0.95, f"mean {mu:.2f}\nSD {sd:.2f}", transform=ax.transAxes,
                va="top", ha="right", fontsize=9,
                bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.9))
    fig.suptitle(f"Per-enzyme grounded thermal parameters across {n} enzymes "
                 "(held fixed as the model baseline)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150); plt.close(fig)
    return p


def plot_example_kcatT(ec, temps_C, out_dir: str, n_examples: int = 10,
                       highlight=None, fname: str = "example_enzyme_kcatT.png"):
    """Per-enzyme temperature responses for a spread of representative enzymes:
    the relative turnover $\\widehat k(T)=k_\\mathrm{cat}(T)/k_\\mathrm{cat}(T_\\mathrm{opt})$
    (left, peaks at each enzyme's Topt) and the native (folded) fraction $f_N(T)$
    (right, collapses at each enzyme's Tm). ``highlight`` is an optional list of
    (enzyme_id, label) to draw thicker/labelled (e.g. top thermal-control enzymes)."""
    plt = _mpl()
    from . import unfolding as U
    T = np.asarray(temps_C, float)
    Tk = T + 273.15
    topt = np.asarray(ec._Topt, float)
    tm = np.asarray(ec._Tm, float)
    hth, sts, cpu, cpt = ec._uHTH, ec._uSTS, ec._uCpu, ec._uCpt
    finite = np.where(np.isfinite(topt) & np.isfinite(tm))[0]
    # exclude extreme outliers (they distort the panel and the colour gradient):
    # keep the central 2nd-98th percentile of Topt before spreading examples across it
    lo, hi = np.percentile(topt[finite], [2, 98])
    core = finite[(topt[finite] >= lo) & (topt[finite] <= hi)]
    order = core[np.argsort(topt[core])]
    # ~n_examples enzymes evenly spaced across the (clipped) Topt range
    qs = np.linspace(0, len(order) - 1, n_examples).round().astype(int)
    idx = list(dict.fromkeys(order[qs].tolist()))
    # add highlighted control enzymes if present
    hi_idx = {}
    if highlight:
        ents = ec.table.entries
        id_to_i = {getattr(e, "enzyme_id", None): i for i, e in enumerate(ents)}
        for eid, lab in highlight:
            i = id_to_i.get(eid)
            if i is not None:
                hi_idx[i] = lab
                if i not in idx:
                    idx.append(i)

    cmap = plt.get_cmap("viridis")
    trange = (topt[idx].min(), topt[idx].max())
    fig, (axk, axn) = plt.subplots(1, 2, figsize=(13, 5))
    for i in idx:
        rk = U.rel_kcat(Tk, hth[i], sts[i], cpu[i], cpt[i], topt[i])
        fn = U.native_fraction(Tk, hth[i], sts[i], cpu[i])
        c = cmap((topt[i] - trange[0]) / (trange[1] - trange[0] + 1e-9))
        hl = i in hi_idx
        lw = 2.6 if hl else 1.3
        alpha = 1.0 if hl else 0.8
        lab = None
        if hl:
            lab = f"{hi_idx[i]} ($T_\\mathrm{{opt}}${topt[i]-273.15:.0f}, $T_m${tm[i]-273.15:.0f} °C)"
        axk.plot(T, rk, color=c, lw=lw, alpha=alpha, label=lab)
        axn.plot(T, fn, color=c, lw=lw, alpha=alpha, label=lab)
    axk.set_xlabel("Temperature (°C)"); axk.set_ylabel("relative turnover $\\widehat k(T)$")
    axk.set_title("Per-enzyme turnover (anchored at each $T_\\mathrm{opt}$)")
    axk.set_ylim(0, 1.05)
    axn.set_xlabel("Temperature (°C)"); axn.set_ylabel("native fraction $f_N(T)$")
    axn.set_title("Per-enzyme folded fraction (collapses at each $T_m$)")
    axn.set_ylim(0, 1.05)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(trange[0] - 273.15, trange[1] - 273.15))
    cb = fig.colorbar(sm, ax=[axk, axn], fraction=0.03, pad=0.02)
    cb.set_label("enzyme $T_\\mathrm{opt}$ (°C)")
    if any(i in hi_idx for i in idx):
        axk.legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle(f"Example per-enzyme temperature responses ({len(idx)} enzymes "
                 "spanning the optimum range)", fontsize=11)
    p = os.path.join(out_dir, fname)
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    return p
