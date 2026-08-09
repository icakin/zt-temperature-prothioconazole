"""Unified command-line interface for etcgem.

Two tiers of commands over the defaults/strain/experiment config model
(see etcgem.config):

Strain-only (no experiment needed -- a strain is runnable on its own):
    etcgem build --strain NAME                 # build provider; print+save summary
    etcgem tpc   --strain NAME [--fits [PATH]]  # nominal TPC + descriptors + plot
    etcgem fba   --strain NAME --temp C [--fits [PATH]]   # single solve at C
    etcgem dltkcat prep|parse --strain NAME ... # DLTKcat tooling (strain-aware)

Strain + experiment (method overlay from configs/experiments/EXP.yaml):
    etcgem sweep     --strain NAME --experiment EXP [--fits [PATH]] [--resume] [--seconds N] [--no-plots]
    etcgem decompose --strain NAME --experiment EXP [--no-plots]
    etcgem sweep     --config PATH ...          # ad-hoc self-contained config escape hatch

Every run writes into strains/NAME/outputs/<tag>/ and dumps the exact merged
config there as resolved_config.yaml. No scientific/numerical code changes --
provider construction (config.build_provider) is byte-for-byte unchanged.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

from . import config
from .config import (build_provider, temperature_grid, resolve, dump_resolved,
                     load_config, strain_dir, calibrate_dCp_to_Ea)
from .enzyme_cost import Perturbation
from .sensitivity import (SensitivityResult, run_sensitivity, _lhs,
                          _make_perturbation, _spearman_matrix)
from .tpc import TPC, compute_tpc

_FITS_DEFAULT = "\0default"   # `--fits` given without an explicit path


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------
def _out_dir(strain, tag):
    return os.path.join(strain_dir(strain), "outputs", tag)


def _run_tag(command, experiment):
    """Output-dir tag ``<command>_<experiment>``, de-doubled when the experiment name
    already repeats the command (e.g. decompose + 'decomposition_sectors' ->
    decompose_sectors, not decompose_decomposition_sectors). Affects NEW runs only."""
    if not experiment:
        return command
    fam = {"decompose": ("decompose", "decomposition"), "control": ("control",),
           "elasticity": ("elasticity",), "sweep": ("sweep",)}.get(command, (command,))
    exp = experiment
    for f in fam:
        if exp == f:
            return command
        if exp.startswith(f + "_"):
            exp = exp[len(f) + 1:]
            break
    return f"{command}_{exp}"


def _fits_path(strain, fits_arg):
    """Resolve the --fits value: None -> no fits; sentinel -> strain default."""
    if fits_arg is None:
        return None
    if fits_arg == _FITS_DEFAULT:
        return os.path.join(strain_dir(strain), "dltkcat", "fits.csv")
    return fits_arg


def _build_pm(cfg, fits_path=None, key="rxn_id"):
    """Build the provider, cap the LP solver, and optionally apply DLTKcat fits
    (Topt/dCp only; kcat/base costs untouched -- the old run_with_fits behaviour)."""
    pm = build_provider(cfg)
    try:
        pm.ec.model.solver.configuration.timeout = int(cfg.get("solver_timeout", 10))
    except Exception:
        pass
    if fits_path:
        from .dltkcat import apply_fits_to_provider
        apply_fits_to_provider(pm, pd.read_csv(fits_path), key=key)
    return pm


def _plot_curve(temps_C, growth, out_dir, title, fname="nominal_tpc.png"):
    try:
        from .plotting import _mpl
        plt = _mpl()
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(temps_C, growth, "k-", lw=2)
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Growth rate (1/h)")
        ax.set_title(title)
        fig.tight_layout()
        p = os.path.join(out_dir, fname)
        fig.savefig(p, dpi=150)
        plt.close(fig)
        return p
    except Exception as e:
        print(f"[tpc] plotting skipped ({e})")
        return None


# ---------------------------------------------------------------------------
# strain-only commands
# ---------------------------------------------------------------------------
def cmd_build(args):
    cfg = resolve(args.strain)
    out_dir = _out_dir(args.strain, "build")
    os.makedirs(out_dir, exist_ok=True)
    pm = _build_pm(cfg)
    temps = temperature_grid(cfg)
    summary = {
        "strain": args.strain,
        "model": pm.name,
        "provider_type": cfg["provider"]["type"],
        "n_enzymes": len(pm.ec.table),
        "default_budget": pm.ec.default_budget,
        "T0_C": cfg.get("T0_C", 30.0),
        "temperature_grid_C": [float(temps[0]), float(temps[-1]), int(len(temps))],
    }
    with open(os.path.join(out_dir, "model_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    dump_resolved(cfg, out_dir)
    print(f"[build] strain={args.strain} model={pm.name}")
    print(f"[build] enzymes={len(pm.ec.table)} budget={pm.ec.default_budget:.4g} "
          f"T0={summary['T0_C']}C grid={temps[0]:.0f}-{temps[-1]:.0f}C n={len(temps)}")
    print(f"[build] wrote {out_dir}/model_summary.json")
    return out_dir


def cmd_tpc(args):
    cfg = resolve(args.strain)
    out_dir = _out_dir(args.strain, "tpc")
    os.makedirs(out_dir, exist_ok=True)
    fits = _fits_path(args.strain, args.fits)
    pm = _build_pm(cfg, fits)
    temps = temperature_grid(cfg)
    tpc = compute_tpc(pm, temps, Perturbation())
    desc = tpc.descriptors(cfg.get("crit_frac", 0.05))
    pd.DataFrame({"temp_C": temps, "growth": tpc.growth}).to_csv(
        os.path.join(out_dir, "nominal_tpc.csv"), index=False)
    with open(os.path.join(out_dir, "descriptors.json"), "w") as fh:
        json.dump(desc.as_dict(), fh, indent=2)
    dump_resolved(cfg, out_dir)
    print(f"[tpc] strain={args.strain}" + (f" +fits {os.path.basename(fits)}" if fits else ""))
    print(f"[tpc] Topt={desc.Topt_C:.1f}°C rmax={desc.rmax:.3f}/h "
          f"CTmax={desc.CTmax_C:.1f}°C niche={desc.niche_width_C:.1f}°C Ea={desc.Ea_eV:.2f}eV")
    if not args.no_plots:
        _plot_curve(temps, tpc.growth, out_dir, f"Nominal TPC — {args.strain}")
    print(f"[tpc] wrote {out_dir}")
    return out_dir


def cmd_fba(args):
    cfg = resolve(args.strain)
    out_dir = _out_dir(args.strain, "fba")
    os.makedirs(out_dir, exist_ok=True)
    fits = _fits_path(args.strain, args.fits)
    pm = _build_pm(cfg, fits)
    tpc = compute_tpc(pm, [args.temp], Perturbation())
    growth = float(tpc.growth[0])
    result = {"strain": args.strain, "temp_C": args.temp, "growth": growth,
              "fits": os.path.basename(fits) if fits else None}
    with open(os.path.join(out_dir, "fba_result.json"), "w") as fh:
        json.dump(result, fh, indent=2)
    dump_resolved(cfg, out_dir)
    print(f"[fba] strain={args.strain} T={args.temp}°C -> growth={growth:.5f} /h")
    print(f"[fba] wrote {out_dir}/fba_result.json")
    return out_dir


# ---------------------------------------------------------------------------
# calibrate-dcp (strain-only): choose provider.default_dCp for a target Ea
# ---------------------------------------------------------------------------
def _set_yaml_scalar(text, key, value):
    """Replace the first `  key: <number>` line's value, preserving indentation
    and any trailing comment. Returns (new_text, n_replaced)."""
    import re
    pat = re.compile(rf"(?m)^(\s*{re.escape(key)}:\s*)(-?\d+(?:\.\d+)?)(.*)$")
    n = [0]
    def _sub(m):
        n[0] += 1
        return f"{m.group(1)}{value}{m.group(3)}"
    return pat.sub(_sub, text, count=1), n[0]


def cmd_calibrate_dcp(args):
    import copy
    cfg = resolve(args.strain)
    if cfg["provider"].get("type") == "toy":
        print("[calibrate] note: toy provider uses dCp_mean, not default_dCp; "
              "calibration/write-back will not affect the toy TPC")
    print(f"[calibrate] strain={args.strain} target Ea={args.target_ea} eV "
          f"(configurable; set to your measured bacterial growth-TPC Ea)")

    g = cfg["temperature_grid"]
    start_C, stop_C, n_grid = float(g["start_C"]), float(g["stop_C"]), int(g["n"])
    spacing = (stop_C - start_C) / (n_grid - 1)

    # A shallow dCp broadens the TPC so CT_max can run well past the configured
    # grid stop. Calibrate on an EXTENDED (moderately coarse) grid so Ea and CT_max
    # are resolved throughout the bisection instead of clipping at the grid edge.
    cal_cfg = copy.deepcopy(cfg)
    cal_stop = max(stop_C, 100.0)
    cal_cfg["temperature_grid"] = {"start_C": start_C, "stop_C": cal_stop,
                                   "n": int(round((cal_stop - start_C) / 1.5)) + 1}
    res = calibrate_dCp_to_Ea(cal_cfg, args.target_ea, lo=args.lo, hi=args.hi, tol=args.tol)
    print(f"[calibrate] calibrated default_dCp={res['dCp']:.3f}  "
          f"achieved Ea={res['Ea']:.3f} eV  resolved CTmax={res['CTmax']:.1f}°C"
          + ("" if res["bracketed"] else "  (target not bracketed; clamped)"))

    yaml_path = os.path.join(strain_dir(args.strain), "strain.yaml")
    with open(yaml_path) as fh:
        text = fh.read()
    new_dcp = round(res["dCp"], 3)
    text, n = _set_yaml_scalar(text, "default_dCp", new_dcp)
    if n != 1:
        raise SystemExit(f"[calibrate] could not update default_dCp in {yaml_path} "
                         f"(matched {n} lines); aborting write")

    # Extend the strain's own grid if the resolved CT_max is near/past its stop.
    grid_changed = False
    if res["CTmax"] > stop_C - 2.0:
        new_stop = float(np.ceil(res["CTmax"] + 5.0))
        new_n = int(round((new_stop - start_C) / spacing)) + 1     # keep spacing
        text, ns = _set_yaml_scalar(text, "stop_C", int(new_stop))
        text, nn = _set_yaml_scalar(text, "n", int(new_n))
        if ns == 1 and nn == 1:
            grid_changed = True
            print(f"[calibrate] CTmax {res['CTmax']:.1f}°C exceeds grid stop {stop_C:.0f}-2°C -> "
                  f"raised temperature_grid to stop_C={int(new_stop)}, n={int(new_n)}")
        else:
            print(f"[calibrate] WARNING: wanted to extend grid but matched "
                  f"stop_C={ns}, n={nn} lines; leaving grid unchanged")

    with open(yaml_path, "w") as fh:
        fh.write(text)
    print(f"[calibrate] wrote provider.default_dCp={new_dcp} to {yaml_path}")

    if grid_changed:
        cfg2 = resolve(args.strain)
        from .config import _nominal_ea_ctmax
        ea2, ctmax2 = _nominal_ea_ctmax(cfg2, new_dcp)
        print(f"[calibrate] on the new strain grid: Ea={ea2:.3f} eV  CTmax={ctmax2:.1f}°C")
    return yaml_path


# ---------------------------------------------------------------------------
# sweep (strain + experiment, or --config escape hatch)
# ---------------------------------------------------------------------------
def _sweep_cfg_out(args):
    """Return (cfg, out_dir) for a sweep, from either --config or strain+experiment."""
    if args.config:
        cfg = load_config(args.config)
        cfg.setdefault("output_dir", "outputs/run")
        return cfg, cfg["output_dir"]
    if not (args.strain and args.experiment):
        raise SystemExit("sweep needs --strain NAME --experiment EXP (or --config PATH)")
    cfg = resolve(args.strain, args.experiment)
    if "sensitivity" not in cfg:
        raise SystemExit(f"experiment '{args.experiment}' defines no `sensitivity` block")
    return cfg, _out_dir(args.strain, _run_tag("sweep", args.experiment))


def _finalize_sweep(cfg, out_dir, pm, result, no_plots, tag="run"):
    result.save(out_dir)
    nom = result.nominal.descriptors()
    summary = {
        "model": pm.name,
        "n_enzymes": len(pm.ec.table),
        "default_budget": pm.ec.default_budget,
        "nominal": nom.as_dict(),
        "descriptor_medians": result.descriptors.median(numeric_only=True).to_dict(),
        "descriptor_iqr": (result.descriptors.quantile(0.75)
                           - result.descriptors.quantile(0.25)).to_dict(),
    }
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    dump_resolved(cfg, out_dir)
    print(f"[{tag}] nominal TPC: Topt={nom.Topt_C:.1f}°C rmax={nom.rmax:.3f}/h "
          f"CTmax={nom.CTmax_C:.1f}°C Ea={nom.Ea_eV:.2f}eV -> {out_dir}")
    if not no_plots:
        try:
            from .plotting import plot_all
            plot_all(result, out_dir)
            print(f"[{tag}] figures written")
        except Exception as e:
            print(f"[{tag}] plotting skipped ({e})")


def _run_oneshot(cfg, out_dir, fits_path, key, no_plots):
    print(f"[run] building provider: {cfg['provider']['type']}"
          + (f"  (+fits {os.path.basename(fits_path)})" if fits_path else ""))
    pm = _build_pm(cfg, fits_path, key)
    temps = temperature_grid(cfg)
    print(f"[run] model={pm.name}  enzymes={len(pm.ec.table)}  "
          f"budget={pm.ec.default_budget:.4g}  T grid={temps[0]:.0f}-{temps[-1]:.0f}°C")
    s = cfg["sensitivity"]
    from .decomposition import build_envelope_samples
    env_samples = build_envelope_samples(pm, cfg, s.get("n_samples", 200), s.get("seed", 1))
    if env_samples is not None:
        print(f"[run] envelope sampling: {cfg['envelope_sampling'].get('mode')} "
              f"(rho={cfg['envelope_sampling'].get('shared_fraction', 0.7)})")
    result = run_sensitivity(
        pm, temps, {k: tuple(v) for k, v in s["parameters"].items()},
        n_samples=s.get("n_samples", 200), seed=s.get("seed", 1),
        group_names=s.get("groups", []), crit_frac=cfg.get("crit_frac", 0.05),
        envelope_samples=env_samples,
    )
    _finalize_sweep(cfg, out_dir, pm, result, no_plots, "run")
    return out_dir


def _run_resume(cfg, out_dir, seconds, fits_path, key, no_plots):
    ckpt = os.path.join(out_dir, "_checkpoint.npz")
    s = cfg["sensitivity"]
    ranges = {k: tuple(v) for k, v in s["parameters"].items()}
    n, seed = int(s.get("n_samples", 200)), int(s.get("seed", 1))
    groups = s.get("groups", [])
    names = list(ranges.keys())
    temps = temperature_grid(cfg)
    U = _lhs(n, len(names), seed)
    sampled = {nm: ranges[nm][0] + U[:, j] * (ranges[nm][1] - ranges[nm][0])
               for j, nm in enumerate(names)}
    samples_df = pd.DataFrame(sampled)

    if os.path.exists(ckpt):
        z = np.load(ckpt)
        curves, done = z["curves"], z["done"]
        if curves.shape != (n, len(temps)):
            curves = np.full((n, len(temps)), np.nan); done = np.zeros(n, bool)
    else:
        curves = np.full((n, len(temps)), np.nan); done = np.zeros(n, bool)

    todo = np.where(~done)[0]
    if len(todo):
        print(f"[resume] building provider ({cfg['provider']['type']}) ...")
        pm = _build_pm(cfg, fits_path, key)
        default_budget = pm.ec.default_budget
        t0, n_run = time.time(), 0
        for i in todo:
            row = {k: float(samples_df.iloc[i][k]) for k in names}
            pert = _make_perturbation(row, default_budget, groups)
            curves[i] = compute_tpc(pm, temps, pert).growth
            done[i] = True; n_run += 1
            if time.time() - t0 > seconds:
                break
        np.savez(ckpt, curves=curves, done=done)
        n_done = int(done.sum())
        print(f"[resume] ran {n_run} this call; {n_done}/{n} done")
        if n_done < n:
            print(f"[resume] PARTIAL {n_done}/{n} -- call again")
            return out_dir
    else:
        print(f"[resume] all {n} samples already done -> finalizing")

    desc_df = pd.DataFrame(
        [TPC(temps, curves[i]).descriptors(cfg.get("crit_frac", 0.05)).as_dict()
         for i in range(curves.shape[0])])
    sens = _spearman_matrix(samples_df, desc_df)
    pm = _build_pm(cfg, fits_path, key)
    nominal = compute_tpc(pm, temps, Perturbation())
    res = SensitivityResult(np.asarray(temps), curves, samples_df, desc_df, sens, nominal)
    _finalize_sweep(cfg, out_dir, pm, res, no_plots, "resume")
    if os.path.exists(ckpt):
        os.remove(ckpt)
    print("[resume] ALL DONE")
    return out_dir


def cmd_sweep(args):
    cfg, out_dir = _sweep_cfg_out(args)
    os.makedirs(out_dir, exist_ok=True)
    fits = None if args.config else _fits_path(args.strain, args.fits)
    if args.config and args.fits not in (None, _FITS_DEFAULT):
        fits = args.fits
    if args.resume:
        return _run_resume(cfg, out_dir, args.seconds, fits, args.key, args.no_plots)
    return _run_oneshot(cfg, out_dir, fits, args.key, args.no_plots)


# ---------------------------------------------------------------------------
# decompose (strain + experiment) -- dispatches to a decomposition module
# ---------------------------------------------------------------------------
def cmd_decompose(args):
    if not (args.strain and args.experiment):
        raise SystemExit("decompose needs --strain NAME --experiment EXP")
    cfg = resolve(args.strain, args.experiment)
    out_dir = _out_dir(args.strain, _run_tag("decompose", args.experiment))
    try:
        from . import decomposition  # noqa: F401
    except Exception:
        print("decomposition module not installed yet")
        return sys.exit(1)
    os.makedirs(out_dir, exist_ok=True)
    dump_resolved(cfg, out_dir)
    # recast mode: range-fair, nominal-centred, magnitude-aware, +Tm (H1.3)
    if cfg.get("recast") or (cfg.get("decomposition", {}) or {}).get("mode") == "recast":
        return decomposition.run_recast(cfg, out_dir, no_plots=args.no_plots)
    return decomposition.run(cfg, out_dir, no_plots=args.no_plots)


# ---------------------------------------------------------------------------
# control (strain-level diagnostic; experiment optional for the control block)
# ---------------------------------------------------------------------------
def cmd_control(args):
    if not args.strain:
        raise SystemExit("control needs --strain NAME")
    cfg = resolve(args.strain, args.experiment)
    tag = _run_tag("control", args.experiment)
    out_dir = _out_dir(args.strain, tag)
    os.makedirs(out_dir, exist_ok=True)
    dump_resolved(cfg, out_dir)
    from . import control
    return control.run(cfg, out_dir, no_plots=args.no_plots)


# ---------------------------------------------------------------------------
# elasticity: equal-perturbation (standardised) sensitivity
# ---------------------------------------------------------------------------
def cmd_elasticity(args):
    if not (args.strain and args.experiment):
        raise SystemExit("elasticity needs --strain NAME --experiment EXP")
    cfg = resolve(args.strain, args.experiment)
    out_dir = _out_dir(args.strain, _run_tag("elasticity", args.experiment))
    os.makedirs(out_dir, exist_ok=True)
    pm = _build_pm(cfg)
    temps = temperature_grid(cfg)
    s = cfg.get("sensitivity", {}) or {}
    h = float(args.h if args.h is not None else s.get("h", 0.10))
    inputs = s.get("inputs") or list(s.get("parameters", {}).keys()) or \
        ["dTopt", "topt_scale", "dCp_scale", "budget_scale"]
    from .sensitivity import run_elasticity
    print(f"[elasticity] model={pm.name} h={h} inputs={inputs}")
    res = run_elasticity(pm, temps, inputs, h=h, crit_frac=cfg.get("crit_frac", 0.05),
                         group_names=s.get("groups", []))
    res.save(out_dir)
    dump_resolved(cfg, out_dir)
    if not args.no_plots:
        try:
            from .plotting import plot_elasticity_heatmap, plot_elasticity_tornado
            plot_elasticity_heatmap(res, out_dir)
            plot_elasticity_tornado(res, out_dir)
        except Exception as e:
            print(f"[elasticity] plotting skipped ({e})")
    print(f"[elasticity] h={h}  dTopt reference scale Δref={res.reference_scales['dTopt_reference_scale_K']} K")
    for D in [d for d in ("rmax", "Ea_eV", "CTmax_C", "Topt_C") if d in res.elasticity.columns]:
        rank = res.elasticity[D].astype(float).abs().sort_values(ascending=False)
        top = rank.index[0]
        print(f"[elasticity] {D:12s} top input = {top} (E={res.elasticity.loc[top, D]:+.2f}); "
              f"ranking: " + ", ".join(f"{i}={res.elasticity.loc[i, D]:+.2f}" for i in rank.index))
    print(f"[elasticity] wrote {out_dir}")
    return out_dir


# ---------------------------------------------------------------------------
# dissect: tuned-model analyses (elasticity + decomposition + identifiability)
# ---------------------------------------------------------------------------
def cmd_dissect(args):
    if not args.strain:
        raise SystemExit("dissect needs --strain NAME")
    from . import dissect
    v3 = _out_dir(args.strain, "calibration_vanderlinden_v3")
    if not os.path.exists(os.path.join(v3, "summary.json")):
        raise SystemExit(f"[dissect] no P2 v3 posterior at {v3}; run `calibrate --vdl` first")
    out_base = os.path.dirname(v3)   # strains/<strain>/outputs
    allow = bool(os.environ.get("ALLOW_GLPK"))
    res = dissect.run(args.strain, v3, out_base,
                      n_proc=int(args.procs) if args.procs is not None else 0,
                      n_draws=int(args.draws) if args.draws is not None else 150,
                      allow_glpk=allow)
    ctrl_dir = dissect.run_control_tuned(args.strain, v3,
                                         _out_dir(args.strain, "control_tuned"), allow_glpk=allow)
    bd = res["baseline_descriptors"]
    print("\n" + "=" * 72)
    print("DISSECT (tuned model, rich BHI) — elasticity + decomposition + identifiability")
    print("=" * 72)
    print(f"  tuned baseline: rmax={bd['rmax']:.3f}  Topt={bd['Topt_C']:.1f}  "
          f"CTmax={bd['CTmax_C']:.1f}  Ea={bd['Ea_eV']:.3f}")
    print(f"  elasticity -> {res['elasticity_dir']}")
    print(f"  decompose  -> {res['decompose_dir']}")
    print(f"  control    -> {ctrl_dir}")
    E = res["elasticity"]
    for D in ("rmax", "CTmax_C", "Topt_C"):
        if D in E.columns:
            rank = E[D].astype(float).abs().sort_values(ascending=False)
            print(f"  {D:8s} top levers: " + ", ".join(f"{i}={E.loc[i, D]:+.2f}" for i in rank.index[:4]))
    return out_base


# ---------------------------------------------------------------------------
# anatomy: reference-operating-point description (curve + enzyme distributions)
# ---------------------------------------------------------------------------
def cmd_anatomy(args):
    """Model-anatomy figures at the reference operating point: the reference TPC
    with descriptors, the per-enzyme thermal parameter densities, and an example
    per-enzyme kcat(T)/f_N(T) panel. The reference is the rich BHI medium by
    default (--medium to override, e.g. glucose_minimal)."""
    if not args.strain:
        raise SystemExit("anatomy needs --strain NAME")
    cfg = resolve(args.strain, args.experiment)
    out_dir = _out_dir(args.strain, "anatomy")
    os.makedirs(out_dir, exist_ok=True)
    pm = _build_pm(cfg)
    temps = temperature_grid(cfg)
    dump_resolved(cfg, out_dir)
    # reference operating point medium (rich BHI by default)
    medium = args.medium or "BHI"
    op_label = _anatomy_set_medium(pm, medium, args.strain)
    # optional highlight: top thermal-control enzymes, if a control run exists
    highlight = []
    ctrl = os.path.join(os.path.dirname(out_dir), "control_control", "thermal_control.csv")
    if os.path.exists(ctrl):
        try:
            import pandas as pd
            c = pd.read_csv(ctrl)
            idc = "enzyme_id" if "enzyme_id" in c.columns else None
            if idc:
                for _, r in c.head(2).iterrows():
                    highlight.append((str(r[idc]),
                                      str(r.get("rxn_id", r[idc]))))
        except Exception:
            pass
    from .plotting import (plot_reference_tpc, plot_enzyme_param_densities,
                           plot_example_kcatT)
    print(f"[anatomy] model={pm.name} grid={temps[0]:.0f}-{temps[-1]:.0f}°C "
          f"medium={medium} enzymes={len(pm.ec.table.entries)}")
    p1 = plot_reference_tpc(pm, temps, out_dir, crit_frac=cfg.get("crit_frac", 0.05),
                            op_label=op_label)
    p2 = plot_enzyme_param_densities(pm.ec, out_dir)
    p3 = plot_example_kcatT(pm.ec, temps, out_dir, highlight=highlight or None)
    import numpy as _np
    print(f"[anatomy] Topt SD={_np.std(pm.ec._Topt):.3f}K  Tm SD={_np.std(pm.ec._Tm):.3f}K "
          f"(reference scales for dTopt/dTm)")
    print(f"[anatomy] wrote {p1}\n          {p2}\n          {p3}")
    return out_dir


def _anatomy_set_medium(pm, medium, strain):
    """Set the anatomy reference operating-point medium and return a display label
    with the effective sector fractions at the curve optimum."""
    from .providers import set_medium
    import os as _os
    if medium == "BHI":
        set_medium(pm, "BHI", bhi_media_csv=_os.path.join("strains", strain, "media", "BHI_media.csv"))
    elif medium == "LB":
        set_medium(pm, "LB", lb_media_csv=_os.path.join("strains", strain, "media", "LB_media.csv"))
    else:
        set_medium(pm, "glucose_minimal", "glc__D", True)
    # effective sector fractions at 37 C (T-dependent when allocation_from_data is on)
    alloc = getattr(pm.ec, "_alloc_from_data", None)
    label = {"BHI": "rich (BHI)", "LB": "rich (LB)"}.get(medium, "glucose-minimal")
    if alloc is not None and hasattr(alloc, "model_alloc"):
        fm, fmaint = alloc.model_alloc(37.0)
        fb = 1.0 - fm - fmaint
        return (f"Reference operating point: {label}\n"
                f"$f_\\mathrm{{metab}}$={fm:.3f}, $f_\\mathrm{{bio}}$={fb:.3f}, "
                f"$f_\\mathrm{{maint}}$={fmaint:.3f} (at 37 °C)")
    return f"Reference operating point: {label}"


# ---------------------------------------------------------------------------
# validate: emergent model vs the trusted MG1655 curve (Van Derlinden, BHI)
# ---------------------------------------------------------------------------
def cmd_validate(args):
    if not args.strain:
        raise SystemExit("validate needs --strain NAME")
    from . import validation as val
    out_dir = _out_dir(args.strain, "validation_trusted")
    res = val.run(args.strain, out_dir, include_secondary=not args.no_secondary)
    print("\n" + "=" * 70)
    print("VALIDATION — emergent model vs strain-matched TPCs (nothing fit to growth)")
    print("=" * 70)
    for k, r in res.items():
        if not isinstance(r, dict) or "abs_R2" not in r:
            continue
        print(f"  {k:18s}[{r.get('role','')[:9]:9s}] absR2={r['abs_R2']:>7} "
              f"RMSE={r['RMSE_per_h']} | rmax obs {r['obs_rmax']} vs pred {r['pred_rmax']} | "
              f"Topt obs {r['obs_Topt_C']} vs pred {r['pred_Topt_C']} | "
              f"Ea obs {r['obs_Ea_eV']} vs pred {r['pred_Ea_eV']}")
    print(f"[validate] wrote {out_dir}")
    return out_dir


# ---------------------------------------------------------------------------
# calibrate: single-curve Bayesian calibration (emcee); prior vs posterior
# ---------------------------------------------------------------------------
def cmd_calibrate(args):
    from . import calibration as cal
    if args.list:
        print(cal.list_defined_curves(args.strain).to_string(index=False))
        return
    if getattr(args, "vdl", False):
        return _cmd_calibrate_vdl(args)
    trusted_noll = bool(args.noll)
    if not args.strain or (not args.curve and not trusted_noll):
        raise SystemExit("calibrate needs --strain NAME --curve CURVE_ID (or --noll / --list)")
    cfg = resolve(args.strain, args.experiment) if args.experiment else {}
    cal_cfg = (cfg.get("calibration") or {}) if isinstance(cfg, dict) else {}
    s = cal_cfg.get("sampler", {})
    priors_cfg = cal_cfg.get("priors", {})
    # trusted Noll fit writes to its own dir (does NOT overwrite the superseded phase-1)
    out_dir = _out_dir(args.strain, "calibration_noll_minimal" if trusted_noll
                       else "calibration_phase1")
    res = cal.run_emcee(
        args.strain, args.curve or "Noll2023_NCM3722", out_dir,
        medium=args.medium or cal_cfg.get("medium", "glucose_minimal"),
        n_walkers=int(args.walkers or s.get("n_walkers", 24)),
        n_steps=int(args.steps or s.get("n_steps", 1500)),
        n_burn=int(args.burn or s.get("n_burn", 500)),
        seed=int(args.seed if args.seed is not None else s.get("seed", 1)),
        n_proc=int(args.procs if args.procs is not None else s.get("n_proc", 0)),
        priors_cfg=priors_cfg, trusted_noll=trusted_noll)
    _calibrate_console_summary(res)
    print(f"[calibrate] wrote {out_dir}")
    return out_dir


def _cmd_calibrate_vdl(args):
    """Unified multi-parameter tuning on Van Derlinden at the rich BHI operating
    point, growth law ON, reconciled single pool with the in-vivo saturation sigma
    freed (P2 v3 -> calibration_vanderlinden_v3/; v1=flat-top/kcat-wasted,
    v2=growth-law-on but redundant pool). Gurobi + pre-flight + warm-start +
    autocorr early stop. ALLOW_GLPK=1 in the env overrides the Gurobi stop-and-flag."""
    import os as _os
    from . import calibration_multi as cm
    out_dir = _out_dir(args.strain, "calibration_vanderlinden_v3")
    res = cm.run(args.strain, out_dir,
                 n_walkers=int(args.walkers or 52),
                 n_steps_max=int(args.steps or 6000),
                 n_burn=int(args.burn or 150),
                 seed=int(args.seed if args.seed is not None else 1),
                 n_proc=int(args.procs if args.procs is not None else 0),
                 allow_glpk=bool(_os.environ.get("ALLOW_GLPK")))
    sm, d, po = res["sampler"], res["descriptors"], res["posterior"]
    print("\n" + "=" * 72)
    print("CALIBRATION v3 — Van Derlinden MG1655/BHI, rich operating point (sigma freed)")
    print("=" * 72)
    print(f"  solver={res.get('solver')}  preflight_single_solve={res.get('preflight_single_solve_ms')} ms")
    print(f"  warm_started={sm['warm_started']}  stop_reason={sm['stop_reason']}")
    print(f"  acceptance={sm['acceptance_fraction']}  autocorr_max={sm['autocorr_time_max']}  "
          f"n_eff={sm['n_eff']}  steps={sm['n_steps']}  wall={sm['wall_time_s']}s")
    print(f"  {'descriptor':10s} {'observed':>10s} {'emergent':>10s} {'posterior':>10s}")
    for k in ("rmax", "Topt_C", "Ea_eV", "CTmax_C"):
        print(f"  {k:10s} {d['observed'][k]:>10} {d['emergent_prior'][k]:>10} {d['posterior_median'][k]:>10}")
    print("  demanded corrections (posterior median [90% CI], constrained?):")
    for k in ("kcat_scale", "kappa_scale", "sigma", "dCp_scale", "dTopt", "dTm", "topt_scale",
              "tm_scale", "f_metab", "f_maint", "ngam_scale", "ngam_steepness", "sigma_disc"):
        if k not in po:
            continue
        p = po[k]
        print(f"    {k:12s} {p['demanded_correction']:>18s}  CI{p['posterior_90CI']}  "
              f"constrained={p['constrained_by_curve']}")
    print(f"[calibrate] wrote {out_dir}")
    return out_dir


def _calibrate_console_summary(res):
    d = res["descriptors"]; sm = res["sampler"]; po = res["posterior"]
    print("\n" + "=" * 68)
    print(f"CALIBRATION SUMMARY — {res['curve']['curve_id']} "
          f"({res['curve']['study']}, {res['curve']['medium_class']})")
    print("=" * 68)
    print(f"  acceptance={sm['acceptance_fraction']}  autocorr_max={sm['autocorr_time_max']}  "
          f"n_eff={sm['n_eff']}  wall={sm['wall_time_s']}s")
    print(f"  {'descriptor':10s} {'observed':>10s} {'prior':>10s} {'posterior':>10s}")
    for k in ("rmax", "Topt_C", "Ea_eV", "CTmax_C"):
        print(f"  {k:10s} {d['observed'][k]:>10} {d['prior'][k]:>10} {d['posterior_median'][k]:>10}")
    print("  demanded corrections (posterior median [90% CI], curve-constrained?):")
    for k in ("kappa_scale", "dCp_scale", "dTopt", "dTm"):
        p = po[k]
        print(f"    {k:12s} {p['demanded_correction']:>9s}  CI{p['posterior_90CI']}  "
              f"constrained={p['constrained_by_curve']}")


# ---------------------------------------------------------------------------
# proteome-sectors: empirical temperature-dependent allocation + validation
# ---------------------------------------------------------------------------
def cmd_proteome_sectors(args):
    from . import proteome_alloc as pa
    cfg = resolve(args.strain, args.experiment)
    data = args.data or os.path.join(strain_dir(args.strain), "proteomics", "tem_proteomic.csv")
    if not os.path.exists(data):
        raise SystemExit(f"proteomics file not found: {data}")
    out_dir = _out_dir(args.strain, "proteome_sectors")
    os.makedirs(out_dir, exist_ok=True)

    # Force sectors on + temperature-dependent allocation for this diagnostic run.
    cfg.setdefault("proteome_sectors", {})
    cfg["proteome_sectors"]["enabled"] = True
    cfg["allocation_from_data"] = data
    pm = _build_pm(cfg)

    # -- PART 1: measured sector fractions vs temperature --
    b2u = pa.build_b_to_uniprot(pm.ec.model, cfg["provider"].get("prot_prefix", "prot_"))
    df = pa.load_temperature_proteome(data, b2u)
    sf = pa.sector_fractions_vs_T(df)
    sf.to_csv(os.path.join(out_dir, "sector_fractions_vs_T.csv"))
    if not args.no_plots:
        pa.plot_sector_fractions(sf, os.path.join(out_dir, "sector_fractions_vs_T.png"))
    ch = sf["f_chaperone"]
    print(f"[proteome] measured sector fractions vs T -> {out_dir}/sector_fractions_vs_T.csv")
    print(f"[proteome] chaperone fraction: {ch.loc[30]:.3f}@30°C -> {ch.loc[43]:.3f}@43°C "
          f"({ch.loc[43]/ch.loc[30]:.1f}x ramp)")

    # -- PART 2: mapping coverage --
    model_enz = set(e.enzyme_id for e in pm.ec.table.entries if e.enzyme_id)
    measured = set(df.dropna(subset=["UniProt"])["UniProt"])
    cov = model_enz & measured
    print(f"[proteome] mapping coverage: {len(cov)}/{len(model_enz)} model enzymes "
          f"({100*len(cov)/max(1,len(model_enz)):.1f}%) have a measured temperature profile")

    # -- PART 4: predicted vs measured usage + sector fractions --
    per_df, corr_df, pred_sec, meas_sec = pa.validate(pm, df)
    per_df.to_csv(os.path.join(out_dir, "validation_enzyme_usage.csv"), index=False)
    corr_df.to_csv(os.path.join(out_dir, "validation_correlations.csv"), index=False)
    pred_sec.to_csv(os.path.join(out_dir, "sector_fractions_predicted.csv"))
    meas_sec.to_csv(os.path.join(out_dir, "sector_fractions_measured_matched.csv"))
    if not args.no_plots:
        pa.plot_usage_pred_vs_meas(per_df, corr_df, os.path.join(out_dir, "usage_pred_vs_meas.png"))
        pa.plot_sector_pred_vs_meas(pred_sec, meas_sec, os.path.join(out_dir, "sector_pred_vs_meas.png"))
    print("[proteome] predicted-vs-measured per-enzyme usage:")
    for _, r in corr_df.iterrows():
        print(f"    {int(r.temp_C):>2}°C  n={int(r.n)}  Spearman ρ={r.spearman:.2f}  log-Pearson R²={r.log_pearson_r2:.2f}")
    dump_resolved(cfg, out_dir)
    print(f"[proteome] wrote {out_dir}")
    return out_dir


# ---------------------------------------------------------------------------
# dltkcat tooling
# ---------------------------------------------------------------------------
def cmd_dltkcat(args):
    from . import dltkcat as dk
    if args.dcmd == "prep":
        if args.strain:
            cfg = resolve(args.strain)
            model = cfg["provider"]["model_path"]
            out = args.out or os.path.join(strain_dir(args.strain), "dltkcat", "input.csv")
            t0 = args.t0 if args.t0 is not None else cfg.get("T0_C", 30.0)
        else:
            model, out, t0 = args.model, args.out, (args.t0 or 30.0)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        dk.build_dltkcat_input(model, out, args.tmin, args.tmax, args.n, t0)

    elif args.dcmd == "parse":
        if args.strain:
            d = os.path.join(strain_dir(args.strain), "dltkcat")
            pred = args.pred or os.path.join(d, "output.csv")
            out = args.out or os.path.join(d, "fits.csv")
        else:
            pred, out = args.pred, args.out
        fits = dk.parse_dltkcat_output(pred, kcat_col=args.kcat_col, key=args.key,
                                       temp_col=args.temp_col, T0_C=args.t0 or 30.0,
                                       log10=not args.no_log10)
        fits.to_csv(out, index=False)
        print(f"[dltkcat] fitted {int(fits['ok'].sum())}/{len(fits)} usable -> {out}")

    elif args.dcmd == "fit":
        fits = dk.fit_predictions(pd.read_csv(args.pred), key=args.key, T0_C=args.t0 or 30.0)
        fits.to_csv(args.out, index=False)
        print(f"[dltkcat] fitted {int(fits['ok'].sum())}/{len(fits)} usable -> {args.out}")

    elif args.dcmd == "csv":
        fits = dk.fit_predictions(pd.read_csv(args.pred), key=args.key, T0_C=args.t0 or 30.0)
        dk.write_csv_table(fits, args.model, args.out, key=args.key, T0_C=args.t0 or 30.0)

    elif args.dcmd == "targets":
        dk.export_targets(args.model, args.out, args.tmin, args.tmax, args.n, args.t0 or 30.0)


# ---------------------------------------------------------------------------
# argument parser
# ---------------------------------------------------------------------------
def build_parser():
    ap = argparse.ArgumentParser(
        prog="etcgem",
        description="Enzyme/temperature-constrained GEM thermal-performance tools.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # --- strain-only ---
    b = sub.add_parser("build", help="build a strain's provider; print+save summary")
    b.add_argument("--strain", required=True)
    b.set_defaults(func=cmd_build)

    t = sub.add_parser("tpc", help="nominal TPC + descriptors for a strain")
    t.add_argument("--strain", required=True)
    t.add_argument("--fits", nargs="?", const=_FITS_DEFAULT, default=None)
    t.add_argument("--key", default="rxn_id", choices=["rxn_id", "enzyme_id"])
    t.add_argument("--no-plots", action="store_true")
    t.set_defaults(func=cmd_tpc)

    fb = sub.add_parser("fba", help="single enzyme-constrained solve at one temperature")
    fb.add_argument("--strain", required=True)
    fb.add_argument("--temp", type=float, required=True, help="temperature in °C")
    fb.add_argument("--fits", nargs="?", const=_FITS_DEFAULT, default=None)
    fb.add_argument("--key", default="rxn_id", choices=["rxn_id", "enzyme_id"])
    fb.set_defaults(func=cmd_fba)

    cd = sub.add_parser("calibrate-dcp",
                        help="[DEPRECATED — not used by the emergent model] choose "
                             "provider.default_dCp so nominal Ea hits a target. The "
                             "emergent model does NOT fit dCp/Ea to the growth curve; "
                             "it grounds per-enzyme dCp in a literature prior + DLTKcat "
                             "and lets Ea emerge.")
    cd.add_argument("--strain", required=True)
    cd.add_argument("--target-ea", dest="target_ea", type=float, default=0.65,
                    help="target rising-limb Ea in eV (default 0.65, metabolic-theory value)")
    cd.add_argument("--lo", type=float, default=-20.0, help="dCp bisection lower bound")
    cd.add_argument("--hi", type=float, default=-3.0, help="dCp bisection upper bound")
    cd.add_argument("--tol", type=float, default=0.02, help="Ea tolerance (eV)")
    cd.set_defaults(func=cmd_calibrate_dcp)

    # --- strain + experiment ---
    sw = sub.add_parser("sweep", help="TPC sensitivity sweep (strain + experiment)")
    sw.add_argument("--strain")
    sw.add_argument("--experiment")
    sw.add_argument("--config", help="ad-hoc self-contained config (escape hatch)")
    sw.add_argument("--fits", nargs="?", const=_FITS_DEFAULT, default=None)
    sw.add_argument("--key", default="rxn_id", choices=["rxn_id", "enzyme_id"])
    sw.add_argument("--resume", action="store_true")
    sw.add_argument("--seconds", type=float, default=35.0)
    sw.add_argument("--no-plots", action="store_true")
    sw.set_defaults(func=cmd_sweep)

    dc = sub.add_parser("decompose", help="variance decomposition (strain + experiment)")
    dc.add_argument("--strain")
    dc.add_argument("--experiment")
    dc.add_argument("--no-plots", action="store_true")
    dc.set_defaults(func=cmd_decompose)

    ct = sub.add_parser("control", help="per-enzyme thermal control + identifiability (strain)")
    ct.add_argument("--strain")
    ct.add_argument("--experiment", help="experiment carrying the control: block (optional)")
    ct.add_argument("--no-plots", action="store_true")
    ct.set_defaults(func=cmd_control)

    el = sub.add_parser("elasticity",
                        help="equal-perturbation (standardised) sensitivity: elasticities "
                             "E[D,p] moving every input by the same step h")
    el.add_argument("--strain")
    el.add_argument("--experiment")
    el.add_argument("--h", type=float, default=None, help="standardised step (default from experiment / 0.10)")
    el.add_argument("--no-plots", action="store_true")
    el.set_defaults(func=cmd_elasticity)

    di = sub.add_parser("dissect",
                        help="dissect the TUNED model (P2 v3 posterior) at the rich BHI "
                             "operating point: posterior-propagated elasticity + envelope-vs-"
                             "magnitude decomposition + per-enzyme identifiability")
    di.add_argument("--strain")
    di.add_argument("--procs", type=int, default=None, help="parallel workers (0=auto)")
    di.add_argument("--draws", type=int, default=None, help="posterior draws for elasticity bands (default 150)")
    di.set_defaults(func=cmd_dissect)

    an = sub.add_parser("anatomy",
                        help="model-anatomy figures at the reference operating point "
                             "(reference TPC + per-enzyme parameter densities + example kcat(T))")
    an.add_argument("--strain")
    an.add_argument("--experiment")
    an.add_argument("--medium", default=None,
                    help="reference operating-point medium (default BHI; e.g. glucose_minimal, LB)")
    an.set_defaults(func=cmd_anatomy)

    va = sub.add_parser("validate",
                        help="emergent model vs trusted strain-matched TPCs "
                             "(Noll glucose-minimal + Erdos LB); absolute rates")
    va.add_argument("--strain")
    va.add_argument("--no-secondary", action="store_true",
                    help="skip the optional Erdos LB secondary cross-check")
    va.set_defaults(func=cmd_validate)

    ca = sub.add_parser("calibrate",
                        help="single-curve Bayesian calibration (emcee): prior vs "
                             "posterior on one measured glucose-minimal TPC")
    ca.add_argument("--strain")
    ca.add_argument("--curve", help="curve id from thermal/ecoli_tpc_curves.csv (legacy)")
    ca.add_argument("--noll", action="store_true",
                    help="fit the trusted Noll glucose-minimal curve with per-point-SD "
                         "likelihood (writes to calibration_noll_minimal/)")
    ca.add_argument("--vdl", action="store_true",
                    help="unified multi-parameter tuning on Van Derlinden (MG1655, BHI) at "
                         "the rich operating point (writes to calibration_vanderlinden/)")
    ca.add_argument("--experiment", help="experiment carrying a calibration: block (optional)")
    ca.add_argument("--medium", default=None)
    ca.add_argument("--walkers", type=int, default=None)
    ca.add_argument("--steps", type=int, default=None)
    ca.add_argument("--burn", type=int, default=None)
    ca.add_argument("--seed", type=int, default=None)
    ca.add_argument("--procs", type=int, default=None, help="parallel workers (0=auto)")
    ca.add_argument("--list", action="store_true", help="list defined curves and exit")
    ca.set_defaults(func=cmd_calibrate)

    pspar = sub.add_parser("proteome-sectors",
                           help="empirical temperature-dependent proteome allocation + "
                                "predicted-vs-measured validation")
    pspar.add_argument("--strain", required=True)
    pspar.add_argument("--experiment")
    pspar.add_argument("--data", help="temperature proteomics CSV "
                       "(default strains/NAME/proteomics/tem_proteomic.csv)")
    pspar.add_argument("--no-plots", action="store_true")
    pspar.set_defaults(func=cmd_proteome_sectors)

    # --- dltkcat tooling ---
    dl = sub.add_parser("dltkcat", help="DLTKcat -> MMRT (Topt, dCp) tooling")
    dsub = dl.add_subparsers(dest="dcmd", required=True)
    p = dsub.add_parser("prep", help="write DLTKcat convert_input CSV (enz, sub, temp)")
    p.add_argument("--strain"); p.add_argument("--model"); p.add_argument("--out")
    p.add_argument("--tmin", type=float, default=5.0)
    p.add_argument("--tmax", type=float, default=55.0)
    p.add_argument("--n", type=int, default=11)
    p.add_argument("--t0", type=float, default=None)
    pr = dsub.add_parser("parse", help="fit MMRT from a DLTKcat prediction table")
    pr.add_argument("--strain"); pr.add_argument("--pred"); pr.add_argument("--out")
    pr.add_argument("--kcat-col", default="kcat"); pr.add_argument("--key", default="rxn_id")
    pr.add_argument("--temp-col", default="Temp_C"); pr.add_argument("--t0", type=float, default=None)
    pr.add_argument("--no-log10", action="store_true")
    f = dsub.add_parser("fit", help="fit MMRT to a filled prediction table")
    f.add_argument("--pred", required=True); f.add_argument("--out", required=True)
    f.add_argument("--key", default="rxn_id"); f.add_argument("--t0", type=float, default=None)
    c = dsub.add_parser("csv", help="fit + write a from_kcat_csv provider table")
    c.add_argument("--pred", required=True); c.add_argument("--model", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--key", default="rxn_id"); c.add_argument("--t0", type=float, default=None)
    tg = dsub.add_parser("targets", help="export a prediction template")
    tg.add_argument("--model", required=True); tg.add_argument("--out", required=True)
    tg.add_argument("--tmin", type=float, default=5.0); tg.add_argument("--tmax", type=float, default=55.0)
    tg.add_argument("--n", type=int, default=11); tg.add_argument("--t0", type=float, default=None)
    dl.set_defaults(func=cmd_dltkcat)

    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    main()
