"""Config loading and provider dispatch for the CLI.

Three-tier config model, merged as ``defaults <- strain <- experiment``:

* ``configs/defaults.yaml``             -- universal METHOD defaults only
  (solver_timeout, crit_frac, a fallback temperature_grid).
* ``strains/<name>/strain.yaml``        -- ORGANISM only (provider block, T0_C,
  temperature_grid); self-sufficient to build + run the model.
* ``configs/experiments/<name>.yaml``   -- OPTIONAL method overlays (kind +
  sensitivity / decomposition blocks).

``resolve(strain, experiment=None)`` returns the merged dict with
``provider.model_path`` injected; ``dump_resolved`` records the exact merged
config into each run's output folder.
"""
from __future__ import annotations

import copy
import json
import os
from typing import Any, Dict, Optional

import numpy as np

from . import providers

# Project-root-relative locations (commands run from the project root). Config now
# lives under configs/; the legacy root locations are kept as a fallback so an old
# checkout still resolves.
DEFAULTS_FILE = "configs/defaults.yaml"
DEFAULTS_FILE_LEGACY = "defaults.yaml"
STRAINS_DIR = "strains"
EXPERIMENTS_DIR = "configs/experiments"
EXPERIMENTS_DIR_LEGACY = "experiments"


def load_config(path: str) -> Dict[str, Any]:
    with open(path) as fh:
        text = fh.read()
    if path.lower().endswith((".yaml", ".yml")):
        try:
            import yaml
            return yaml.safe_load(text)
        except ImportError as e:
            raise ImportError("Install pyyaml, or use a .json config.") from e
    return json.loads(text)


# ---------------------------------------------------------------------------
# three-tier config model
# ---------------------------------------------------------------------------
def _deep_merge(base: Dict[str, Any], over: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge ``over`` onto a copy of ``base`` (dicts merge, scalars/
    lists replace)."""
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def strain_dir(name: str) -> str:
    return os.path.abspath(os.path.join(STRAINS_DIR, name))


def load_defaults() -> Dict[str, Any]:
    for p in (DEFAULTS_FILE, DEFAULTS_FILE_LEGACY):
        if os.path.exists(p):
            return load_config(p)
    return {}


def load_strain(name: str) -> Dict[str, Any]:
    return load_config(os.path.join(strain_dir(name), "strain.yaml"))


def load_experiment(name: str) -> Dict[str, Any]:
    for d in (EXPERIMENTS_DIR, EXPERIMENTS_DIR_LEGACY):
        p = os.path.join(d, f"{name}.yaml")
        if os.path.exists(p):
            return load_config(p)
    return load_config(os.path.join(EXPERIMENTS_DIR, f"{name}.yaml"))  # raise on the new path


def resolve(strain: str, experiment: Optional[str] = None) -> Dict[str, Any]:
    """Compose defaults <- strain <- experiment and inject provider.model_path.

    Organism keys come from the strain; method keys from defaults, overridden by
    the experiment. The absolute model path is derived from the strain folder +
    ``provider.model_file``; ``output_dir`` is left for the CLI to set per run.
    """
    cfg = _deep_merge(load_defaults(), load_strain(strain))
    if experiment:
        cfg = _deep_merge(cfg, load_experiment(experiment))
    model_file = cfg.get("provider", {}).get("model_file")
    if model_file:
        cfg["provider"]["model_path"] = os.path.join(strain_dir(strain), "model", model_file)
    cfg.setdefault("_strain", strain)
    if experiment:
        cfg["_experiment"] = experiment
    return cfg


def dump_resolved(cfg: Dict[str, Any], out_dir: str) -> str:
    """Write the exact merged config used for a run into its output folder."""
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "resolved_config.yaml")
    try:
        import yaml
        with open(path, "w") as fh:
            yaml.safe_dump(cfg, fh, sort_keys=False)
    except ImportError:
        with open(path.replace(".yaml", ".json"), "w") as fh:
            json.dump(cfg, fh, indent=2)
        path = path.replace(".yaml", ".json")
    return path


def build_provider(cfg: Dict[str, Any]):
    p = cfg["provider"]
    T0 = (cfg.get("T0_C", 30.0)) + 273.15
    kind = p["type"]
    if kind == "toy":
        pm = providers.toy_ecoli_core(
            T0=T0, seed=p.get("seed", 0),
            topt_mean_offset=p.get("topt_mean_offset", 8.0),
            topt_sd=p.get("topt_sd", 6.0),
            dCp_mean=p.get("dCp_mean", -8.0),
            target_fraction=p.get("target_fraction", 0.6),
        )
    elif kind == "gecko":
        enzyme_params = p.get("enzyme_params")
        if enzyme_params and not os.path.isabs(enzyme_params) and cfg.get("_strain"):
            cand = os.path.join(strain_dir(cfg["_strain"]), enzyme_params)
            if os.path.exists(cand):
                enzyme_params = cand
        # EMERGENT magnitude: budget = P_total x f_metab x sigma from independent
        # data/literature (not the growth-calibrated GECKO bound). Active when
        # provider.p_total is set; f_metab from the (measured) sector nominal.
        p_total = p.get("p_total")
        sigma = p.get("sigma", 0.45)
        budget_override = None
        if p_total is not None:
            f_metab = (cfg.get("proteome_sectors") or {}).get("f_metab", 0.285)
            budget_override = float(p_total) * float(sigma) * float(f_metab)
        pm = providers.from_gecko(
            model_path=p["model_path"], T0=T0,
            default_Topt_offset=p.get("default_Topt_offset", 7.0),
            default_dCp=p.get("default_dCp", -12.0),
            prot_prefix=p.get("prot_prefix", "prot_"),
            pool_id=p.get("pool_id", "prot_pool"),
            biomass_rxn=p.get("biomass_rxn"),
            target_fraction=p.get("target_fraction"),
            pool_scale=p.get("pool_scale", 1.0),
            thermal_model=p.get("thermal_model", "mmrt"),
            ngam_temperature=p.get("ngam_temperature", False),
            ngam_rxn=p.get("ngam_reaction"),
            enzyme_params=enzyme_params,
            budget_override=budget_override,
            enzyme_params_use_dCpt=(p.get("dcp_from", "table") != "prior"),
            dcp_prior_kJ=p.get("dcp_prior_kJ", -4.0),
            close_free_o2_sinks=cfg.get("close_free_o2_sinks", True),
        )
        if budget_override is not None:
            print(f"[emergent] pool budget = P_total({p_total}) x f_metab({f_metab}) "
                  f"x sigma({sigma}) = {budget_override:.4g} g/gDW (not growth-calibrated)")
        # record exactly which reactions the closure touched, for the resolved config
        cfg["closed_free_o2_sinks"] = list(pm.closed_free_o2_sinks)
    elif kind == "csv":
        pm = providers.from_kcat_csv(
            model_path=p["model_path"], csv_path=p["csv_path"], T0=T0,
            default_Topt_offset=p.get("default_Topt_offset", 8.0),
            default_dCp=p.get("default_dCp", -8.0),
            biomass_rxn=p.get("biomass_rxn"),
            target_fraction=p.get("target_fraction", 0.6),
        )
    else:
        raise ValueError(f"Unknown provider type: {kind}")

    # Opt-in DLTKcat per-enzyme Topt/dCp overlay (baked into the strain). Applied
    # before sectors so the sector auto-calibration sees the final envelope. In
    # unfolding mode only Topt is used; enzymes without a fit keep their grounded
    # (or default) parameters.
    dl = p.get("dltkcat_fits")
    if dl:
        if not os.path.isabs(dl) and cfg.get("_strain"):
            cand = os.path.join(strain_dir(cfg["_strain"]), dl)
            if os.path.exists(cand):
                dl = cand
        import pandas as pd
        from .dltkcat import apply_fits_to_provider
        apply_fits_to_provider(pm, pd.read_csv(dl), key=p.get("dltkcat_key", "rxn_id"))

    # Opt-in proteome sectors (Basan/Scott). Absent or disabled -> untouched.
    ps = cfg.get("proteome_sectors")
    if ps and ps.get("enabled"):
        from .sectors import add_proteome_sectors
        # thread the emergent in-vivo saturation (budget = P_total x f_metab x sigma)
        # so a free sigma_sat perturbation can scale both sector caps by sigma/sigma_nom
        ps.setdefault("sigma_nom", float(p.get("sigma", 0.45)))
        add_proteome_sectors(pm, ps)

        # Opt-in temperature-dependent allocation from measured proteomics. Only
        # active when sectors are enabled AND a data file is configured; otherwise
        # sectors stay temperature-independent (unchanged behaviour).
        alloc_data = cfg.get("allocation_from_data")
        if alloc_data:
            if not os.path.isabs(alloc_data) and cfg.get("_strain"):
                cand = os.path.join(strain_dir(cfg["_strain"]), alloc_data)
                if os.path.exists(cand):
                    alloc_data = cand
            from . import proteome_alloc as pa
            b2u = pa.build_b_to_uniprot(pm.ec.model, p.get("prot_prefix", "prot_"))
            df = pa.load_temperature_proteome(alloc_data, b2u)
            # MEDIUM-matched, temperature-dependent sector fractions (per LB/Glucose/
            # Glycerol series). The default (Glucose) drives the nominal; set_medium
            # switches the active medium so LB uses the LB fractions (growth-law).
            sfm = pa.sector_fractions_by_medium(df)
            default_med = cfg.get("default_medium_proteome", "Glucose")
            alloc = pa.TemperatureAllocation.from_medium_fractions(sfm, default_medium=default_med)
            alloc.set_active_medium("glucose_minimal")
            pm.ec._alloc_from_data = alloc
            fb = {m: round(float(sfm[(sfm.medium == m) & (sfm.temp_C == 37)]["f_bio"].mean()), 3)
                  for m in sfm.medium.unique() if ((sfm.medium == m) & (sfm.temp_C == 37)).any()}
            print(f"[alloc] MEDIUM-matched sector allocation from {os.path.basename(str(alloc_data))} "
                  f"(f_bio@37C by medium: {fb})")
    return pm


def temperature_grid(cfg: Dict[str, Any]) -> np.ndarray:
    g = cfg["temperature_grid"]
    return np.linspace(g["start_C"], g["stop_C"], int(g["n"]))


# ---------------------------------------------------------------------------
# dCp -> Ea calibration
# ---------------------------------------------------------------------------
def _nominal_ea_ctmax(cfg: Dict[str, Any], dCp: float):
    """Rebuild the provider with provider.default_dCp = dCp and return the
    nominal (unperturbed) rising-limb Ea (eV) and CT_max (°C)."""
    from .tpc import compute_tpc
    from .enzyme_cost import Perturbation
    c = copy.deepcopy(cfg)
    c["provider"]["default_dCp"] = float(dCp)
    pm = build_provider(c)
    try:
        pm.ec.model.solver.configuration.timeout = int(c.get("solver_timeout", 10))
    except Exception:
        pass
    d = compute_tpc(pm, temperature_grid(c), Perturbation()).descriptors(
        c.get("crit_frac", 0.05))
    return float(d.Ea_eV), float(d.CTmax_C)


def calibrate_dCp_to_Ea(cfg: Dict[str, Any], target_Ea_eV: float,
                        lo: float = -20.0, hi: float = -3.0, tol: float = 0.02,
                        max_iter: int = 40, verbose: bool = True) -> Dict[str, Any]:
    """Choose provider.default_dCp so the nominal rising-limb Ea hits a target.

    The nominal Ea increases monotonically with |dCp| (a more negative dCp gives
    a steeper, narrower TPC), so this bisects dCp in ``[lo, hi]``, rebuilding the
    provider and recomputing the nominal TPC each step. Returns a dict with the
    calibrated ``dCp``, achieved ``Ea`` and ``CTmax`` (on the config's grid),
    ``iters`` and whether the target was ``bracketed``.
    """
    ea_lo, _ = _nominal_ea_ctmax(cfg, lo)    # steepest -> highest Ea
    ea_hi, _ = _nominal_ea_ctmax(cfg, hi)    # shallowest -> lowest Ea
    if verbose:
        print(f"[calibrate] Ea(dCp={lo})={ea_lo:.3f}eV  Ea(dCp={hi})={ea_hi:.3f}eV  "
              f"target={target_Ea_eV:.3f}eV")
    bracketed = min(ea_lo, ea_hi) <= target_Ea_eV <= max(ea_lo, ea_hi)
    if not bracketed:
        # clamp to the nearest endpoint
        pick = lo if abs(ea_lo - target_Ea_eV) < abs(ea_hi - target_Ea_eV) else hi
        ea, ctmax = _nominal_ea_ctmax(cfg, pick)
        if verbose:
            print(f"[calibrate] WARNING: target {target_Ea_eV:.3f}eV not bracketed "
                  f"by [{min(ea_lo, ea_hi):.3f}, {max(ea_lo, ea_hi):.3f}]; "
                  f"clamping to dCp={pick}")
        return {"dCp": float(pick), "Ea": ea, "CTmax": ctmax,
                "iters": 0, "bracketed": False}

    a, b = lo, hi
    best_dcp, best_ea, best_ctmax = hi, ea_hi, None
    it = 0
    for it in range(1, max_iter + 1):
        mid = 0.5 * (a + b)
        ea, ctmax = _nominal_ea_ctmax(cfg, mid)
        best_dcp, best_ea, best_ctmax = mid, ea, ctmax
        if verbose:
            print(f"[calibrate] iter {it:2d}  dCp={mid:7.3f}  Ea={ea:.3f}eV  CTmax={ctmax:.1f}°C")
        if abs(ea - target_Ea_eV) <= tol:
            break
        # Ea decreasing in dCp: too-high Ea -> need shallower (larger) dCp
        if ea > target_Ea_eV:
            a = mid
        else:
            b = mid
    return {"dCp": float(best_dcp), "Ea": float(best_ea), "CTmax": float(best_ctmax),
            "iters": it, "bracketed": True}
