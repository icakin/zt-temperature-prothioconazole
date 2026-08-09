"""Bayesian calibration of the etc-GEM against a measured growth TPC (Phase 0+1).

The emergent model is the *a-priori* prediction -- the PRIOR. This module asks the
inverse question on ONE well-characterised glucose-minimal curve: what would the
borrowed/uncertain parameters have to be to reproduce the data, and how well does
the curve constrain them?

Method (fixed): the etc-GEM is a DETERMINISTIC simulator (parameters -> TPC) with a
Gaussian residual-noise model, so we do EXACT likelihood-based inference with a
gradient-free ensemble sampler (emcee). NOT ABC, NOT Stan/brms (the FBA simulator
is non-differentiable).

Phase-1 free set (only knobs the elasticity says can move the fit targets):
    kappa_scale  -- in-vivo translation efficiency; the MAGNITUDE / r_max lever
    dCp_scale    -- MMRT curvature; the rising-limb E_a lever
    dTopt        -- uniform optimum shift
    dTm          -- uniform melting-temperature shift; the upper-limit lever
plus a nuisance noise scale sigma_noise. The metabolic saturation sigma and the
measured sector fractions are held FIXED at their emergent values (sigma is
degenerate with kappa_scale for magnitude; a later phase explores it).

Positive-only parameters are sampled in log space so the emcee walkers live in an
unconstrained R^n; priors are centred on the emergent values with provenance-based
widths (see build_priors / the emitted priors table).
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .enzyme_cost import Perturbation
from .tpc import compute_tpc

# free-parameter vector layout (fixed order); positive params carried in log space
PARAM_NAMES = ["dTopt", "dTm", "log_dCp_scale", "log_kappa_scale", "log_sigma"]
NAT_NAMES = ["dTopt", "dTm", "dCp_scale", "kappa_scale", "sigma_noise"]
NDIM = len(PARAM_NAMES)


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def curve_path(strain: str) -> str:
    return os.path.join("strains", strain, "thermal", "ecoli_tpc_curves.csv")


def load_curve(strain: str, curve_id: str) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Return (temps_C, rate_per_h, meta) for one curve id (raw absolute rate)."""
    df = pd.read_csv(curve_path(strain))
    sub = df[df["curve_id"].astype(str) == str(curve_id)].copy()
    if sub.empty:
        raise SystemExit(f"curve id {curve_id!r} not found in {curve_path(strain)}")
    sub = sub.sort_values("temp_C")
    temps = sub["temp_C"].to_numpy(float)
    rates = sub["rate_per_h"].to_numpy(float)
    meta = {"curve_id": str(curve_id),
            "study": str(sub["study"].iloc[0]),
            "medium_class": str(sub["medium_class"].iloc[0]),
            "defined": bool(sub["defined"].iloc[0]),
            "n": int(len(sub)),
            "temp_min_C": float(temps.min()), "temp_max_C": float(temps.max()),
            "obs_rmax": float(rates.max()),
            "obs_Topt_C": float(temps[int(np.argmax(rates))]),
            "units": "1/h"}
    return temps, rates, meta


def list_defined_curves(strain: str) -> pd.DataFrame:
    df = pd.read_csv(curve_path(strain))
    g = (df[df["defined"]].groupby(["curve_id", "study", "medium_class"])
         .agg(n=("temp_C", "size"), Tmin=("temp_C", "min"), Tmax=("temp_C", "max"),
              rmax=("rate_per_h", "max")).reset_index())
    return g


# ---------------------------------------------------------------------------
# priors
# ---------------------------------------------------------------------------
def build_priors(obs_rates: np.ndarray, cfg: Optional[dict] = None) -> Dict:
    """Priors centred on the emergent values, widths from provenance. Documented
    so every choice is auditable; overridable via the `calibration.priors` cfg block."""
    cfg = cfg or {}
    # weakly-informative noise scale from the data scatter (HalfNormal scale)
    sigma_scale = float(cfg.get("sigma_halfnormal_scale",
                                max(0.02, 0.5 * float(np.std(obs_rates)))))
    pri = {
        # additive temperature shifts: centred 0, scale = per-enzyme Topt/Tm SD
        "dTopt":  {"dist": "normal", "loc": 0.0, "scale": float(cfg.get("dTopt_sd", 4.25)),
                   "note": "per-enzyme Topt SD (4.25 K); optimum-prediction spread"},
        "dTm":    {"dist": "normal", "loc": 0.0, "scale": float(cfg.get("dTm_sd", 6.13)),
                   "note": "per-enzyme Tm SD (6.13 K); melting-proteome prediction error"},
        # multiplicative shape/magnitude: LogNormal about 1 (mean of log = 0)
        "dCp_scale":   {"dist": "lognormal", "log_loc": 0.0,
                        "log_scale": float(cfg.get("log_dCp_scale_sd", 0.30)),
                        "note": "MMRT curvature multiplier ~ LogNormal about 1 (sd 0.30 on log)"},
        "kappa_scale": {"dist": "lognormal", "log_loc": 0.0,
                        "log_scale": float(cfg.get("log_kappa_scale_sd", 0.70)),
                        "note": "in-vivo translation efficiency ~ LogNormal about 1 "
                                "(sd 0.70 on log ~ up to ~4x either way); broad/uncertain"},
        # nuisance noise
        "sigma_noise": {"dist": "halfnormal", "scale": sigma_scale,
                        "note": "Gaussian residual sd on absolute rate (1/h); "
                                "weakly-informative from data scatter"},
        # hard support bounds: keep the sampler in a physical, numerically
        # well-behaved region (GLPK can stall on the degenerate LPs that extreme
        # cost distributions produce). These comfortably contain the posterior.
        "bounds": {"dTopt": [-12.0, 12.0], "dTm": [-15.0, 15.0],
                   "dCp_scale": [0.25, 4.0], "kappa_scale": [0.05, 20.0],
                   "sigma_noise": [1e-4, 2.0]},
    }
    return pri


def theta_to_natural(theta: Sequence[float]) -> Dict[str, float]:
    dTopt, dTm, log_dCp, log_kappa, log_sigma = theta
    return {"dTopt": float(dTopt), "dTm": float(dTm),
            "dCp_scale": float(np.exp(log_dCp)),
            "kappa_scale": float(np.exp(log_kappa)),
            "sigma_noise": float(np.exp(log_sigma))}


def theta_to_pert(theta: Sequence[float]) -> Perturbation:
    n = theta_to_natural(theta)
    return Perturbation(dTopt=n["dTopt"], dTm=n["dTm"],
                        dCp_scale=n["dCp_scale"], kappa_scale=n["kappa_scale"])


def log_prior(theta: Sequence[float], pri: Dict) -> float:
    """log prior density in the sampled (log-for-positives) parameterisation,
    including the change-of-variable Jacobians for the log-transformed params."""
    nat = theta_to_natural(theta)
    b = pri["bounds"]
    for k, (lo, hi) in b.items():
        if not (lo <= nat[k] <= hi):
            return -np.inf
    lp = 0.0
    # dTopt, dTm: Normal
    for k in ("dTopt", "dTm"):
        s = pri[k]["scale"]
        lp += -0.5 * (nat[k] / s) ** 2 - np.log(s * np.sqrt(2 * np.pi))
    # dCp_scale, kappa_scale: LogNormal about 1 -> Normal on log with Jacobian +log(x)
    log_dCp, log_kappa, log_sigma = theta[2], theta[3], theta[4]
    for k, u in (("dCp_scale", log_dCp), ("kappa_scale", log_kappa)):
        s = pri[k]["log_scale"]
        # density of log-param u ~ Normal(0, s) (this IS the LogNormal in u-space)
        lp += -0.5 * (u / s) ** 2 - np.log(s * np.sqrt(2 * np.pi))
    # sigma_noise: HalfNormal(scale) on sigma, sampled as log_sigma -> +log(sigma) Jacobian
    s = pri["sigma_noise"]["scale"]
    sig = nat["sigma_noise"]
    lp += (-0.5 * (sig / s) ** 2 + np.log(np.sqrt(2 / np.pi) / s)) + log_sigma
    return float(lp)


# ---------------------------------------------------------------------------
# simulator + likelihood
# ---------------------------------------------------------------------------
@dataclass
class Simulator:
    pm: object
    temps_C: np.ndarray

    def predict(self, theta: Sequence[float]) -> np.ndarray:
        """Predicted absolute growth rate (1/h) at the curve temperatures."""
        pert = theta_to_pert(theta)
        return compute_tpc(self.pm, self.temps_C, pert).growth


def log_likelihood(theta, sim: Simulator, obs: np.ndarray) -> float:
    sigma = float(np.exp(theta[4]))
    pred = sim.predict(theta)
    if not np.all(np.isfinite(pred)):
        return -np.inf
    r = obs - pred
    n = len(obs)
    return float(-0.5 * np.sum((r / sigma) ** 2) - n * np.log(sigma * np.sqrt(2 * np.pi)))


def log_likelihood_sd(theta, sim: Simulator, obs: np.ndarray, obs_sd: np.ndarray) -> float:
    """Gaussian likelihood with a per-point measurement variance plus a free model-
    discrepancy variance: var_i = obs_sd_i^2 + sigma_disc^2 (theta[4]=log sigma_disc).
    Weights each temperature by its measured precision while absorbing structural
    model error."""
    sigma_disc = float(np.exp(theta[4]))
    pred = sim.predict(theta)
    if not np.all(np.isfinite(pred)):
        return -np.inf
    var = np.asarray(obs_sd, float) ** 2 + sigma_disc ** 2
    r = obs - pred
    return float(-0.5 * np.sum(r ** 2 / var + np.log(2.0 * np.pi * var)))


# --- module globals for multiprocessing workers (non-picklable cobra model) ---
_SIM: Optional[Simulator] = None
_OBS: Optional[np.ndarray] = None
_OBS_SD: Optional[np.ndarray] = None
_PRI: Optional[Dict] = None


def load_noll_curve(strain: str):
    """Load the trusted Noll glucose-minimal curve (temps 28-45 C with rate + SD)
    via the shared validation loader, so calibration and validation read the SAME
    data. Returns (temps_C, rate, sd, meta)."""
    from . import validation as V
    spec = next(s for s in V.CURVES if s.key == "noll_minimal")
    cv = V.load_curve(strain, spec)
    sd = cv["sd"]
    if sd is None:
        sd = np.full_like(cv["rate"], 0.05)
    sd = np.where(np.isfinite(sd) & (sd > 0), sd, np.nanmedian(sd[sd > 0]) if np.any(sd > 0) else 0.05)
    meta = {"curve_id": "Noll2023_NCM3722", "study": cv["study"], "strain": cv["strain"],
            "medium_class": "glucose_minimal", "medium_detail": cv["medium_detail"],
            "n": int(len(cv["temps_C"])),
            "temp_min_C": float(cv["temps_C"].min()), "temp_max_C": float(cv["temps_C"].max()),
            "obs_rmax": float(cv["obs_rmax"]), "obs_Topt_C": float(cv["obs_Topt_C"]),
            "units": "1/h", "has_sd": True}
    return cv["temps_C"], cv["rate"], sd, meta


def _set_solver_timeout(pm, seconds):
    """Cap each LP solve so a slow/degenerate perturbation cannot stall the run
    (a timed-out solve returns no growth, which compute_tpc reads as rate 0)."""
    try:
        pm.ec.model.solver.configuration.timeout = int(seconds)
    except Exception:
        pass


def _worker_init(strain, curve_id, medium, priors, trusted_noll=False):
    global _SIM, _OBS, _OBS_SD, _PRI
    from .config import resolve, build_provider
    from .providers import set_medium
    pm = build_provider(resolve(strain))
    _set_solver_timeout(pm, 2)   # cap slow/degenerate LP solves so no eval stalls
    set_medium(pm, medium, "glc__D", True)
    if trusted_noll:
        temps, rates, sd, _ = load_noll_curve(strain)
        _OBS_SD = sd
    else:
        temps, rates, _ = load_curve(strain, curve_id)
        _OBS_SD = None
    _SIM = Simulator(pm, temps)
    _OBS = rates
    _PRI = priors


def _log_prob_worker(theta):
    lp = log_prior(theta, _PRI)
    if not np.isfinite(lp):
        return -np.inf
    if _OBS_SD is not None:
        return lp + log_likelihood_sd(theta, _SIM, _OBS, _OBS_SD)
    return lp + log_likelihood(theta, _SIM, _OBS)


def make_log_prob(sim: Simulator, obs: np.ndarray, pri: Dict, obs_sd=None):
    """Single-process log-prob closure (used when no pool)."""
    def _lp(theta):
        lp = log_prior(theta, pri)
        if not np.isfinite(lp):
            return -np.inf
        if obs_sd is not None:
            return lp + log_likelihood_sd(theta, sim, obs, obs_sd)
        return lp + log_likelihood(theta, sim, obs)
    return _lp


# ---------------------------------------------------------------------------
# runner
# ---------------------------------------------------------------------------
def run_emcee(strain: str, curve_id: str, out_dir: str, *, medium="glucose_minimal",
              n_walkers: int = 24, n_steps: int = 1500, n_burn: int = 500,
              seed: int = 1, n_proc: int = 0, priors_cfg: Optional[dict] = None,
              trusted_noll: bool = False, progress: bool = True) -> Dict:
    """Run the single-curve emcee calibration and write all outputs.

    ``trusted_noll``: fit the trusted Noll glucose-minimal curve (loaded via the
    validation loader) with a per-point-SD + model-discrepancy likelihood; the 5th
    free parameter is then sigma_disc (structural model error) rather than a single
    homoscedastic sigma_noise."""
    import emcee

    os.makedirs(out_dir, exist_ok=True)
    np.random.seed(seed)

    # build the main-process provider (for pred TPCs + priors) and the curve
    from .config import resolve, build_provider
    from .providers import set_medium
    pm = build_provider(resolve(strain))
    _set_solver_timeout(pm, 2)
    set_medium(pm, medium, "glc__D", True)
    if trusted_noll:
        temps, obs, obs_sd, meta = load_noll_curve(strain)
        # model-discrepancy sigma needs room to absorb structural misfit (the
        # emergent model under-predicts by ~0.5 1/h), so widen its HalfNormal prior
        priors_cfg = dict(priors_cfg or {})
        priors_cfg.setdefault("sigma_halfnormal_scale", 0.3)
    else:
        temps, obs, meta = load_curve(strain, curve_id)
        obs_sd = None
    sim = Simulator(pm, temps)
    pri = build_priors(obs, priors_cfg)

    # initialise walkers around the emergent point, spread over a fraction of the
    # prior widths so the ensemble covers the space but no walker starts in an
    # extreme, numerically degenerate region (keeps GLPK well-behaved at burn-in).
    rng0 = np.random.default_rng(seed)
    b = pri["bounds"]
    p0 = np.column_stack([
        np.clip(rng0.normal(0.0, 0.6 * pri["dTopt"]["scale"], n_walkers), *b["dTopt"]),
        np.clip(rng0.normal(0.0, 0.6 * pri["dTm"]["scale"], n_walkers), *b["dTm"]),
        rng0.normal(0.0, 0.6 * pri["dCp_scale"]["log_scale"], n_walkers),
        rng0.normal(0.0, 0.6 * pri["kappa_scale"]["log_scale"], n_walkers),
        np.log(pri["sigma_noise"]["scale"]) + 0.2 * rng0.standard_normal(n_walkers),
    ])

    n_proc = n_proc or max(1, (os.cpu_count() or 2) - 2)
    t0 = time.time()
    if n_proc > 1:
        from multiprocessing import Pool
        pool = Pool(processes=n_proc, initializer=_worker_init,
                    initargs=(strain, curve_id, medium, pri, trusted_noll))
        try:
            sampler = emcee.EnsembleSampler(n_walkers, NDIM, _log_prob_worker, pool=pool)
            _run_chain(sampler, p0, n_steps, progress)
        finally:
            pool.close(); pool.join()
    else:
        lp = make_log_prob(sim, obs, pri, obs_sd)
        sampler = emcee.EnsembleSampler(n_walkers, NDIM, lp)
        _run_chain(sampler, p0, n_steps, progress)
    wall = time.time() - t0

    # diagnostics
    try:
        tau = sampler.get_autocorr_time(tol=0)
    except Exception:
        tau = np.full(NDIM, np.nan)
    tau_max = float(np.nanmax(tau)) if np.isfinite(tau).any() else float("nan")
    burn = int(n_burn if np.isnan(tau_max) else max(n_burn, 2 * tau_max))
    burn = min(burn, n_steps - 10)
    thin = max(1, int(tau_max / 2)) if np.isfinite(tau_max) else 1
    flat = sampler.get_chain(discard=burn, thin=thin, flat=True)
    accept = float(np.mean(sampler.acceptance_fraction))
    n_eff = flat.shape[0] if not np.isfinite(tau_max) else float(
        n_walkers * (n_steps - burn) / tau_max)

    result = {
        "strain": strain, "curve": meta, "medium": medium,
        "sampler": {"n_walkers": n_walkers, "n_steps": n_steps, "burn": burn,
                    "thin": thin, "acceptance_fraction": round(accept, 3),
                    "autocorr_time": [None if not np.isfinite(x) else round(float(x), 1) for x in tau],
                    "autocorr_time_max": None if not np.isfinite(tau_max) else round(tau_max, 1),
                    "n_eff": round(float(n_eff), 1), "wall_time_s": round(wall, 1),
                    "n_proc": n_proc, "seed": seed},
        "priors": pri,
    }
    _finalise_outputs(out_dir, flat, sim, obs, pri, meta, result, pm, obs_sd=obs_sd)
    return result


def _run_chain(sampler, p0, n_steps, progress):
    sampler.run_mcmc(p0, n_steps, progress=progress)


# ---------------------------------------------------------------------------
# outputs: posterior summary, prior-vs-posterior TPC, corrections, corner
# ---------------------------------------------------------------------------
def _posterior_natural(flat: np.ndarray) -> Dict[str, np.ndarray]:
    return {"dTopt": flat[:, 0], "dTm": flat[:, 1],
            "dCp_scale": np.exp(flat[:, 2]), "kappa_scale": np.exp(flat[:, 3]),
            "sigma_noise": np.exp(flat[:, 4])}


def _ci(x, q=(5, 50, 95)):
    return [float(v) for v in np.percentile(x, q)]


def _finalise_outputs(out_dir, flat, sim, obs, pri, meta, result, pm, obs_sd=None):
    np.save(os.path.join(out_dir, "chain_flat.npy"), flat)
    post = _posterior_natural(flat)
    noise_key = "sigma_disc" if obs_sd is not None else "sigma_noise"
    post[noise_key] = post.pop("sigma_noise")

    # --- per-parameter posterior summary + demanded corrections ---
    prior_sd = {"dTopt": pri["dTopt"]["scale"], "dTm": pri["dTm"]["scale"],
                "dCp_scale": pri["dCp_scale"]["log_scale"],
                "kappa_scale": pri["kappa_scale"]["log_scale"]}
    rows, summary = [], {}
    emergent = {"dTopt": 0.0, "dTm": 0.0, "dCp_scale": 1.0, "kappa_scale": 1.0, noise_key: None}
    for k in ("dTopt", "dTm", "dCp_scale", "kappa_scale", noise_key):
        lo, med, hi = _ci(post[k])
        summary[k] = {"prior_center": emergent[k],
                      "posterior_median": round(med, 4),
                      "posterior_90CI": [round(lo, 4), round(hi, 4)]}
        # is the posterior ~ the prior (curve does not constrain it)?
        constrained = None
        if k in ("dTopt", "dTm"):
            constrained = (hi - lo) < 1.6 * prior_sd[k]      # tighter than the prior
            corr = f"{med:+.1f} K"
        elif k in ("dCp_scale", "kappa_scale"):
            post_logsd = float(np.std(np.log(post[k])))
            constrained = post_logsd < 0.8 * prior_sd[k]
            corr = f"x{med:.2f}"
        else:
            corr = f"{med:.3f} 1/h"
        summary[k]["constrained_by_curve"] = None if constrained is None else bool(constrained)
        summary[k]["demanded_correction"] = corr
        rows.append({"parameter": k, "prior_center": emergent[k],
                     "prior_width": prior_sd.get(k, pri["sigma_noise"]["scale"]),
                     "posterior_median": round(med, 4),
                     "ci5": round(lo, 4), "ci95": round(hi, 4),
                     "demanded_correction": corr,
                     "constrained_by_curve": constrained})
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "demanded_corrections.csv"), index=False)
    result["posterior"] = summary

    # --- descriptors: prior vs posterior-median vs observed ---
    from .tpc import TPC
    # ~1 C grid (coarse enough to be fast for the predictive band, fine enough for descriptors)
    t_lo = min(float(sim.temps_C.min()), 5.0)
    t_hi = max(float(sim.temps_C.max()), 55.0)
    # ~1.7 C grid: fine enough to read Topt/CTmax, coarse enough to keep the
    # posterior-predictive band (many curves) affordable
    dense = np.linspace(t_lo, t_hi, int(round((t_hi - t_lo) / 1.7)) + 1)
    prior_curve = compute_tpc(pm, dense, Perturbation()).growth
    med_theta = [np.median(flat[:, j]) for j in range(NDIM)]
    post_med_curve = compute_tpc(pm, dense, theta_to_pert(med_theta)).growth
    d_obs = TPC(sim.temps_C, obs).descriptors(0.05)
    d_prior = TPC(dense, prior_curve).descriptors(0.05)
    d_post = TPC(dense, post_med_curve).descriptors(0.05)
    result["descriptors"] = {
        "observed": {"rmax": round(d_obs.rmax, 3), "Topt_C": round(d_obs.Topt_C, 1),
                     "Ea_eV": round(d_obs.Ea_eV, 3), "CTmax_C": round(d_obs.CTmax_C, 1)},
        "prior": {"rmax": round(d_prior.rmax, 3), "Topt_C": round(d_prior.Topt_C, 1),
                  "Ea_eV": round(d_prior.Ea_eV, 3), "CTmax_C": round(d_prior.CTmax_C, 1)},
        "posterior_median": {"rmax": round(d_post.rmax, 3), "Topt_C": round(d_post.Topt_C, 1),
                             "Ea_eV": round(d_post.Ea_eV, 3), "CTmax_C": round(d_post.CTmax_C, 1)},
    }

    # --- posterior-predictive band (subsample the chain) ---
    rng = np.random.default_rng(0)
    n_pp = min(120, flat.shape[0])
    pick = rng.choice(flat.shape[0], size=n_pp, replace=False)
    pp = np.vstack([compute_tpc(pm, dense, theta_to_pert(flat[i])).growth for i in pick])
    pp_lo, pp_med, pp_hi = np.percentile(pp, [5, 50, 95], axis=0)
    np.savez(os.path.join(out_dir, "posterior_predictive.npz"),
             temps_C=dense, lo=pp_lo, med=pp_med, hi=pp_hi,
             prior=prior_curve, obs_T=sim.temps_C, obs=obs,
             obs_sd=(obs_sd if obs_sd is not None else np.zeros_like(obs)))

    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(result, fh, indent=2)

    _plot_prior_vs_posterior(out_dir, dense, prior_curve, pp_lo, pp_med, pp_hi,
                             sim.temps_C, obs, meta, obs_sd=obs_sd)
    _plot_corner(out_dir, flat, noise_key=noise_key)


def _plot_prior_vs_posterior(out_dir, T, prior, lo, med, hi, obsT, obs, meta, obs_sd=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(T, lo, hi, color="tab:orange", alpha=0.25,
                    label="posterior predictive 90% CI")
    ax.plot(T, med, color="tab:orange", lw=2, label="posterior median")
    ax.plot(T, prior, color="tab:blue", lw=2, ls="--", label="emergent prior")
    if obs_sd is not None:
        ax.errorbar(obsT, obs, yerr=obs_sd, fmt="o", color="k", ms=5, capsize=3,
                    label=f"data ({meta['study']}, ±SD)")
    else:
        ax.plot(obsT, obs, "o", color="k", ms=6, label=f"data ({meta['study']})")
    ax.set_xlabel("Temperature (°C)"); ax.set_ylabel("Growth rate (1/h)")
    ax.set_title(f"Prior vs posterior — {meta['curve_id']} ({meta['medium_class']}, "
                 f"n={meta['n']})")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(out_dir, "prior_vs_posterior_tpc.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    return p


def _plot_corner(out_dir, flat, noise_key="sigma_noise"):
    import matplotlib
    matplotlib.use("Agg")
    import corner
    import matplotlib.pyplot as plt
    # show natural units for interpretability
    nat = np.column_stack([flat[:, 0], flat[:, 1], np.exp(flat[:, 2]),
                           np.exp(flat[:, 3]), np.exp(flat[:, 4])])
    labels = ["dTopt (K)", "dTm (K)", "dCp_scale", "kappa_scale", noise_key]
    fig = corner.corner(nat, labels=labels, show_titles=True,
                        title_fmt=".2f", quantiles=[0.05, 0.5, 0.95],
                        truths=[0, 0, 1, 1, None])
    p = os.path.join(out_dir, "corner.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p
