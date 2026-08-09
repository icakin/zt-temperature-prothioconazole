# Claude Code prompt — Bayesian calibration Phase 0 + 1: setup + single-curve proof of concept (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). This is an ADDITIVE inverse analysis: it
does NOT change the emergent model and does NOT touch the report (report write-up is deferred
until all phases are checked). Growth only.

NOTE TO USER: launch in an auto-approving mode so it runs unattended. This is a multi-hour run
(many model evaluations); parallelise and use sensible chain lengths.

GOAL: the emergent model is the a-priori prediction (the PRIOR). Here we ask the inverse
question on ONE well-characterised curve: what would the (borrowed/uncertain) parameters have to
be to reproduce the data, and how well does the curve constrain them? Deliverable = the first
"prior vs posterior" figure + a demanded-corrections table + a degeneracy (corner) plot. Nothing
here is relabelled as a prediction.

METHOD CHOICE (fixed): the etc-GEM is a DETERMINISTIC simulator (parameters -> TPC) and we can
write a Gaussian residual-noise model, so we do EXACT likelihood-based inference with a
GRADIENT-FREE sampler (emcee ensemble MCMC). NOT ABC (reserved for the high-dim per-enzyme Phase
4), NOT Stan/brms (the FBA simulator is non-differentiable).

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
src/etcgem/{enzyme_cost.py (Perturbation; the biosynthesis/translation cap kappa*v_bio <=
f_bio*P_total; _costs/_costs_unfolding),sectors.py,tpc.py (compute_tpc, TPC, descriptors),
providers.py,config.py,cli.py}, strains/eciML1515/strain.yaml, the per-curve validation code +
outputs under strains/eciML1515/outputs/percurve_validation/ (to find the encodable
glucose-minimal curves and their absolute rates), and the elasticity outputs
(elasticity_elasticity/) to confirm which knobs move which descriptor. Do NOT modify the report
or the emergent-model defaults.

PART A - free-parameter set + magnitude knob (Phase 0)
- The elasticity determines the free set (free only knobs that can move the fit targets):
  * MAGNITUDE / peak (rmax): rmax is translation-cap limited, so the metabolic pool (budget/
    sigma) can NOT move it — the lever is the biosynthesis/translation cap. Add a NEW knob
    `kappa_scale` to Perturbation: a multiplier on the translation coefficient kappa in the
    biosynthesis cap (kappa/kappa_scale * v_bio <= f_bio*P_total), so kappa_scale > 1 raises the
    achievable rmax. Wire it through the cap; verify a small +kappa_scale raises rmax and leaves
    the cold-side shape ~unchanged. (This is the in-vivo translation-efficiency parameter — a
    genuinely uncertain/borrowed quantity, the honest magnitude lever.)
  * SHAPE: dCp_scale (rising-limb curvature / Ea), dTopt (optimum), dTm (upper limit).
- Phase-1 free set = {kappa_scale, dCp_scale, dTopt, dTm} + a noise scale sigma_noise (nuisance).
  Keep the metabolic saturation sigma and the measured sector fractions FIXED at their emergent
  values in Phase 1 (sigma is degenerate with kappa_scale for magnitude; note this, explore in a
  later phase). Everything else stays at the emergent/measured value.
- Provide a thin SIMULATOR wrapper: theta -> Perturbation -> compute_tpc(pm, curve_temps, pert)
  -> predicted absolute growth rate (1/h) at the curve's measured temperatures. Evaluate ONLY at
  the curve's measured temperatures (not a dense grid) for speed.

PART B - priors + likelihood (Phase 0)
- Priors centred on the emergent values, widths from provenance (document every choice in a
  priors table / json):
  * dTopt ~ Normal(0, ~4 K)  (scale ~ the per-enzyme Topt SD, 4.25 K)
  * dTm   ~ Normal(0, ~6 K)  (scale ~ the per-enzyme Tm SD, 6.13 K; also ~ Tm-prediction error)
  * dCp_scale ~ LogNormal about 1 (e.g. sd ~0.3 on log), truncated positive
  * kappa_scale ~ LogNormal about 1, broad (translation efficiency uncertain, ~ up to ~2x)
  * sigma_noise ~ HalfNormal (or Exponential), weakly-informative from the data scatter
- Likelihood: Gaussian residuals on ABSOLUTE rate — obs_rate_i ~ Normal(model_rate(theta, T_i),
  sigma_noise). log-posterior = log-prior + sum_i log N(obs_i | model_i, sigma_noise).

PART C - the fit curve (Phase 0)
- From the per-curve validation, select ONE encodable glucose-minimal curve on RAW ABSOLUTE
  rates (1/h) with good temperature coverage (e.g. Bennett-Lenski 2001 or Mohr-Krawiec 1980 —
  pick the better-covered one; report which and why). Fit under the glucose-minimal medium at
  that curve's measured temperatures. Confirm the units are 1/h (convert doubling times if
  needed) before fitting.

PART D - inference (Phase 1)
- Use emcee (gradient-free ensemble MCMC). Sensible settings: e.g. >= 2*ndim walkers (>= ~16),
  a burn-in + sampling run long enough to converge, initialise walkers in a small ball around
  the emergent point (theta = 0 shifts, scales = 1). PARALLELISE the log-prob over walkers
  (multiprocessing / emcee pool). Record acceptance fraction and integrated autocorrelation time;
  thin/burn accordingly. If emcee is unavailable, a likelihood-tempered SMC is an acceptable
  substitute — but a gradient-free LIKELIHOOD sampler, not ABC.
- Keep the run bounded (report wall-time); this is a proof of concept, not the final long run.

PART E - outputs (NO report edits)
- Save under strains/eciML1515/outputs/calibration_phase1/:
  * posterior samples (chain) + a summary json (per-parameter posterior mean/median/90% CI,
    acceptance, autocorr time, n_eff, the chosen curve, wall-time).
  * PRIOR-vs-POSTERIOR TPC figure: the emergent/prior predicted TPC, the posterior-predictive
    band (median + 90% credible interval) from the posterior, and the empirical data points, all
    on RAW ABSOLUTE rate (1/h). This is the headline figure (triptych panels 1-2).
  * DEMANDED-CORRECTIONS table (csv): per free parameter, the prior (centre/width), the posterior
    (median + 90% CI), and the implied correction in interpretable units (e.g. dTm = +X K,
    kappa_scale = xY). Flag any parameter whose posterior ~ its prior (curve does not constrain
    it).
  * CORNER plot of the posterior (marginals + pairwise), to expose degeneracies (expect a
    kappa_scale<->shape or dCp<->dTopt trade-off).
  * A short console SUMMARY: does the posterior close the rmax / high-T gaps, what corrections it
    demands, and which parameters are well- vs poorly-constrained.
- Do NOT edit reports/etcgem/*. Report write-up is a later, separate step after all phases.

VERIFY (report all)
1. kappa_scale works (raising it raises rmax; cold-side ~unchanged); the free set and priors used.
2. The fitted curve id, its temperatures, unit check (1/h).
3. Posterior: per-parameter median + 90% CI and the demanded corrections; which parameters the
   curve constrains vs not; the main degeneracy seen.
4. Sampler convergence: acceptance fraction, autocorrelation time, effective sample size.
5. The prior-vs-posterior TPC figure and corner plot exist; the emergent model and the report are
   UNCHANGED.

CONSTRAINTS
- Additive inverse analysis only: do NOT change the emergent-model defaults, the core MMRT/
  unfolding/enzyme-cost math (beyond adding the kappa_scale knob), or the report.
- Exact likelihood + gradient-free sampler (emcee). Not ABC, not Stan/brms.
- Priors centred on the emergent values with documented, provenance-based widths; nothing tuned
  to make the fit look good.
- Fit RAW ABSOLUTE rates (1/h); single glucose-minimal curve only (per-curve and joint fits are
  later phases).
- New module (e.g. src/etcgem/calibration.py) + a CLI subcommand (e.g. `etcgem calibrate
  --strain eciML1515 --curve <id>`); config-driven priors/sampler settings. pip deps (emcee,
  corner) installed with --break-system-packages if needed.
- Autonomous; commit in parts: "add kappa_scale (translation-efficiency) knob",
  "calibration module: simulator wrapper + priors + Gaussian likelihood",
  "emcee single-curve calibration + CLI", "phase-1 outputs: prior-vs-posterior TPC, corrections
  table, corner plot".
```
