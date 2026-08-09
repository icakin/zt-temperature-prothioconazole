# Claude Code prompt — re-run the Bayesian calibration on the TRUSTED Noll NCM3722 minimal curve (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). This REDOES the Phase 0+1 Bayesian
calibration (previously run on the mis-scaled Cooper curve) on the trusted, strain-matched
Noll glucose-minimal curve. Additive inverse analysis; the emergent model is unchanged; NO
report edits (write-up deferred). Growth only.

NOTE TO USER: launch in an auto-approving mode so it runs unattended. Multi-hour (emcee).

WHY: the earlier calibration fit curve 39668, which we found was mis-scaled ~8x (peak 0.10 vs
true ~1.0 h-1), so its result (kappa_scale x0.18, "5x too efficient") was an artefact. We now
fit the trusted Noll (Katipoglu-Yazan et al. 2023) K-12 NCM3722 defined glucose-minimal curve,
which has per-temperature growth rate with SD and n wells. EXPECT the magnitude correction to
FLIP: the emergent peak (~0.55 h-1) UNDER-predicts Noll (~1.0-1.13 h-1), so the fit should now
demand kappa_scale > 1 (translation efficiency higher, not lower).

Trusted curve: strains/eciML1515/thermal/sources/noll2023_ncm3722/noll2023_ncm3722_tpc.csv
(rate_per_h + rate_sd_per_h + n_wells; 28-45 C have rates; 27 C is yield-only, no rate).

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
src/etcgem/{calibration.py (the existing emcee calibration: PARAM_NAMES, load_curve,
build_priors, the simulator wrapper, likelihood),validation.py (how the trusted Noll curve +
its SD are loaded, incl. the medium and the noll_minimal_AAopen variant),cli.py (calibrate +
validate commands),providers.py (set_medium),tpc.py}, and the Noll CSV above. Do NOT change the
emergent-model defaults or the report.

PART A - point the calibration at the trusted Noll curve (with its SD)
- Make the calibration load the Noll curve from the SOURCE CSV (reuse validation.py's trusted-
  curve loader / CurveSpec so calibration and validation read the SAME data), NOT from the
  retired ecoli_tpc_curves.csv. Carry temp_C, rate_per_h AND rate_sd_per_h (and n_wells).
- Fit temperatures = the rows with a rate (28-45 C; drop the 27 C yield-only row).
- Medium = glucose_minimal (Noll = defined glucose minimal + 6 trace amino acids). Predict under
  glucose_minimal for the headline; ALSO report a variant with the 6 amino-acid uptakes opened
  (reuse the existing noll_minimal_AAopen handling) as a sensitivity check.

PART B - use the measured SD in the likelihood (improvement over the single sigma_noise)
- Replace the single free sigma_noise with a per-point measurement variance plus a free MODEL-
  DISCREPANCY term: for each temperature i, likelihood variance = rate_sd_i^2 + sigma_disc^2,
  where sigma_disc >= 0 is a free nuisance (structural model error), prior HalfNormal. This
  weights each point by its measured precision and still absorbs model misfit. Keep the exact
  Gaussian likelihood on raw ABSOLUTE rate (1/h).
- Free parameter set (unchanged physical knobs): {kappa_scale, dCp_scale, dTopt, dTm} + sigma_disc.
  Keep the metabolic saturation sigma and the measured sector fractions FIXED at emergent values
  (sigma is degenerate with kappa_scale for magnitude; note this). Priors as already defined
  (centred on emergent values; dTopt ~ N(0,4.25 K), dTm ~ N(0,6.13 K), dCp/kappa LogNormal).

PART C - run emcee (longer than the proof-of-concept)
- Gradient-free ensemble MCMC, parallelised, initialised near the emergent point. Run LONGER
  than Phase-1 (which gave n_eff~157): aim for n_eff >~ 400 per parameter (more walkers and/or
  steps); report acceptance, autocorrelation time, n_eff. Keep the per-solve timeout / finite
  cost guards added earlier. Report wall-time.

PART D - outputs (NO report edits)
- Save under strains/eciML1515/outputs/calibration_noll_minimal/ (do NOT overwrite the stale
  calibration_phase1/; note it is superseded):
  * chain + summary.json (per-parameter posterior median/90% CI, acceptance, autocorr, n_eff,
    curve = Noll, medium, wall-time).
  * PRIOR-vs-POSTERIOR TPC on raw absolute rate: emergent/prior TPC, posterior-predictive band
    (median + 90% CI), and the Noll data points WITH SD error bars.
  * demanded_corrections.csv: per free parameter, prior vs posterior (median + 90% CI) and the
    interpretable correction (expect kappa_scale > 1 now; dTopt maybe + toward the ~40-42 C
    optimum; dTm to place CTmax; dCp poorly constrained). Flag any parameter the curve doesn't
    constrain.
  * corner plot (degeneracies; expect kappa<->envelope and dTopt<->dCp ridges).
  * console SUMMARY: does the posterior now reach the Noll peak? the demanded kappa_scale
    (and whether it hits a ceiling — if the metabolic pool binds before the peak is reached,
    that is a finding: note it and that freeing sigma/pool would be the next lever). State that
    the cold-side rising limb is UNCONSTRAINED (Noll starts at 27 C), so dCp/dTopt are set by
    the plateau + falling limb only.
- Do NOT edit reports/etcgem/*.

VERIFY (report all)
1. Calibration reads the trusted Noll curve (with SD); fit temps 28-45 C; medium glucose_minimal.
2. Likelihood uses rate_sd_i^2 + sigma_disc^2; free set {kappa_scale, dCp_scale, dTopt, dTm, sigma_disc}.
3. Posterior: per-parameter median + 90% CI and demanded corrections; confirm kappa_scale > 1
   (magnitude flipped vs the old artefact) or, if it hit the pool ceiling, say so.
4. Convergence: acceptance, autocorr, n_eff (>~400 target). Prior-vs-posterior TPC + corner exist.
5. Emergent model and report UNCHANGED.

CONSTRAINTS
- Additive inverse analysis only: do NOT change emergent-model defaults, core math, or the report.
- Exact likelihood + gradient-free sampler (emcee). Not ABC, not Stan/brms.
- Priors centred on emergent values; nothing tuned to flatter the fit; magnitude stays a
  posterior on kappa_scale, not a hand-set number.
- Fit RAW ABSOLUTE rates (1/h), single Noll minimal curve (Erdos/rich and joint/per-curve fits
  are later phases).
- Autonomous; commit in parts: "calibration: load trusted Noll curve + per-point SD likelihood",
  "emcee re-run on Noll minimal (longer chains)",
  "calibration_noll_minimal outputs: prior-vs-posterior TPC, corrections, corner".
```
