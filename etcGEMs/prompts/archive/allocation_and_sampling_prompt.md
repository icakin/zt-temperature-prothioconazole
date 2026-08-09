# Claude Code prompt — proteome-sector allocation + correlated/calibrated thermal sampling (M1.2)

Self-contained. Run from the project root (`.../MICROADAPT/etcGEMs`). Two related
additions that make the allocation SPACE and the parameter SAMPLING realistic, so
the sweep/decomposition become calibrated uncertainty rather than a nominal scan.
Both must be BACKWARD-COMPATIBLE (default off -> current behaviour and numbers
unchanged). Implement as TWO separate commits: (1) proteome sectors, (2) thermal
sampling.

The project is already on the `src/etcgem` layout with the defaults/strain/
experiment config model and a two-tier CLI. Concretely, the current state you are
extending is:
- src/etcgem/enzyme_cost.py : `Perturbation` (fields dTopt, topt_scale, dCp_scale,
  budget, group_alloc) and `EnzymeConstrainedModel` (pool constraint, set_budget,
  refresh_params, per-entry Topt/dCp/base_cost).
- src/etcgem/sensitivity.py : `_lhs`, `_make_perturbation`, `run_sensitivity`.
- src/etcgem/decomposition.py : `run_decomposition(pm, temps_C, allocation_ranges,
  envelope_ranges, ...)` — allocation params today are budget_scale + `alloc_<grp>`;
  envelope params are dTopt/topt_scale/dCp_scale, applied via `_make_perturbation`.
- src/etcgem/dltkcat.py : `fit_mmrt`, `fit_predictions`, `apply_fits_to_provider`
  (the mutate-single-entries-then-refresh_params pattern — reuse it).
- src/etcgem/cli.py : two tiers — build/tpc/fba/dltkcat (strain only), sweep/
  decompose/control (strain + experiment). config.py resolves defaults <- strain
  <- experiment and injects provider.model_path + output_dir.
- Config homes: organism biophysics -> strains/<name>/strain.yaml ; method/ranges
  -> experiments/<name>.yaml ; universal method defaults -> defaults.yaml.
Read these before changing anything. Do NOT alter the MMRT / core cost math or the
ANOVA/Shapley logic in decomposition.py.

=====================================================================
PART 1 - PROTEOME-SECTOR ALLOCATION (mechanistic allocation axis)
=====================================================================
Replace the single scalar pool (+ coarse group multipliers) with a coarse-grained
proteome partition after Basan 2015 / Scott 2010 (proposal refs 6-7): total
proteome mass fraction P_total splits into three sectors summing to 1:
  f_metab : metabolic enzymes  -> the existing pool budget = f_metab * P_total
  f_bio   : biosynthesis/ribosomes -> a translation cap on growth
  f_maint : maintenance/housekeeping -> proteome overhead + maintenance ATP (NGAM)
This is the natural coordinate for WP2's "maintenance -> biosynthesis" reallocation
and makes B0 (H1.2) an explicit allocation trade-off with an interior optimum.

Model layer (new file src/etcgem/sectors.py, or extend enzyme_cost/providers):
- add_proteome_sectors(pm, cfg): given P_total and nominal fractions, wire on the
  existing cobra model:
  * metabolic pool: set the existing pool constraint ub = f_metab * P_total.
  * biosynthesis cap: add a constraint  translation_coeff * v_biomass <=
    f_bio * P_total  (Scott/Basan growth law). AUTO-CALIBRATE translation_coeff at
    build time so that at the nominal split and T0 the translation cap and the
    metabolic pool are co-limiting (both just bind); report the calibrated value.
  * maintenance: available proteome for metab+bio = (1 - f_maint) * P_total; and if
    a maintenance-ATP reaction exists (detect ATPM / NGAM), set its lower bound
    proportional to f_maint / f_maint_nominal. Higher maintenance -> less proteome
    and more ATP burned -> lower growth, more respiration, lower CUE.
- EnzymeConstrainedModel.set_allocation(f_metab, f_maint): update the three bounds
  in place (f_bio = 1 - f_metab - f_maint); validate the simplex (>=0, sum<=1).
  Parallels set_budget.
- Defaults: back out P_total from the model's existing pool bound given a nominal
  f_metab; ship f_metab~0.5, f_maint~0.15, f_bio~0.35 and a growth-law kappa ~4-5/h
  fallback if auto-calibration is off. Document these are calibratable vs growth-
  rate proteomics (refs 6-7).

Config + integration:
- strain.yaml (organism): add
    proteome_sectors: {enabled: false, P_total: null, f_metab: 0.5, f_maint: 0.15,
      atpm_reaction: null, translation_coeff: auto}
  When enabled, config.build_provider (or the CLI) calls add_proteome_sectors after
  building the provider.
- Perturbation: add optional fields f_metab, f_maint (default None). In compute_tpc,
  when they are set, call ec.set_allocation(...) INSTEAD of set_budget; when None,
  behaviour is exactly as today.
- Extend sensitivity._make_perturbation so an allocation sample may contain f_metab /
  f_maint / a single `maint_to_bio` scalar (shift mass maintenance->biosynthesis at
  fixed f_metab) as well as the existing budget_scale / alloc_<grp>. run_sensitivity
  and run_decomposition then accept these as allocation params with NO other change
  (the ANOVA/Shapley code is untouched). Example experiment allocation block when
  sectors are on:  allocation_params: {f_metab: [0.4, 0.6], f_maint: [0.05, 0.25]}

=====================================================================
PART 2 - CORRELATED / CALIBRATED THERMAL SAMPLING (envelope uncertainty)
=====================================================================
Today the envelope is perturbed by two global knobs (dTopt shift, topt_scale). Add
a per-enzyme thermal-parameter SAMPLER with a realistic correlation structure and,
optionally, real DLTKcat posteriors -> calibrated uncertainty (deliverable M1.2).

One-factor sampler (shared organism thermal regime + per-enzyme idiosyncrasy). For
each ensemble member draw shared Z, Z2 ~ N(0,1); per enzyme eps_i, eps2_i ~ N(0,1):
  Topt_i = mean_Topt_i + sd_Topt_i * ( sqrt(rho)*Z  + sqrt(1-rho)*eps_i )
  dCp_i  = mean_dCp_i  + sd_dCp_i  * ( sqrt(rho)*Z2 + sqrt(1-rho)*eps2_i )
rho = shared_fraction in [0,1]: rho=1 -> a coherent whole-proteome shift (like the
current global dTopt); rho=0 -> independent per-enzyme draws (unrealistic; explores
mismatched-optima space). Captures that thermostability co-varies across a genome.
- Modes: "knobs" (current behaviour, default), "correlated" (mean = nominal per-
  enzyme Topt/dCp; sd from config), "posterior" (mean+sd per enzyme from DLTKcat
  fits; enzymes without a fit fall back to correlated defaults).
- New file src/etcgem/thermal_sampling.py:
    sample_thermal(pm, n, mode, rho, sd_cfg/posterior, seed) -> yields per-sample
    per-enzyme (Topt array, dCp array);
    apply_thermal_sample(pm, sample) mutates the entries and calls refresh_params
    (reuse the dltkcat.apply_fits_to_provider pattern).
- Integration: give run_sensitivity and run_decomposition an optional envelope
  sampler. When an `envelope_sampling` block is present, the ENVELOPE samples are
  drawn per-enzyme by sample_thermal and applied via apply_thermal_sample before
  compute_tpc (instead of the dTopt/topt_scale/dCp_scale knobs); the crossed design,
  descriptor extraction and ANOVA/Shapley math are UNCHANGED (each envelope sample j
  is now a per-enzyme vector rather than a 3-knob point). When absent, identical to
  now. Restore nominal params after each evaluation.

Predictor posteriors from DLTKcat (make the fit uncertainty real):
- fit_mmrt is linear least squares in (dH, dS, dCp), so it has a parameter
  covariance Cov = sigma^2 * (X^T X)^-1 (sigma^2 = residual variance). Extend
  fit_mmrt/fit_predictions to also return Cov and per-enzyme sd_Topt, sd_dCp got by
  sampling (dH,dS,dCp) ~ N(coef, Cov) and computing Topt/dCp per draw (Topt is
  nonlinear in dH,dCp -> sample, don't linearise). FLOOR sigma using DLTKcat's
  global skill (R^2 ~ 0.6, log10 RMSE ~ 0.9) so optimistic residuals don't
  understate uncertainty. Write Topt_sd, dCp_sd columns into fits.csv. The
  "posterior" sampler reads strains/<name>/dltkcat/fits.csv.

Config (method overlay -> experiment yaml, e.g. added to default.yaml /
decomposition.yaml):
  envelope_sampling:
    mode: correlated            # knobs | correlated | posterior
    shared_fraction: 0.7        # rho
    topt_sd_K: 4.0              # default per-enzyme sd (correlated / fallback)
    dcp_sd_frac: 0.3
    posterior_from: dltkcat     # posterior mode -> strains/<name>/dltkcat/fits.csv

VERIFY (do these; report)
1. HARD GATE - backward compatibility: with proteome_sectors.enabled=false and no
   envelope_sampling block, a toy and an eciML1515 smoke run of `sweep` and
   `decompose` reproduce the current descriptors within 1e-9 of a run on the
   pre-change code. Do not proceed if this fails.
2. Sectors (toy + eciML1515): sweeping f_metab<->f_bio at fixed f_maint shows an
   INTERIOR growth optimum (co-limitation); raising f_maint lowers growth. Confirm
   set_allocation validates the simplex and add_proteome_sectors reports the
   auto-calibrated translation_coeff.
3. Thermal sampler: rho=1 (correlated) reproduces a coherent global Topt shift
   (organismal Topt tracks the shared factor ~ dTopt); rho=0 decorrelates per-enzyme
   optima. Report how the organismal-TPC ensemble and the enzyme-parameter space
   differ between rho=0 and rho=1.
4. Posterior mode: fit_mmrt returns Topt_sd/dCp_sd; on synthetic noisy kcat(T) the
   sds shrink with more points and grow with noise, and the R^2 floor is applied.
   Run a small posterior-mode envelope ensemble on eciML1515 using the existing
   strains/eciML1515/dltkcat/fits.csv and confirm it writes.
5. Update README + docs/RUNBOOK with an "Allocation sectors" and a "Calibrated
   thermal sampling (M1.2)" section; add example experiment yamls
   (experiments/sectors.yaml, experiments/calibrated.yaml) demonstrating each.

CAVEATS TO DOCUMENT
- Sector growth-law/maintenance couplings (translation_coeff, maintenance ATP,
  P_total) are calibratable vs growth-rate proteomics; defaults are order-of-
  magnitude and auto-calibration co-limits at the nominal point only.
- The one-factor correlation is a deliberately simple stand-in for the true
  (phylogenetic/structural) covariance of thermostability; rho is the single knob.
- Posterior sd from a linear-LS fit is a local Gaussian approximation floored by
  DLTKcat's global skill, not a full Bayesian posterior.

CONSTRAINTS
- Backward compatible: disabled -> identical numbers (item 1 is a gate).
- Reuse Perturbation/set_budget/set_allocation/refresh_params, compute_tpc,
  TPC.descriptors, sensitivity._lhs/_make_perturbation, the decomposition ANOVA
  code, and the apply_fits_to_provider mutate-then-refresh pattern.
- Two commits: "add proteome-sector allocation (Basan/Scott)" and "add correlated +
  DLTKcat-posterior thermal sampling (M1.2)".
```
