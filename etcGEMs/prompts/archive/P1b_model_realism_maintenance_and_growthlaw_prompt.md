# Claude Code prompt — P1b: model-realism fixes — restore T-dependent maintenance under sectors + growth-rate-dependent biosynthesis cap + new envelope/maintenance knobs (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). This fixes two diagnosed problems and adds
the knobs the re-tune will use. It does NOT re-run the Bayesian (P2) or the P3 analyses — it
prepares the model, verifies on the EMERGENT (untuned) curve, and stops. Growth only.

NOTE TO USER: launch in an auto-approving mode.

DIAGNOSIS (from comparing our model to the Li yeast + MRes E. coli etcpy on GitHub):
1. FLAT TOP: our T-dependent maintenance (NGAM(T)) is switched OFF whenever the sector model is
   wired (`... and self._sectors is None` in enzyme_cost.set_temperature). Li/MRes call
   `set_NGAMT` at EVERY temperature — maintenance always rises with T, which rounds the peak.
   Restore that under the sector model.
2. MAGNITUDE SHORTFALL: Li/MRes use a SINGLE undivided GECKO pool with NO translation/biosynthesis
   cap. We added the proteome sector model, whose biosynthesis cap (translation_coeff*v_bio <=
   f_bio*P_total) BINDS at the rich peak. Our f_bio is ~static (measured), but the bacterial growth
   law says the ribosome/biosynthesis fraction rises with growth rate. Make f_bio growth-rate-
   dependent so the cap is self-consistent and can reach fast rich rates.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
src/etcgem/{enzyme_cost.py (Perturbation; set_temperature — the sector-gated NGAM block ~L328;
_costs_unfolding; kcat_scale/kappa_scale/dTm),sectors.py (the biosynthesis cap
translation_coeff*v_bio <= f_bio*P_total; f_bio=1-f_metab-f_maint(-f_chaperone); ATPM/maintenance
handling),unfolding.py (ngam_T, native_fraction),providers.py,config.py,cli.py (anatomy +
validation)}, strains/eciML1515/strain.yaml, and the proteome/sector data under
strains/eciML1515/{proteomics,outputs/proteome_sectors}. Do NOT run the Bayesian calibration or
the sensitivity/decomposition/identifiability.

PART A - restore T-dependent maintenance UNDER the sector model (round the peak)
- In set_temperature (or the sector maintenance path), apply the ngam_T(T) temperature factor even
  when sectors are wired: scale the sector-owned ATPM lower bound by ngam_T(T)/ngam_T(T0) so
  maintenance rises with temperature under sectors (combined with the existing f_maint sector
  scaling). This mirrors Li/MRes `set_NGAMT` being called at every T. Remove the
  `and self._sectors is None` gate (or add the sector branch).
- Add two maintenance knobs to Perturbation (no-ops at 1, for the unified set / re-tune):
  * ngam_scale: multiplier on the maintenance amplitude (the ~8.5 baseline);
  * ngam_steepness: multiplier on the temperature-dependence (the activation term in ngam_T),
    i.e. how fast maintenance rises with T (the peak-rounding lever).
  Centre both at the ported Li getNGAMT values.
- VERIFY on the EMERGENT model at the rich (BHI) operating point: the peak is now ROUNDED — a
  proper interior optimum with a gradual pre-collapse decline — not the flat plateau.

PART B - growth-rate-dependent proteome partition (the COUPLED f_bio/f_metab growth law)
- The caps are currently translation_coeff (kappa)*v_bio <= f_bio*P_total (biosynthesis) and
  Sum cost_i*|v_i| <= f_metab*B (metabolic pool; B = the current pool RHS scaling, P_total*sigma or
  P_total — use whatever the code uses), with f_bio/f_metab ~fixed. Add a CONFIG TOGGLE (e.g.
  biosynthesis_growth_law: true/false; default keeps the current FIXED partition for comparison)
  that makes the partition growth-rate-dependent AND PROTEOME-CONSERVING (Scott/Hwa):
        f_bio(mu)   = f_bio_0   + slope*mu
        f_metab(mu) = f_metab_0 - slope*mu     (the SAME slope: proteome moves metab -> ribosomes)
  so f_bio + f_metab is conserved. This is the CORRECT growth law; do NOT implement the
  f_metab-fixed version (it violates proteome conservation and over-relaxes).
- Both caps stay LINEAR under v_bio = mu (no iteration):
        biosynthesis:  (kappa - slope*P_total) * v_bio         <= f_bio_0   * P_total
        metabolic:     Sum cost_i*|v_i| + (slope*B) * v_bio     <= f_metab_0 * B
  The extra `+ slope*B*v_bio` term in the metabolic constraint is the coupling: as growth rises,
  proteome is reserved for ribosomes and the metabolic budget tightens. Implement by adjusting the
  v_bio coefficient in the biosynthesis constraint and adding the v_bio term to the metabolic pool
  constraint. GUARD: kappa - slope*P_total > 0 (slope < kappa/P_total); clamp/warn otherwise.
- This self-limits growth correctly: at high mu the shrinking metabolic sector binds and sets the
  maximal rate — the ribosome<->metabolism trade-off, which is the point.
- Ground f_bio_0, f_metab_0 and slope in data: fit the biosynthesis and metabolic proteome
  fractions vs growth rate across the DeyuWang medium/temperature conditions (and/or the E. coli
  growth-law literature). Use a SINGLE consistent slope (f_bio slope = -f_metab slope). Store as
  config values; cite the source.
- VERIFY on the EMERGENT model at the rich operating point: with the growth-law toggle ON, the
  achievable rmax RISES toward the observed ~2.4 (report vs the fixed-partition ~1.03 and observed
  2.40), and report WHICH constraint now binds at the peak (expect the coupled metabolic sector) —
  i.e. that the trade-off, not a static cap, sets the maximum.

PART C - tm_scale knob (Tm spread; the falling-shoulder lever)
- Add tm_scale to Perturbation: multiplies each enzyme's (Tm - mean_Tm) about the distribution
  mean (mirroring topt_scale for optima), controlling how synchronised the unfolding collapse is
  (narrow => sharp cliff; broad => gradual shoulder). No-op at 1.

PART D - regenerate emergent diagnostics (no re-tune)
- Regenerate the emergent reference TPC and the emergent-vs-Van-Derlinden validation showing
  (i) the rounded peak (maintenance restored) and (ii) the growth-law f_bio lifting the rich
  rmax. Report the new emergent descriptors (Topt, rmax, CTmax, Ea) under both f_bio modes.
- Keep NGAM-under-sectors ON by default; keep the growth-law f_bio as a toggle (report both).
  Do NOT re-run the Bayesian, sensitivity, decomposition, or identifiability.

VERIFY (report all)
1. NGAM(T) active under sectors; the emergent rich peak is rounded (not flat). ngam_scale,
   ngam_steepness, tm_scale added (no-ops at 1).
2. Coupled growth-law toggle: proteome-conserving f_bio(mu)=f_bio_0+slope*mu, f_metab(mu)=
   f_metab_0-slope*mu; both linear (biosynthesis coeff kappa-slope*P_total guarded >0; metabolic
   pool gains +slope*B*v_bio). f_bio_0/f_metab_0/slope values + source. With it ON: emergent rich
   rmax vs fixed (~1.03) vs observed (2.40), and which constraint binds at the peak (expect the
   coupled metabolic sector — the trade-off setting the max).
3. The unified free set ready for the re-tune is: {dTopt, topt_scale, dCp_scale, dTm, tm_scale,
   kcat_scale, kappa_scale, f_metab, f_maint, ngam_scale, ngam_steepness} — plus the growth-law
   slope available as a config value (optionally a broad-prior knob later; NOTE it is degenerate
   with kappa_scale via kappa_eff, so free only one).
4. Emergent diagnostics regenerated; Bayesian/sensitivity/decomposition NOT run.

CONSTRAINTS
- Model-realism fixes + new knobs only. Do NOT re-run the Bayesian, sensitivity, decomposition, or
  identifiability (those are the P2 re-run + P3, separate).
- Keep toggles: NGAM-under-sectors (default ON) and the coupled growth-law partition (default
  OFF/comparison) so we can compare against the current fixed-partition behaviour.
- The growth law MUST be proteome-conserving (f_bio up <=> f_metab down, one shared slope); do NOT
  ship the f_metab-fixed version.
- Match the maintenance restoration to the Li/MRes `set_NGAMT`/`getNGAMT` behaviour we ported.
- Ground the growth-law slope/offsets in data or literature and cite it.
- Autonomous; commit in parts: "restore T-dependent NGAM under sectors + ngam_scale/ngam_steepness
  knobs", "coupled growth-law partition (f_bio up / f_metab down, linear, toggle)",
  "add tm_scale knob", "regenerate emergent diagnostics (rounded peak + growth-law magnitude)".
```
