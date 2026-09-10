# v2 preregistration: transcript-informed, enzyme- and temperature-constrained GEM of *Z. tritici*

**Status: FROZEN — approved by I. Cakin, 17 August 2026, after one external
review round (revision 2). Amendments, if ever needed, are separate dated
addenda that never relax a threshold retroactively. The v1 manuscript and all
its outputs remain untouched.**

Successor to `docs/v2_model_design_charter.md`, instantiated on the recorded
Gate 0 outcome (`tables/revision/v2gate0/gate0_results.json`, branches **2 + 5**).
Where this document and the charter differ, this document governs.

Revision 2 (pre-freeze review): N1 redefined as a maximising null under μ ≤ μ_R
(equality-clamping prohibited); N0 redefined as a fold-fitted log-linear
growth→respiration map fed M1T's own growth predictions, defined at all
temperatures and independent of N1-cal; carbon uptake removed from the
sharpness metric; exchange bounds frozen with an audit and exchange-cap-driven
classification (§3.1); P_min biological-plausibility failure rule added.

---

## 0. Evidential framing (fixed)

> Because the untreated thermal performance curve and its optimum were
> characterised and published in v1, v2 analyses of untreated growth are
> design-informed consistency tests rather than genuinely held-out predictions.
> The primary predictive tests concern forward respiration and energetic state,
> respiratory branch allocation, transfer across temperatures without
> transcriptome measurements, and dose–temperature generalisation under
> parameters calibrated only in pre-specified reference conditions.

All Gate 0 quantities (sector shares, r, P_min, coverage) are architecture-design
evidence and are never presented as independent tests. A principled failure is an
acceptable primary outcome; failures are attributed only to *classes* of missing
mechanism.

## 1. Architecture as fixed by Gate 0

**Retained (branch 2):** the reference-normalised translation-capacity ceiling,

μ_c ≤ μ_R,c = μ_obs(21) · [φ_R,c / φ_R,21] · [k_R(T_c) / k_R(21)],

active **only at the transcript-measured temperatures 15, 21 and 27 °C**, using
the measured mass-weighted φ_R values recorded at Gate 0
(0.1244, 0.1337, 0.0816; r = 0.610 [0.469, 0.773]) with replicate uncertainty
propagated. k_R Q10 = 2.0 primary, {1.5, 3.0} scenario band, never fitted.

**Rejected (branch 5):** linear interpolation of φ_R. Consequently **no
translation ceiling exists at 18, 24, 26 or 28 °C** — those temperatures are
predicted by the network alone (pool, kcat(T), maintenance, stoichiometry). No
nonlinear φ_R trajectory may be fitted as rescue.

φ_met is measured flat across 15/21/27 °C (0.1484, 0.1489, 0.1500); at
unmeasured temperatures it is carried as the linear interpolation/extension of
these values (in practice near-constant ≈ 0.149). This is permitted because the
measured trajectory is flat; it is a declared assumption, listed in section 9.

**Pre-stated structural expectations (written before any model exists, so that
observing them is not reinterpreted later):**

1. The primary model **cannot generate the 24 °C peak** — a direct consequence
   of branch 5 plus the flat DLTKcat ensemble and flat φ_met. Its predicted
   untreated TPC will be: suppressed at 15 °C (ceiling), roughly flat across
   18–26 °C, suppressed at 27 °C (ceiling), and roughly flat again at 28 °C.
2. It will therefore **overpredict growth at 28 °C** (observed 0.0396 h⁻¹
   against a plateau near the 21 °C level), and this overprediction is the
   pre-registered signature of a missing supra-optimal mechanism outside the
   measured sectors and predicted kinetics.
3. It will **underpredict respiration and overpredict CUE** at supra-optimal
   temperatures (growth-optimal allocation carries no wasteful ATP demand).

If any of these expectations is *not* observed, that is a reportable surprise in
its own right.

## 2. Model ladder (unchanged from charter)

M0 (thermal ecModel, shared pool) → **M1T (primary: + sector allocation +
ceiling at measured temperatures)** → M2 (+ drug) → M3 (+ selected relative
transcript bounds, β ∈ {0, 0.5, 1}, AOX always excluded) → M4 (+ turnover cost,
identifiability test only).

**The translation relationship is an upper bound everywhere it appears
(μ ≤ μ_R); it is never imposed as an equality or growth target.** The only
equality on growth anywhere in this design is N1-cal's clamp to *observed*
training growth during maintenance fitting. This wording governs all
implementation.

Nulls:

- **N1** — loose-pool stoichiometric null that *maximises* growth subject to
  μ ≤ μ_R at 15, 21 and 27 °C, and maximises unconstrained network growth at
  temperatures without a ceiling; same frozen hierarchy and FVA as M1T. (An
  equality-clamped N1 is prohibited: it would force the ceiling to bind
  whenever feasible and make the slackness falsification test impossible by
  construction.)
- **N0** — empirical growth-only respiration map: log R = α + γ·log μ, fitted
  on observed untreated training pairs within the same outer series-block
  folds. **Functional form frozen now: log-linear.** At each held-out
  condition and each P_eff, N0 is fed **M1T's own growth prediction**:
  R̂_N0,c(P_eff) = f_train(μ̂_M1T,c(P_eff)). N0 therefore exists at every
  temperature, including the four network-only ones, and shares M1T's growth
  residuals exactly — so the E_phys difference between them isolates whether
  the metabolic network predicts *respiration* better than a bare empirical
  growth–respiration relationship. N0 is a purely empirical comparator and
  does not use the N1-cal maintenance law.

Each rung must beat the rung below on held-out conditions or is not retained.

## 3. Calibration (all of it)

- **P_eff:** profiled, never point-fitted, over [P_min, 0.25 g gDW⁻¹]. P_min is
  recomputed at build time under a frozen rule: the smallest pool supporting
  observed 21 °C growth with DLTKcat kcats evaluated at 21 °C (the Gate 0 value
  0.0206 pool units, computed under uniform baseline kcats against a built
  bound of 0.09, is an order-of-magnitude screen only). All results reported as
  functions of P_eff across the full profile. **Failure rule, fixed now:** if
  the recomputed P_min exceeds 0.25 g gDW⁻¹, or observed 21 °C growth is
  infeasible throughout the declared literature range, the metabolic-pool
  architecture fails its biological-plausibility gate and is reported as such.
  The range is not expanded after inspection.
- **Maintenance (N1-cal):** v_maint = a_GAM·μ + b_NGAM, both ≥ 0, no
  temperature or drug dependence. Fitted within outer temperature-block folds,
  untreated conditions only, through N1-cal — the loose-pool stoichiometric
  network with **observed** training growth imposed as equality, no ceiling
  involved — using paired series-block bootstrap draws of (μ, R). The frozen
  pair transfers unchanged to N0, N1, every P_eff point of M1T, and M2–M4.
  Pre-registered adequacy diagnostic: correlation of training residuals with
  temperature; a monotonic pattern is reported into the structural-necessity
  map and never repaired with a temperature term. Secondary directional test:
  fit on sub-optimal untreated conditions only, predict supra-optimal
  respiration/CUE without refitting.
- **Drug (K_app):** one apparent-potency constant, occupancy 1/(1 + D/K_app),
  Hill coefficient 1, calibrated on the **normalised 24 °C dose–response,
  growth only** (v1 EC₅₀ there: 1.40 mg L⁻¹). Explicitly phenomenological.
  CYP51 kcat(T): three fixed scenarios — constant (primary), DLTKcat
  trajectory, generic Arrhenius from the broad enzyme set — reported jointly,
  never selected by fit.
- Frozen and never fitted: v1 sector gene lists (map03010 ribosomal; ecZtGEM
  draw-enzyme metabolic), mass weighting, Q10 band, biomass composition,
  respiratory stoichiometry, the lexicographic hierarchy (max growth → min
  enzyme use → FVA), optimizer tolerances.

### 3.1 Exchange bounds (frozen; audit mandatory)

At network-only temperatures, M1T growth must emerge from the pool, kcat(T),
maintenance and stoichiometry — and N1's "loose pool" is interpretable only if
its growth is finite for declared reasons. Inherited exchange caps could
silently generate either. Therefore:

- The medium is the v1 sucrose-minimal set. **Every condition-dependent
  exchange bound in force is enumerated in the build log at first execution,
  with its numeric value and provenance** (v1 model default, medium
  definition, or physical requirement). Substrate uptake was never measured,
  so any finite carbon-uptake cap is a **model prior**, and is declared as
  such.
- **No bound derived from measured respiration, measured growth, or any
  inferred carbon uptake enters any forward prediction.** (v1's
  condition-constrained oxygen bounds are specifically prohibited here.)
  Oxygen, water, ions and phosphate are unconstrained except by
  stoichiometry.
- **Dual/activity values are reported for every potentially growth-limiting
  exchange, per condition, per P_eff**, alongside the pool and ceiling duals
  in the constraint-activity map.
- **Classification rule:** any prediction in which a finite exchange cap
  carries a nonzero dual at the growth optimum is classified
  **exchange-cap-driven**, reported separately, and excluded from claims
  about pool-, kinetics- or ceiling-generated behaviour. As a diagnostic
  (non-voting), the carbon cap is halved and doubled to show which
  conclusions are cap-sensitive.
- If N1's unconstrained-growth solution at any temperature is unbounded or
  determined solely by a default cap, N1 is reported as non-neutral at that
  condition and the surplus comparison there falls back to N0 alone.

## 4. Metrics (frozen)

- **E_μ:** log-scale RMSE on untreated growth, 21 °C excluded, reported in two
  pre-registered subsets that are never pooled: ceiling temperatures {15, 27}
  and network-only temperatures {18, 24, 26, 28}. Hard admissibility gates:
  correct sign of suppression at 15 °C and at 27 °C relative to 21 °C; no
  observed growth above the ceiling at a measured temperature (violation =
  charter Outcome C for the proxy). The 24 °C peak is **reported against the
  pre-stated expectation of failure** (section 1) and is not an admissibility
  gate.
- **E_phys:** stacked residual vector (Δlog μ, Δlog R) over held-out
  conditions, scored with the full series-block bootstrap covariance
  (μ and R share a DO trace; cross-condition dependence preserved), shrinkage
  regularisation and dimension normalisation frozen at implementation and
  recorded in the build log before any real scoring. CUE is derived only:
  report trajectories, the signed supra-optimal bias B_CUE,hot, and
  dose–temperature bias surfaces. Effective P/O likewise diagnostic only.
- **ΔE_net(P_eff):** E_phys(best null) − E_phys(M1T) across the entire
  profile, series-block bootstrap, Stage-1 materiality convention (0.25 × the
  bootstrap noise scale of the target quantity; win threshold 0.90). A surplus
  claim requires materiality over a nontrivial portion of the profile, not
  only near P_min.
- **Sharpness W_y:** median FVA-interval width / observed cross-condition
  range, for **respiration, CUE and AOX share only** — carbon uptake was never
  measured, so it has no observed range to normalise against. Carbon uptake is
  reported as **unvalidated latent-output FVA** in physical units (or relative
  to its own predicted midpoint) and is never called sharpness against
  observations. The CUE feasible interval is a transformed joint interval of
  the underlying growth and respiration fluxes, not an independently validated
  measurement interval, and is labelled so. A condition whose feasible
  interval spans most of the observed range is "compatible", never
  "predicted". Report W_y for M1T against N1; a strong claim requires
  reduction.
- **Constraint-activity map:** utilisation and duals for the ceiling, pool,
  maintenance and key respiratory constraints, per condition, per P_eff — a
  primary result, since a mechanism that never binds has explained nothing.

## 5. Validation layers

1. **Untreated thermal physiology** (design-informed): E_μ subsets, gates,
   E_phys, B_CUE,hot, activity map. 21 °C growth excluded as anchor.
2. **Six RNA conditions:** ceilings and binding at measured temperatures;
   energetic-state errors; surplus vs N0/N1; the two drug runs at 2 mg L⁻¹
   (direct inhibition + control transcriptome vs + drug transcriptome), whose
   difference is the model's estimate of transcript-associated adaptation; AOX
   transcript rank-agreement as secondary check (AOX never constrained).
3. **Network-only temperatures {18, 24, 26, 28}:** the strongest genuine test
   surface. Reported separately, never pooled with ceiling temperatures.
   Includes the pre-registered 28 °C overprediction check.
4. **Drug generalisation:** fixed K_app and untreated-trained maintenance
   transferred to all other doses and temperatures; growth, respiration, CUE.
   Interpretation bound by the v1 shift-versus-scale C2 result: transfer
   failure does not discriminate exposure from amplitude mechanisms; transfer
   success shows sufficiency of apparent potency, not uniqueness of CYP51.
5. **Negative controls:** transcript permutation among abundance-matched
   enzymes; temperature-label shuffling; random constrained-enzyme sets
   (M3 only); constant-sector replacement; temperature-only empirical
   interpolation of the physiology; M3 vs M2.

## 6. Pre-registered falsification

- **Ceiling:** observed growth above the ceiling (with uncertainty) at
  15/21/27 °C → the transcript-to-active-ribosome proxy fails. Ceiling slack
  (utilisation < 0.95 or zero dual) at 15 °C → translation is not the
  cold-limb mechanism; at 27 °C → the falling limb reverts to "unexplained"
  (φ_met is flat, so no sector fallback exists — pre-acknowledged).
- **Network surplus:** μ̂ ≈ μ_R at ceiling temperatures AND pool/enzyme
  constraints immaterial to growth AND no material ΔE_net over the profile →
  the ecModel adds nothing beyond the ceiling; report as such.
- **M3:** constraints rarely bind, improvement only at β = 1, or permuted
  transcripts match → transcript level adds nothing beyond sectors.
- **M4:** expected non-identifiable versus the matched generic maintenance
  null; reported as the class-level statement "at least one unmeasured
  temperature-dependent energetic process is required within the tested model
  class."
- **Whole framework:** feasible ranges spanning the observed range
  (permissive), or point accuracy contingent on the secondary objective.
- Retention always requires: beating the immediate nested null on held-out
  conditions, correct qualitative responses, acceptable uncertainty
  calibration, and stability across the fixed scenario ensembles. Small RMSE
  drops retain nothing. No threshold, gate, subset or scenario may change
  after real scoring begins.

## 7. Scenario discipline

One primary configuration: ceiling at measured temperatures, Q10 = 2, constant
CYP51 kcat, β = 0, v1 sector lists, mass-weighted shares, K_app from normalised
24 °C. Factorised challenges only (thermal Q10 band; sector definitions and
weighting; three CYP51 scenarios; β ensemble after M3 activation; CYP51
scenario × temperature). No full grid.

## 8. Endpoints (all publishable, fixed now)

**Positive:** the ceiling explains the measured-temperature thermal shape and
the network adds material held-out energetic and branch predictions across the
profile. **Partial:** the ceiling explains growth at measured temperatures; the
network adds no surplus. **Structural:** the pre-stated expectations of section
1 hold and the energetic collapse and/or drug response demonstrably require
unmeasured processes — delivered as the structural-necessity map naming the
missing measurement classes (intracellular drug/desthio, proteomics, direct
turnover assays, a denser thermal transcriptome around the optimum). Gate 0
already guarantees the last item on that list is non-empty: three RNA
temperatures cannot resolve the flat-then-cliff ribosomal trajectory that the
24 °C peak requires.

## 9. Declared assumptions (non-exhaustive list frozen here)

Transcript share proxies proteome sector share (sector level only);
mass-weighting via longest-isoform protein MW; k_R temperature dependence is a
generic biosynthetic Q10; φ_met near-constant extension; uniform baseline kcat
where DLTKcat is missing (21.4% of metabolic-sector mass — median-curve
fallback); nominal external dose proxies exposure; biomass composition
temperature-invariant; single isolate.

## 10. Build plan and record

Scripts 58+ in the repo, outputs under `tables/revision/v2build/`. Input hashes
recorded at first execution; this document is committed before that execution
and not edited afterwards (amendments, if ever needed, are separate dated
addenda that never relax a threshold retroactively). Estimated effort:
N1-cal + M0/N1/M1T with profile and folds ≈ 1–2 weeks; M2 and challenges ≈ 1
week; M3/M4 conditional on earlier rungs.

## 11. Sign-off

Approved to freeze and begin the build: ______________________  Date: __________
