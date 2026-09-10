# v2 preregistration — Amendment A1 (dated addendum)

**Status: draft for sign-off (I. Cakin). Executes nothing until signed. This
amendment supplements `docs/v2_preregistration.md` (FROZEN); it does not edit
that document, and it relaxes no threshold retroactively — the failed
plausibility gate stands as a committed result. Where this amendment and the
prereg conflict, this amendment governs v2.1 analyses only.**

Date: 17 August 2026. Motivated by the pre-registered plausibility-gate failure
recorded in `tables/revision/v2build/` (scripts 58–59), reviewed adversarially
by ChatGPT; this text implements that review's requirements.

---

## A0. What v2.1 is

The pre-registered absolute-proteome architecture failed its
biological-plausibility gate: at the maximum declared protein-pool bound, the
raw DLTKcat ecModel attained ~31% of observed 21 °C growth (P_min = 3.47 g gDW⁻¹
vs the 0.25 bound). Because a global turnover multiplier and the pool size are
algebraically non-identifiable (only their product matters), v2.1 **abandons
absolute pool inference** and calibrates **one effective-capacity scalar** at
untreated 21 °C. v2.1 is a **conditional effective-capacity transfer
experiment**, not a repaired absolute ecModel:

> Given one globally calibrated effective catalytic-capacity scale at 21 °C, do
> the predicted relative temperature responses, allocations, respiration, CUE
> and drug effects transfer?

## A1. Effective capacity C_eff★

- k_i,v2.1(T) = κ★ · k_i,DLTKcat(T) for **every enzyme-constrained catalytic
  arm**: DLTKcat-covered and median-fallback arms alike, including transport
  enzymes, complexes (all subunits, limiting-subunit rule as built), both
  directions of reversible arm pairs, and CYP51. κ★ preserves every enzyme's
  temperature ratio exactly. κ★ is **not** applied to ribosomal k_R(T) and not
  to non-enzymatic reactions. κ★ is fixed once and never varies across
  temperature, transcript condition, or dose.
- **Calibration rule (frozen):** reporting convention P_ref = 0.25 g gDW⁻¹
  (pool bound = P_ref × φ_met(T)); κ★ = the minimum κ such that network
  μ_max(21 °C) ≥ the condition-mean observed growth μ̄21 (computed from the 21
  retained control vials, the adopted target of §A3), with **no translation
  ceiling active at 21 °C** and the N1-cal full-data maintenance (b = 6.0640)
  imposed. Found by bisection on κ ∈ [1, 100], terminating at relative interval
  width < 10⁻³; solver tolerances as frozen in the prereg. Expected value from
  the linear pre-check: κ★ ≈ 13.9. It is reported as a **13.9-fold aggregate
  effective-capacity correction**, never as a claim that individual enzymes are
  13.9-fold faster in vivo; it may absorb predictor bias, assignment errors,
  missing isoenzymes, stoichiometry errors, biomass composition and pool
  bookkeeping indistinguishably.
- The P_eff profile of the frozen prereg is **replaced**: primary results are
  reported at C_eff★ alone, ΔE_net(C_eff★). Secondary structural sensitivity:
  the relaxation path C_eff = λ·C_eff★ for λ ∈ {1, 1.5, 2, 4, 8, 16},
  approaching N1 as λ → ∞. This path has **no** plausibility interpretation; it
  only shows how fast the network surplus dissolves as capacity is relaxed.

## A2. Translation ceiling placement

The 21 °C ceiling is **removed as an active constraint** (it would impose
translation–metabolism co-limitation at the calibration anchor by
construction — the double-use error already rejected once). 21 °C remains the
normalisation denominator of the relative translation trajectory only.
Ceilings remain active, as upper bounds only, at **15 and 27 °C**.

## A3. Target definition, resolved once

All v2.1 calibration, ceiling comparison and scoring use **condition means of
the per-series fitted rates** (μ̄_c, R̄_c over the retained control vials), with
covariance from the same series-block bootstrap draws (means recomputed within
each draw). Gate 0 is immutable design evidence; its median-based (Hill r₀)
cold-limb check is labelled an initial feasibility screen that did not
guarantee compatibility with the final mean-based target. **Consequence
accepted and carried forward: against the adopted target, the primary
Q10 = 2 translation ceiling is already falsified at 15 °C.** Q10 = 1.5 remains
a labelled sensitivity scenario and is **not** promoted to primary.

## A4. What does not change

- **Maintenance:** the N1-cal foldwise estimates transfer unchanged. The
  capacity correction does not repair or recalibrate the maintenance law; its
  boundary solution (a_GAM = 0) and monotonic temperature residual are
  committed structural findings.
- **Objective:** min-total-enzyme-usage remains primary. κ rescales the
  enzyme-cost objective by 1/κ without changing its argmin, so the loose-pool
  route pathology is *not* expected to be resolved by κ (a binding pool may
  even favour enzyme-cheap wasteful routes more strongly). pFBA-style
  min-total-flux runs as **objective sensitivity only**; any respiration, CUE
  or branch claim that reverses between objectives is objective-dependent and
  not mechanistically identified. Switching the primary objective now would be
  result-chasing and is prohibited.
- Sector lists, weighting, φ_met treatment, Q10 band, K_app rule, CYP51
  scenario ensemble, exchange-bound audit (§3.1), negative controls, and all
  prereg falsification logic not superseded above.

## A5. Carried-forward findings (already observed; never v2.1 validation)

P_min = 3.47 g gDW⁻¹; κ★ ≈ 13.9 under the stated convention; maintenance
degenerates to pure NGAM; monotonic maintenance residual (ρ = −0.93); 15 °C
primary-ceiling violation; 27 °C ceiling binding (+16%); N1 low-CUE
route selection under min-enzyme with loose pool.

## A6. Pre-stated expectations for v2.1 (frozen before execution)

1. Untreated 21 °C growth is matched by construction and is never scored.
2. The 15 °C Q10 = 2 ceiling failure is unchanged by κ.
3. The 27 °C ceiling result is unchanged by κ.
4. A global multiplier preserves all DLTKcat temperature ratios and therefore
   **cannot** generate a 24 °C peak or a 28 °C decline.
5. The monotonic maintenance inadequacy remains.
6. **No directional expectation** is assigned to respiration or CUE after
   capacity calibration.
7. Any apparent repair of respiration/CUE must beat N0/N1 and survive the
   objective sensitivity to count.
8. Absolute proteome-plausibility claims are permanently withdrawn.

## A7. Withdrawn claims (permanent)

Absolute metabolic protein mass; absolute enzyme saturation; physical
plausibility of pool utilisation; absolute enzyme-abundance requirements; any
use of κ★ as evidence for systematically higher in-vivo turnover.

## A8. Evidential strength ordering (corrected)

κ★ consumes one growth anchor, but untreated respiration still informs the
foldwise maintenance law, so the untreated dataset is cross-predicted rather
than held out collectively. The stronger transfer tests, in order: drug-dose
respiration and CUE (maintenance never saw any drug condition); fold-held-out
temperatures under foldwise maintenance; supra-optimal extrapolation from the
sub-optimal-only maintenance fit.

## A9. Rejected alternatives (recorded)

Per-enzyme kcat priors/sampling — uncertainty propagation only, secondary;
two-point temperature calibration — rejected (spends the thermal signal under
test); DLTKcat-covered-subnetwork restriction — rejected as primary (topology
selection bias), retained as a burden-decomposition diagnostic; GECKO-style
iterative enzyme-specific tuning — rejected (adaptively edits relative costs
using the calibration phenotype); alternative kcat predictors — optional
factorised sensitivity, never selected for reducing the deficit.

## A10. Execution

Script 60 implements κ★ calibration and the v2.1 Layer-1/2 scoring under this
amendment; outputs → `tables/revision/v2build/`. The v2 paper structure
becomes two clean stages: (1) the pre-registered structural failure (14-fold
aggregate capacity deficit; temperature-dependent baseline energetic burden;
primary cold-ceiling failure; transcript-temperature sparsity), then (2) the
amended conditional transfer experiment under one effective-capacity anchor.

## A11. Sign-off

Approved to freeze Amendment A1 and execute script 60: __________  Date: ______
