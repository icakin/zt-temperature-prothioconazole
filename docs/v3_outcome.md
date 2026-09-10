# v3 outcome — ETC-GEM one-anchor dose–temperature transfer

**Status: committed result, 17 August 2026. Executed under the frozen
`docs/v3_preregistration.md`. Scripts: `etc_build.py`, `63_calibrate.py`,
`64_thermal_gate0.py`. Outputs: `tables/revision/v3build/`. The v1 manuscript and
all its outputs are untouched.**

## Bottom line

**Outcome = FAILURE at Gate A (prereg §12).** A genome-scale stoichiometric
ETC-GEM cannot mechanistically represent prothioconazole's graded, dose- and
temperature-dependent inhibition of growth, because the drug target — and every
enzyme in the sterol-biosynthesis pathway — carries large excess metabolic
capacity. Per the frozen decision rules, no drug-transfer model is built, and no
sectors / RNA constraints / extra parameters are introduced afterward.

One substantive positive survives, but it is calibration, not validation (§6.3):
the untreated **growth thermal optimum is reproduced mechanistically** by
cold-Arrhenius enzyme activation plus a supra-optimal maintenance term, without
requiring the (non-credible) TemStaPro differential-damage layer.

## Calibration (steps 1–2, `63_calibration.json`)

- Maintenance NGAM(T) = NGAM_0 + α·max(0,T−24)²; NGAM_0 = 2.072 (pinned to
  sub-optimal T≤24 respiration), α = 0.773 (single fitted supra-optimal coeff).
- κ★ = **5.70** (effective pool 0.513 g/gDW at 21 °C); growth(21) = 0.0525 vs
  observed 0.0524. Reported as a 5.7-fold aggregate *effective-capacity*
  correction, never a per-enzyme biochemical multiplier. Sign note: with the
  translation ceiling of v2 removed, DLTKcat capacity under-supplies (κ★ > 1),
  the same direction as — but milder than — v2's ~14× deficit.

## Thermal layer (step 3, `64_thermal_gate0.json`)

- Selected Ea = 0.55 eV, ΔT_m = −4 K (tie rule not needed; min RMSE).
- Predicted untreated μ(T): 0.030 / 0.040 / 0.052 / 0.066 / 0.064 / 0.052 / 0.034
  across 15–28 °C vs observed 0.042 / 0.047 / 0.052 / 0.065 / 0.054 / 0.041 /
  0.036. Predicted peak 24 °C = observed peak 24 °C. Untreated log-μ RMSE 0.191.
- The peak arises from cold-Arrhenius (f_cold rising with T) × supra-optimal
  maintenance (rising above 24 °C); a modest hot term (ΔT_m = −4 K) improves the
  26–28 °C limb but is not the source of the optimum.

## Gate 0 (§4)

| Gate | Result | Meaning |
|------|--------|---------|
| A — CYP51 coupling | **FAIL** | g(a) flat (1.000→0.991) for a∈[0.01,1]; cliff to 0 only at a=0. Binary, not graded; >99% of the decline is one step. Fails monotone-graded, ≥5-distinct, and ≤50%-single-step rules. Identical at 15/21/24/27 °C. |
| B — TemStaPro credibility | PASS (|ΔT_m|=4 K ≤ 10 K) | But Tm IQR = 2.07 K → ranked-vs-uniform (M_R vs M_U) underpowered by construction; the hot layer is near-off at the fit, so no enzyme-specific thermal-stability claim is made regardless. |
| C — enzyme-layer activity | PASS | Pool binds (util 1.00, positive shadow price) at 6–7/7 temps. |
| D — T-potency capacity | **FAIL** | Δ_T(a) = 0 for all a: g(a,T) is temperature-invariant. The network is structurally incapable of temperature-dependent potency. Consequence of Gate A, not independent. |

## Why the drug axis is dead (diagnostic, not a prereg change)

Capacity-reduction sweeps on the whole sterol pathway (drug target unchanged):

| reaction | g(a) behaviour |
|----------|----------------|
| r_0317 CYP51 14α-demethylase | flat → cliff to 0 at a=0 |
| r_0238 C-4 methyl sterol oxidase | flat → 0.971 at a=0 |
| r_0231 C-14 sterol reductase | flat → 0.931 at a=0 |
| r_0698 lanosterol synthase | flat → cliff to 0.020 |
| r_0300 squalene monooxygenase | flat (alt route) → 0.984 |

Every essential sterol enzyme is flat-then-cliff. In a constraint-based model
without enzyme saturation kinetics, flux through a low-demand essential pathway is
binary — enzyme in excess (cut it, nothing happens) or below the stoichiometric
requirement (growth collapses). Sterol demand (∝ biomass ∝ growth) is metabolically
tiny, so no sterol enzyme is capacity-limiting. Graded dose-response is a *kinetic*
phenomenon the stoichiometric layer cannot produce. This is the same reason the v1
CYP51 safety-margin model was demoted and the v2 κ-corrected transfer failed, now
reached mechanistically and definitively.

## Scientific conclusion for the paper

The temperature-dependent potency of prothioconazole that the assays document is
**not derivable from genome-scale stoichiometry plus enzyme–temperature
constraints**. It requires the phenomenological / kinetic description (the
internally cross-validated empirical response surface). The GEM's legitimate,
bounded contribution is that untreated thermal-optimum growth can be reconstructed
from cold-side enzyme kinetics and a supra-optimal maintenance burden — a
mechanistic account of *where the optimum comes from*, with no predictive claim on
the drug response. No mechanistic drug-transfer claim is made.

## What is permanently withdrawn / not claimed

Any mechanistic CYP51 dose-transfer prediction; any temperature-dependent drug
potency from network structure; differential enzyme thermal damage (TemStaPro
ranking, Tm IQR 2.07 K); absolute proteome allocation. Consistent with prereg §12
failure branch: no post-hoc rescue.
