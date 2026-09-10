# v2 model design charter: transcript-informed, enzyme- and temperature-constrained GEM of *Z. tritici*

**Status: design charter — nothing executed, nothing frozen. The full preregistration
is written and frozen only after Gate 0 runs and selects its branch. The v1
manuscript, its models and its figures are untouched by everything in this
document.**

Product of a three-round design review (Claude ↔ ChatGPT, mediated by I. Cakin),
August 2026. Project 12129.

---

## 0. Evidential framing (goes verbatim at the top of the preregistration)

> Because the untreated thermal performance curve and its optimum were
> characterised and published in v1, v2 analyses of untreated growth are
> design-informed consistency tests rather than genuinely held-out predictions.
> The primary predictive tests concern forward respiration and energetic state,
> respiratory branch allocation, transfer across temperatures without
> transcriptome measurements, and dose–temperature generalisation under
> parameters calibrated only in pre-specified reference conditions.

The scientific question v2 asks (a deliberate reversal of v1's logic):
**conditional on the observed transcriptional state, is that state
quantitatively sufficient — through the metabolic network — to generate the
measured physiology?** v1 imposed physiology and read fluxes; v2 imposes
expression-derived capacities and predicts physiology.

A principled failure is an acceptable primary outcome and is designed for from
the start: the deliverable is a structural-necessity map — which parts of the
thermal and drug response the measured sectors, predicted kinetics and
stoichiometry can generate, and which parts demonstrably require unmeasured
processes (protein turnover, intracellular exposure, membrane costs, or others
within the tested model class; failures are attributed only to *classes* of
missing mechanism, never uniquely).

---

## 1. Architecture ladder

| id | model | role |
|---|---|---|
| M0 | thermal ecModel, shared protein pool | baseline |
| **M1T** | M0 + transcript-derived two-sector allocation + **translation-capacity ceiling** | **primary model** |
| M2 | M1T + direct CYP51 inhibition (drug module) | drug extension |
| M3 | M2 + selected enzyme-level relative transcript bounds (A′) | nested extension |
| M4 | M3 + proteostasis turnover cost (C) | pre-registered identifiability test only |
| N0 | empirical null: respiration ≈ a + b·μ_R (fitted within training folds), matched energetic law | growth-proxy null |
| N1 | growth clamped to translation ceiling, loose pool, same hierarchy and FVA | growth-clamped stoichiometric null |

Genome-wide absolute E-Flux (the original architecture A) is **rejected**: one
global transcript-to-protein scale silently fixes hundreds of gene-specific
efficiencies η_i to 1, and 743 hard bounds can manufacture bottlenecks that
mimic predictive success. A′ survives only as M3: *relative* fold-change bounds
E_max_i,c = E_max_i,ref · FC^β on a pre-registered enzyme subset (unambiguous
GPR, counts above floor, plausibly abundance-regulated, constraint active
somewhere, no unmodelled bypass), with β ∈ {0, 0.5, 1} fixed, never fitted.
**AOX is excluded from all transcript constraints** so branch-allocation
agreement with AOX transcripts remains a check, not a construction.

Each rung must beat the rung below on held-out conditions against
matched-complexity nulls, or it is not retained.

## 2. The translation-capacity ceiling (primary thermal mechanism)

Reference-normalised — a conditional thermal ceiling, not an absolute
first-principles prediction:

μ_R,c = μ_obs(21) · [φ_R,c / φ_R,21] · [k_R^lit(T_c) / k_R^lit(21)],  with μ_c ≤ μ_R,c

- φ_R: ribosomal-sector transcript share (mass-weighted, published v1 module
  lists — frozen before this design existed).
- k_R^lit(T): ribosomal elongation-rate temperature dependence fixed from
  literature; never fitted.
- The absolute conversion (λ_R: transcript share → assembled active ribosomes)
  is unidentifiable and is absorbed by the 21 °C anchoring; the claim is the
  *relative thermal trajectory* only.

Why this is the primary mechanism: the DLTKcat kcat ensemble is essentially
flat over 15–28 °C (median apparent Ea ≈ 0.001 eV, run-verified), so the
metabolic layer carries no rising-limb signal; ribosome kinetics can. Empirical
priors: the cold limb ratio μ(15)/μ(21) = 0.706 falls inside the Q10 ∈
[1.5, 2] ribosome-kinetics band (0.66–0.78); the ribosome module score tracks
growth at r = +0.93 across the six RNA-seq conditions.

**Pre-registered structural result (the 24 °C incompatibility).** Under linear
interpolation of φ_R between 21 and 27 °C, no value of r = φ_R(27)/φ_R(21)
both clears the observed 24 °C peak (needs r ≥ 0.955 at Q10 = 2) and makes the
ceiling bind at 27 °C (needs r ≤ 0.494). The intervals do not overlap at any
Q10. Translation-only + linear interpolation is therefore rejected *in
advance*, and three outcomes are pre-registered:

- **A (mild decline, ceiling slack at 27):** translation does not generate the
  falling limb; test metabolic-sector limitation; otherwise classify the hot
  limb as missing mechanism.
- **B (24 °C ceiling violation):** rejects the proxy + interpolation *jointly*
  (φ_R(24) is interpolated); no nonlinear interpolation may be fitted as
  rescue; report that an unmeasured nonlinear sector trajectory could
  reconcile it.
- **C (violation at a transcript-measured temperature, esp. 27):** direct
  falsification of the transcript-to-active-ribosome proxy and/or kinetic
  transfer; interpolation cannot be blamed.

Ceiling-violation tests are primary only at transcript-measured temperatures
(15/21/27); at interpolated temperatures they are joint mechanism +
interpolation tests and are labelled so.

**Attribution discipline:** assigning the falling limb to the metabolic sector
requires, jointly: translation slack there; the pool or a specific enzyme
constraint binding with nonzero dual; causal relaxation removing the decline;
and survival across the whole P_eff profile. Otherwise: "hot limb structurally
unexplained."

## 3. Metabolic pool: P_eff as a profile, not a parameter

P_eff (= pool × saturation, not separately identifiable) is **profiled** over
[P_min, P_max^lit], where P_min is the smallest value supporting observed 21 °C
growth and P_max^lit comes from a pre-specified filamentous-fungal protein
content × metabolic-fraction range. The profile is a structural sweep whose
endpoints mean something: near P_min the pool is maximally restrictive; near
P_max^lit the pool is slack and M1T degenerates continuously into N1. All
surplus results are reported as ΔE_net(P_eff) across the whole profile — a
claim that holds only in a tiny neighbourhood of P_min is a tuned pool, not a
finding. At every profile point report: predicted respiration, pool and
translation utilisation + duals, maintenance contribution, feasible respiration
interval.

## 4. Maintenance: shared nuisance, foldwise, untreated-only

v_maint,c = a_GAM·μ_c + b_NGAM, both ≥ 0, **no temperature or drug dependence**
(temperature-dependent excess demand is exactly what the structural-necessity
analysis must detect; letting maintenance absorb it would kill M4 before
testing it). v1's maintenance values are NOT imported (they were fitted to the
same control physiology Layer 1 scores — a leakage channel now closed).

- Fitted within outer temperature-block folds on untreated conditions only,
  through the frozen common calibration vehicle **N1-cal**: the N1
  stoichiometric network and loose pool, same fixed biomass and respiratory
  coupling assumptions, same secondary hierarchy, with **observed
  training-condition growth imposed as an equality** (μ_c = μ_c^obs) and **no
  transcript-derived translation ceiling anywhere in the nuisance fit**. It
  cannot be M1T (maintenance would depend on P_eff), and it cannot be
  ceiling-clamped N1 (translation-model error would masquerade as GAM/NGAM).
  The fit uses paired series-block bootstrap draws of (μ, R), propagating
  their shared-trace covariance rather than treating growth as error-free.
- The frozen pair (a_GAM, b_NGAM) is passed unchanged, per fold, to N0's
  energetic law, predictive N1 (growth set by the translation ceiling), every
  P_eff point of M1T, and later M2–M4. No model gets its own maintenance fit.
  This separates two questions cleanly: (1) given observed biomass production,
  what temperature-invariant affine ATP demand accounts for baseline
  respiration; (2) can the architecture predict that biomass production and
  its respiration under the frozen energetic law.
- **Adequacy diagnostic (pre-registered):** report the correlation of the
  N1-cal training residuals (R^obs − R^N1-cal) with temperature. A strong
  monotonic pattern means the temperature-invariant baseline is already
  structurally inadequate; it is reported as such and **never repaired by
  refitting a temperature term** — that inadequacy is itself input to the
  structural-necessity map.
- Never fitted to any drug condition: the untreated-trained law transfers
  unchanged across the entire dose surface (Layer 4 transfer test).
- Secondary directional test: fit maintenance on sub-optimal untreated
  conditions only, predict supra-optimal respiration/CUE without refitting —
  the direct test of whether baseline energetics can explain the collapse.

## 5. Drug module (M2)

v_CYP51,c ≤ kcat_CYP51(T_c) · E_CYP51,c · 1/(1 + D/K_app), Hill coefficient
fixed at 1. K_app is explicitly phenomenological (absorbs uptake, efflux,
partitioning, binding — D is extracellular), calibrated on the **normalised
24 °C dose–response only** (EC₅₀ there is 1.40 mg L⁻¹; expect K_app of order
1), growth only; respiration and CUE stay as drug-response validation.

CYP51 kcat(T) is unresolved by design — three fixed scenarios, never selected
by fit: (1) **constant (primary)**; (2) DLTKcat trajectory (decreasing);
(3) generic monotonic Arrhenius with coefficient fixed from the broad enzyme
set. Report what survives all three.

At 2 mg L⁻¹, run both: direct inhibition + control transcriptome, and direct
inhibition + drug transcriptome; the difference is the model's estimate of the
transcript-associated adaptation component. Drug transcriptomes never replace
the molecular perturbation.

**C2 constraint on interpretation (from the v1 shift-versus-scale result):**
failure of one K_app to transfer across temperatures falsifies the
temperature-invariant apparent-potency module but does **not** discriminate
temperature-dependent exposure from amplitude changes, sterol toxicity, or
other drug–temperature interactions; successful transfer shows sufficiency of
a simple apparent-potency model, not uniqueness of CYP51 inhibition. No
toxic-intermediate sink unless CYP51 limitation fails reproducibly, and then
only against a matched one-parameter generic drug-stress null.

## 6. RNA-seq usage split

- **Input:** control-transcriptome sector shares (mass-weighted, v1 module
  lists); drug transcriptomes in the Layer-2 conditional comparisons.
- **Optional input (M3 only):** selected metabolic enzyme fold changes.
- **Held out:** proteostasis/chaperone/protease transcripts (reserved as the
  independent qualitative check on M4's predicted turnover burden); AOX
  transcripts (branch-allocation check).
- **Descriptive only:** efflux and drug-response genes that cannot affect
  intracellular drug without a transport model.
- Within-gene fold changes only for enzyme-level use; never absolute cross-gene
  abundance. Expression floor — low counts never become knockouts. GPR rules:
  complexes = limiting subunit with stoichiometry; isozymes additive;
  multifunctional = one shared pool; ambiguous = unconstrained.
- Replicate uncertainty: bootstrap the 7–8 RNA replicates; predictions are
  distributions over RNA variation × kcat scenario × alternative optima.

## 7. Optimisation and sharpness

Frozen lexicographic hierarchy: (1) maximise growth (subject to ceilings);
(2) minimise total enzyme usage; (3) FVA on respiration, carbon uptake, AOX and
key exchanges at/near the optimum; report CUE intervals. Sharpness gate: W_y =
median feasible-interval width / observed cross-condition range, with a
pre-registered ceiling (or required reduction vs N1). A point estimate inside
an enormous feasible interval is "compatible", never "predicted". Qualitative
changes under small near-optimal relaxations disqualify strong mechanism
claims. A mechanism that never binds (zero dual) has explained nothing — the
constraint-activity map is a primary result.

## 8. Metrics (pre-registered, three + gates)

1. **E_μ** — log-scale RMSE on untreated growth, excluding 21 °C (anchor), with
   hard shape gates: correct limb signs; predicted optimum within one
   experimental temperature interval of 24 °C; no observed growth above the
   translation ceiling. Non-RNA temperatures reported separately.
2. **E_phys** — stacked residual vector (Δlog μ, Δlog R) across held-out
   conditions, scored with the **full series-block bootstrap covariance**
   (μ and R come from one DO trace; cross-condition dependence preserved),
   with pre-frozen shrinkage/regularisation and dimension normalisation. CUE
   is derived, never a third likelihood term; report the signed supra-optimal
   CUE bias B_CUE,hot and dose–temperature bias surfaces. Same for effective
   P/O.
3. **ΔE_net(P_eff)** — E_phys(best null) − E_phys(M1T), across the full
   profile, with series-block bootstrap and the Stage-1 materiality
   convention. This is the metric that prevents republishing r = 0.93 in
   mechanistic notation.

## 9. Validation layers and evidence hierarchy

Layers: (1) untreated thermal physiology (design-informed consistency; 21 °C
growth excluded); (2) six transcript-informed conditions (growth not
independent here; ceilings, binding, energetic errors, surplus vs N0/N1; AOX
agreement secondary); (3) temperatures without RNA-seq — the most important
for the translation architecture; interpolated (18, 24, 26) and any
extrapolated conditions reported separately, never pooled; (4) drug
generalisation (held-out doses and temperatures under fixed K_app and
untreated-trained maintenance); (5) negative controls — transcript permutation
among abundance-matched enzymes, temperature-label shuffling, random
constrained-enzyme sets, constant-sector replacement, temperature-only
empirical interpolation, and M3-vs-M2.

Evidence ranking (weakest → strongest): architecture consistency (untreated
TPC shape) → cross-predicted untreated respiration → supra-optimal
extrapolation of energetic state → network surplus across the P_eff profile →
branch prediction with sharpness → drug transfer → structural-necessity
mapping.

## 10. Gate 0: registered feasibility audit (the only next step)

Runs on existing files only (vsd matrix, v1 module lists, published r₀ values,
DLTKcat table). Script, thresholds, tolerances, input hashes, output schema and
the complete decision tree are written and frozen **before** it runs; the
branch taken is recorded without modifying the script. Branching uses the
nominated primary configuration ONLY (v1 gene lists, mass-weighted shares,
primary k_R curve, constant CYP51 kcat, β = 0); alternative preprocessing
diagnoses robustness but does not vote. All Gate-0 quantities (φ_R, φ_met at
15/21/27; cold-limb bracket; r-inequalities; P_min; DLTKcat coverage) are
thereafter **design evidence**, never presented as independent tests. This is
a registered feasibility audit preceding architecture freeze — not a blinded
preregistration.

Pre-written terminal branches: (1) cold-limb bracket fails → translation is
not the cold-limb mechanism; fall back to the fixed global-Ea architecture
(Ea ∈ {0.3, 0.55, 0.8} eV, scenario-fixed, never fitted). (2) cold limb
feasible + measured-temperature ceilings valid → retain translation ceiling.
(3) 27 °C ceiling slack → falling limb not attributed to translation; test
sector limitation, else "unexplained". (4) 27 °C hard violation → reject the
translation proxy. (5) 24 °C violation only → reject proxy + linear
interpolation jointly. (6) P_min outside the literature pool range → metabolic
pool module flagged implausible. (7) poor DLTKcat sector coverage → restrict
claims to covered subnetwork or apply the pre-written missing-enzyme fallback.

## 11. Scenario discipline

No fully crossed grid. One primary configuration (translation variant,
constant CYP51 kcat, β = 0, v1 sector lists, mass-weighted, K_app from
normalised 24 °C) plus factorised challenges: thermal-generator
(translation vs three global-Ea); sector (lists/weighting); drug (three CYP51
scenarios, sector model fixed); M3 (β ensemble, only after activation);
selected interactions only where the interaction is the scientific question
(CYP51 scenario × temperature).

## 12. Pre-registered falsification (summary)

- **Translation ceiling:** observed growth (with uncertainty) above the
  ceiling at a transcript-measured temperature; or slack on the limb it is
  meant to generate (binding = utilisation ≥ 0.95 + nonzero dual, pre-fixed).
- **Network surplus:** μ̂ ≈ μ_R everywhere AND removing pool/enzyme constraints
  barely changes growth AND M1T fails to materially beat N0/N1 on held-out
  respiration/energetics → the ecModel adds nothing beyond the ceiling.
- **M3:** transcript constraints rarely bind; improvement only at β = 1;
  permuted transcripts perform similarly.
- **M4 (C):** turnover coefficient at zero/boundary; growth improves but not
  respiration+CUE jointly; matched generic maintenance null performs equally;
  activation requires a fitted Tm shift. Expected outcome: non-identifiable —
  reported as such, with failure attributed to a *class* of missing
  temperature-dependent energetic processes, not to turnover uniquely.
- **Drug module:** per section 5, with the C2 ambiguity stated in both
  directions.
- **Whole framework:** feasible ranges span most of the observed range
  (permissive, not predictive); or point accuracy depends on an arbitrary
  secondary objective.

Retention additionally requires: improvement over the immediate nested null on
held-out conditions, correct qualitative temperature and dose responses,
acceptable uncertainty calibration, and stability under the pre-specified
sensitivity ensembles. A small RMSE drop retains nothing.

## 13. Three publishable endpoints (fixed before any result exists)

- **Positive:** translation capacity explains the TPC shape; network
  constraints add material held-out energetic and branch predictions across
  the P_eff profile.
- **Partial:** translation capacity explains the growth shape; the network
  adds no surplus.
- **Structural failure:** even the ceiling cannot reproduce the thermal
  profile, or energetic/drug behaviour requires unidentified latent processes
  — reported as the structural-necessity map naming the missing measurement
  classes (intracellular drug/desthio, proteomics, direct turnover assays).

All three are publishable without redefining success after seeing results.

## 14. Open items before Gate 0 can be written

1. ~~Maintenance-calibration vehicle~~ **Settled: N1-cal** — loose-pool
   stoichiometric model with observed growth clamped within untreated training
   folds; translation ceiling excluded from nuisance fitting; fitted
   (a_GAM, b_NGAM) transferred unchanged to N0, predictive N1, and the full
   M1T P_eff profile (section 4).
2. Fix the literature values: k_R(T) source and Q10 band; P_max^lit range;
   pre-specify both.
3. Freeze the Gate-0 script + decision tree (one document, hashes recorded).
4. Sign-off by I. Cakin (and G. Yvon-Durocher if this becomes the next
   project) before anything runs.

*Estimated effort after sign-off: Gate 0, one afternoon. M0/N1/M1T build with
maintenance folds and profile, roughly one to two weeks. Drug module and
challenges, another week. Entirely separable from, and without risk to, the
submitted v1 manuscript.*
