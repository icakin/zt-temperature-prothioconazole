# Preregistration v4 — transcriptome-driven growth prediction (FROZEN)

**Status: FROZEN before any model run. Amendments after freezing are logged in
§7 and invalidate affected gates. Failure of any gate → the analysis is archived
in `tables/revision/v4/` and does NOT enter the manuscript. No silent re-runs
with altered mappings.**

Data (final; no new experiments): DO-derived growth/respiration posteriors
(7 T × 8 doses); RNA-seq for 6 conditions (15/21/27 °C × 0/2 mg L⁻¹, n = 7–8);
ztGEM v03 (plain) and ecZtGEM (sensitivity only). Settled and not revisited:
drug-from-stoichiometry (gate-failed v3), sign audit, AOX interpretation
(condition-constrained + caveats stands regardless of v4 outcomes).

---

## Approach B — control-trained ribosome forecast (runs FIRST, runs regardless)

**Model.** Ribosome-only transcript-allocation proxy. Per-sample ribosome score
s = mean log₂(normalized count + 1) over the 117-gene ribosome module (identical
set to script 81). Condition score = mean over samples. Fit
log(growth) = a + b·s on the THREE CONTROL conditions only (n = 3). Predict
growth of the three drug conditions from their scores, blind to their measured
growth. Measured growth = posterior median of the condition's biomass-specific
growth rate; its 95% CrI is the comparison interval.

**Also run (secondary):** leave-one-temperature-pair-out — drop both conditions
at one temperature, fit a, b on the remaining four, predict the held-out pair;
report all six predictions and residuals in C C⁻¹ h⁻¹ and on the log scale.

**Gates (frozen):**
- **B1:** each of the three blind drug-condition predictions falls inside the
  measured growth 95% CrI.
- **B2:** predicted log-relative drug effect satisfies δ̂₂₇ > δ̂₁₅ where
  δ̂_T = log(ĝ_{T,2}/g_{T,0}) (attenuation reproduced).
- **Reporting:** residuals reported in full whatever the outcome. No correlation
  statistics on 3 predictions.

**Placement if passed:** supplementary figure + ≤2 main-text sentences framed as
"a coarse-grained transcript-allocation proxy"; explicitly state-based, not
causal (post-treatment RNA; reverse causation not excluded). If failed:
archived, one honesty sentence MAY be added to Discussion at author discretion.

**Uncertainty:** bootstrap RNA replicates within condition (500 draws) for the
score; joint with growth-posterior draws for gate B1/B2 robustness
(report P(gate holds) over draws; gate passes on the point estimate, the
probability is reported alongside).

---

## Approach A — expression-constrained FBA (E-Flux, plain ztGEM)

**Frozen specification (no element may be tuned against growth):**
- Model: ztGEM v03, sucrose minimal medium, identical exchange bounds, identical
  NGAM, standard biomass objective — same in all six conditions. ecZtGEM version
  (expression scales relative enzyme draws; pool 0.09 g/gDW and kcat(T) frozen)
  is a prespecified sensitivity analysis only, labelled "transcript-proxy
  enzyme-allocation".
- Expression input: DESeq2 size-factor–normalized counts (one global
  normalization across all 45 samples), condition mean per gene, pseudocount +1.
  Nonnegative abundances; no VST z-scores; no per-condition rescaling.
- GPR mapping: AND = min, OR = sum.
- Bound mapping: v_max(r, cond) = V · e(r, cond) / max_cond' e(r, cond') where
  e is the GPR-mapped expression and V is one global scale set so that the
  UNCONSTRAINED control-15 °C model is not expression-limited at its optimum
  (V chosen once, before any growth comparison, from feasibility alone).
  Reactions with no measured GPR gene: unconstrained. Reversible reactions:
  symmetric bounds.
- Prediction: maximize biomass per condition. Output = predicted relative growth
  (6 values, arbitrary common scale).

**Gates (frozen):**
- **A1 (ordering):** Kendall τ ≥ 0.867 against measured growth over the six
  conditions, evaluated on RESOLVABLE pairs only — a pair is resolvable when the
  two conditions' measured growth posteriors differ with P ≥ 0.95. Equivalent
  statement: at most one discordant resolvable pair.
- **A2 (sign):** predicted δ_T = log(ĝ_{T,2}/ĝ_{T,0}) < 0 at all three
  temperatures.
- **A3 (attenuation):** predicted δ₂₇ > δ₁₅.
- **A4 (network specificity — required for any manuscript use):** the E-Flux
  prediction must beat BOTH (i) the ribosome-only predictor of Approach B
  (lower total squared log-residual over the six conditions under
  leave-one-temperature-pair-out) and (ii) null mappings: condition-label
  permutation of expression profiles and gene-to-GPR scrambling
  (degree-preserving; 200 scrambles) — observed τ must exceed the 95th
  percentile of each null distribution.
- **Robustness:** RNA-replicate bootstrap (500) for τ CI; reported, not gated.

**Placement if all gates pass:** an added prediction layer, NOT a replacement of
the condition-constrained analysis, and NOT support for the AOX causal claim.
Frozen wording ceiling: "transcript-informed capacity constraints reproduced the
condition ordering and attenuation of growth, indicating the regulatory state
was sufficient to reconstruct the observed phenotype within this experiment."
If any gate fails: archived in `tables/revision/v4/` with the frozen protocol
and outputs; not in the manuscript.

---

## Explicitly excluded (frozen)
Method leaderboards (GIMME/iMAT selection after results); threshold/percentile
integration methods; parameter tuning of medium, NGAM, pool, GPR rules, floors
or caps against measured growth; multi-sector fitted allocation models; ML on
six outcomes; mediation claims; RNA samples treated as independent growth
outcomes (effective n = 6 throughout); ΔFBA; dynamic FBA; any modification of
the v3 conclusions or the AOX text.

## Order of execution
1. Approach B (script 84). 2. Approach A (script 85) with nulls. 3. Single
verdict report `tables/revision/v4/v4_verdict.md`.

## §7 Amendment log
1. (Clarification, before gates evaluated) "Measured growth posterior" =
   Bayesian bootstrap over replicate vials at the exact T × dose cells — the
   same source the condition-constrained flux analysis uses. Chosen over
   evaluating the Sharpe–Schoolfield posterior because the RNA conditions
   coincide with measured cells.
2. (Before Approach A ran; after B's first run) Ribosome module unified to the
   117 KEGG map03010 genes (matching Fig. 5C and the module-trait analysis);
   the first B run used a broader 153-gene annotation-OR set mislabelled
   "117 genes" in script 81. B was run under BOTH definitions: verdict
   identical (B1 FAIL, B2 pass) — no gate outcome affected.

## OUTCOME (post-run)
Approach B: FAIL (B1). Approach A: FAIL (A1, A2, A4). Both archived in
`tables/revision/v4/`; see `v4_verdict.md`. Neither enters the manuscript.
