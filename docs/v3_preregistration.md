# v3 preregistration — ETC-GEM one-anchor dose–temperature transfer

**Status: DRAFT for freeze. Once approved by I. Cakin this document is FROZEN;
any change is a dated amendment appended below, never an edit. Design agreed
across four adversarial ChatGPT rounds (17 August 2026) plus four Claude
counter-points, all resolved; this document instantiates that agreement with the
data facts and remaining frozen constants. The v1 manuscript and all its outputs
remain untouched by everything here. Outputs → `tables/revision/v3build/`;
scripts 63–66.**

---

## 0. The question

Does a genome-scale stoichiometric network, combined with calibrated effective
catalytic capacity, a low-dimensional thermal activation/damage layer,
temperature-dependent maintenance and **one** CYP51 binding constant, **predict
temperature-dependent prothioconazole inhibition of growth** at temperatures
where no drug data were used in fitting — better than empirical dose–temperature
models given the same outcome information?

This is a **one-anchor dose–temperature transfer test**, not leave-temperature-out
prediction: untreated growth at every target temperature *is* used to fit the
thermal baseline. The only genuinely held-out information is the drug response at
the six non-anchor temperatures.

## 1. What is dead and not resurrected (committed, from v2)

1. Absolute-proteome interpretation of the ecModel (P_min = 3.47 g gDW⁻¹ vs 0.25
   bound; ~14× aggregate capacity deficit). **Permanently withdrawn.**
2. κ★ as a literal biochemical multiplier; a single κ★ transferring respiration/CUE
   (lost to a log-linear null in 100 % of draws). κ★ returns here **only** as an
   effective catalytic-resource scalar, never interpreted physically.
3. DLTKcat kcat(T) thermal *shape* as the thermal layer (median apparent
   Ea ≈ 0.001 eV; flat 15–28 °C; cannot generate the 24 °C peak). DLTKcat is used
   here **only** as the reference enzyme ranking at 21 °C.
4. RNA-seq as a foundational constraint (3 temperatures, too sparse to interpolate
   a 24 °C peak). Not in the primary model at all.

## 2. Data facts (frozen)

- Grid: T ∈ {15, 18, 21, 24, 26, 27, 28} °C × dose ∈ {0, 0.06, 0.12, 0.25, 0.5,
  1.0, 2.0, 4.0} mg L⁻¹ prothioconazole, 3 replicate series nominal (some cells
  1–2; e.g. 28 °C/0.06 = 1). Fitted μ and R per series from `tables/physiology/
  derived_N0_R_results_with_carbon.csv`.
- **Seven** nonzero doses. Transfer set = 7 nonzero doses × 6 target temperatures
  {15, 18, 21, 26, 27, 28} = **42 held-out drug cells**. 24 °C is the drug anchor.
- Base network `models/gem/ztGEM_v03.xml` (validated stoichiometry, biomass,
  growth reaction r_2111). ecModel machinery `models/gem/ecZtGEM_full.xml`.
- DLTKcat kcat ensemble + TemStaPro Tm per enzyme (existing).
- N_inoc/viability: the 0.5-viability question is **out of scope here**; v3 uses
  the derived rates as committed. Any future N_inoc correction is a separate,
  later amendment and does not gate this build.

## 3. Frozen architecture — primary ETC-GEM

Enzyme-constrained flux bounds on plain ztGEM_v03:

    v_i(T) ≤ κ★ · k_i,21^DLT · f_cold(T; Ea) · f_hot,i(T; θ_H) · E_i
    Σ_i MW_i E_i ≤ P_eff            (P_eff an EFFECTIVE resource, not measurable grams)

Maintenance:

    NGAM(T) = NGAM_0 + α · max(0, T − 24)²

CYP51 inhibition — active CYP51 capacity scaled by:

    a_CYP51(D) = 1 / (1 + D/K_app)

Primary RNA integration: **none**.

### 3.1 Thermal functions (frozen)

Cold side (global, not per-enzyme):

    f_cold(T; Ea) = exp[ −(Ea/R) · (1/T − 1/294.15) ]           T in K, 294.15 = 21 °C

Discrete candidates Ea ∈ {0.3, 0.55, 0.8} eV. Selecting among them on untreated
growth counts as **one fitted discrete parameter**.

Hot side (per-enzyme, TemStaPro-ranked):

    f_hot,i(T) = 1 / (1 + exp[(T − T_m,i − ΔT_m) / s])

**s = 2 K, FROZEN before execution.** Fit only ΔT_m. (If s were also fitted the
inventory would rise from five to six — it is not.)

Uniform-damage null M_U: every T_m,i replaced by the **median** TemStaPro T_m,
same fixed s, same fitted ΔT_m.

## 4. Gate 0 — feasibility, no fitting beyond the specified diagnostic sweeps

Gates run on the untreated-calibrated model (α, κ★, Ea-scenario, ΔT_m, pool
convention installed via §6 steps 1–3). No drug-response observation is used in
any gate. **All gate statistics are recorded for all three Ea scenarios; the
selected primary Ea must pass.**

### Gate A — CYP51 coupling (at 24 °C)
Sweep residual CYP51 activity a ∈ {1, 0.9, …, 0.1, 0.05, 0.02, 0.01, 0}; compute
g(a) = μ(a)/μ(1). Pass iff ALL hold:
- g(a) monotone nondecreasing in a within solver tolerance;
- full inhibition (a=0) reduces growth by enough to span the normalised 24 °C
  response used for K_app;
- ≥ 5 activity values give predictions differing by > 1 % of untreated growth;
- no single adjacent step accounts for > 50 % of total predicted decline.

Fail ⇒ response flat / binary / nonmonotone / too small ⇒ **no drug model is
built**.

### Gate B — TemStaPro credibility
Record full T_m distribution **before** fitting. Credibility boundary
**|ΔT_m| ≤ 10 K, FROZEN.** If the minimum offset to place appreciable damage
inside 15–28 °C exceeds 10 K: TemStaPro loses physical interpretation; the ranked
model is not primary; the uniform hot-side model may still proceed as an
explicitly phenomenological null; **no claim about enzyme-specific thermal
stability is allowed.** Passing does not validate TemStaPro — the ranked-vs-uniform
comparison (§9) remains mandatory. Also record the T_m spread (IQR); if the spread
is negligible, note that M_R vs M_U is underpowered by construction.

### Gate C — enzyme-layer activity
After κ★, Ea, hot-side and pool convention installed, at each of the seven
untreated temperatures require EITHER global pool utilisation ≥ 0.95 with a
positive growth shadow price, OR ≥ 1 temperature-dependent enzyme-capacity
constraint with a positive growth shadow price. If neither: the layer is inert at
that T. Require activity at **all seven** for the full ETC claim. Any failing T
cannot support an enzyme-mechanism claim there; **failure at ≥ 2 temperatures
terminates the predictive ETC build.**

### Gate D — temperature-dependent potency *capacity* (direction-blind, Claude point 1)
Before using any drug-response observation, on the frozen untreated thermal model,
for residual CYP51 activities a ∈ {0.1, 0.2, …, 0.9} compute g(a,T) = μ(a,T)/μ(1,T)
and

    Δ_T(a) = max_T log g(a,T) − min_T log g(a,T).

**Gate D passes iff Δ_T(a) ≥ 0.05 for at least two ADJACENT interior activity
levels a ∈ {0.2, …, 0.8}.** Additionally: repeat with solver tolerances tightened
tenfold, and each qualifying Δ_T(a) must change by < 10 %. Evaluate the range
**direction-blind** (hotter may raise or lower inhibition); **do not compare the
predicted direction with observed z until after Gate D is recorded.** The selected
primary Ea must pass; record Δ_T(a) for all three Ea scenarios. If the primary
fails, **stop the ETC predictive build** — the network is structurally incapable of
producing a temperature interaction, and this is reported as the finding. Passing
establishes only structural capacity, not correctness.

## 5. Primary target and scoring

Mechanistic response (temperature-dependent inhibition):

    z(T,D) = log[ μ(T,D) / μ(T,0) ]

Primary score (each temperature equally weighted, replicate-count-agnostic):

    RMSE_z = sqrt( (1/6) Σ_T (1/n_D) Σ_{D>0} [ ẑ(T,D) − z_obs(T,D) ]² )

Secondary: RMSE_logμ on absolute growth at the same held-out cells (tests the full
predictor incl. calibrated baseline). Never score extrapolated EC₅₀; predictions
scored only at doses actually observed at that temperature.

Paired series-block bootstrap; **every draw repeats the full calibration sequence**
including Ea selection and K_app fitting.

## 6. Calibration order (frozen — no parameter trading)

1. **α** from untreated respiration, T_c = 24 fixed, conditional on observed
   untreated growth. Not fitted to predicted growth; never refit after seeing drug
   predictions. (Uses the N1-cal r_4046-equality maintenance machinery.)
2. **κ★** to untreated 21 °C growth only. Pool is an effective resource; no
   absolute-proteome claim.
3. **Thermal layer:** for each fixed Ea candidate, fit the single ΔT_m on untreated
   growth at the remaining temperatures with α, κ★ frozen. Select Ea by **minimum
   untreated log-growth RMSE**; tie rule = **lower Ea**. The untreated TPC is
   calibration only — zero validation evidence.
4. **K_app** from the normalised 24 °C dose response μ(24,D)/μ(24,0) only. No other
   drug-treated temperature enters any ETC fit.

## 7. Maintenance misspecification — explicit (Claude point 3)

The selected maintenance function represents **supra-optimal energetic burden
only**. It cannot represent the already-observed monotonic maintenance residual
across the full temperature range (Spearman −0.93, cold included). Consequently,
**cold-side energetic inadequacy remains deliberately unresolved**, and its effects
on fitted growth capacity may be absorbed by κ★, Ea and the hot-side calibration.
No claim is made that the model identifies temperature-dependent maintenance below
24 °C. This is a known structural misspecification accepted for parsimony, not a
successful energetic model. No cold-side maintenance parameter is added.

## 8. Ea-scenario sign disclosure — mandatory (Claude point 2)

For each Ea ∈ {0.3, 0.55, 0.8}, repeat the untreated thermal calibration and fit
that scenario's single K_app at 24 °C. On the fixed observed-dose grid define
z̄_Ea(T) = mean over nonzero doses of ẑ_Ea(T,D), and the temperature–potency
gradient γ_Ea = OLS slope of z̄_Ea(T) ~ T over the six transfer temperatures.
Report all three γ_Ea, all three full z̄_Ea(T) trajectories, and Gate D ranges for
all scenarios. **If sign(γ_Ea) differs across scenarios, any mechanistic
temperature–potency claim is downgraded to "conditional on the untreated-growth-
selected thermal scenario," and a bootstrap win does not override the downgrade.**

## 9. Frozen comparators

- **H1 (matched one-anchor Hill):** same fitted untreated baseline as ETC; fit
  K_H and h on 24 °C data only; transfer both unchanged to the six targets. One
  more continuous parameter than the ETC drug module — acceptable; ETC must beat a
  strong dose-shape ablation. **A Hill curve is never refitted per target
  temperature.**
- **H0 (one-parameter Hill):** fix h = 1, fit only K_H at 24 °C.
- **G1 (same-information tensor GAM):** trained on all seven untreated controls +
  all 24 °C dose observations + **no** nonzero-dose target-temperature data; predict
  the identical transfer set. Basis, smoothing, RE **frozen before execution**
  (tensor te(T, log10(D+0.03)), cr bases, REML, series RE; matched to the v2
  empirical-surface conventions). Fair matched-information statistical comparator.
- **G2 (full-information temperature-block GAM):** secondary benchmark only; per
  target T, fit on all other temperature surfaces and predict the omitted one.
  Receives far more drug information than ETC — the achievable statistical ceiling,
  **not** a matched ablation. Never merged with G1 in the headline.

### 9.1 Null-family redundancy — honesty (Claude point 4)
Report pairwise correlation of out-of-fold z predictions, pairwise correlation of
cellwise losses, and pairwise RMSE differences among H0/H1/G1. **Collapse rule
(pre-frozen):** if two nulls have prediction correlation ≥ 0.98 AND absolute RMSE
difference < 0.10 ν_z, treat them as one effective flat-transfer null family when
interpreting evidence. All numerical comparisons kept; beating highly correlated
nulls is not described as repeated independent validation.

## 10. Noise scale ν_z (frozen definition)

For each target cell j = (T,D), D>0: z_j = log μ̄_{T,D} − log μ̄_{T,0}. In
series-block bootstrap draw b, z_j^{(b)} = log μ̄_{T,D}^{(b)} − log μ̄_{T,0}^{(b)},
the **same** resampled 0-dose observations used for every dose at a given T
(preserving ratio covariance). A draw is valid for j only when both treated
numerator and same-T control denominator are represented.

**Eligible cells** for ν_z: treated condition ≥ 2 independent series; same-T
control ≥ 2 independent series; ≥ 80 % of draws yield a valid z_j^{(b)}. Singletons
stay in all predictive scores but do **not** determine the materiality floor.

For each eligible cell: s_j = 1.4826 · median_b | z_j^{(b)} − median_{b'} z_j^{(b')} |
(bootstrap MAD). Global scale:

    ν_z = median_{j∈E} s_j             (E = eligible transfer cells)

Cell-balanced, robust, covariance-propagating, not deflated by singletons.

**Adequacy rule:** require ≥ 12 eligible target cells spanning ≥ 4 of 6 target
temperatures. If this fails, the 0.25 ν_z materiality criterion is declared **not
estimable** — report bootstrap win probabilities and raw RMSE differences but **no
formal materiality win** (no post-hoc pooled fallback). Singleton-cell prediction
intervals remain reportable from the fitted predictive distribution but do not
enter ν_z.

## 11. Win conditions (all must hold for the strong claim)

For the claim *"network stoichiometry and differential enzyme thermal capacity,
combined with one drug-binding constant, predict temperature-dependent
prothioconazole inhibition":*

**Overall transfer** — against BOTH H1 and G1:

    P( RMSE_z^ETC < RMSE_z^null ) ≥ 0.90
    median( RMSE_z^null − RMSE_z^ETC ) > 0.25 ν_z

**Genome-scale hot-side contribution** — same two conditions for M_R vs M_U.
Otherwise TemStaPro ranking is decoration.

**Block robustness** — ETC must: beat each principal null at ≥ 4 of 6 target
temperatures; be worse than a principal null by more than the materiality floor at
≤ 1 temperature; preserve correct monotonic dose direction at every target
temperature.

**Absolute-growth guardrail** — RMSE_logμ^ETC ≤ 1.10 × RMSE_logμ^G1.

**Uncertainty** — report 50/80/95 % interval coverage and width on transfer cells;
ETC proper interval score no worse than G1. Coverage is not itself a win gate.

## 12. Decision outcomes (pre-stated)

- **Full success:** ETC-ranked beats H1, G1 and M_U under every primary rule ⇒
  "A genome-scale, differentially temperature-constrained network predicts
  temperature-dependent drug inhibition from one potency calibration." (Downgraded
  to scenario-conditional if §8 sign disagrees.)
- **Partial:** uniform ETC beats H1 and G1 but M_R does not beat M_U ⇒ "The ETC-GEM
  transfers the dose response; enzyme-specific TemStaPro ranking adds no predictive
  information." No differential-damage claim.
- **Statistical tie:** ETC within materiality of G1 with fewer fitted interaction
  terms ⇒ "Semi-mechanistic compression with comparable predictive performance."
  Not superiority.
- **Failure:** ETC fails to beat H1 or G1, or fails a Gate ⇒ "The genome-scale
  structure does not add predictive information beyond empirical dose–temperature
  models under the available data." **No** added sectors, RNA constraints, extra
  temperature parameters, or temperature-dependent K_app introduced afterward.

## 13. Final parameter inventory (exactly five fitted)

1. κ★ — 21 °C untreated growth.
2. Ea — discrete selection among three, on untreated growth.
3. ΔT_m — untreated thermal growth, s fixed at 2 K.
4. α — untreated respiration, conditional on observed growth.
5. K_app — normalised 24 °C dose response.

H1 null has K_H, h (two, same anchor). ETC must beat it despite one fewer drug
parameter.

## 14. What the model can and cannot claim

CAN (if it wins): a genome-scale stoichiometric network + calibrated effective
capacity + low-dimensional thermal activation/damage + maintenance + one CYP51
binding constant predicts condition-specific temperature-dependent inhibition.

CANNOT, ever: absolute proteome allocation; an independently predicted thermal
optimum; enzyme-specific thermal kinetics from DLTKcat; unique identification of
heat-damage mechanism; intracellular drug exposure; temperature-dependent
maintenance below 24 °C.

## 15. Build order (scripts 63–66)

1. Reconstruct clean enzyme-constrained ztGEM_v03; verify mass balance; reproduce
   unconstrained FBA growth. **(no fitting)**
2. Add reference-temperature enzyme constraints; calibrate α, κ★ (§6.1–6.2).
3. Add cold Arrhenius + TemStaPro-ranked hot side; fit ΔT_m per Ea; select Ea
   (§6.3).
4. Run Gate 0 (A/B/C/D) on the untreated-calibrated model. Hard stops.
5. If Gates pass: fit K_app (§6.4); run transfer; score vs H0/H1/G1/G2 and M_U with
   the full bootstrap; §8 scenario disclosure; §9.1 redundancy; §11 win rules.
6. Every added layer cross-validated against the simpler model beneath it.

## 16. Sign-off

Approved to freeze v3 preregistration and execute the build: __________  Date: ______
