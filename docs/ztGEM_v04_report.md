# ztGEM v0.4 — the model pushed as far as the data allow

Project 12129 · 2026-08-07 · supersedes the v0.3 analysis section; model file
unchanged (`ztGEM_v03.xml`), all new analysis layers in `advance_ztGEM_v04.py`.

## What was added over v0.3

1. **Maintenance fitted to your own data.** Using the four sub-optimal control
   temperatures (15/18/21/24 °C, where AOX transcript is at baseline), surplus
   ATP scales with growth: ATP_maintenance ≈ 79.2 × µ mmol gDW⁻¹ h⁻¹ (NNLS,
   NGAM ≈ 0, r = 0.83). The yeast default (0.7 fixed) is replaced by this fit.
2. **Required-AOX recomputed under fitted maintenance**, with FVA showing the
   AOX flux is *unique* given the constraints (min = max), i.e. the network has
   no alternative dissipation route at the parsimonious optimum.
3. **E-Flux transcriptome-constrained models** (expression sets reaction
   capacity; calibrated at the single condition 21/0, all others held out).
4. **Full flux readouts per condition**: effective P/O, ATP budget, CYP51 flux,
   subsystem flux sums.

## What the model now says (Figure G2)

**A–B. Required AOX, refined.** With fitted maintenance the required
alternative-respiration shares are 0 / 14 / 0 / 25 / 4 / 42 % of O₂ uptake
(15/0 … 27/2), unique by FVA, and the required AOX flux tracks measured AOX
transcript at ρ = 0.99 (p < 0.001). Maintenance-assumption whiskers span the
earlier NGAM scan; the ordering never changes. The falsifiable respirometry
prediction stands, now with an uncertainty band: drug conditions and
supra-optimal heat require alternative respiration, cool controls require none.

**C. A continuous uncoupling axis.** The effective P/O ratio (ATP synthase
flux per 2e⁻ to O₂) collapses monotonically from 1.28 (15/0) to 0.07 (27/2).
This is the model's restatement of the CUE collapse in bioenergetic units, and
a second experimentally checkable quantity (ADP/O respirometry).

**D. An unprompted convergence.** To satisfy the measured physiology the
parsimonious flux state spontaneously engages pentose-phosphate + glutathione
cycling in the stressed conditions (glutathione-pathway flux 0 → 8
mmol gDW⁻¹ h⁻¹), and this model-required redox cycling correlates with the
measured glutathione module score (ρ = 0.94, p = 0.005). Nothing in the
constraints mentions glutathione; the network chose it. This independently
corroborates the transcriptome's glutathione arm as part of the energetic
stress state.

**E. An informative failure.** E-Flux models, where the transcriptome alone
sets reaction capacities, predict growth across the held-out conditions with
r ≈ −0.2 to −0.3 (worse than useless). Conclusion: the drug's growth effect is
not encoded in transcript abundances — it is post-transcriptional (direct
CYP51 enzymatic inhibition) and the uncoupled state is regulated rather than
flux-optimal. This is a quantitative argument for the paper's mechanistic
framing (enzyme-level inhibition + programmed stress response) and against
purely transcriptional readouts of drug action.

**F. Carbon budget closes** in all six conditions (consistency check under
RQ = 1, as before).

## Remaining limits (and why they are structural, not fixable here)

Plain stoichiometric FBA cannot make equal CYP51 inhibition temperature-
dependent (linear response at any temperature) and cannot generate the
uncoupled state from optimality (it must be imposed via measured physiology or
regulation). Both point the same way: the temperature dependence lives in
enzyme kinetics and regulation, i.e. the etc-GEM layer (temperature-dependent
kcat/stability on top of this GEM). The scaffold is ready for that: Zt-native
GPRs, AOX explicit, maintenance fitted, condition constraint sets defined.
Empirical anchor for a temperature-dependent CYP51 capacity term:
doi:10.1128/aem.03965-14 (M. graminicola CYP51 is temperature-sensitive).

## Files

- `advance_ztGEM_v04.py` — all v0.4 analyses (reproducible)
- `ztGEM_v04_required_AOX_fitted_maintenance.csv`
- `ztGEM_v04_condition_detail.csv` (P/O, ATP budget, CYP51 flux, per condition)
- `ztGEM_v04_subsystem_fluxsums.csv`
- `ztGEM_v04_eflux_predictions.csv` (the honest negative)
- `ztGEM_v04_maintenance_fit.json`, `maintenance_fit_points.csv`
- `FigureG2_ztGEM_v04_main.png` — main-figure candidate

---

# Deep dive into the E-Flux negative (point 6) — added same day

## (a) The failure is robust
Across E-Flux design choices (capping percentile 90/95/99, cap clipping,
linear expression), held-out prediction of growth stays at r = −0.10 to −0.29.
Using log-scale expression flips it to a still-useless r = +0.21. The negative
is not an artefact of one parameterisation (`eflux_robustness.csv`).

## (b–c) Why it fails, precisely
E-Flux overpredicts growth everywhere except the calibration condition, and
worst exactly where the drug acts (5.9× at 15/2). Diagnosing which
expression caps actually bind at the optimum: essentially ONE reaction ever
binds — ATP synthase (r_0226; its GPR contains mitochondrially-encoded
subunits, which sit at low apparent expression in a poly(A) RNA-seq library,
an artefact worth noting). Every other transcript-derived cap is slack: the
network routes around them. Two conclusions: (i) transcript abundance carries
almost no binding information about metabolic capacity in this dataset;
(ii) the processes that actually limit growth under drug — ribosome/translation
shutdown and direct enzymatic CYP51 inhibition — are respectively OUTSIDE
metabolic-model scope and INVISIBLE to transcript-based capacities. The
transcriptome describes the response; it does not set the limits
(`eflux_error_decomposition.csv`).

## (d) Hybrid rescue: one post-transcriptional fact fixes the cool conditions —
## and exposes an excess effect at 27 °C
Injecting the single enzymatic fact (CYP51 activity retained at 2 mg/L, from
the Hill fits: 43% / 50% / 86% at 15/21/27 °C) as a cap on the CYP51 reaction
in the drug-condition E-Flux models predicts drug/control growth ratios of:

  15 °C: predicted 0.431 vs measured 0.440   (near-exact)
  21 °C: predicted 0.504 vs measured 0.459   (close)
  27 °C: predicted 0.725 vs measured 0.545   (overpredicts retention)

Reading: at 15–21 °C, target-level inhibition alone quantitatively accounts
for the entire growth effect. At 27 °C, measured suppression exceeds what
target engagement explains — an excess consistent with the combined
heat+drug energetic burden (the uncoupled stress state) rather than CYP51
inhibition per se. Caveat: the 27 °C Hill parameters are wide (EC50 CrI
3.2–11.2), so the input retention spans roughly 0.72–0.97 and the gap is
suggestive rather than conclusive. Also flagged: the pipeline's own two
estimates of the 27 °C drug effect at 2 mg/L disagree (Hill-curve retention
0.86 vs effect-size-model retention 0.55) — worth reconciling in the
physiology pipeline independently of any modelling
(`eflux_hybrid_cyp51.csv`).

## Literature status for Z. tritici GEMs (full sweep)
No published, dedicated Z. tritici GEM exists. The only species-specific
reconstruction indexed anywhere is the 2014 CoReCo scaffold we already
assessed (BioModels MODEL1302010031). The CarveFungi preprint
(bioRxiv 2023.08.23.554328; deep-learning localisation + CarveMe carving)
generated models for 834 fungi — species list not stated in the paper; the
2.6 GB Zenodo archive (10.5281/zenodo.7413265) may contain one and can be
checked on the Mac (Zenodo unreachable from this sandbox). The CarveFungi
pipeline itself could also be run directly on IPO323 using the eggNOG
annotation you already have, as an independent second reconstruction to
cross-validate ztGEM. Bottom line: ztGEM v0.3/0.4 is, to the best of a
thorough search, the first working curated GEM for this pathogen.

---

# Breaking the ceiling: a CYP51 safety-margin temperature layer (v0.5 analysis)

The "plain FBA cannot make CYP51 inhibition temperature-dependent" limit was
broken with ONE added concept: an ABSOLUTE, temperature-dependent CYP51
capacity, instead of proportional flux caps. The mechanism ("safety margin"):

  - CYP51 capacity follows enzyme kinetics: Arrhenius, fitted Ea = 0.47 eV
    (textbook enzymatic range). Empirical anchor for heat lability of this
    species' CYP51: aem.03965-14 (active at 22 C, no activity at 37 C).
  - Ergosterol demand scales with measured control growth mu(T).
  - Potency is set by the margin sigma(T) = capacity/demand:
    cold slows the enzyme faster than demand falls -> small margin -> the drug
    bites hard (synergy); supra-optimal heat cuts demand faster than capacity
    -> large margin -> the drug is tolerated (antagonism).

Model comparison on the seven measured EC50s (equal 2-parameter complexity):
margin model AICc 11.05 < van't Hoff binding-only 13.08; statistically tied
with the no-temperature-dependence null (10.85) given n = 7, so this is a
mechanistic best-candidate, not proof. Crucially it predicts the right ORDER
(most sensitive cold, least at 27 C) with a plausible Ea.

Implemented inside ztGEM (absolute Vmax(T) on r_0317, competitive inhibition
1/(1+D/Ki), growth ceiling = measured control mu(T)), the network now returns
temperature-dependent EC50s:

  15 C: GEM 1.93 vs observed 1.12
  21 C: GEM 2.05 vs observed 2.04   (calibration temperature)
  27 C: GEM 4.91 vs observed 5.00

The GEM's 27/15 potency ratio is 2.5x vs the observed 4.5x: enzyme kinetics +
demand explains roughly 60% of the interaction in log units, with the residual
(mostly the extra cold-side synergy) left for binding thermodynamics and the
stress-state protection documented by the transcriptome. Together with the
hybrid E-Flux result (target inhibition fully explains the drug effect at
15-21 C; excess suppression at 27 C), the model now provides a quantitative
decomposition of the temperature x fungicide interaction into enzyme-kinetic,
demand-side, and stress-state components — the question the review posed in
section 6.3.

Files: cyp51_margin_model.py, cyp51_margin_fit.csv, cyp51_margin_params.json,
cyp51_margin_reduced_fit.json, ztGEM_cyp51_temperature_ec50.csv,
FigureG3_cyp51_margin.png.

---

# Hardening pass (v0.4.1): Monte-Carlo uncertainty + full GECKO expansion

## Monte-Carlo uncertainty propagation
120 draws per condition sampling every measured and assumed input (mu and
respiration posteriors, drug-effect CIs, biomass carbon fraction 0.48 +/- 0.03,
maintenance 79.2 +/- 12). All 720 draws feasible. The headline claims survive
intact (`mc_uncertainty_quantiles.csv`; Figure G2 now carries these 95% CIs):

  - Required AOX share: cool controls are 0% in EVERY draw; 27/2 requires
    30-47% in every draw; 21/2 requires 2-38%; 15/2 spans 0-35% (median 13%,
    the one condition whose CI touches zero — honest caveat).
  - Effective P/O: ordered collapse preserved (15/0 median 1.28 -> 27/2
    median 0.07) with non-overlapping cool-control vs hot-drug bands.
  - Glutathione-pathway flux rises monotonically across the stress gradient
    in the medians; model CUE bands overlap measured CUE in all conditions.

## Full GECKO expansion (ecZtGEM_full.xml)
The one-representative-enzyme simplification is gone: 2,418 enzymatic
reactions expanded into 3,504 isozyme/direction arms with complex (AND)
subunit costs, 742 enzymes with real MWs, per-gene draw reactions. Verified
end-to-end through `etcgem build` + `etcgem tpc`: same emergent shape
(Topt 28 C, rising-limb Ea 0.59 eV) and the uncalibrated magnitude moved from
~9x low (lite) to ~4x high (full) — bracketing the measured rmax, as expected
for placeholder kcats, and exactly what DLTKcat + calibration will pin down.
The strain folder now ships the full model (`strain.yaml` updated; lite kept
alongside for speed).

## Where this leaves the modelling
Discovery, robustness, and pipeline-integration are complete. Remaining steps
are external by nature: TemStaPro Tm table, DLTKcat kcats, Gab's calibration,
and (optionally) the CarveFungi cross-check — all documented in
zt_ipo323/STATUS.md.

---

# Loose-end adjudication: the 27 C drug-effect discrepancy (resolved to a diagnosis)

Raw dose-response data (fig1_extra_dose_response_by_temp.csv, growth, n=3 per
cell) give drug retention at 2 mg/L of 0.36 (15 C), 0.51 (21 C), 0.80 (27 C).
Comparison:

  15 C: raw 0.36 | Hill 0.43 | effect-size model 0.44   (all agree)
  21 C: raw 0.51 | Hill 0.50 | effect-size model 0.46   (all agree)
  27 C: raw 0.80 | Hill 0.86 | effect-size model 0.55   (effect-size diverges)

Diagnosis: the Hill fit is consistent with the raw rates it was fitted to; the
effect-size model (built on the biomass-corrected, carbon-specific rates from
the TPC pipeline stage) is the outlier at 27 C only. Either the biomass
correction legitimately changes the 27 C drug effect (biomass differences
under drug at supra-optimal temperature), or a processing step overstates it.
THIS IS THE ITEM TO RECONCILE IN THE PHYSIOLOGY PIPELINE (same family as the
156-vs-157-row dataset question).

Implication for the modelling conclusions, stated honestly: the safety-margin
GEM predicted 27 C retention of 0.73 and we previously called the gap to the
effect-size value (0.55) "excess suppression beyond target engagement". If the
raw/Hill value (0.80-0.86) is the right benchmark, that gap largely closes and
target-level engagement + the kinetic margin explains the 27 C drug effect
almost fully. Which benchmark is right depends on the biomass-correction
question above. The required-AOX analysis used the effect-size condition
values; its cool-control (0%) and 27/2-highest ordering is insensitive to
this, but the exact 27/2 percentage would shift if the effect-size layer is
revised. Re-run advance_ztGEM_v04.py after any pipeline fix — everything
regenerates from the tables.

---

# RESOLUTION: the "Hill vs effect-size" 27 C question (audit rerun, complete)

1. REPRODUCED. The effect-size value (0.545 retention at 27 C) is exactly what
   the per-dose Sharpe-Schoolfield curves give when evaluated at 27 C
   (recomputed from the posterior draws: 0.547). Not a bug — a definition.
2. BIOMASS CORRECTION EXONERATED. At replicate level, raw and
   biomass-corrected growth retentions are IDENTICAL at every temperature
   (the correction cancels in the dose ratio). The 156-vs-157-row question is
   unrelated to this discrepancy.
3. ROOT CAUSE. The SS effect-size estimator pools ACROSS temperatures within
   each dose (global lnB0/E/Eh/Th per dose), so it structurally smooths away
   any temperature-localised change in the drug effect — its retention ramps
   gently 0.44 -> 0.62 across 15-28 C by construction. The local 27 C data
   (0.80, n=3/3) and the Hill fit (0.86, which pools across doses WITHIN
   27 C) agree with each other and are the right estimators for
   temperature-specific drug effects.
4. RECOMMENDATION (pipeline): keep SS effect sizes for global dose trends
   only; derive temperature-specific fungicide effects from the Hill model or
   local contrasts. The manuscript's temperature-specific claims are already
   Hill-based, so no published number changes. No refit required.
5. GEM SENSITIVITY. Recomputing required AOX with LOCAL drug effects instead
   of SS-based ones: drug conditions need 20/22/23% of O2 through AOX at
   15/21/27 C (vs 14/25/42% under SS values); cool controls remain 0% under
   both. The categorical prediction (drug and heat require alternative
   respiration, unstressed cool controls do not) is estimator-robust; the
   steepness of the gradient across drug conditions depends on the estimator
   and should be quoted with that caveat.

---

# LAYER-1 GROUNDING COMPLETE: TemStaPro + DLTKcat (2026-08-08)

## TemStaPro melting temperatures -> unfolding mode

ProtT5 embeddings + TemStaPro thresholds for all 742 model enzymes (7.5 h CPU,
IPO323 MacBook). Pseudo-Tm range 35.1-64.3 C, mean 37.6 C; converted to
Topt/Tm/dCpt via temstapro_to_bestparams_zt.py (center 24 C, slope 0.5; no
Z. tritici meltome exists for homology calibration, flagged in TmTag).
Switched strain.yaml to thermal_model: unfolding with these params:

- Emergent Topt: 28.0 -> 24.0 C — exactly the measured growth optimum,
  with NOTHING fitted to the growth curve.
- CTmax: off the 35 C grid edge to 36.8 C on a 5-45 C grid — a real
  emergent collapse (growth cliff at 35.5-37 C) driven by the most
  thermolabile essential enzymes.
- Remaining shape gap: measured growth loses ~37% between 24 and 27 C,
  the model only ~4% — the falling limb starts too late. This is the main
  calibration target (and/or a Tm offset: TemStaPro pseudo-Tm is an
  in-vitro-like stability proxy; in-vivo functional loss starts earlier).

## DLTKcat temperature-dependent kcats

DLTKcat cannot resolve Mycgr3 genes (UniProt lookup) or Yeast9 metabolite
names (PubChem lookup), so the run was done fully offline
(dltkcat/run_dltkcat_zt.py): sequences from the IPO323 proteome; SMILES via
ModelSEED compounds (BiGG alias -> MetaNetX alias -> KEGG annotation ->
exact-name matching), with per-reaction fallback to an alternative real
substrate when the chosen one had no structure. Coverage: 2,356/3,330 arm
reactions (71%); 25,916 predictions on a 5-45 C grid; median kcat ~3.7 s^-1;
62% of reactions predicted kcat rising with T (median +0.08 log10 over 40 C —
DLTKcat's known weak temperature signal). CYP51 (r_0317) kcat declines with
temperature, consistent with the safety-margin picture.

etcgem dltkcat parse: 77/2,356 MMRT fits passed the built-in ok filter
(peaked, r2>0.8). Of these, 17 had fitted Topt ABOVE the enzyme's TemStaPro
Tm — thermodynamically impossible, and in unfolding mode they poison
rel_kcat (dHt blows up when Topt > Tm): the whole-model TPC collapsed to
rmax 0.000 (growth pinned at ~7e-5 /h at all temperatures, Ea = nan).
Filter applied (Topt in [2,45] C and Topt <= Tm - 2 C): 60 clean fits,
TPC restored and unchanged (Topt 24.0 / rmax 0.216 / CTmax 36.8).

>> ACTION FOR GAB: apply_fits_to_provider() should validate fitted Topt
   against the enzyme's Tm before overriding unfolding-mode params — one
   unphysical fit on an essential reaction silently kills the whole curve.

Note on expectations: apply_fits deliberately preserves kcat magnitudes
(only Topt/dCp/dCpt are overridden), so DLTKcat does NOT move rmax — the
~3.4x rmax excess vs the measured 0.063 /h is the job of `etcgem calibrate`.

## Strain state

zt_ipo323 is now fully Layer-1 grounded: real Tm for 3,504/3,504 arms (100%),
DLTKcat thermal fits where defensible (60 arms), full GECKO expansion with
real MWs, measured TPC as validation target, fitted maintenance. Emergent-
curve scorecard vs measurement: Topt exact (24 C), rising-limb Ea 0.52 vs
measured 0.36 [0.13, 0.58] eV (consistent), CTmax 36.8 vs literature ~33-34 C
(a few degrees warm), rmax 3.4x high (uncalibrated by design), falling limb
too flat (main gap). Ready for the calibration stage — Gab's call on settings.

---

# CarveFungi CHECK RESOLVED (2026-08-08): drafts EXIST — claim reworded

The Zenodo archive was grepped end-to-end on the Mac (first attempt truncated
at 83% — disk full — and stopped alphabetically at P, i.e. just before the
Z section; a re-run with a streaming pipe settled it). CarveFungi DOES
contain Z. tritici: six strain assemblies (IPO323 as both GCF_000219625.1
and the MG2/Ensembl annotation, ST99CH-1A5/-1E4/-3D1/-3D7) plus Z. brevis,
~11 ensemble variants each (74 SBMLs, now in gem/carvefungi_zt/).

Head-to-head, CarveFungi IPO323 (main ensemble model) vs ztGEM v0.3:

| | CarveFungi IPO323 | ztGEM v0.3 |
|---|---|---|
| reactions / metabolites | 2,191 / 1,695 | 3,959 / 2,722 |
| genes | 553 | 1,890 (3.4x; only 342 of theirs overlap ours) |
| growth prediction | 1.08 /h default medium, 3.25 /h + sucrose (absurd; uncalibrated) | matches measured 0.052 /h after condition constraints |
| CYP51 (prothioconazole target) | ABSENT (only a 14-demethyllanosterol transport; no sterol 14a-demethylase reaction) | curated, native GPR (Mycgr3G110231 + CPR) |
| alternative oxidase | ABSENT | r_AOX, transcriptome-supported (rho=0.99) |
| enzyme constraints / temperature | none | GECKO + TemStaPro Tm + DLTKcat kcat(T) |
| curation / validation | automated CarveMe carving, none | gap-filled, validated vs measured TPC |

The decisive point for the manuscript: the CarveFungi drafts cannot
represent EITHER half of the temperature x prothioconazole mechanism —
no drug target, no AOX. They are genome-coverage drafts, not analysis-ready
models, and nothing about Z. tritici was ever analysed or published from
them.

REVISED NOVELTY CLAIM (use this wording): "the first curated,
enzyme- and temperature-constrained genome-scale metabolic model of
Z. tritici, validated against experimental physiology. Automated draft
reconstructions of Z. tritici exist within the 834-species CarveFungi
collection (Castillo et al., bioRxiv 2023.08.23.554328; Zenodo 7413265)
but are uncurated, carry no enzyme or thermal constraints, and lack both
the azole target CYP51 and the alternative oxidase."
Citing and comparing against them strengthens, not weakens, the paper.

---

# PHASE-1 CALIBRATION COMPLETE (2026-08-08, emcee on the IPO323 MacBook)

`etcgem calibrate --strain zt_ipo323 --curve zt12129_control` — 24 walkers x
1000 steps (300 burn), 4 global parameters (dTopt, dTm, dCp_scale,
kappa_scale) + noise, fitted to the 7-point measured control TPC.
acceptance 0.31, n_eff 185. A second, independent 500-step chain run in the
cloud sandbox (checkpointed) agreed in direction on all parameters.

Posterior medians [90% CI]: kappa_scale x1.10 [0.32, 3.16];
dTm -4.7 K [-11.6, 9.0]; dTopt +11.3 K [-11.6, 11.8];
dCp_scale x1.11 [0.62, 2.23]. ALL flagged constrained=False — the single
7-point growth curve under-determines the four global knobs (expected for
phase 1; the parameters trade off along a ridge).

Descriptors (observed | emergent prior | posterior median):
rmax 0.063 | 0.216 | 0.123 /h — the kappa/pool ridge halves the gap.
CTmax <=28 (censored: data end at 28 C) | 37.4 | 32.5 C — falling limb
sharpens, driven by the negative dTm.
Topt 24.0 | 24.0 | 29.1 C — degrades as the cost of the sharper limb;
Ea 0.43 | 0.53 | 1.28 eV — cold limb steepens correspondingly.

Honest reading: phase-1 calibration confirms the two headline gaps
(magnitude and falling-limb sharpness) are correctable within physically
plausible parameter shifts (a ~1.1x in-vivo efficiency factor and a ~5 K
Tm offset), but a single aggregate growth curve cannot constrain all four
directions at once. The right next step is Gab's multi-parameter stage
(the Van-Derlinden-style unified tuning with sigma/f_metab freed), ideally
adding the drug conditions as additional curves — six conditions would
break the ridge that one curve cannot.

Figure 6 of the manuscript set now shows: (A) emergent vs calibrated TPC
with the 90% posterior band over the measured points; (B) required AOX
share per condition (Monte-Carlo 95% bands, Spearman rho = 0.99 vs AOX
transcript); (C) EC50(T) — observed Hill vs the CYP51 safety-margin model
and the ecFBA implementation; (D) the predicted effective P/O collapse
under drug + heat. Scripts: figures_final/ (fig_style.R/.py + per-figure
scripts, R and matplotlib, one shared style system).

## Final calibration (tight dTopt prior, Mac run, 2026-08-08)

Settings: 24 walkers x 1,000 steps (300 burn), dTopt prior sigma = 1 K
(strain.yaml calibration.priors), acceptance 0.447, max autocorr 73.8,
n_eff 227.8, wall 3 h 45 m.

| parameter | posterior median | 90% CI | curve-constrained |
|---|---|---|---|
| dTopt | +1.1 K | [-0.8, +2.6] | no (prior-dominated; ridge artifact gone) |
| dTm | **-14.1 K** | [-15.0, -11.0] | **yes** (abuts -15 K bound) |
| dCp_scale | x1.50 | [0.86, 2.52] | no |
| kappa_scale | x0.92 | [0.31, 2.99] | no |
| sigma_noise | 0.058 /h | [0.042, 0.077] | — |

Posterior-median descriptors: rmax 0.185 (obs 0.063), Topt 20.5 C (obs 24.0),
CTmax 23.3 C (obs 28.0), Ea 0.92 eV (obs 0.43).

Verdict: the calibration is an identifiability analysis, not a better model.
With the optimum pinned, the only curve-constrained parameter is dTm, driven
to its bound (in-vivo thermal failure ~14 K below sequence-predicted Tm —
biologically sensible), but the posterior overcorrects the falling limb
(median collapse at ~22-23 C, below the data) and cannot close the ~3x rmax
gap because kappa is prior-bound. Decision taken in the manuscript: Figure 6A
shows the uncalibrated emergent model only (its parameter-free Topt match is
the result); the calibration is shown honestly as Figure S11D; all downstream
analyses (condition-constrained pFBA, EC50 scan) use the uncalibrated model
and are calibration-independent. Superseded phase-1 numbers (dTopt +11.3 K
ridge) are retained in the outputs history but no longer cited anywhere.
