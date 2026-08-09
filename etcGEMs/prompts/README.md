# Prompts — etcGEM project

Each file here is a prompt to paste into **Claude Code** (run from the project root
`.../MICROADAPT/etcGEMs`). Most are autonomous and modify code, run analyses, and
commit; launch Claude Code in an auto-approving mode (accept edits + allow commands,
e.g. `--dangerously-skip-permissions`) to run them unattended. They build on each
other, so run them roughly in the order below.

## Layout

- `prompts/repo_restructure_prompt.md` — the current / most-recently run task (top level).
- `prompts/archive/` — **all executed historical prompts**, kept for provenance. The index
  below refers to files now under `archive/` (run roughly in that order to reproduce the build).

Status legend: ✅ run · ▶ in progress · ⏳ pending  (all archived prompts have been executed;
the ⏳ markers below are historical and no longer accurate)

## Recommended run order

**Structure & packaging**
1. `reorg_prompt.md` — ✅ reorganise to a `src/etcgem` package, per-strain folders.
   *(Superseded by #2 for the config model; kept for history.)*
2. `project_restructure_prompt.md` — ✅ config model = `defaults / strain / experiment`
   + two-tier CLI (`build/tpc/fba` strain-only; `sweep/decompose/control` + experiment).

**Core analyses**
3. `decomposition_prompt.md` — ✅ allocation-vs-envelope Shapley variance decomposition (H1.3).
4. `control_identifiability_prompt.md` — ✅ per-enzyme thermal control coefficients + identifiability (H1.1/H1.2).
5. `allocation_and_sampling_prompt.md` — ✅ proteome-sector allocation + correlated / DLTKcat-posterior thermal sampling (M1.2).

**Data (DLTKcat) + figures + runs**
6. `dltkcat_sweep_prompt.md` — ✅ run DLTKcat → per-enzyme MMRT params → eciML1515 sweep.
7. `native_figures_and_runs_prompt.md` — ✅ native descriptor-interval & sector-tradeoff figures; run the `calibrated` & `sectors` experiments with plots.

**Fixes & model realism**
8. `fixes_identifiability_ea_sectordecomp_prompt.md` — ✅ proteome-wide identifiability, dCp→Ea calibration, sector-based H1.3 decomposition.
9. `per_enzyme_dcp_decouple_prompt.md` — ✅ per-enzyme thermal heterogeneity (dCp/Topt) to decouple breadth & Ea.
   *(Interim fix; #12 is the more principled, data-grounded route.)*

**Reporting**
10. `quarto_report_prompt.md` — ✅ set up the Quarto report (PDF/docx), `reports/etcgem/`.
11. `flesh_out_report_prompt.md` — ✅ full write-up: Background → Model → Design → Results → Next steps.
12. `report_addendum_structural_caveat_prompt.md` — ✅ add the "clean separation is partly
    structural / growth-fit can't resolve it / temperature-dependent allocation" caveat
    + refs (Mairet2021, Wang2026). *(Targeted edit to the fleshed-out report.)*

**Data-grounded model**
13. `align_with_mres_unfolding_model_prompt.md` — ✅ two-state **Tm-based unfolding** thermal
    model from the MRes/Li work; grounds the *envelope* in per-enzyme Topt/Tm (reuses the
    MRes repo's parameter tables — no BRENDA/TOMER needed for E. coli).
14. `ingest_temperature_proteomics_prompt.md` — ✅ ingest E. coli temperature proteomics
    (DeyuWang `tem_proteomic.csv`) to ground **temperature-dependent proteome allocation**
    and validate predicted vs measured enzyme usage. (Measured chaperone sector ramps
    2.7× by 43 °C; wired in as an opt-in run + validation.)

**Consolidation & report**
15. `complete_model_consolidation_and_report_prompt.md` — ⏳ assemble ONE **complete**
    data-grounded model (unfolding envelope + measured temperature-dependent sector
    allocation + DLTKcat), validate it against the empirical E. coli TPC and the proteome,
    re-run sensitivity / decomposition / control **on the complete model**, and
    **restructure the report** (complete model → validation → results; incremental
    build-up moved to a supplementary construction section). **Run after #13 and #14.**

**Fully emergent, per-curve validation (current frontier)**
16. `emergent_tpc_per_curve_validation_prompt.md` — ⏳ make the TPC a genuine PREDICTION:
    remove the growth-fit knobs (`pool_scale`, `calibrate-dcp`); make **magnitude** emergent
    (pool = P_total × f_metab × σ, from proteome + literature σ) and **Ea** emergent (grounded
    per-enzyme dCp), set the **medium per curve** from the Smith 2019 database
    (`smithtp/hotterbetterprokaryotes`), and validate against **individual** E. coli curves
    (many per-curve fits) + the ~0.85 eV growth-Ea benchmark. **Growth only** at this stage.
    **Run after #15.**

17. `rich_medium_and_percurve_medium_prompt.md` — ⏳ add an **LB / rich medium**
    (availability, from the MRes `media.csv` / Machado 2018) so the 24 rich-broth E. coli
    curves are predicted under the correct condition; re-validate per-curve on **RAW
    ABSOLUTE rates (1/h)** — magnitude is emergent, so drop peak-normalisation as the primary
    metric (shape kept as secondary); **highlight minimal vs LB panels**; switch the ablation
    plot to absolute too; document the medium treatment. (Only 2/26 curves were
    defined-minimal before; the rest were shape-only.) **Run after #16.**
18. `medium_matched_sector_allocation_prompt.md` — ⏳ make the proteome-sector fractions
    **medium-matched** (LB vs Glucose vs Glycerol, all from the DeyuWang proteomics) so the
    ribosome/translation cap is medium-dependent and LB grows faster than minimal (the
    growth-law effect, data-grounded not fitted). Fixes the #17 finding that rmax was stuck
    at 0.34/h on every medium because `f_bio` was fixed; re-runs the per-curve absolute
    validation. **Also upgrades the report** to a reproducible, equation-backed write-up:
    model **equations** (MMRT, unfolding, sector caps, decomposition), a **parameter-
    provenance table** (source + coverage for every parameter), and detailed in-silico-
    experiment and validation descriptions. **Run after #17 finishes.**

**Analysis upgrades**
19. `equal_perturbation_sensitivity_prompt.md` — ⏳ replace the sensitivity analysis's
    hand-set, unequal parameter ranges + rank correlation with an **equal-perturbation
    (standardised elasticity)** analysis: move every input by the same standardised step and
    report a dimensionless, comparable magnitude (elasticity) per TPC descriptor, so the
    ranking reflects the **model's structural leverage** ("which parts of the model drive the
    TPC") rather than our range choices. Handles the additive-vs-multiplicative units issue
    (reference scale for the `dTopt` shift). The complementary **uncertainty-weighted** version
    (ranges = real measurement uncertainty, to prioritise what to *measure*) is a deliberate
    follow-up, not in this prompt. Rewrites the report's sensitivity section. **Run after #18.**
20. `recast_decomposition_with_tm_prompt.md` — ⏳ reconcile the two conflicting analyses. Adds a
    **Tm perturbation knob** (`dTm`) so the unfolding/high-T side is finally tested — Tm was
    never perturbed, which is why `CTmax` looked "insensitive to everything". Recasts the
    Shapley **decomposition** to match the elasticity: equal, standardised, **nominal-centred**
    ranges (fixes the envelope-favouring wide ranges + the off-nominal `f_metab` sampling),
    reports **variance shares AND magnitude**, keeps the crossed design for the interaction
    term, and adds Tm to the envelope (split kinetic vs stability). Centerpiece: **median+IQR
    absolute-TPC band plots** per axis (allocation→height, kinetic envelope→cold side/Ea,
    Tm→hot side/CTmax). Rewrites both report sections into **one consistent narrative**,
    correcting the old "envelope dominates the rate" claim. **Run after #19.** ✅ run.

**Report cleanup**
21. `report_clarity_single_operating_point_prompt.md` — ⏳ report-quality pass after #20:
    enforce **one operating point** across the whole report (the sector-trade-off subsection
    still cites the old `f_metab≈0.285`; move everything to the strain nominal `≈0.483` and
    audit for other stale numbers); make **Tm appear in the sensitivity figures** (it's in the
    table but not the heatmap/tornado); **remove the incremental, self-referential tone** that
    narrates corrections of our own earlier analyses; and do a **whole-report clarity rewrite**
    — define every term/acronym in full, lead each section with a plain-language summary, write
    self-contained prose for a reader who hasn't followed the work (per the collaborator's
    feedback). No model-math changes. **Run after #20.**
22. `report_model_anatomy_and_global_knobs_prompt.md` — ⏳ add a **"how the model works"**
    section right after the mathematical model description (before validation): the **reference
    TPC** at the glucose-minimal
    strain nominal, **density plots** of the per-enzyme parameter distributions (Topt/Tm/dCp
    across ~2,560 enzymes; their SDs are the elasticity reference scales), and an **example
    per-enzyme kcat(T) panel** (~10 enzymes) so the reader gets a feel for the model. Also makes
    the methodology explicit: the sensitivity/decomposition knobs are **global scalars** that
    move the whole enzyme distribution (not per-enzyme free parameters), and those analyses run
    on the **single glucose-minimal reference TPC**, not per empirical curve. Plotting + report
    only. **Run after #21.**

**Bayesian calibration (emergent vs what-the-data-demand)**
23. `bayesian_calibration_phase0_1_prompt.md` — ⏳ Phase 0+1 of the Bayesian analysis: the
    emergent model is the **prior**; ask the inverse question on **one** glucose-minimal curve —
    what would the borrowed/uncertain parameters have to be to fit the data? Adds a
    `kappa_scale` (translation-efficiency) magnitude knob (the elasticity-correct lever for
    rmax), frees {kappa_scale, dCp_scale, dTopt, dTm} with provenance-based priors, and does
    **exact** likelihood inference with a gradient-free sampler (**emcee**; not ABC, not Stan/
    brms — the FBA simulator is non-differentiable). Outputs the **prior-vs-posterior TPC**, a
    **demanded-corrections** table, and a **corner/degeneracy** plot. No report changes (write-up
    deferred until all phases checked). Later phases: 1b per-curve-separate, 2 cross-check vs
    independent data, 3 posterior→sensitivity, 4 joint/per-enzyme ABC. **Run after the reporting
    prompts settle.**

**Trusted-data rebuild**
24. `redo_report_with_trusted_curves_prompt.md` — ⏳ rebuild all report analyses on the two
    trusted, **strain-matched** curves and retire the low-confidence Smith compilation:
    **Noll 2023** (K-12 NCM3722, defined glucose-minimal, with per-T SD) as the **minimal**
    curve and **Erdos 2026** (MG1655 wt, LB) as the **rich** curve. Re-runs the **absolute-rate
    validation** (minimal vs rich magnitude contrast), regenerates the model-intrinsic analyses
    (anatomy/elasticity/decomposition/control/ablation, unchanged results), and rewrites the
    report's validation + data-provenance sections (documenting the ~8× Cooper error and the
    digitization caveats). **Excludes the Bayesian tuning** (redone later). **Run after the data
    corrections + reporting prompts.**

25. `bayesian_calibration_noll_minimal_prompt.md` — ⏳ redo the Bayesian calibration (Phase 0+1)
    on the **trusted Noll NCM3722 glucose-minimal curve** instead of the mis-scaled Cooper data.
    Repoints the calibration at the `sources/` Noll CSV, upgrades the likelihood to use the
    **measured per-temperature SD** (`rate_sd_i² + sigma_disc²`), and runs **longer emcee chains**.
    Expect the magnitude correction to **flip** (kappa_scale > 1: emergent ~0.55 under-predicts
    Noll ~1.0 h⁻¹). Outputs to `calibration_noll_minimal/`; no report edits. **Run after #24.**

**MG1655 rebuild (Van Derlinden, rich medium) — supersedes the Noll/Erdos line**
P1. `P1_rich_bhi_reference_and_vanderlinden_validation_prompt.md` — ⏳ switch the reference
    operating point to a **rich medium (BHI)** and validate the **emergent** model against the
    one exact-strain curve, **Van Derlinden 2012** (MG1655, BHI, 7–46 °C). Wires `"BHI"` into
    `set_medium` (availability, curated `BHI_media.csv`, no pinned uptakes), sets the rich
    reference + regenerates anatomy, drops Noll from validation. No Bayesian, no sensitivity/
    decomposition re-run (those are P3, on the tuned model). Documents the BHI-as-LB/rich-proxy
    (even the Rothia GEM didn't encode BHI in silico) and the rich-is-easier caveat. **First of
    the P1→P2→P3 sequence.**
P1b. `P1b_model_realism_maintenance_and_growthlaw_prompt.md` — ⏳ **model-realism fixes** before
    re-tuning, from the Li/MRes structural comparison. (1) Restores **T-dependent maintenance under
    the sector model** (it was gated off — the cause of the flat top; Li/MRes run NGAM(T) always)
    → rounds the peak; adds `ngam_scale`/`ngam_steepness` knobs. (2) Adds a **growth-rate-dependent
    proteome partition** — proteome-conserving `f_bio(μ)=f_bio_0+slope·μ` / `f_metab(μ)=f_metab_0−slope·μ`,
    both linear (biosynthesis `κ_eff=κ−slope·P_total`; metabolic pool gains `+slope·B·v_bio`) — the
    fix for the ~1.5× magnitude shortfall (Li/MRes have no biosynthesis cap; ours binds at the rich
    peak). The ribosome↔metabolism trade-off self-limits growth correctly. (3) Adds `tm_scale`. Verifies on the **emergent** curve (rounded peak +
    lifted rich rmax); **no re-tune**. **Run before re-running P2.**
P1c. `P1c_reconcile_proteome_pools_prompt.md` — ⏳ **fixes the ceiling P2c named**. The ~1.55 cap was
    NOT structural — it's a redundant enzyme-mass double-accounting: our sMOMENT sector pool *and* the
    leftover GECKO base `prot_pool_exchange` (0.091 g/gDW) both bind, and `kcat_scale` only touches the
    former. Reconciles them to a **single proteome budget** (the sector/growth-law pool) so the
    magnitude levers act on the true binding pool. Reports how much of the gap is the fixable
    redundancy vs an above-literature σ (~0.77). No re-tune, no report rewrite. **Run before the P2
    re-run.** (Then re-run P2 on the reconciled model.)
P2. `P2_bayesian_tuning_vanderlinden_prompt.md` — ⏳ Bayesian tuning on Van Derlinden at the rich
    (BHI) operating point. Adds the **`kcat_scale`** knob (global metabolic-turnover magnitude
    lever, replaces the degenerate `budget_scale`), fits the **unified param set** (envelope +
    {kcat_scale, kappa_scale} + allocation) with **provenance priors** (measured allocation tight;
    magnitude/envelope broad) and a free discrepancy term (Van Derlinden has no SD). Van Derlinden's
    7–46 °C span constrains the envelope Noll couldn't. Headline: does it reach ~2.4, is the demanded
    `kcat_scale` physically plausible (pool-vs-kcat), does `dTm` fix the ~4 °C-too-hot CTmax. No
    report edits. **Inspect after this.** **Run after P1.**
P2b. `P2b_report_restructure_prompt.md` — ⏳ **report-only** restructure (runs concurrently with P2;
    report files only). Consolidates to **ONE supplement** (`supplementary.qmd`): folds the report's
    back-matter section in, and moves "How the model works" (anatomy) and the "Measured proteome"
    block there. Cross-document figure links are turned into **prose pointers** (SDs restated inline).
    **Removes** the moribund views: Spearman rank-consistency (Fig 9 / Table 5, + the S1 full table),
    calibrated-uncertainty (Fig 10 / Table 6), and the sector-sweep rank figure (Fig 14), cleaning
    all dangling refs. **Run before P3.**
P2c. `P2c_ceiling_diagnostic_prompt.md` — ⏳ small diagnostic: **name the ~1.55 structural growth
    ceiling**. The P2 re-run showed kcat unlocked (×2.56, plausible) but rmax capped at ~1.55 by a
    non-proteome constraint. Relaxes the proteome limits at the rich optimum and reads **shadow
    prices / active bounds** to identify what binds (carbon/O₂ uptake, ATPM/NGAM, GAM, or a biomass
    precursor), classifies it **fixable vs structural**, confirms by relaxing the culprit, and writes
    a NOTE for P3. **Run before P3.**
P3. `P3_tuned_model_analyses_and_report_prompt.md` — ⏳ dissect the **tuned (calibrated) model** at
    the rich BHI operating point (P2 posterior medians): rebuild the **elasticity** (unified set:
    kcat_scale in, budget_scale out, kappa_scale added; **posterior-propagated** uncertainty bands,
    no Spearman), the **decomposition** (envelope vs magnitude[allocation+catalysis]; median+IQR
    reachable-curve figure), and **identifiability**. Adds a **calibration section** to the report
    (prior-vs-posterior, demanded corrections, pool-vs-kcat finding), rebuilds the analysis sections
    on the tuned model, writes the interpretation (folding in the P2c named ceiling), and keeps P2b's
    single-supplement structure. **Run after P2 + P2b + P2c.**

**Report — model description (runs alongside the Bayesian; report-only)**
- `report_model_equations_comprehensive_rewrite_prompt.md` — ⏳ comprehensive, publication-grade
  rewrite of the **"Model equations, parameters and methods"** section in the style of **Li et al.
  2021** (Nat Commun 12:190): every governing equation written out in full with all terms defined
  (ec constraints, MMRT kcat(T), two-state native fraction, sector partition + biosynthesis/metabolic
  caps, the growth-law coupling, NGAM(T), descriptors, correction knobs), the model **hierarchy**,
  the **design choices**, a clear **"departures from Li et al."** subsection (single-pool→sector
  partition, growth law, grounded params, likelihood vs ABC, E. coli), and the updated **in-silico
  experimental design**. Report-only; safe to run concurrently with P2. ✅ run.
- `report_four_layers_update_prompt.md` — ⏳ follow-up (the equations rewrite already ran): updates
  the stale conceptual **"four layers"** overview to match — Layer 3 gains the **growth-rate-dependent
  coupled growth law** + restored **NGAM(T)** (distinct from the static medium-matched allocation),
  Layer 4 reframed from "nothing is fit" to **emergent-prior + Bayesian-inverse**, and stale
  glucose-minimal numbers (r_max 0.55) replaced with the **rich (BHI) reference**. Report-only.

_Deferred — respiration & CUE:_ once the growth model is trusted, add respiration
(CO₂/O₂ flux) and CUE TPC outputs and validate the Smith benchmarks (E_μ > E_R; CUE rises/
flat with warming; max CUE ≈ 0.4–0.6). Full library-scale validation against the 29 Smith
2021 environmental strains needs draft GEMs for those strains.

## The through-line

The scientific thread these encode: build an enzyme- & temperature-constrained GEM of
E. coli (eciML1515), decompose its TPC into a genome-set **envelope** and an
allocation-set **magnitude** (H1.1–H1.3), then progressively replace assumed pieces
with data — DLTKcat kcat(T), the MRes/Li Tm-based unfolding envelope (#13), and
measured temperature-dependent proteome allocation (#14). The recurring caveat
(captured in #12) is that the clean envelope/allocation split is partly a consequence
of model structure (peak-normalisation + temperature-independent allocation) and is
only resolvable with proteome/flux data — which #13 and #14 begin to supply.

## The decisive test (done by #15)

The recurring open question is whether the clean allocation-vs-envelope "division of
labour" is real or an artefact of two assumptions (peak-normalisation +
temperature-independent allocation). #13 removes the first (Tm-based unfolding) and #14
removes the second (measured temperature-dependent allocation). #15 assembles both into
one complete model, validates it against the empirical TPC, and re-runs the
decomposition on it — so whatever division of labour survives on the canonical,
data-grounded model is the real, defensible result (expected: more coupled, with `rmax`
gaining envelope share and larger interaction terms than the assumed-model 100/0 split).
