# Revision edit package — figures 1–5 audit + inferential upgrades

**For: Gabriel Yvon-Durocher sign-off. Prepared from the full figure audit
(all five main figures data-verified; two text–data mismatches found and fixed)
plus the inferential-upgrade round. Main figures are FINAL as committed in
`tables/revision/fig1G|fig2|fig3|fig4|fig5/`; everything below is captions and
manuscript text. Scripts 75–82 reproduce every figure and statistic.**

---

## 1. Figure swaps (files ready)

| Manuscript figure | Replace with | Key changes |
|---|---|---|
| Figure 1 | `fig1G/Figure1_withG` | + panel G (activation energies) and ΔE_a strip |
| Figure 2 | `fig2/Figure2_redesign` | 5 panels; old C/D → one %-effect panel; heatmap full-width |
| Figure 3 | `fig3/Figure3_redesign` | fixed gene labels (CYP51 "ns" everywhere), DEG totals, grouped heatmap, panel D points on log2 |
| Figure 4 | `fig4/Figure4_redesign` | same conventions; labels suppressed in 15→21 facet; shared ±1.3 heatmap scale |
| Figure 5 | `fig5/Figure5_redesign` | 3 panels; B faceted; AOX–CUE → supplementary |
| NEW supp | `fig2/FigureS_resp_dose_multiples` | respiration fits diagnostic |
| NEW supp | `fig5/FigureS_AOX_CUE` (+ `fig5_loo_sensitivity.csv`) | demoted AOX–CUE map + LOO |
| NEW supp | `inference/FigureS_interaction_gsea` | pathway-level temperature × drug interaction |

## 2. Final captions
(Full text agreed in review; supplied per figure in the working thread. Headlines:)
- **Fig 1G:** E_a from the Sharpe–Schoolfield fits; respiration higher at every
  dose; across-dose contrast ΔĒ_a = 0.34 eV [0.21–0.48] vs |E_CUE| = 0.37 [0.26–0.49]
  — cross-model concordance, not an independent test; per-dose Pr in caption;
  "no systematic dose dependence was resolved" (backed by `inference/Ea_dose_contrasts.csv`).
- **Fig 2C:** % respiration change per mg/L; largest at 24 °C (+27.9 [15.4–42.4]);
  in 86% of draws the largest of seven; 95% CrI excluded zero at 15/21/24 only.
- **Fig 3:** fixed annotation set (CYP51/ERG11, AOX, CDR1, RPL7) labelled in all
  facets irrespective of significance; CYP51 ns at all temperatures (padj 0.83/
  0.056/0.94); heatmap = mean member-gene log2FC (not NES); panel D on log2(x+1),
  CIs across biological samples.
- **Fig 4:** title → "…limited toward the thermal optimum and extensive above it"
  ("regimes", not "phases"); fixed-set points unlabelled in 15→21 facet (caption
  explains); heatmap shares the ±1.3 scale with Fig 3; rows ordered by 15→27 effect.
- **Fig 5:** pooled = arithmetic mean of separately estimated NES; "most
  SIGNIFICANTLY ENRICHED pathways" in the lower-left (68% of 28; 32% of all 105);
  per-temperature sensitivity r = 0.57/0.68; B slopes = "gene-level alignment
  coefficients" on the 8,137 genes tested in all contrasts; C descriptive r,
  n = 6 condition means.

## 3. Results insertions

**(a) Formal interaction (new core inference; place with the Fig 3/5 text):**
> A factorial model with a temperature × prothioconazole interaction confirmed
> that the transcriptional drug response is temperature-dependent (likelihood-
> ratio test: 1,076 genes; 27 °C × drug contrast: 896 genes). At the pathway
> level, the interaction was dominated by attenuation of the drug response at
> 27 °C: repression of the gene-expression machinery weakened (ribosome
> interaction NES +3.44, FDR < 10⁻⁴; spliceosome +2.72; RNA polymerase +2.38),
> and induction of xenobiotic/P450 detoxification also weakened (NES −2.3),
> whereas ABC-efflux induction showed no significant interaction — the drug
> response that persists at supra-optimal temperature (Fig. S_interaction).
> At the sample level, the ribosomal-module response to the drug was 0.84 ±
> 0.33 log₂ units weaker at 27 °C than at 15 °C (temperature × drug interaction
> F(2,39) = 14.2, P = 2.4 × 10⁻⁵).

**(b) AOX guard sentence (with Fig 3D):**
> The apparent attenuation of AOX induction at 27 °C was not individually
> resolved (temperature × drug interaction P = 0.052); the contraction of the
> drug response is established at the module and transcriptome level.

**(c) Ribosome/translation observation (marker-plus-hypothesis framing):**
> Drug-induced repression of the ribosomal regulon was strongly temperature-
> dependent (sample-level interaction above), such that at 27 °C — where heat
> alone had already down-shifted the translational machinery (Fig. 4) — the
> drug removed substantially less residual ribosomal expression. Across the six
> conditions, ribosomal-module expression covaried with growth rate (r = 0.93,
> descriptive). These observations are consistent with a translation-capacity
> hypothesis for the loss of potency at supra-optimal temperature, in which the
> drug's growth effect is limited by how much translational capacity remains to
> be repressed; testing this causally requires perturbation experiments.

**(d) Target-axis sign audit (narrowed claim):**
> The observed loss of potency with warming was not explained by the measured
> CYP51-axis quantities: CYP51 and CPR transcript abundance, predicted CYP51
> turnover, and efflux-transporter expression all declined or remained flat with
> warming, individually and in combination predicting stable or increased —
> rather than decreased — potency (Supplementary sign audit). Within the
> enzyme-constrained metabolic model, reducing CYP51 capacity produced a
> temperature-invariant, threshold-like growth response, further indicating that
> the temperature dependence of potency arises above or downstream of
> single-target capacity.

**(e) State-convergence (honest uncertainty; NOT in abstract):**
> Projecting the untreated 27 °C transcriptome onto the fixed drug-response
> direction indicated a displacement of the control state toward the drug-
> responsive state (cos = 0.45, ≈20% of the total heat shift; sample bootstrap
> 95% CI 0.27–0.56). This displacement was directionally consistent but did not
> exceed a correlation-preserving sample-label permutation null at α = 0.05
> (P = 0.096) and is reported descriptively.

**(f) GEM/AOX softening (replaces current claims; per agreed wording):**
> Across the six condition means, AOX expression covaried inversely with CUE
> (descriptive r = −0.98, n = 6; leave-one-out range −0.99 to −0.98, Fig. S…).
> This identifies AOX induction as a marker of the low-CUE state but does not
> establish that AOX activity caused the decline in CUE. Under the minimum-flux
> objective the condition-constrained model admitted solutions with increased
> AOX flux; this allocation was not retained under a minimum-protein objective,
> and AOX transcript abundance correlated strongly with the simple respiration-
> to-growth ratio (ρ = 0.89) as well as with the model-derived AOX share
> (ρ = 0.99) across the same six conditions. The model therefore demonstrates
> compatibility with increased AOX use but provides limited evidence beyond the
> imposed physiological ratios for a uniquely AOX-mediated mechanism.

**(g) CDR1 counterexample sentence (with Fig 4):**
> Notably, CDR1 responded in opposite directions to the two stressors — induced
> by prothioconazole (Fig. 3) but repressed by warming (Fig. 4) — indicating
> that the drug–heat convergence is a programme-level phenomenon rather than
> gene-by-gene identity.

## 4. Discussion adjustments
- AOX paragraph → "candidate contributor and marker of the stress state rather
  than an established causal mediator"; add that the drug–warming overlap was
  transcriptome-wide (AOX ≈ 0.1% of the convergence signal; needs one
  supplementary-methods paragraph for the projection analysis).
- "near-silence" (15→21 °C) → "limited transcriptional response (49 genes)".
- Replace any "AOX module" phrasing (it is a single gene): global search.
- AOX naming: single alternative-oxidase gene (Mycgr3G72918; annotated AOX2 by
  orthology, hereafter AOX) — one Methods parenthetical; "AOX" everywhere else.

## 5. Explicitly NOT added (agreed discipline)
Gene-level Wilcoxon tests (pseudoreplication); post-hoc power analyses;
mediation/SEM on n = 6; new coexpression modules; dose-anchoring K-vs-IC₅₀ as
mechanistic evidence (supplement footnote at most); shift-vs-scale result
(stays out per prereg); any smooth T-model fitted to 3 RNA temperatures;
26 °C ad-hoc explanations (influence diagnostic queued instead).

## 6. Remaining queue (supplement-grade, non-blocking)
Posterior-predictive-check package for the Bayesian models; RNA-seq batch/
influence audit; GEM flux-variability analysis (AOX identifiability); 26 °C
leverage check.
