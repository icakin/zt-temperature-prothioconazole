# Claude Code prompt — P3: dissect the TUNED (calibrated) model — sensitivity, decomposition, identifiability — and finish the report (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Run AFTER P2 v3 (Bayesian tuning on Van
Derlinden, outputs in strains/eciML1515/outputs/calibration_vanderlinden_v3/). The report already
has a single supplementary.qmd, the refreshed emergent validation, and the updated four-layer
overview committed — build on that; do NOT recreate a back-matter supplement. Growth only. Multi-
hour run (decomposition grid + elasticity + control + render). Use Gurobi (PART 0).

NOTE TO USER: launch in an auto-approving mode, from the venv that has Gurobi.

SCOPE (decided): the sensitivity, decomposition and identifiability analyses run on the TUNED
(calibrated) model ONLY, at the rich (BHI) operating point, with the SAME unified parameter set the
calibration used and uncertainty carried by the P2 v3 posterior. The emergent model appears only in
the validation section (a-priori check). Report arc: validation (emergent) -> calibration (tuned) ->
dissect the tuned model (sensitivity -> decomposition -> identifiability) -> interpret.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
strains/eciML1515/outputs/calibration_vanderlinden_v3/ (P2 v3 posterior: summary.json,
demanded_corrections.csv, corner, prior-vs-posterior TPC),
strains/eciML1515/outputs/{pool_reconciliation,ceiling_diagnostic}/NOTE.md (pool-reconciliation
story for the interpretation), src/etcgem/{enzyme_cost.py (Perturbation incl. kcat_scale, kappa_scale,
sigma),calibration.py (PARAM_NAMES, the v3 free set),sensitivity.py (elasticity),decomposition.py,
control.py,providers.py (set_medium BHI, reconciled single pool),sectors.py (growth law),cli.py,
plotting.py}, reports/etcgem/report.qmd and reports/etcgem/supplementary.qmd. READ P2 v3's ACTUAL
numbers and write the report from them — do not assume outcomes.

PART 0 - solver: use Gurobi (same result, faster/cleaner)
- Switch provider/worker models to Gurobi (`pm.ec.model.solver = "gurobi"`) with a GLPK fallback and
  a one-solve pre-flight; PRINT which solver is in use. Gurobi is installed + licensed in this venv.
  If it is NOT active, STOP with a clear message (do not grind the decomposition grid on GLPK) unless
  the user set ALLOW_GLPK.

PART A - define the TUNED model (rich BHI operating point + P2 v3 posterior)
- Build a helper that returns the tuned Perturbation from P2 v3's POSTERIOR MEDIANS across the FULL
  unified set (dTopt, topt_scale, dCp_scale, dTm, tm_scale, kcat_scale, kappa_scale, sigma, f_metab,
  f_maint, ngam_scale, ngam_steepness) and sets the RICH (BHI) operating point (set_medium(pm,"BHI"),
  LB rich sector allocation, growth law ON, reconciled single pool). This "tuned model" is the
  baseline for ALL P3 analyses. Load a sample of posterior DRAWS (e.g. 100-200) for uncertainty
  propagation. SANITY: confirm it reproduces the v3 headline (rmax ~2.14 [1.98, 2.29], Topt ~38.4,
  CTmax ~45.9, Ea ~0.92) before proceeding.

PART B - unified parameter set for the analyses (match the v3 calibration)
- Set the sensitivity (elasticity) and decomposition INPUT set to the SAME unified set as P2 v3:
  envelope/stability {dTopt, topt_scale, dCp_scale, dTm, tm_scale} + magnitude/catalysis
  {kcat_scale, kappa_scale, sigma} + allocation {f_metab, f_maint} + maintenance {ngam_scale,
  ngam_steepness}. (sigma_disc is a noise term, NOT a model lever — exclude it. budget_scale stays
  dropped, degenerate with kcat_scale/sigma.) All perturbations are equal standardised steps about
  the TUNED operating point.

PART B1 - OUTPUT TAGGING (do NOT overwrite the stale untagged dirs)
- Write the tuned analyses to explicitly-tagged dirs so they do not clobber the pre-reconciliation
  ones and the report references only the current results:
  * strains/eciML1515/outputs/{elasticity_tuned, decompose_tuned, control_tuned}/
- The existing untagged dirs (elasticity_elasticity, decompose_decomposition_*, control_control) are
  pre-reconciliation / pre-growth-law and STALE — leave in place but do NOT reuse or reference them.

PART C - tuned-model sensitivity (elasticity), uncertainty-aware (no Spearman)
- Run the equal-perturbation elasticity around the TUNED operating point (rich BHI) with the unified
  set. Report the normalised elasticity per descriptor (tornado + heatmap). PROPAGATE the posterior:
  recompute the elasticities over the posterior DRAWS to get an uncertainty band / IQR on each (this
  replaces the deleted calibrated-uncertainty view). Do NOT reintroduce the Spearman rank view.
  NOTE (changed by the reconciliation): on the single pool the magnitude levers (sigma, kcat_scale)
  now genuinely MOVE r_max at the rich peak — report their elasticity explicitly (this contrasts with
  the pre-reconciliation behaviour where the metabolic pool had slack; say so).

PART D - tuned-model decomposition
- Recast allocation-vs-envelope decomposition on the TUNED model at the rich operating point, same
  unified set. Grouping: envelope {dTopt, topt_scale, dCp_scale, dTm, tm_scale} vs magnitude
  {allocation f_metab,f_maint + catalysis kcat_scale,kappa_scale,sigma + maintenance ngam_scale,
  ngam_steepness}; also report the finer split (kinetic / stability / allocation / catalysis /
  maintenance). Keep the median+IQR reachable-TPC band figure (per lever) as the central figure, and
  report variance shares AND magnitude (IQR) side by side. The crossed grid is expensive: run it at
  the posterior-median tuned point for the headline; if affordable, propagate a few posterior draws
  for a share-uncertainty note, else state the point-estimate caveat and let the elasticity bands
  carry the uncertainty.

PART E - enzyme identifiability (tuned model)
- Recompute the per-enzyme thermal control coefficients + identifiability on the TUNED model at the
  rich operating point (control_tuned).

PART F - report: add the calibration result, rebuild the tuned analyses, write the HONEST interpretation
- Add a "Calibration: what the data demand" section (from calibration_vanderlinden_v3/): the prior-
  vs-posterior TPC, the demanded-corrections table, the corner/degeneracy plot, and the headline from
  the ACTUAL v3 posterior:
  * MAGNITUDE MOSTLY CLOSED: emergent r_max 1.04 -> tuned 2.14 [1.98, 2.29] vs observed 2.40 (~89% of
    the gap closed).
  * sigma = 0.867 [0.671, 0.986] railed toward its physical 1.0 ceiling (well above the 0.4-0.5
    literature range); kcat_scale a modest 1.25 [0.86, 2.64]; the two share the magnitude on an
    anti-correlated ridge.
  * SHAPE captured: CTmax 51.8 -> 45.9 (matches observed) via dTm -5.6 K [-7.2, -2.8]; Topt 38.4 vs
    40; Ea held near the ~0.85 bacterial benchmark (posterior 0.92) with sigma_disc ~0.20 absorbing
    the possibly-flattened digitized 0.64 rather than chasing it; sectors at measured (f_metab 0.281,
    f_maint 0.341); kappa/ngam/topt_scale/dCp/tm_scale near prior.
- FRAME THE RESIDUAL HONESTLY — and CORRECT THE RECORD (this OVERTURNS the earlier "structural ~1.55
  ceiling" conclusion):
  * The apparent ~1.55 ceiling was a PROTEOME DOUBLE-ACCOUNTING artefact — two redundant enzyme-mass
    pools (the sMOMENT sector pool + a leftover GECKO base `prot_pool_exchange` bound), so the
    magnitude levers saturated against a pool they could not touch. The ceiling diagnostic identified
    this and the pool reconciliation (single sector/growth-law budget) fixed it. Cite the
    pool_reconciliation / ceiling_diagnostic NOTEs. Do NOT repeat the old "model structurally under-
    represents the rich/minimal growth-law boost" claim.
  * The remaining 2.14 -> 2.40 residual is the sigma <= 1 PHYSICAL CEILING: sigma pins near its bound,
    so the rich peak demands near-MAXIMAL in-vivo enzyme saturation. Honest reading: E. coli operates
    at higher in-vivo capacity than the yeast-derived 0.4-0.5 average — the fast-grower point,
    consistent with "departures from Li" (Li's yeast sat comfortably at sigma=0.5 and was never
    stress-tested against a 3-5x faster organism; 2.40 h-1 ~ a 17 min doubling is a real, high-end
    E. coli rich-medium rate, and the observed value is a figure-digitized high-end point).
  * ADD THE SOFT CAVEAT: a parameter pinned to its bound signals the metabolic budget (P_total x
    f_metab) is tight at the peak — and the coupled growth law deliberately shifts proteome to
    ribosomes (shrinks f_metab) exactly at fast growth, so part of the residual is the growth-law
    allocation squeeze at the rich peak. Frame as a concrete, quantified modelling frontier, NOT an
    unexplained failure.
- Rebuild the sensitivity (elasticity, tuned, uncertainty-aware), decomposition (tuned), and
  identifiability (tuned) sections from the *_tuned dirs, replacing the old glucose-minimal-nominal
  versions. assemble.py must reference the *_tuned dirs, not the stale untagged ones. REMOVE the
  stale placeholders in report.qmd (the "rebuilt at the rich reference ... on the calibrated model in
  a later phase" line, and the "single reference growth curve (glucose-minimal)" framing) — those
  analyses are now DONE on the tuned rich model.
- DE-FOREGROUND RESULTS FROM THE SETUP SECTIONS. The "complete model" / four-layers overview AND the
  equations/methods section must describe the model STRUCTURALLY and CONCEPTUALLY only — NO
  quantitative results ahead of the sections that present them. Strip from those setup sections any:
  tuned/emergent r_max values (1.04, 2.14, 2.40), sigma/kcat/dTm posterior values, magnitude-closure
  or "gap" statements, validation R^2/CTmax/Ea/Topt numbers, and any "the model reaches / demands /
  under-predicts ..." claim. Replace inline numbers with forward cross-refs in prose ("quantified in
  the Calibration section", "see Validation"). Layer 3 may still SAY it adds the growth law + NGAM(T);
  Layer 4 may still SAY the emergent model is a prior refined by Bayesian inversion — WITHOUT the
  rich-reference numbers. The ONLY place results appear up front is the ABSTRACT. Sweep Background/
  objectives + "complete model" for the same issue.
- Section order (main body): Background/objectives -> complete model -> equations/methods ->
  validation (Van Derlinden, emergent) -> calibration (tuned) -> sensitivity (tuned) ->
  decomposition (tuned) -> control/identifiability (tuned) -> interpretation -> next steps. Single
  supplement unchanged (anatomy, measured proteome, ablations, surviving tables).
- Update the abstract, objectives and interpretation to the tuned-model narrative (shape grounded +
  magnitude mostly closed by calibration + residual = near-maximal sigma / growth-law budget
  squeeze). Remove any remaining "later phase" placeholders.
- Update assemble.py to collect the new tuned figures/tables + the calibration figures; re-run
  assemble.py; quarto render BOTH report and supplementary; confirm no unresolved crossrefs.

VERIFY (report all)
0. Solver = gurobi (or aborted, not silently GLPK).
1. The tuned model = rich BHI + reconciled pool + growth law ON + P2 v3 posterior medians (full
   unified set incl. sigma, tm_scale, ngam_*); it reproduces the v3 headline (rmax ~2.14, Topt ~38.4,
   CTmax ~45.9, Ea ~0.92). Analyses' input set = the unified set (sigma in; sigma_disc excluded;
   budget_scale out). Outputs in *_tuned dirs; stale untagged dirs not used.
2. Tuned elasticity: which levers drive which descriptors, with posterior uncertainty bands; sigma
   and kcat_scale shown to move r_max at the rich peak; no Spearman view anywhere.
3. Tuned decomposition: envelope-vs-magnitude shares + finer split + median+IQR reachable-curve
   figure.
4. Identifiability on the tuned model.
5. Calibration section added with the v3 numbers; residual framed as sigma-ceiling + growth-law
   squeeze (NOT structural), with the pool-reconciliation record correction; interpretation updated.
6. Setup sections DE-FOREGROUNDED: four-layers/complete-model overview + equations/methods carry NO
   downstream results — only structure + forward cross-refs (results up front only in the abstract);
   confirm by grepping those sections for stray numeric results.
7. Section order as above; one supplement; both PDFs build; no placeholders or dangling refs.

CONSTRAINTS
- All downstream analyses on the TUNED model at the rich (BHI) operating point; do NOT fall back to
  the glucose-minimal default. The emergent model appears only in validation.
- Unified parameter set (sigma in, budget_scale out, sigma_disc excluded from levers). Uncertainty
  via the P2 v3 posterior (propagated through the elasticity), NOT a re-introduced calibrated-
  uncertainty ensemble.
- Keep the single-supplement structure; do NOT re-create a back-matter supplement or the removed
  Spearman/calibrated/sector-rank views.
- Write the report from P2 v3's ACTUAL numbers; keep emergent vs tuned distinct; CORRECT the old
  "structural ceiling" narrative per PART F.
- Autonomous; commit in parts: "tuned-model helper (v3 posterior medians + rich BHI + reconciled pool)",
  "analyses: unified param set (sigma in, budget_scale out) + *_tuned output tagging",
  "tuned-model elasticity (posterior-propagated) + decomposition + identifiability",
  "report: add calibration section (v3) + rebuild tuned analyses + de-foreground setup + honest sigma-ceiling interpretation".
```
