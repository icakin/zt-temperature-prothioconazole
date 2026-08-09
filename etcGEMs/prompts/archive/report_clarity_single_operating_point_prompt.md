# Claude Code prompt — report cleanup: one operating point, Tm in the figures, remove the incremental tone, and a whole-report clarity rewrite (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Run AFTER #20 (which already landed:
dTm knob, Tm in elasticity + decomposition, recast decomposition, IQR-band TPC figure). This
prompt is report-quality work — fixing inconsistencies #20 left and making the whole document
readable by someone who has NOT followed the development. It re-runs one analysis and
regenerates figures, but does NOT change model math.

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

CONTEXT (observed problems in the current report/outputs):
- The proteome-sector-trade-off subsection still cites the OLD sectors run's optimum
  (f_metab ~ 0.285), while the rest of the report is centred on the strain nominal
  (f_metab ~ 0.483, glucose-minimal). The whole report must use ONE operating point.
- dTm is present in elasticity_table.csv but does NOT appear in the FIGURES of the
  "What generates TPC variation (global sensitivity)" section (heatmap / tornado).
- The report still carries an incremental, self-referential tone that narrates corrections of
  our own earlier analyses (e.g. "An earlier version of this section did diverge from the
  elasticity; that was an artefact of the analysis, corrected below."). We want findings
  presented as they stand, not as a diff against previous runs.
- A collaborator found the prose confusing: "I also found some of the prose and explanation
  pretty confusing at times ... might be useful for how to frame reports in future." The report
  should assume less prior knowledge and be written out completely and self-containedly.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
reports/etcgem/report.qmd (WHOLE file), reports/etcgem/assemble.py, src/etcgem/plotting.py,
src/etcgem/{sectors.py,sensitivity.py,cli.py}, the sector-trade-off experiment/config and its
outputs under strains/eciML1515/outputs/ (the sectors run), and the elasticity outputs under
strains/eciML1515/outputs/elasticity_elasticity/. Do NOT change core model math.

PART A - ONE operating point across the whole report (strain nominal)
- The strain nominal (glucose-minimal) is f_metab ~ 0.483 (f_bio ~ 0.191). Make the ENTIRE
  report use this single operating point.
- The proteome-sector-trade-off analysis/subsection currently reports an optimum from an OLD
  sectors run (f_metab ~ 0.285). Re-run (or re-reference) that analysis so it is centred on /
  consistent with the strain nominal, and regenerate its figure(s)/table(s)/numbers. If the
  sector-trade-off legitimately reports a model-optimal f_metab that differs from the measured
  nominal, state BOTH explicitly and label which is which (measured operating point vs
  model-optimal) rather than silently mixing them.
- AUDIT the whole document: every figure, table and inline number must refer to the same
  operating point and the same current model. Find and fix any other stale outputs/numbers from
  superseded runs. Report what was stale and what it now is.

PART B - Tm must appear in the "What generates TPC variation" FIGURES
- dTm is in the elasticity table but not in the section's figures. Regenerate the elasticity
  HEATMAP (dTm as a row) and the TORNADO(s) (dTm as a bar) and make sure THESE regenerated
  figures are the ones embedded in the report (via assemble.py). Verify by opening/greping the
  rendered assets that dTm/Tm is actually present in the figures, not just the table.
- FIX the recast-decomposition figure (the "variance SHARE vs MAGNITUDE" two-panel figure —
  Figure 10 in the current report): panel (a), the "3-group variance share (range-fair)" bar
  chart, has its LEGEND OVERLAPPING the bars (the rmax bar fills to 1.0 and the legend sits on
  top of it). Move the legend OUT of the plotting area so it never overlaps the bars — e.g. a
  single shared legend placed above or below both panels, or anchored outside the axes
  (bbox_to_anchor). Since panel (b) already carries the same allocation / kinetic / Tm-stability
  legend, prefer one shared legend for the figure rather than one per panel. Regenerate the
  figure and confirm no overlap in the rendered PDF.

PART C - remove the incremental / self-referential tone (whole report)
- Delete all language that frames results as corrections of earlier versions of our own
  analysis, or that narrates the development history. Concrete instances to remove (search for
  these and similar): "An earlier version of this section did diverge ... artefact ... corrected
  below"; "was largely a structural artifact"; "overturns the assumed picture"; "a far cry from
  the ~100% allocation share ... older peak-normalised / scalar-pool convention"; "prior to the
  final emergent-magnitude refinement"; "the older ... convention"; "corrected below"; and any
  "old/earlier/previous run/version" comparisons. Present the model and every result AS IT
  STANDS NOW, as one finished piece of work.
- KEEP genuine, current scientific caveats and limitations: what is measured vs borrowed, the
  Topt/Ea coupling, the fact that which axis "dominates" is measure-dependent, and the rmax /
  high-temperature gaps. These are present-tense findings, not incremental narrative — state
  them plainly without referencing past versions of the analysis.

PART D - whole-report CLARITY rewrite (assume less knowledge; write it out completely)
- Rewrite the report so it is self-contained and understandable by a reader who has NOT
  followed the development. Assume little prior knowledge.
- Define every term and acronym IN FULL on first use, and explain what it is (not just expand
  the acronym): TPC (thermal performance curve), etc-GEM, FBA / flux-balance analysis,
  enzyme-constrained model, kcat and kcat(T), MMRT (macromolecular rate theory), two-state
  unfolding / native fraction, Tm (melting temperature), proteome sector, allocation vs
  envelope, elasticity, Shapley / variance decomposition, interaction term, in-vivo vs in-vitro
  kcat, sigma (saturation), NGAM (maintenance), CUE, and each descriptor (Topt, rmax,
  CTmin/CTmax, Ea, niche width).
- Lead each major section with a short plain-language summary (what this section does and why)
  BEFORE any equations. Then give the detail. Write complete, self-contained explanations
  rather than terse allusions. Use prose and paragraphs for explanation (not dense bullet
  dumps). Keep all equations, the parameter-provenance table and the rigour; make the
  surrounding prose readable and explanatory. Keep a plain, direct tone; do not over-format.
- Ensure consistent terminology and notation throughout, and that the document reads as one
  coherent, current, self-contained report.

PART E - render + check
- Re-run assemble.py; `quarto render`; confirm the PDF (and docx if configured) builds with no
  unresolved crossrefs or citations. Skim the rendered PDF to confirm: Tm appears in the
  sensitivity figures; the operating point is consistent (no f_metab ~ 0.285 anywhere unless
  explicitly labelled as model-optimal vs measured); no incremental-comparison language
  remains; terms are defined; each section opens with a plain-language summary.

VERIFY (report all)
1. Single operating point: confirm the sector-trade-off subsection and all report numbers use
   the strain nominal (f_metab ~ 0.483); list anything that was stale (e.g. 0.285) and is now
   fixed or explicitly relabelled.
2. Tm visible in the sensitivity heatmap AND tornado (not just the table); and the recast-
   decomposition figure (Fig 10) panel (a) legend no longer overlaps the bars.
3. Incremental/self-referential language removed (grep the qmd for "earlier", "artifact/
   artefact", "overturns", "older", "previous", "corrected below", "convention" and confirm
   none remain in a self-referential sense).
4. Terms/acronyms defined on first use; each major section opens with a plain-language summary.
5. PDF builds cleanly.

CONSTRAINTS
- Report-quality + figure/one-analysis regeneration only. Do NOT change core MMRT/unfolding/
  enzyme-cost/decomposition math.
- All numbers in the prose must match the regenerated outputs (read them; do not carry over
  stale figures).
- Do not remove genuine scientific caveats/limitations; only remove development-history and
  self-referential-correction framing.
- Autonomous; commit in parts: "report: single operating point (sector trade-off re-centred)",
  "report: Tm in sensitivity figures", "report: remove incremental/self-referential tone",
  "report: whole-report clarity rewrite (define terms, plain-language section leads)".
```
