# Claude Code prompt — rebuild all report analyses on the two TRUSTED, strain-matched curves (Noll NCM3722 minimal + Erdos MG1655 LB rich); retire the old compilation; NO Bayesian tuning (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Growth only. This REPLACES the validation
dataset and re-runs the pipeline + report against two trusted curves. It does NOT touch the
Bayesian calibration (that is redone separately, later).

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

WHY: the previously-used Smith-derived compilation (strains/eciML1515/thermal/ExpGrowth.csv,
ecoli_tpc_curves.csv, 26 curves) is no longer trusted — one curve (Bennett-Lenski/Cooper) was
found to be mis-scaled by ~8x, and for most curves the exact strain, medium and data robustness
are uncertain. We now have two high-quality, MODEL-MATCHED curves and will validate against
those instead:
- MINIMAL: Noll (Katipoglu-Yazan et al. 2023), E. coli K-12 NCM3722, defined glucose minimal
  (modified MOPS + 0.5 g/L glucose + 6 trace amino acids), 27-45 C, mean +/- SD per T, n wells.
  File: strains/eciML1515/thermal/sources/noll2023_ncm3722/noll2023_ncm3722_tpc.csv
- RICH: Erdos et al. 2026 (unpublished), E. coli K-12 MG1655 wt, LB, 0 kanamycin, ~16-45 C.
  File: strains/eciML1515/thermal/sources/erdos2026_unpubl/erdos2026_mg1655_wt_nokan_tpc.csv
  (figure-digitized by eye, +/- ~0.1 h-1; treat as lower-precision).
An additional MG1655 rich curve (Van Derlinden 2012, BHI, 7-46 C) is available under
sources/vanderlinden2012_intfoodmicro/ as an OPTIONAL independent rich cross-check.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first (take stock of the
CURRENT state before assuming anything): reports/etcgem/report.qmd (whole file) +
reports/etcgem/assemble.py, the per-curve validation code and its outputs, src/etcgem/
{providers.py (set_medium),config.py,tpc.py,sectors.py,cli.py}, strains/eciML1515/strain.yaml,
and the three trusted-curve CSVs under strains/eciML1515/thermal/sources/. Do NOT run or modify
the Bayesian calibration (src/etcgem/calibration.py, the `calibrate` command, calibration_phase1
outputs) in this prompt.

PART A - make the two trusted curves the canonical validation set; retire the old compilation
- Load the validation curves DIRECTLY from the source CSVs: Noll (glucose_minimal) and Erdos
  (LB rich). Carry their strain, medium_detail, and (for Noll) the per-temperature SD and n.
- RETIRE the Smith-derived 26-curve compilation from the MAIN analyses and report. Do not delete
  the files, but stop using ExpGrowth.csv / ecoli_tpc_curves.csv as the validation source.
  If retained at all, only as a clearly-labelled SUPPLEMENTARY "legacy compilation (low
  confidence)" note — not in the main validation.
- Medium handling (availability, via set_medium): Noll -> glucose_minimal; Erdos -> LB.
  For Noll, note the 6 trace amino acids; treat as glucose_minimal for the headline, and (if
  cheap) also report a variant with those 6 amino-acid uptakes opened as a sensitivity check.

PART B - re-run VALIDATION on raw ABSOLUTE rates: minimal (Noll) vs rich (Erdos)
- Predict the emergent model TPC under glucose_minimal and compare to Noll; predict under LB and
  compare to Erdos. Compare on RAW ABSOLUTE growth rate (1/h). Where Noll SDs exist, weight/report
  them (and use them as the natural error band).
- Report, per curve: absolute R^2 / RMSE, predicted vs observed Topt, rmax, CTmax, and the
  emergent vs observed rising-limb Ea (compare to the ~0.85 eV bacterial-growth benchmark).
- HEADLINE contrast: the minimal-vs-rich magnitude. Does the medium-matched sector allocation
  reproduce the observed rich>minimal peak (Erdos ~2.35 vs Noll ~1.0 h-1, ~2.3x)? Report the
  model's LB/minimal rmax ratio against the observed ratio, and read any residual as a finding
  (kcat scale / sigma / medium composition), not a tuning failure (no tuning here).
- Figures: model-vs-data overlay for each curve on absolute rates (Noll with SD error bars);
  a minimal-vs-rich panel. Replace the old 26-curve small-multiples/So distribution figures.

PART C - re-run the model-intrinsic analyses for consistency (regenerate all figures/tables)
- Re-run the reference TPC / anatomy, the equal-perturbation elasticity, the recast decomposition,
  the per-enzyme control/identifiability, and the ablations, at the glucose-minimal reference
  operating point. These are model-intrinsic (they do NOT depend on the empirical curves), so the
  numbers should be essentially unchanged — re-run only to regenerate outputs consistently and to
  confirm nothing shifted. Note this explicitly in the summary.
- OPTIONAL (only if cheap): also compute the reference TPC (and, if quick, the elasticity) at the
  LB operating point, so the anatomy/validation can show both media. Keep scope bounded.

PART D - update the WHOLE report
- Rewrite the validation section around the two trusted curves: lead with "validated against two
  high-quality, strain-matched curves - a defined glucose-minimal curve (K-12 NCM3722, Noll 2023,
  with per-temperature SD) and a rich LB curve (MG1655 wt, Erdos 2026)", present the absolute-rate
  fits, the descriptors, and the minimal-vs-rich magnitude contrast as the headline validation.
- Add a short, factual DATA-PROVENANCE paragraph explaining why the validation set changed: the
  older compilation had a documented ~8x scaling error in one curve and uncertain strains/media/
  robustness, so we moved to strain-matched primary sources. State the digitization caveats
  honestly (Noll = transcribed from the paper's Table 1 with SD; Erdos = figure-digitized,
  unpublished, lower precision; Van Derlinden available as an optional independent rich cross-check).
- Update the reference-TPC/anatomy comparison to use Noll (minimal). Keep the model-intrinsic
  sensitivity/decomposition/control/ablation sections (unchanged results). Update the abstract,
  results text, and next-steps to reflect the new validation basis. Keep the writing clear,
  self-contained, and free of incremental "this overturns our earlier..." tone (consistent with the
  clarity pass).
- Re-run assemble.py; quarto render; confirm the PDF builds with no unresolved crossrefs.

VERIFY (report all)
1. Validation now uses ONLY Noll (minimal) + Erdos (rich); the old 26-curve compilation is retired
   from the main report.
2. Absolute-rate fits (R^2/RMSE), predicted vs observed Topt/rmax/CTmax/Ea for each curve; the
   model's minimal-vs-rich rmax ratio vs the observed ~2.3x.
3. Model-intrinsic analyses re-ran and are essentially unchanged (state any differences).
4. Report re-renders; validation + provenance sections read clearly; digitization caveats stated.

CONSTRAINTS
- NO Bayesian tuning in this prompt (do not run/modify calibration).
- Nothing is fit to the growth curves (emergent model); medium is availability, not pinned uptake.
- Absolute rates are primary; peak-normalised shape only as a secondary diagnostic if useful.
- Be honest about the Erdos digitization uncertainty; prefer Noll (with SD) as the quantitative
  anchor and Erdos as the rich-medium comparator.
- Match all numbers in the prose to the regenerated outputs.
- Autonomous; commit in parts: "retire legacy compilation; adopt Noll + Erdos as validation set",
  "validation: absolute-rate minimal (Noll) vs rich (Erdos) + medium contrast",
  "re-run model-intrinsic analyses (anatomy/elasticity/decomposition/control/ablation)",
  "report: rebuild validation + data-provenance around trusted strain-matched curves".
```
