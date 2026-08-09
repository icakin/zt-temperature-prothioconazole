# Claude Code prompt — consolidate the COMPLETE etcGEM model, validate it, re-run analyses on it, and restructure the report (autonomous, comprehensive)

Run from the project root (`.../MICROADAPT/etcGEMs`). This is a large, multi-part
change. Run AFTER the unfolding (#13) and proteome (#14) prompts have committed.

NOTE TO USER: launch Claude Code in an auto-approving mode so it runs unattended
(this does model config changes, several sweeps, and a full report rewrite). Everything
is committed in git and revertible.

---

```
Work AUTONOMOUSLY end to end; commit in logical parts; print a summary at the end.
GOAL: stop presenting the report as an incremental build-up from an assumed nominal
model. Instead: (1) assemble ONE complete, data-grounded model; (2) validate it
against the empirical E. coli TPC and the measured proteome; (3) run the sensitivity /
decomposition / control analyses ON THAT complete model; (4) restructure the report so
the complete validated model is the headline and the incremental "what each layer adds"
becomes a model-construction / supplementary section. Report results honestly even if
the clean allocation-vs-envelope division of labour weakens.

TAKE STOCK FIRST (do not assume). Read: strains/eciML1515/strain.yaml (thermal_model,
proteome_sectors, allocation_from_data), src/etcgem/{config.py,providers.py,sectors.py,
enzyme_cost.py,tpc.py,decomposition.py,control.py,sensitivity.py,plotting.py},
experiments/*.yaml, reports/etcgem/{report.qmd,assemble.py}, and the existing output
dirs under strains/eciML1515/outputs/. Note the current measured 30 C sector fractions
(f_metab~0.285, f_bio~0.341, f_chaperone~0.077, f_other~0.298) and the proteomics file
strains/eciML1515/proteomics/tem_proteomic.csv. The MRes empirical growth TPC is in the
cloned MResProject repo (data/ExpGrowth.csv) — use it for validation.

============================================================
PART A - assemble the COMPLETE model (one configuration)
============================================================
Make the complete model the strain's model (edit strains/eciML1515/strain.yaml):
- thermal_model: unfolding   (per-enzyme Topt/Tm grounded from the MRes tables /
  Leuenberger; the Tm-based falling limb).
- proteome_sectors: enabled: true, with:
  * nominal sector fractions SET FROM the measured 30 C proteome (not the hand-set
    0.5/0.15). Reconcile the data's 4 COG sectors (metabolic, biosynthesis/ribosome,
    chaperone, other) with the model's sectors; represent the chaperone/stress sector's
    temperature dependence explicitly if the machinery allows, else fold chaperone+other
    into maintenance — DOCUMENT the mapping.
  * allocation_from_data: <path to sector_fractions_vs_T.csv> so the sector fractions
    are the MEASURED temperature-dependent curve (interpolated per temperature),
    anchored at 30 C.
- Apply DLTKcat per-enzyme Topt/dCp where good fits exist (use_dltkcat_fits / the fits
  csv), falling back to the grounded/unfolding defaults.
- Build it (`etcgem build --strain eciML1515`) and confirm it grows and the sectors +
  unfolding are active in the resolved config. Report enzyme coverage for Topt/Tm
  (grounded vs default) and for the proteome sector mapping.
Commit: "assemble complete data-grounded etcGEM (unfolding + measured T-dependent sectors)".

============================================================
PART B - VALIDATE the complete model
============================================================
- Empirical TPC fit: compute the complete model's nominal growth TPC over the
  empirical temperature range and score it against the empirical E. coli growth data
  (MResProject data/ExpGrowth.csv, and/or the Smith 2019 curve used by the MRes).
  Report R^2 / RMSE and produce a figure `complete_vs_empirical_tpc.png` (model line +
  empirical points). State Topt, rmax, CT_max and how they compare to observed.
- Proteome validation: reuse / refresh the predicted-vs-measured per-enzyme allocation
  and the measured vs modelled sector fractions across 16-43 C (the proteome_* figures/
  tables already exist — regenerate against the complete model).
- Write a short validation summary (validation_summary.json).
Commit: "validate complete model vs empirical TPC and measured proteome".

============================================================
PART C - re-run the HEADLINE analyses on the complete model
============================================================
Run WITH figures (no --no-plots), into clearly named complete-model output dirs:
- Global sensitivity sweep (envelope + allocation axes appropriate to the complete
  model: envelope = per-enzyme Topt/Tm perturbations; allocation = the data-anchored
  sector fractions / total budget). 
- Allocation-vs-envelope decomposition (the Shapley split) — allocation axis perturbs
  the MEASURED, temperature-dependent sector allocation; envelope axis perturbs the
  grounded Topt/Tm. This is the decisive, defensible version of the H1.3 result.
- Per-enzyme control + proteome-wide identifiability.
- Calibrated (correlated / DLTKcat-posterior) uncertainty ensemble.
Report, for the decomposition, phi_A/phi_E and the interaction term per descriptor, and
COMPARE to the old assumed-model numbers (rmax was ~100% allocation, interactions tiny):
state plainly whether, on the complete model, rmax picks up envelope share and the
interactions grow (i.e. whether the clean division of labour was partly structural).
Commit: "re-run sensitivity / decomposition / control / calibrated on the complete model".

============================================================
PART D - ablations for the model-construction / supplementary section
============================================================
Produce a small set of ablation runs (nominal TPC + key descriptors + fit-to-empirical
R^2 for each) to show the incremental contribution of each layer, into clearly named
dirs (e.g. ablation_mmrt, ablation_sectors_off, ablation_sectors_Tindep, ablation_no_dltkcat):
- mmrt (peak-normalised, no unfolding);
- unfolding but scalar pool (sectors off);
- unfolding + temperature-INDEPENDENT sectors (allocation_from_data off);
- complete model.
This gives an honest "what each ingredient adds" comparison. SANITY: the mmrt+scalar
ablation should reproduce the pre-unfolding nominal numbers.
Commit: "ablation runs quantifying each model layer's contribution".

============================================================
PART E - RESTRUCTURE the report around the complete model
============================================================
Rewrite reports/etcgem/report.qmd to this structure (keep the writing quality; update
the abstract to match):
1. Background and objectives (largely unchanged).
2. The complete model — describe the final data-grounded etcGEM (unfolding/Tm envelope,
   DLTKcat kcat(T), measured temperature-dependent sector allocation), with provenance
   (iML1515/GECKO/cobrapy, Hobbs MMRT, Li2021/Pettersen2023, Qiu2024 DLTKcat,
   Basan2015/Scott2010, Mairet2021/Wang2026 for temperature-dependent allocation,
   Madkaikar2023 + Leuenberger2017 for grounded Topt/Tm).
3. Validation — fit to the empirical E. coli TPC (@fig complete_vs_empirical_tpc) and
   the proteome (predicted-vs-measured allocation, chaperone ramp). Lead with "it
   reproduces the data".
4. What generates TPC variation (global sensitivity) — on the complete model.
5. Allocation vs envelope decomposition — on the complete model; report the numbers and
   explicitly revisit the "is it forced or emergent?" question with the honest answer
   (likely more coupled than the assumed-model 100/0 split).
6. Per-enzyme thermal control and identifiability.
7. Interpretation and caveats (update: keep the peak-normalisation/growth-fit/
   temperature-dependent-allocation caveat, now largely ADDRESSED by the complete model;
   note remaining limits — Ea/breadth per-enzyme dCp, chaperone-sector fidelity,
   E. coli-specific, first-order identifiability).
8. Next steps.
9. Supplementary: model construction — the ablation series (Part D) showing each layer's
   incremental contribution and the fit improvement; the DLTKcat pipeline; parameter
   sources and coverage.
- Update reports/etcgem/assemble.py so all figure/table sources point at the new
  complete-model run dirs (and the ablation dirs for the supplementary).
- `python reports/etcgem/assemble.py` then `cd reports/etcgem && quarto render`; confirm
  _output/report.pdf builds with no unresolved crossrefs/citations. Add any missing refs.
Commit: "restructure report: complete validated model -> validation -> sensitivity ->
decomposition; incremental construction moved to supplementary".

VERIFY (report all)
1. Complete-model build: coverage of grounded Topt/Tm and proteome sector mapping;
   sectors + unfolding + allocation_from_data all active in the resolved config.
2. Validation: R^2/RMSE vs the empirical TPC; nominal Topt/rmax/CT_max; proteome
   correlations.
3. Decomposition on the complete model: phi_A/phi_E + interaction per descriptor, and
   the head-to-head vs the assumed-model numbers (did rmax's allocation share fall and
   interactions rise?).
4. Report re-renders to _output/report.pdf with the new structure.

CONSTRAINTS
- Keep numerical/scientific code correct; you may add config/plumbing but do not alter
  MMRT/unfolding/cost math or the decomposition ANOVA/Shapley logic.
- Everything data-grounded must be honestly reported, including where the complete model
  fits WORSE or where the division of labour weakens — that is the intended outcome.
- Preserve the existing incremental outputs where useful; do not delete committed
  results (the ablations regenerate what the supplementary needs).
- Autonomous; commit in the parts above; leave a clean working tree.
```
