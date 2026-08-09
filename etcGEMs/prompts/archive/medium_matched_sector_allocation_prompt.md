# Claude Code prompt — medium-matched proteome-sector allocation (make the ribosome cap medium-dependent) + re-validate (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Run AFTER #17 has FINISHED and
committed (do not run concurrently). Growth only. Data-grounded, not fitted.

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read: src/etcgem/
{sectors.py,providers.py (set_medium),config.py,tpc.py}, the per-curve validation code and
strains/eciML1515/outputs/percurve_validation/, and reports/etcgem/report.qmd.

CONTEXT / FINDING TO FIX (from #17): rmax is identical on LB vs glucose-minimal (0.34/h)
because the biosynthesis/translation-sector cap binds and it is FIXED — the model cannot
grow the ribosome sector for rich media, so it misses the ~2-3x LB>minimal growth boost
(the growth-law effect, phi_ribosome ~ mu). The DeyuWang proteomics we already use contains
per-MEDIUM columns (LB16/25/30/37/43, Glucose25/30/37, Glycerol25/30/37), so the sector
fractions can be made MEDIUM-matched from measured data — no fitting.

PART A - compute MEDIUM-matched sector fractions from the proteomics
- From strains/eciML1515/proteomics/tem_proteomic.csv, compute mass-weighted sector
  fractions (metabolic/biosynthesis/chaperone/other via COG, as already done) SEPARATELY
  for each medium series: LB (16-43 C), Glucose (25-37 C), Glycerol (25-37 C). Output
  sector_fractions_by_medium.csv (rows: medium,temp_C; cols: f_metab,f_bio,f_chaperone,
  f_other). Report the biosynthesis fraction f_bio for LB vs Glucose at ~37 C (expect LB
  higher — faster growth, more ribosomes).

PART B - make the sector allocation medium-dependent
- Extend the sector layer so the sector fractions (and hence the translation cap and the
  pool budget = P_total*f_metab*sigma) are taken from the medium-matched curve: when the
  medium is LB use the LB fractions; when glucose_minimal use the Glucose fractions
  (Glycerol available too). Interpolate in temperature within each medium's measured range
  and clamp outside it (Glucose covers only 25-37 C). Wire this to the medium selected by
  set_medium so a run under LB uses LB f_bio and a run under minimal uses Glucose f_bio.
  Keep backward-compatible: if no per-medium data, fall back to the single measured curve.
- This should make the ribosome/translation cap higher under LB than minimal, so rmax under
  LB > rmax under minimal, straight from measured allocation (not fitted).

PART C - re-run the per-curve ABSOLUTE validation with medium-matched allocation
- Re-run the per-curve validation (#17): each curve predicted under ITS medium AND its
  medium-matched sector fractions. Recompute the PRIMARY absolute R^2/RMSE (raw 1/h) and the
  secondary shape R^2 and emergent Ea, split by medium. Report: does LB rmax now exceed
  minimal rmax and move toward the observed LB range (~1-2/h)? By how much does the absolute
  gap close? Update percurve_table.csv / percurve_summary.json and the plots (raw absolute,
  minimal vs LB highlighted, as in #17).

PART D - report the medium-matched change
- Update the model + validation sections to describe the medium-matched sector allocation:
  the ribosome/translation cap and pool budget now use MEDIUM-matched measured fractions
  (LB vs Glucose), the growth-law mechanism it captures, the LB-vs-minimal rmax difference,
  the closed/residual absolute gap, and the medium-dependent emergent Ea. State clearly it
  is data-grounded, not fitted.

PART E - make the report's MODEL / METHODS / VALIDATION description comprehensive
(a general rigor upgrade for the whole report; do it thoroughly, with equations)
1. MODEL EQUATIONS - add a proper model/methods section giving the governing equations as
   Quarto LaTeX (numbered, cross-referenced), covering at least:
   - MMRT kcat(T):  ln k(T) = ln(kB T/h) + (-dH0 - dCp(T-T0))/(R T) + (dS0 + dCp ln(T/T0))/R;
   - two-state unfolding native fraction  f_N(T) = 1/(1+exp(-dGu(T)/(R T)))  with dGu(T)
     (Tm, dHu, dCp,u);
   - effective per-flux enzyme cost  c_i(T) = base_cost_i / (rel_kcat_i(T) * f_N_i(T));
   - enzyme pool / sector constraints:  sum_i c_i(T)|v_i| <= f_metab * P_total ;
     translation cap  translation_coeff * v_biomass <= f_bio * P_total ; maintenance NGAM(T);
     pool budget  P_metab = P_total * f_metab * sigma ; and the MEDIUM-matched, temperature-
     dependent sector fractions f_sector(medium, T);
   - TPC descriptor definitions (Topt, rmax, CTmin/CTmax at the crit fraction, thermal
     breadth, Boltzmann-Arrhenius Ea from the rising limb);
   - the allocation-vs-envelope variance decomposition (two-group functional ANOVA + Shapley,
     phi_A = S_A + S_AE/2), per-enzyme thermal control coefficients, and the identifiability
     proxy.
2. DATA / PARAMETERISATION - add a provenance TABLE: each parameter (symbol), its SOURCE, its
   COVERAGE, and whether grounded-from-independent-data or a stated literature value (never
   fitted to growth): kcat_ref & MW (GECKO ecModel of iML1515 / BRENDA); Topt (Li-Engqvist
   2019); Tm (Leuenberger 2017 / Meltome); dCp (DLTKcat + MMRT-literature prior); sector
   fractions f_metab/f_bio/f_chaperone (DeyuWang proteomics, per medium & temperature);
   P_total (~0.5 g/gDW, literature); sigma (~0.4-0.5, Davidi 2016 / Heckmann 2020); medium
   (Smith 2019 metadata; LB from MRes media.csv / Machado 2018). Give the numeric values/
   ranges used and their references.
3. IN-SILICO EXPERIMENTS - describe each analysis in enough detail to reproduce: the global
   sensitivity sweep (inputs, ranges, LHS sample size), the crossed decomposition design
   (n_allocation x n_envelope, Shapley split), the two-stage per-enzyme control +
   identifiability, the calibrated-uncertainty sampling, and the per-curve validation design.
   Map each to its objective (O1-O6).
4. VALIDATION - describe the validation data and procedure: the Smith 2019 empirical E. coli
   TPCs (per-curve source studies + media), the RAW ABSOLUTE (1/h) per-curve comparison, the
   independent Ea benchmark (~0.85 eV), and the proteome (predicted-vs-measured) validation;
   state the metrics (absolute R^2/RMSE primary, shape secondary) and how undefined-medium /
   non-1/h curves are handled.
- Ensure all equations render and all crossrefs/citations resolve; add any missing refs.
  Re-run assemble.py; quarto render; confirm the PDF builds.

VERIFY (report all)
1. f_bio LB vs Glucose at 37 C (LB should be higher); the resulting rmax under LB vs minimal.
2. Per-curve absolute fit before (fixed-f_bio, #17 baseline) vs after (medium-matched):
   median abs_R2 split by medium; how much the LB magnitude gap closed. Emergent Ea per medium.
3. Report now contains: the model equations (rendered), the parameter-provenance table, the
   detailed in-silico-experiment descriptions, and the detailed validation section. PDF builds
   with no unresolved crossrefs/citations.

CONSTRAINTS
- Data-grounded, not fitted: fractions come from the measured per-medium proteome, never
  from the growth curve. Magnitude stays emergent (P_total*f_metab*sigma), Ea emergent.
- Medium is availability; do not pin uptakes. Reuse sectors.py/set_medium/compute_tpc.
- The MODEL EQUATIONS / provenance / experiment / validation write-up must match the ACTUAL
  code and parameters (read the modules; do not invent equations or values).
- Autonomous; commit in parts: "compute medium-matched sector fractions",
  "medium-dependent sector allocation (ribosome cap)", "re-validate per-curve",
  "report: model equations, parameter provenance, detailed experiments + validation".
```
