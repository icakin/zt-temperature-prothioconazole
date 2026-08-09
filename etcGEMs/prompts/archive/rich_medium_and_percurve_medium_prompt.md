# Claude Code prompt — add a rich (LB) medium, validate per-curve with the correct medium, and document it (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Run AFTER #16 (emergent per-curve
validation) has committed. Growth only. The medium is an experimental INPUT (availability),
never fit to the growth curve.

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
src/etcgem/providers.py (set_medium — currently toggles carbon sources + O2 only),
the per-curve validation code and its outputs under
strains/eciML1515/outputs/percurve_validation/ (percurve_table.csv,
percurve_smallmultiples.png, percurve_R2_distribution.png, percurve_summary.json),
and reports/etcgem/report.qmd (validation section). The current state: of 26 E. coli
curves only 2 are defined (glucose_minimal); the other 24 are rich_broth and were run on
the DEFAULT minimal-glucose model (identical prediction for all, shape-only, abs_R2 blank).
Fix: give the model an LB/rich medium so the rich curves are predicted under the RIGHT
condition and their magnitude enters validation; distinguish minimal vs LB in the plots;
document the whole medium treatment.

PART A - add an LB / rich medium (availability, not pinned uptakes)
- Extend set_medium (or add set_medium(..., medium="glucose_minimal" | "LB" | list)) so it
  can open a RICH medium: the LB component uptakes (20 amino acids, nucleosides/bases, and
  common ions/vitamins) plus a carbon source and O2. Encode as AVAILABILITY: open the
  `EX_<met>_e_REV` uptake for each LB component present in the model, close other carbon
  sources' non-LB uptakes as appropriate, and let the enzyme-constrained model determine
  actual uptake (do NOT pin uptake rates).
- SOURCE the LB composition from the MResProject repo `data/media.csv` (it already has a
  GEM-compatible LB definition; reuse it) or the Machado et al. 2018 LB composition. Map
  the component metabolite ids to this model's BiGG exchanges (report which LB components
  were found vs missing). Document the LB definition and its source.
- Keep glucose_minimal as the other medium. Report the emergent rmax under LB vs minimal.

PART B - re-run per-curve validation on RAW ABSOLUTE rates (magnitude is emergent)
- Since magnitude is now emergent, ABSOLUTE fit is the PRIMARY validation, not peak-
  normalised shape. Compare the model's predicted specific growth rate mu(T) in 1/h
  directly against the empirical rates in 1/h — do NOT normalise to peak for the headline
  metric.
- UNIT CHECK first: confirm the empirical `r` in ExpGrowth.csv is specific growth rate in
  1/h (values ~0.1–4 are consistent with 1/h; convert doubling-time-based entries via
  mu = ln(2)/td if needed). FLAG and DROP any curve whose rate is in relative/arbitrary
  units (not convertible to 1/h) — those can't be validated on magnitude; record which
  were dropped and why.
- Map each curve's medium field to a model medium: glucose_minimal -> glucose minimal;
  rich_broth -> LB. Predict each curve's TPC under ITS medium (predictions now DIFFER by
  medium: LB should give higher rmax than minimal). For each curve compute:
  * PRIMARY: absolute R^2 / RMSE on raw mu(T) in 1/h;
  * SECONDARY (diagnostic): peak-normalised shape R^2, and the rising-limb Ea (scale-
    invariant; compare to the ~0.85 eV benchmark).
- Update percurve_table.csv (per-curve medium, abs_R2 [primary], shape_R2 [secondary],
  obs/pred rmax and Topt) and percurve_summary.json (median ABS_R2 and shape_R2 split by
  medium: minimal vs LB; number of curves dropped for units). Report whether LB brings the
  rich curves' absolute magnitude into range and whether the emergent Ea (~0.95 eV) holds.

PART C - plot RAW absolute rates; distinguish minimal vs LB; fix the ablation plot too
- percurve_smallmultiples.png must plot RAW ABSOLUTE rates (y-axis = growth rate, 1/h):
  the model mu(T) line in 1/h against the empirical points in 1/h, each panel on its own
  y-scale so magnitude is visible (NOT normalised to peak). Visually flag each panel by
  medium type — coloured border/title (minimal vs LB), a per-panel medium label, and a
  legend; group the minimal-medium panels together so they stand out. Show the per-panel
  absolute R^2. (A peak-normalised shape version may be kept as a SECONDARY supplementary
  figure, clearly labelled, but the primary is absolute.)
- percurve_R2_distribution.png: show the distribution of ABSOLUTE R^2 (primary), split/
  coloured by medium (minimal vs LB); shape-R^2 distribution can be a secondary panel.
- ALSO switch the ablation comparison figure (strains/eciML1515/outputs/ablation_comparison.png)
  from "norm to peak" to RAW ABSOLUTE growth rate (1/h) vs the empirical points, so the
  ablations are compared on magnitude too.

PART D - document the medium treatment in the report
- Add a clear paragraph/subsection to the validation section stating: the medium is set as
  AVAILABILITY (open the medium's component uptakes, minimal salts default, uptake not
  pinned — the enzyme constraints set it); how many curves fall in each medium class
  (2 glucose-minimal, 24 rich/LB), where the LB definition comes from (MRes media.csv /
  Machado 2018), how curves without a documented/encodable medium (or non-1/h units) are
  handled, and that validation is now on RAW ABSOLUTE rates (1/h), with peak-normalised
  shape kept only as a secondary diagnostic. Report results SPLIT BY MEDIUM: the ABSOLUTE
  (magnitude) fit and the emergent Ea (~0.85 eV benchmark) for the minimal curves vs the LB
  curves, and read any residual magnitude gap as a finding (kcat / sigma / medium-composition
  uncertainty), not a tuning failure. Update the abstract to reflect the absolute (not
  normalised) validation and the emergent magnitude story.
- Re-run assemble.py so the report picks up the updated figures/tables; `quarto render`.

VERIFY (report all)
1. LB medium: number of LB components opened vs missing; emergent rmax under LB vs minimal.
2. Units: number of curves dropped for non-1/h/relative units.
3. Per-curve ABSOLUTE fit (primary): median abs_R2 split by medium (minimal vs LB); does LB
   bring the rich curves' absolute magnitude into a realistic range? Shape_R2 (secondary)
   and the emergent Ea value.
4. Plots show RAW absolute rates and clearly distinguish minimal vs LB panels; the ablation
   figure is now absolute; report re-renders.

CONSTRAINTS
- Medium is availability, never a fit; do not pin uptake rates and do not tune anything to
  the growth curve. Keep magnitude emergent (pool = P_total x f_metab x sigma) and Ea
  emergent (grounded dCp) as in #16.
- Reuse set_medium, the per-curve validation code, compute_tpc; add plumbing only.
- Autonomous; commit in parts: "add LB/rich medium (availability)",
  "per-curve validation with correct medium per curve + minimal/LB split",
  "highlight minimal vs LB in per-curve plots + report medium documentation".
```
