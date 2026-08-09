# Claude Code prompt — fully EMERGENT TPC (shape + magnitude) validated PER CURVE (autonomous, comprehensive)

Run from the project root (`.../MICROADAPT/etcGEMs`). Run AFTER #15 (complete-model
consolidation) has committed. This changes the model philosophy: nothing is fit to the
growth TPC — both the envelope (Ea) and the magnitude (B0/rmax) EMERGE from
independently-grounded inputs and the specified medium; the empirical curves are used
ONLY to validate, per individual curve.

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

---

```
Work AUTONOMOUSLY end to end; commit in logical parts; print a summary. PRINCIPLE:
the TPC must be a genuine a-priori prediction. Do NOT tune any parameter to the growth
curve — not the pool budget, not dCp/Ea. Ground every input in INDEPENDENT data (or the
experimental medium), let the whole TPC emerge, and validate against INDIVIDUAL
empirical curves (many per-curve fits, not one pooled R²). Read src/etcgem/{config.py,
providers.py,enzyme_cost.py,sectors.py,tpc.py,unfolding.py,cli.py,dltkcat.py},
strains/eciML1515/strain.yaml, and reports/etcgem/report.qmd before changing anything.
Report honestly, including where the emergent prediction fits worse.

PART A - remove the tuning knobs (both magnitude and envelope)
- Magnitude: stop using pool_scale and the growth-calibrated GECKO pool bound as the
  budget. Remove/neutralise pool_scale (set to 1.0 and deprecate).
- Envelope: stop using calibrate-dcp-to-Ea. The strain must NOT carry a default_dCp
  chosen to hit a target Ea. Keep the CLI command but mark it deprecated / not used for
  the emergent model.

PART B - EMERGENT MAGNITUDE: compute the pool budget from measured proteome + literature
- Set the metabolic-enzyme pool budget as  P_total * f_metab * sigma  where:
  * P_total = total protein per gDW (literature ~0.5 g/gDW for E. coli; a config value);
  * f_metab = metabolic-sector proteome fraction FROM the measured proteomics
    (strains/eciML1515/proteomics/), not a hand-set number;
  * sigma = average in-vivo enzyme saturation, an INDEPENDENT literature value
    (~0.4–0.5, Davidi 2016 / Heckmann 2020) carried as a config value with a stated
    range — NOT fitted to the growth curve.
- Correct the model kcats for the in-vivo/in-vitro gap consistently with sigma (do not
  double-count). Document the budget derivation. rmax then EMERGES; report it.

PART C - EMERGENT ENVELOPE: ground per-enzyme dCp; let Ea emerge
- Per-enzyme dCp from INDEPENDENT sources: DLTKcat kcat(T) fits where available, and for
  the remaining enzymes an empirical dCp prior (mean/spread from the MMRT literature —
  measured enzyme dCp values, NOT from any growth TPC). Topt/Tm stay grounded
  (Li–Engqvist / Leuenberger, as in the unfolding model). No dCp value is chosen to hit
  an Ea target. The organismal Ea then EMERGES; report it and compare to the per-curve
  rising limbs AND to the independent benchmark below as a TEST, not a target.
- BENCHMARK (independent, do not fit to it): the mean short-term activation energy for
  mesophilic BACTERIAL GROWTH is ~0.84-0.88 eV (Smith et al. 2019, Nat Commun; median
  0.84) — well above the 0.65 eV MTE value we previously (wrongly) calibrated to, which
  is why the model over-predicts the low-T rising limb. The emergent Ea should land near
  ~0.85 if the grounded per-enzyme curvatures are right; report the gap as a result.

PART D - MEDIUM as an INPUT (per curve), determined from provenance
- Obtain the Smith et al. 2019 TPC database + metadata from
  github.com/smithtp/hotterbetterprokaryotes (paper: Nat Commun 2019,
  10.1038/s41467-019-13109-1). The compiled database records per-curve growth conditions
  "where possible", so use its medium/condition field; where blank, trace the curve's
  source study. Map each ID in strains/eciML1515/thermal/ExpGrowth.csv to its source
  study and its growth medium / carbon source / aerobic-anaerobic status.
- SELECT the validation curves you can encode: prefer E. coli K-12 curves grown in a
  DOCUMENTED, DEFINED (ideally minimal) medium. For each selected curve, set the medium
  as AVAILABILITY (open the exchange reactions for the medium's components — carbon
  source, minimal salts, O2 if aerobic — close the rest); do NOT pin uptake rates (let
  the enzyme-constrained model determine uptake). Keep the medium constant across the
  curve's temperatures. Record which curves were usable and which were dropped (no
  documented/encodable medium, non-K-12, etc.) and why.

PART E - PER-CURVE validation (many fits, not one pooled R^2)
- For EACH selected curve: predict the emergent TPC (same grounded model; only the
  medium differs) at the curve's measured temperatures, and score the fit for THAT curve:
  * shape R^2 / RMSE (each normalised to its own peak), and
  * ABSOLUTE R^2 / RMSE on raw growth rates where the medium is defined enough to expect
    the right magnitude (the emergent-magnitude test).
  Also record per-curve emergent Topt, rmax, CT_max vs the curve's observed values.
- OUTPUTS: a per-curve table (curve id, medium, n, shape_R2, abs_R2, pred vs obs Topt/
  rmax/CT_max) and figures: small-multiples of model-vs-curve for the selected curves,
  and a distribution/summary of per-curve R^2 (so the result is "the model predicts N
  individual curves with this distribution of skill", not a single pooled number).

PART F - report update
- Replace the single global TPC-fit with the PER-CURVE validation: lead with "the
  emergent model (nothing fit to growth) predicts individual E. coli TPCs" and show the
  distribution of per-curve skill + a few example overlays. State plainly where it fits
  well vs poorly (e.g. low-T overestimation) and read the discrepancies as findings
  (missing cold-side physiology; kcat / sigma uncertainty), not as tuning failures.
- Update the Model and Validation sections and the abstract to reflect the emergent
  (unfit) magnitude and Ea and the per-curve validation. Add refs as needed
  (Davidi2016, Heckmann2020 for sigma; Smith2019 for the TPC database). Re-run
  assemble.py and `quarto render`.

VERIFY (report all)
1. No knob fits the growth curve: pool_scale removed, calibrate-dcp unused, budget =
   P_total*f_metab*sigma with all three from data/literature, dCp grounded.
2. Emergent rmax and Ea (values), and how they compare to the selected curves.
3. Per-curve validation: number of curves used, and the distribution of shape_R2 /
   abs_R2. Confirm the comparison is per-curve, not pooled.
4. Report re-renders with the per-curve validation.

CONSTRAINTS
- Nothing may be fit to the growth TPC. Independent-data grounding and the experimental
  medium are the only allowed inputs; carry sigma / P_total as stated literature values
  with ranges, not hidden knobs.
- Do not change core MMRT/unfolding/cost math; add config/plumbing only.
- Autonomous; commit in parts: "remove growth-fit knobs (pool_scale, calibrate-dcp)",
  "emergent magnitude from measured proteome x sigma", "grounded per-enzyme dCp (emergent Ea)",
  "per-curve medium + validation", "report: emergent per-curve validation".
```
