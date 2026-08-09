# Claude Code prompt — report: add a "model anatomy at the reference operating point" section + make the global-knob / single-operating-point methodology explicit (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Run AFTER #21 (report clarity pass).
Report + plotting work; NO model-math changes.

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

WHY: readers need a feel for HOW the model works before the sensitivity/decomposition results,
and the methodology of those analyses is currently under-explained. Two specific gaps:
- The sensitivity/decomposition "knobs" are GLOBAL scalars that move the WHOLE per-enzyme
  parameter distribution (dTopt = rigid shift of all enzyme optima; topt_scale = stretch/compress
  their spread; dTm = shift all melting temperatures; dCp_scale = scale all curvatures;
  budget_scale / f_metab / f_maint = whole-cell allocation). They are NOT per-enzyme free
  parameters — each enzyme keeps its grounded value; the knobs only move the distribution. This
  is a low-dimensional (~6 knob) projection of the true per-enzyme space (~2,560 enzymes x 3
  thermal params). The report must say this plainly.
- The elasticity/decomposition run on ONE model TPC at ONE operating point (the glucose-minimal
  strain nominal) over a temperature grid; they do NOT iterate over the ~26 empirical curves.
  The empirical curves are used only in the validation section (each under its own medium). The
  report must make this explicit so readers don't think the sensitivity is an average over curves.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
src/etcgem/{enzyme_cost.py (the per-enzyme arrays _Topt/_Tm/_dCp/_base and how the global knobs
map onto them in _costs/_costs_unfolding; Perturbation),unfolding.py,mmrt.py (per-enzyme
rel kcat(T) / native fraction),tpc.py,plotting.py,providers.py,config.py,cli.py},
strains/eciML1515/strain.yaml, the elasticity outputs (elasticity_elasticity/, incl.
reference_scales.json: the Topt/Tm SDs 4.25 / 6.13 K) and the control-coefficient outputs, and
reports/etcgem/report.qmd + reports/etcgem/assemble.py. Do NOT change model math.

PART A - new figures (add to plotting.py; save under the strain outputs and into report assets)
1. REFERENCE TPC: the model growth TPC at the glucose-minimal strain nominal operating point,
   plotted on RAW ABSOLUTE rate (1/h) over the temperature grid, with the key descriptors marked
   (Topt, rmax, CTmax, CTmin, Ea). State the operating point (glucose-minimal; sector fractions
   f_metab~0.483, f_bio~0.191, f_maint~...). Reuse the existing nominal TPC if available; just
   present it clearly as THE reference curve.
2. ENZYME-PARAMETER DENSITY PLOTS: pull the per-enzyme grounded parameter arrays actually in the
   loaded model (per-enzyme Topt, Tm, dCp; and if easy, kcat_ref/MW/base_cost) and plot their
   DISTRIBUTIONS (histogram/KDE, one panel per parameter) across all ~2,560 enzymes. Annotate
   each with the mean and SD; the Topt and Tm SDs should match the reference scales used by the
   elasticity (4.25 / 6.13 K) — note this connection in the caption (the knobs move these very
   distributions; their SD sets the standardised step).
3. EXAMPLE PER-ENZYME kcat(T) PANEL: for ~10 representative enzymes spanning the Topt/Tm range
   (optionally include a few of the highest thermal-control enzymes from the control analysis,
   labelled), plot each enzyme's temperature response — the relative kcat(T) and/or the native
   fraction f_N(T) — on one panel (or a small-multiple), so the reader sees the per-enzyme
   heterogeneity in optima, curvature and denaturation that underlies the organismal TPC. Label
   each curve with its Topt/Tm.

PART B - new report SECTION, placed AFTER "The model" (the mathematical description) and BEFORE
the validation section
- This section is a concrete DESCRIPTION of the model (its inputs and reference behaviour), not a
  result — it makes the equations tangible before validation. Title it clearly (e.g. "How the
  model works: the reference operating point and its enzyme-level parameters"). Lead with a short
  plain-language summary: before testing or dissecting the model, we open it up at a single fixed
  reference condition.
- Content, in plain prose with the figures:
  * Define the REFERENCE OPERATING POINT: the glucose-minimal strain nominal, with its medium and
    sector fractions. Show the reference TPC (fig A1) and read off its descriptors. This is the
    predicted reference curve that the next (validation) section then tests against data.
  * Show the enzyme-parameter DENSITIES (fig A2): the model is ~2,560 enzymes each with its own
    grounded Topt, Tm and dCp; describe the spread (means/SDs). Make the point that these grounded
    per-enzyme values are held FIXED as the model's baseline.
  * Show the example kcat(T) PANEL (fig A3): what an individual enzyme's temperature response
    looks like, and how they differ across enzymes.
  * FORESHADOW (one sentence only): note that the later sensitivity/decomposition analyses do not
    tune enzymes individually — they move these whole distributions with a few global knobs — with
    a cross-reference forward to that section (full explanation lives there, not here).

PART C - make the global-knob / single-operating-point methodology explicit AT the sensitivity +
decomposition sections (where those analyses live)
- In the "What generates TPC variation" (elasticity) section, add a clear paragraph (cross-
  referencing the anatomy figures in Part B): the knobs are GLOBAL scalars that move the whole
  per-enzyme distribution — dTopt slides the Topt density, topt_scale stretches it, dTm slides the
  Tm density, dCp_scale scales the curvatures, and budget_scale / f_metab / f_maint set whole-cell
  allocation. They are NOT per-enzyme free parameters; each enzyme keeps its grounded value. So
  this is a ~6-dimensional summary of a ~7,700-parameter per-enzyme space; the heterogeneity is
  preserved but not individually tuned. Contrast briefly with the per-enzyme control-coefficient
  analysis, which DOES act enzyme-by-enzyme.
- In both the elasticity and decomposition sections, state plainly that these analyses are
  computed on the SINGLE glucose-minimal reference TPC over a temperature grid, NOT per empirical
  curve; the ~26 empirical curves appear only in the validation section, each under its own
  medium. Keep each knob's plain definition where first used.

PART D - assemble + render + verify
- Update assemble.py to collect the new figures; re-run assemble.py; `quarto render`; confirm the
  PDF builds with no unresolved crossrefs; skim to confirm the new section reads clearly and sits
  between validation and sensitivity.

VERIFY (report all)
1. The three new figures exist and are embedded: reference TPC (absolute, descriptors marked),
   enzyme-parameter densities (Topt/Tm/dCp with mean/SD; Topt/Tm SD ~ 4.25 / 6.13 K), example
   per-enzyme kcat(T) panel (~10 enzymes).
2. The new anatomy section is placed AFTER "The model" (the mathematical description) and BEFORE
   the validation section.
3. The report states explicitly, AT the sensitivity/decomposition sections, (a) the knobs are
   global scalars on the whole enzyme distribution, not per-enzyme, and (b) those analyses run on
   the single glucose-minimal reference TPC, not per empirical curve; the anatomy section
   foreshadows this with a forward cross-reference.
4. PDF builds.

CONSTRAINTS
- Plotting + report only; do NOT change core MMRT/unfolding/enzyme-cost/decomposition math.
- Pull the enzyme parameter values and example kcat(T) curves from the ACTUAL loaded model (do
  not fabricate distributions); numbers/means/SDs in the prose must match the figures.
- Keep one operating point (glucose-minimal strain nominal) consistent with the rest of the
  report (per #21).
- Autonomous; commit in parts: "plots: reference TPC + enzyme-parameter densities + example
  kcat(T) panel", "report: model-anatomy section (operating point + enzyme distributions)",
  "report: make global-knob / single-operating-point methodology explicit".
```
