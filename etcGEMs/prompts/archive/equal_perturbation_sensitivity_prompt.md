# Claude Code prompt — equal-perturbation (standardised elasticity) sensitivity analysis + report (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Run on the current complete model
(post #18: unfolding envelope + medium-matched measured sector allocation). Growth only.

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

WHY: the existing sensitivity analysis (see reports/etcgem/report.qmd "What generates TPC
variation") draws a Latin-hypercube sample with HAND-SET, UNEQUAL ranges per input (e.g.
dCp_scale over x0.5-2.0 but budget_scale only x0.7-1.1) and scores influence with a Spearman
RANK correlation. Two problems: (1) a dial we let swing over a wider range can rank higher
just for that reason; (2) rank correlation measures monotonic *consistency*, not the
*magnitude* of the effect, so it cannot say which dial drives the most variation. We want to
answer the structural question "which parts of the model most drive the shape of the TPC" by
moving every dial by an EQUAL, STANDARDISED amount and reporting a dimensionless magnitude
(an elasticity). (The complementary uncertainty-weighted version — ranges set by each
parameter's real measurement uncertainty, to prioritise what to MEASURE — is a deliberate
FOLLOW-UP and is NOT part of this prompt.)

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
src/etcgem/{sensitivity.py,enzyme_cost.py (Perturbation + effective_Topt/dCp),tpc.py
(compute_tpc, TPC.descriptors),cli.py (cmd_sweep / _finalize_sweep),plotting.py},
strains/eciML1515/strain.yaml and strains/eciML1515/experiments/ (the sweep experiment /
sensitivity block), and reports/etcgem/report.qmd (the "What generates TPC variation
(global sensitivity)" section) + reports/etcgem/assemble.py.

GOAL: add an EQUAL-PERTURBATION sensitivity mode that moves every input by the same
standardised step and reports a dimensionless, comparable magnitude (elasticity) per TPC
descriptor, so the ranking reflects the MODEL'S structural leverage, not our range choices.
Keep the existing LHS/Spearman sweep intact (do not delete it); this is an added mode.

PART A - implement standardised (equal) perturbation + elasticities
- Add a function (e.g. run_elasticity in sensitivity.py) that, for the nominal complete
  model, perturbs EACH input by the SAME standardised relative step h (config, default
  h=0.10) around its nominal value, using central finite differences, and computes for each
  TPC descriptor D (Topt, rmax, CTmax, Ea, thermal breadth — whatever TPC.descriptors
  returns) a NORMALISED sensitivity coefficient (elasticity):
        E[D,p] = ( D(p+) - D(p-) ) / ( 2 * h * D_nominal )
  i.e. fractional change in the descriptor per standardised step of the input. Because every
  input uses the SAME h, the E values are directly comparable in magnitude across inputs.
- Handle the units problem explicitly (this is the crux):
  * Multiplicative inputs (nominal 1): topt_scale, dCp_scale, budget_scale, and any sector /
    group-allocation multipliers -> perturb value = 1 +/- h. Relative step = h.
  * Additive-shift input with nominal 0 (dTopt, in K): a percentage of zero is undefined, so
    define a REFERENCE SCALE Δref and perturb by +/- h*Δref. Use the standard deviation of
    the model's per-enzyme Topt distribution (compute it from the loaded enzyme params) as
    Δref; record the numeric Δref used. This makes dTopt move by "h of a natural temperature
    scale", comparable to an h fractional move of a multiplicative dial.
  * Sector fractions if swept as absolute fractions (f_metab etc., in [0,1]): perturb value =
    nominal*(1 +/- h), clamp to valid range, and renormalise the complementary sectors so the
    budget still sums correctly; document the handling.
- Also record the RAW (unnormalised) ΔD alongside each elasticity for context, and the sign.
- Cover the same input set the current sweep uses on the complete model (envelope knobs
  dTopt/topt_scale/dCp_scale, the pool budget_scale, and the sector-allocation knob(s) that
  are active post-#18). If sectors are active, include the medium-matched allocation knob so
  the allocation axis is represented alongside the envelope axis.

PART B - outputs
- elasticity_table.csv: rows = inputs, cols = descriptors, values = normalised sensitivity
  coefficients (plus a companion raw-ΔD table and a signs table, or extra columns).
- reference_scales.json: h, and Δref for dTopt (and any other non-multiplicative input).
- Figures via plotting.py: (i) a tornado/bar chart per key descriptor (rmax and Ea at least,
  ideally all), inputs sorted by |elasticity|, signed colour; (ii) a heatmap of elasticities
  (inputs x descriptors) to sit alongside/replace the old Spearman heatmap. Save nominal
  descriptor values too.
- Wire into the CLI: either a new subcommand (e.g. `etcgem elasticity --strain eciML1515
  --experiment ...`) or a `--mode equal` flag on `sweep`. Keep it config-driven (read h and
  the input list from the experiment's sensitivity block; default h=0.10). Keep backward
  compatibility with the existing sweep.

PART C - run it on the complete model
- Run the equal-perturbation analysis on eciML1515 (current strain.yaml, complete model).
  Report, per descriptor, the ranked inputs by |elasticity|, and specifically: which
  input(s) most drive rmax (expect the allocation/budget axis) and which most drive Ea and
  the thermal limits (expect the envelope knobs dCp/Topt). Note where the elasticity ranking
  DIFFERS from the old Spearman ranking and why (magnitude vs consistency; equal vs unequal
  ranges).

PART D - report
- Rewrite the "What generates TPC variation (global sensitivity)" section:
  * Briefly state the limitation of the previous approach (hand-set unequal ranges; and rank
    correlation captures monotonic consistency, not effect magnitude), so it reads as an
    improvement, not a silent change.
  * Describe the equal-perturbation method plainly and give the equation for E[D,p], the
    standardised step h, and the reference-scale handling for the additive dТopt shift.
  * Present the results: the elasticity heatmap + the per-descriptor tornado(s), and state in
    words which parts of the model drive which features of the curve (allocation -> rmax;
    envelope -> Ea/limits; interactions where seen).
  * State clearly that this answers the STRUCTURAL question ("which parts of the model drive
    the TPC"), and that a complementary UNCERTAINTY-WEIGHTED analysis — varying each input
    across its real measurement uncertainty to prioritise what to measure — is a planned
    follow-up.
  * Keep the old Spearman result only if useful, clearly labelled as the rank-consistency
    view; otherwise move it to the supplementary.
- Update assemble.py to collect the new figures/tables; re-run assemble.py; `quarto render`;
  confirm the PDF builds with no unresolved crossrefs.

VERIFY (report all)
1. h and the dTopt reference scale Δref actually used; confirm every input is perturbed by
   the same standardised step.
2. The elasticity table and, per descriptor, the top-ranked inputs by |elasticity|; the input
   that most drives rmax and the input that most drives Ea.
3. How the ranking compares to the old Spearman ranking (what changed and why).
4. Report re-renders; the new section explains the method, the equation, and the result, and
   flags the uncertainty-weighted version as follow-up.

CONSTRAINTS
- Do NOT change the core MMRT/unfolding/enzyme-cost math; add the analysis + plumbing only.
- Equal-perturbation only. Do NOT implement the uncertainty-weighted ranges here.
- Elasticities are LOCAL (around the nominal); say so. Optionally, as a SECONDARY check only,
  you may add an equal-relative-RANGE global sweep (all inputs +/- the same width) with
  standardised regression coefficients to confirm the local ranking and expose interactions —
  but the finite-difference elasticity is the headline and is what the report leads with.
- The write-up must match the actual implementation (read the code; use the real h, Δref,
  inputs and numbers).
- Autonomous; commit in parts: "equal-perturbation elasticity sensitivity (standardised step)",
  "elasticity outputs + figures + CLI", "report: structural sensitivity via elasticities".
```
