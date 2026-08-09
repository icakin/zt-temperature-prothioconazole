# Claude Code prompt — update the conceptual "four layers" overview to match the current model + the rewritten equations. REPORT-ONLY (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). PURE report writing on report.qmd. No code,
no analyses, no output regeneration. Safe to run concurrently with the Bayesian (report file only).

CONTEXT: the "Model equations, parameters and methods" section was already comprehensively
rewritten (full equations, hierarchy, departures-from-Li, in-silico design). But the CONCEPTUAL
"four stacked layers" overview in the "# The complete model" section is now STALE and out of sync:
- Layer 3 describes only the MEDIUM-MATCHED STATIC sector allocation and (incorrectly) calls that
  "the growth law". The model now also has the GROWTH-RATE-DEPENDENT, proteome-conserving growth
  law (f_bio(mu)=f_bio_0+slope*mu, f_metab(mu)=f_metab_0-slope*mu) and the RESTORED T-dependent
  maintenance NGAM(T) that rounds the peak — neither is reflected.
- Layer 4 ("everything emerges; nothing is fit to the growth curve") predates the Bayesian tuning
  and quotes stale glucose-minimal numbers (e.g. r_max=0.55 at the glucose-minimal operating point).
The reference operating point is now the RICH (BHI) medium.

---

```
Work AUTONOMOUSLY; commit report changes only; print a summary. Read reports/etcgem/report.qmd —
BOTH the "# The complete model" four-layers overview AND the rewritten "Model equations, parameters
and methods" equation blocks (so the overview aligns with them). Edit ONLY report.qmd (and
supplementary.qmd if strictly needed). Do NOT modify code, run analyses, or touch the validation/
calibration outputs.

Update the conceptual "four layers" overview so it is CURRENT and maps cleanly onto the rewritten
equation blocks (keep it as the plain-language, intuition-first companion — do NOT restate the
math, do NOT delete it):

1. LAYER 3 (proteome allocation) — revise so it distinguishes and includes all of:
   - the MEASURED, medium- and temperature-matched sector fractions (the static allocation) — keep;
   - the GROWTH-RATE-DEPENDENT, proteome-conserving growth law (Scott 2010): as growth rises, the
     cell shifts proteome from the metabolic sector into ribosomes (f_bio up, f_metab down, same
     slope), so the ribosome<->metabolism TRADE-OFF — not a static cap — sets the maximal rate.
     Make clear this is the correct bacterial growth law and that the earlier "medium-matched f_bio
     = growth law" wording was only the static part.
   - the biosynthesis (translation) cap and the TEMPERATURE-DEPENDENT maintenance NGAM(T), noting
     the maintenance now rises with temperature (which rounds the peak into a realistic optimum).
   Cross-reference the corresponding equation blocks rather than repeating them.

2. LAYER 4 — reframe from "everything emerges; nothing is fit to the growth curve" to the CURRENT
   two-part picture:
   - the EMERGENT model is the a-priori prediction — every parameter grounded in independent data
     (thermal parameters, measured allocation, literature sigma/P_total), nothing fit to the growth
     curve; this is what is validated against the exact-strain curve;
   - the BAYESIAN tuning is a SEPARATE, labelled INVERSE analysis layered on top that asks what
     corrections the data demand (not "nothing is ever fit").
   Replace the stale glucose-minimal illustrative numbers (e.g. r_max=0.55) with the current RICH
   (BHI) reference values (read them from the current emergent reference-TPC / validation outputs;
   if unsure of a number, state it qualitatively rather than quoting a stale one).

3. Ensure the four-layer terminology (sector names, "biosynthesis cap", "NGAM(T)", "growth law",
   "emergent vs tuned") matches the equation-block terminology exactly — no drift between the
   intuitive overview and the formal equations.

Keep the overview concise and readable. quarto render report; confirm the PDF builds with no
unresolved crossrefs.

VERIFY (report all)
1. Layer 3 now includes the growth-rate-dependent (coupled, proteome-conserving) growth law AND the
   restored T-dependent NGAM(T), distinct from the static medium-matched allocation.
2. Layer 4 reframed to emergent-prior + Bayesian-inverse; stale glucose-minimal numbers replaced
   with current rich (BHI) values (or stated qualitatively).
3. Four-layer terminology matches the equation blocks (no drift); overview does not restate the math.
4. report.pdf builds with no unresolved crossrefs.

CONSTRAINTS
- Report-only: edit report.qmd (supplementary.qmd only if needed). No code, analyses, or output
  regeneration. Safe alongside the Bayesian.
- Keep the four-layer overview as the intuitive companion (do not delete, do not duplicate the
  equations); align it with the already-rewritten equation blocks.
- Autonomous; commit: "report: update four-layer overview (growth law + NGAM(T) in Layer 3;
  emergent-vs-Bayesian Layer 4; rich reference numbers)".
```
