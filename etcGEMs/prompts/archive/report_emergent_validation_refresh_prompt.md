# Claude Code prompt — refresh the emergent-model validation plots on the current (reconciled) model and insert them into the report (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Small refresh: regenerate the EMERGENT
(a-priori, nothing-fit) validation against Van Derlinden on the CURRENT model and update the
report's validation figure/table/numbers. Run BEFORE the P2 v3 Bayesian. No Bayesian, no
sensitivity/decomposition. Growth only.

NOTE TO USER: launch in an auto-approving mode.

CONTEXT: the model has moved since the report's validation section was last generated — P1b
restored T-dependent maintenance (rounded peak) + added the growth law, and P1c reconciled the
redundant enzyme pools into a single budget. The emergent numbers barely changed (rich rmax ~1.04),
but the report should show the CURRENT emergent model, configured as it will be tuned.

---

```
Work AUTONOMOUSLY; commit report + figure changes; print a summary. Read first: src/etcgem/
{validation.py,providers.py (set_medium BHI, reconciled pool),sectors.py (growth-law toggle),
cli.py,plotting.py}, strains/eciML1515/outputs/{validation_trusted,pool_reconciliation}, and
reports/etcgem/report.qmd (validation section + fig-valcurves + tbl-val). Do NOT run the Bayesian,
sensitivity, decomposition, or identifiability; do NOT change model defaults.

PART A - regenerate the emergent validation on the CURRENT model
- Predict the EMERGENT model (nothing fit to growth) at the rich (BHI) operating point, configured
  EXACTLY as it will be tuned in P2 v3: reconciled SINGLE proteome pool (P1c, default on), the
  COUPLED GROWTH LAW ON, T-dependent maintenance restored (rounded peak). Validate against Van
  Derlinden (MG1655, BHI, 7-46 C) on RAW ABSOLUTE rate (1/h).
- Regenerate the emergent-vs-data figure (model line + Van Derlinden points, absolute rate, 7-46 C,
  showing the ROUNDED peak) and the descriptor table (predicted vs observed T_opt, r_max, CT_max,
  E_a, and abs R^2/RMSE). Report the current numbers (expect ~ rich r_max 1.04, T_opt ~37, E_a
  ~0.64 with the growth law, abs R^2 ~0.15).

PART B - insert/update in the report
- Replace the stale validation figure (fig-valcurves) and table (tbl-val) with the refreshed ones,
  and update the surrounding prose to the CURRENT numbers and model description: the emergent model
  reproduces the SHAPE (T_opt, E_a, the rounded peak, CT_max) but under-predicts the absolute rich-
  medium peak (~1.04 vs observed ~2.40) a priori — an honest gap, read as a finding not a failure.
  Keep this framed as the EMERGENT / a-priori result (the Bayesian tuning is a later, separate
  section). Do NOT anticipate the tuned result here.
- Re-run assemble.py so the report picks up the refreshed figure/table; quarto render; confirm the
  PDF builds with no unresolved crossrefs. Keep the single-supplement structure.

VERIFY (report all)
1. The emergent validation was regenerated on the reconciled + growth-law-ON + maintenance-fixed
   model at rich BHI; the current descriptors (T_opt, r_max, CT_max, E_a, abs R^2).
2. The report's fig-valcurves + tbl-val + prose now show the CURRENT emergent numbers and the
   rounded peak; framed as a-priori (not tuned).
3. report.pdf builds; single-supplement structure intact.

CONSTRAINTS
- Emergent (nothing fit) validation + report update only. No Bayesian, sensitivity, decomposition,
  or model-default changes.
- Use the model configured as it will be tuned (reconciled pool, growth law ON, maintenance on).
- Keep emergent vs tuned distinct; do not pre-empt the tuned result.
- Autonomous; commit: "report: refresh emergent validation figure/table on the reconciled +
  growth-law model (rounded peak; a-priori rich r_max ~1.04 vs 2.40)".
```
