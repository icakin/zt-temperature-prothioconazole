# Claude Code prompt — fix the "which enzymes matter" figure/table: rank BOTH by control coefficient, make them consistent, and drop the uninformative identifiability histogram (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). REPORT-ONLY + a cheap figure re-plot from
EXISTING saved control-coefficient data. Do NOT re-run run_control (the control coefficients are
already computed in control_tuned) and do NOT change the model or re-run the Bayesian/sensitivity/
decomposition. Growth only.

NOTE TO USER: launch in an auto-approving mode.

CONTEXT: In the "Which enzymes matter" section, Fig 8 (control_thermal.png, #fig-ctrlthermal) and
Table 7 (tbl-control) are INCONSISTENT. Table 7 ranks enzymes by `thermal_screen` (usage x thermal
sensitivity — only a CANDIDATE SCREEN), while Fig 8 sorts/plots the CONTROL COEFFICIENTS and is out
of sync with the current run, uses raw reaction IDs, and is sorted by an invisible (~0) series. We
want BOTH to use the CONTROL COEFFICIENT (the actual model-computed influence of each enzyme's
thermal parameters on the envelope), ranked identically, with gene/enzyme names. Also Fig 9
(control_identifiability.png, #fig-ctrlident) is a one-giant-bar histogram and should be removed.

Data already present: strains/eciML1515/outputs/control_tuned/thermal_control_annotated.csv and
reports/etcgem/assets/tables/thermal_control_annotated.csv have columns: rank, gene, enzyme_name,
enzyme_id (UniProt), rxn_id, EC (empty in this GECKO build), thermal_screen, group, and the control
coefficients CC[Topt_C,Topt_i], CC[Topt_C,dCp_i], CC[CT_max_C,Topt_i], CC[CT_max_C,dCp_i],
CC[niche_width_C,Topt_i], CC[niche_width_C,dCp_i].

---

```
Work AUTONOMOUSLY; commit in parts; print a summary. Read first: src/etcgem/control.py
(plot_thermal_control ~L395, plot_identifiability_hist ~L438, how thermal_screen and the CC columns
are produced), reports/etcgem/assemble.py (the control figure/table copy rules), reports/etcgem/
report.qmd (the "Which enzymes matter ..." section incl. tbl-control, #fig-ctrlthermal, #fig-ctrlident
and the surrounding prose + Interpretation section), and the annotated CSVs above. Do NOT re-run
run_control; regenerate the figure from the SAVED control-coefficient CSV.

PART A - define ONE control-coefficient ranking measure (from the existing CC columns)
- Add a single scalar per enzyme, the ENVELOPE CONTROL COEFFICIENT, = the MAXIMUM absolute value
  across the envelope control coefficients {CC[Topt_C,*], CC[CT_max_C,*], CC[niche_width_C,*]} over
  both per-enzyme parameters (Topt_i, dCp_i). Call it e.g. `control_coeff`. (Rationale: it is the
  largest influence that enzyme's thermal parameters have on ANY envelope feature; note that for
  individual enzymes CC on Topt is ~0 while CC on CT_max / niche width carry the signal — state this.)
- Re-sort the annotated table by `control_coeff` DESCENDING and re-number `rank`. Write the re-ranked
  annotated table back (control_tuned + assets). Keep thermal_screen as a column but it is NO LONGER
  the ranking.

PART B - regenerate Fig 8 (control_thermal) from the CC ranking, consistent with the table
- Fix plot_thermal_control (or a small re-plot that reads the re-ranked annotated CSV) so the figure:
  * ranks by `control_coeff` and takes the SAME top-K as the table shows (top ~15), so the figure and
    table list the SAME enzymes in the SAME order;
  * plots `control_coeff` as the bar height (single clean series). If a breakdown is wanted, use a
    grouped bar of the two MEANINGFUL features |CC[CT_max_C,Topt_i]| and |CC[niche_width_C,Topt_i]|
    and DROP the ~0 CC[Topt_C,Topt_i] series that is currently invisible;
  * labels the x-axis with GENE names (fall back to enzyme_name, then rxn_id, if gene is missing),
    NOT raw reaction IDs;
  * title/axis reflect "control coefficient (influence on the thermal envelope)".
- Regenerate the PNG from the CURRENT control_tuned data so it is in sync (fix the staleness); make
  sure assemble.py copies the freshly regenerated figure. Verify the figure's enzymes now MATCH the
  table's top rows exactly.

PART C - update Table 7 (tbl-control) to rank by / display the control coefficient
- tbl-control should read the re-ranked annotated CSV and show: rank, gene, enzyme (enzyme_name),
  UniProt (enzyme_id), and the `control_coeff` value (rename the displayed column to e.g. "control
  coefficient"). Drop the thermal_screen column from the display (or keep it clearly labelled as the
  "screen" only). Update the caption: ranked by CONTROL COEFFICIENT (model-computed influence on the
  thermal envelope), with thermal_screen described only as the candidate-screening step; identities
  joined locally; EC not stored in this GECKO build.

PART D - remove Fig 9 (identifiability histogram)
- The identifiability distribution is a single uninformative bar. REMOVE the figure: delete the
  ![...]{#fig-ctrlident} image line and any @fig-ctrlident cross-references in the prose; the
  identifiability RESULT stays, carried by tbl-ident and the text (e.g. ~30/7680 = 0.4% identifiable).
  Optionally stop assemble.py copying control_identifiability.png. Do not remove tbl-ident.

PART E - update the section + interpretation prose to the control-coefficient framing
- In "Which enzymes matter" and the Interpretation section, replace "thermal-screen control score"
  language with "control coefficient". State the current TOP enzymes by control coefficient (read the
  re-ranked table — expect e.g. the UDP-N-acetylmuramoyl... / glmS type steps; report the ACTUAL top
  few by name), whether control is still concentrated in a few enzymes, and the point that individual
  enzymes barely control T_opt but a few control CT_max / niche width. Keep the identifiability
  argument (most per-enzyme parameters not learnable from growth).

PART F - build + verify
- Re-run assemble.py; quarto render BOTH report and supplementary; confirm no unresolved crossrefs
  (esp. no dangling @fig-ctrlident). Keep the single-supplement structure.

VERIFY (report all)
1. `control_coeff` defined (max |CC| over envelope descriptors x {Topt_i,dCp_i}); annotated table
   re-ranked by it (rank renumbered); thermal_screen demoted to a screen column.
2. Fig 8 and Table 7 now show the SAME enzymes in the SAME order, ranked by control coefficient,
   labelled with gene/enzyme names; the figure is regenerated from current control_tuned (not stale);
   the invisible ~0 Topt series is gone.
3. tbl-control displays gene / enzyme / UniProt / control coefficient; caption updated.
4. Fig 9 removed; no @fig-ctrlident dangling refs; tbl-ident + identifiability prose retained.
5. Section + interpretation prose use "control coefficient"; the ACTUAL new top enzymes named.
6. Both PDFs build; one supplement.

CONSTRAINTS
- Report-only + cheap re-plot from saved CC data; do NOT re-run run_control or any heavy analysis; no
  model change. Identities local (no web).
- Fig 8 and Table 7 MUST be consistent (same measure, same enzymes, same order).
- Autonomous; commit in parts: "control: add envelope control-coefficient measure + re-rank annotated table",
  "report: regenerate Fig 8 from control coefficients (gene labels), consistent with Table 7",
  "report: Table 7 ranks/display by control coefficient; thermal_screen demoted to screen",
  "report: remove uninformative identifiability histogram (Fig 9); keep table + prose".
```
