# Claude Code prompt — P2b: report restructure — ONE supplement, move two sections out, remove moribund views. REPORT-ONLY (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). PURE report surgery on the Quarto sources —
no code, no analyses, no output regeneration. Runs BEFORE P3, and may run concurrently with the
long P2 calibration (P2 owns code + calibration_vanderlinden/ outputs; this owns only the report
.qmd files, so there is no file overlap). Edit ONLY reports/etcgem/report.qmd and
reports/etcgem/supplementary.qmd.

DECISION: there must be exactly ONE supplement — the standalone `supplementary.qmd` (which builds
supplementary.pdf). report.qmd must NOT keep a separate back-matter "Supplementary" section.

NOTE TO USER (concurrency): fine to run alongside P2 (report files are disjoint from P2's code +
calibration outputs). Only the git index is shared; if you want zero contention, run in a separate
worktree, otherwise it's low-risk since the changed files don't overlap.

---

```
Work AUTONOMOUSLY; commit report changes only; print a summary. Read reports/etcgem/report.qmd and
reports/etcgem/supplementary.qmd fully first. Edit ONLY those two files. Do NOT modify any
src/etcgem/*.py, validation.py, calibration code, sectors, or analysis outputs; do NOT re-run
analyses; do NOT run assemble.py (figures already exist in reports/etcgem/assets/); do NOT touch
the validation figure/caption or tbl-val (P2/P3 own it — leave the Erdos mention for P3).

GOAL: one supplement (supplementary.qmd). Move two main-body sections into it, fold the report's
existing back-matter section into it, and delete the moribund views. Because the supplement is a
SEPARATE rendered document, any main-body @cross-reference to a figure/table that moves to the
supplement will break — so those references must be converted to plain-text pointers to the
Supplementary Information (NOT Quarto @refs) or removed, and any essential numbers restated inline.

PART A - consolidate to ONE supplement (supplementary.qmd)
- Move the report's existing back-matter section "# Supplementary: model construction and layer
  contributions" (the ablation build-up: @fig-ablation, @tbl-ablation, @fig-dltkcat) OUT of
  report.qmd and INTO supplementary.qmd as a new section. report.qmd should end at Interpretation
  / Next steps with NO separate supplementary section.

PART B - MOVE two sections from the main body into supplementary.qmd
1. The whole section "# How the model works: the reference operating point and its enzyme-level
   parameters" with its figures @fig-reftpc, @fig-enzdens, @fig-enzkcat.
2. The "**Measured proteome.**" block from the Validation section, with @fig-protsec, @fig-protpred,
   @tbl-proteome.
   Give them clear section numbers in the supplement. Keep their internal figure labels so they
   resolve WITHIN the supplement.

PART C - fix the now-cross-document references (they will break)
- Find every main-body reference to a MOVED figure/table (@fig-reftpc, @fig-enzdens, @fig-enzkcat,
  @fig-protsec, @fig-protpred, @tbl-proteome, @fig-ablation, @tbl-ablation, @fig-dltkcat) and
  either (a) rewrite it as plain text pointing to the Supplementary Information (e.g. "the
  reference operating point (Supplementary Information)"), or (b) drop it where it's not needed.
  Do NOT leave a Quarto @ref pointing across documents.
- Where a main-body sentence relied on a moved figure for a NUMBER, restate the number inline. In
  particular the sensitivity section uses the enzyme Topt/Tm standard deviations (4.25 / 6.13 K) as
  its standardised step scales — keep those values written in the sensitivity text itself rather
  than referring to the moved @fig-enzdens.

PART D - REMOVE entirely (content + every dangling reference)
- The "**Rank-consistency view (Spearman).**" subsection + @fig-sens (Fig 9) + @tbl-sens (Table 5).
- The "**Calibrated uncertainty.**" subsection + @fig-calens (Fig 10) + @tbl-calibrated (Table 6).
- @fig-secsens (Fig 14, sectors_sensitivity) — remove ONLY this figure; KEEP the Proteome-sector
  trade-off section and @fig-sectrade.
- In the supplement: if S1 "Full global-sensitivity table" is the Spearman/LHS sensitivity table,
  REMOVE it too (the Spearman view is gone entirely); keep the control/identifiability/provenance
  supplement tables. If S1 is actually the elasticity table, keep it.
- Clean all dangling mentions: objective O6 ("... calibrated uncertainty ..."), the in-silico-design
  mentions of "Spearman indices" and the "calibrated-uncertainty ensemble", and the elasticity-
  section sentence that contrasts with "the rank-correlation ranking (@fig-sens)". Grep to confirm
  NONE of fig-sens, tbl-sens, fig-calens, tbl-calibrated, fig-secsens remain, and no broken @refs.

PART E - render + verify (no assemble.py)
- quarto render BOTH report.qmd and supplementary.qmd. Confirm both PDFs build with NO unresolved
  cross-references or citations. There must be exactly one supplement document.

VERIFY (report all)
1. report.qmd has NO back-matter supplement; supplementary.qmd is the single supplement and now
   contains: the anatomy section (+3 figs), the Measured-proteome block (+2 figs, 1 table), and
   the model-construction/ablation section — plus its existing control/identifiability/provenance
   tables.
2. All previous main-body references to the moved figures are now plain-text pointers or removed
   (no cross-document @refs); the Topt/Tm SD numbers remain stated in the sensitivity text.
3. Removed everywhere: Spearman subsection + fig-sens + tbl-sens; calibrated-uncertainty subsection
   + fig-calens + tbl-calibrated; fig-secsens; and S1 if it was the Spearman table. Grep clean.
4. report.pdf and supplementary.pdf both build with no unresolved crossrefs.

CONSTRAINTS
- Edit ONLY report.qmd and supplementary.qmd. No code, analyses, assemble.py, or output regen.
- Exactly ONE supplement (supplementary.qmd). No back-matter supplement in report.qmd.
- No cross-document @refs; broken links become prose pointers or are removed; essential numbers
  restated inline. Do NOT touch the validation figure/caption/tbl-val.
- Commit only the .qmd changes: "report: consolidate to a single supplement (fold back-matter in)",
  "report: move anatomy + measured-proteome to the supplement; prose-ify cross-refs",
  "report: remove Spearman + calibrated-uncertainty + sector-rank views and clean references".
```
