# Claude Code prompt — repository cleanup & restructure (phased, on a branch, verification-gated) (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Structural cleanup per docs/
RESTRUCTURE_PROPOSAL.md. Do it ON A BRANCH, commit in logical phases, and GATE on the report still
building. Do NOT change any model code behaviour, analysis results, or the report's scientific
content — this is moving/renaming/pruning files and rewiring paths only. Preserve git history with
`git mv`.

NOTE TO USER: launch in an auto-approving mode. Everything happens on a new branch; you review and
merge.

DECISIONS ALREADY MADE (apply these):
- Superseded / quick / diagnostic output runs -> ARCHIVE in-repo under
  strains/eciML1515/outputs/_archive/ (kept on disk, untracked/gitignored — not deleted).
- tpc_pipeline/ (legacy, 0 files tracked) and the root outputs/toy_run/ stray -> REMOVE.
- DLTKcat/ -> KEEP in place (already gitignored).
- Executed historical prompts -> ARCHIVE into prompts/archive/.

---

```
Work AUTONOMOUSLY end to end; commit in phases with clear messages; print a summary. Read first:
docs/RESTRUCTURE_PROPOSAL.md, reports/etcgem/assemble.py (RUNS map + copy rules), reports/etcgem/
report.qmd + supplementary.qmd (asset references), src/etcgem/{config.py,cli.py} (config search
paths + the _out_dir tag logic), strains/eciML1515/{strain.yaml,config.yaml}, and prompts/README.md.

PHASE 0 - SAFETY BASELINE (must pass before any change)
- Create and switch to a branch: `git checkout -b repo-restructure`.
- Establish the known-good baseline: run assemble.py then `quarto render` BOTH report.qmd and
  supplementary.qmd. Copy the two rendered PDFs somewhere outside the tree (e.g. /tmp/baseline_*.pdf)
  to diff against at the end. If the baseline does NOT build, STOP and report — do not restructure a
  broken build.

PHASE 1 - root clutter + scripts
- Delete stray files: .DS_Store, verify_control.log, verify_decomp.log, verify_eci.log, verify_p2.log,
  and any *.bak / *.bak2 (all already gitignored; `git rm --cached` any that are tracked, then delete
  on disk).
- Move collaborator_email_draft.md -> docs/correspondence/ (create the dir).
- Create scripts/ and `git mv push_updates.sh scripts/`.
- Commit: "cleanup: remove stray root files; move helper script to scripts/".

PHASE 2 - configuration consolidation (verify resolution after)
- Consolidate config into configs/: `git mv defaults.yaml configs/defaults.yaml`; create
  configs/examples/ and `git mv configs/example_gecko.yaml configs/toy.yaml configs/examples/`;
  create configs/experiments/ and `git mv experiments/*.yaml configs/experiments/` (then remove the
  now-empty experiments/).
- Update src/etcgem/{config.py,cli.py} config search paths to the new locations (defaults at
  configs/defaults.yaml; experiment overlays at configs/experiments/<name>.yaml). Keep backward
  behaviour identical otherwise.
- Confirm strains/eciML1515/config.yaml is SUPERSEDED by strain.yaml: grep the codebase for any
  reference to "config.yaml"; if nothing loads it (resolution uses strain.yaml + defaults +
  experiments), `git rm` it. If anything DOES depend on it, keep it and note that in the summary.
- VERIFY resolution still works: run one quick CLI command (e.g. `etcgem build --strain eciML1515`
  or the lightest available) and confirm it resolves config + runs without a path error.
- Commit: "config: consolidate defaults/examples/experiments under configs/; drop superseded strain config.yaml".

PHASE 3 - remove legacy
- Remove tpc_pipeline/ (0 files tracked — delete the directory from disk; ensure it is gitignored so
  it cannot reappear) and the root outputs/toy_run/ (`git rm -r outputs/toy_run`; if outputs/ is then
  empty, remove it). Grep the repo first to confirm nothing imports/reads tpc_pipeline or
  outputs/toy_run; if something does, STOP and report.
- Commit: "remove legacy tpc_pipeline scaffold and stray root outputs/toy_run".

PHASE 4 - curate strains/eciML1515/outputs/ (ARCHIVE, do not delete)
- FIRST derive the CANONICAL set: trace which output dirs assemble.py actually copies into
  reports/etcgem/assets AND which of those assets report.qmd / supplementary.qmd actually reference.
  That traced set = KEEP. Everything else in outputs/ = ARCHIVE.
  * Expected KEEP (verify against the trace, do not assume): sweep_default, sweep_dltkcat_ext,
    sweep_calibrated, sweep_sectors, proteome_sectors, anatomy, ablation_* + ablation_*.csv/png,
    validation_trusted, calibration_vanderlinden_v3, elasticity_tuned, decompose_tuned, control_tuned.
  * Expected ARCHIVE: calibration_vanderlinden + _v2, calibration_noll_minimal, calibration_phase1,
    control_control, control_control_quick, decompose_decomposition_quick/recast/sectors,
    elasticity_elasticity, percurve_validation, validation (legacy), sweep_quick,
    sweep_calibrated_quick, model_realism_p1b, ceiling_diagnostic, pool_reconciliation, default,
    build, tpc, fba. (If the trace shows the report still uses "recast" or legacy "validation", KEEP
    those instead — trust the trace over this list.)
- Create strains/eciML1515/outputs/_archive/; move each ARCHIVE dir into it; then `git rm -r --cached`
  the archived paths and add `strains/eciML1515/outputs/_archive/` to .gitignore, so they stay on
  disk (recoverable) but leave the tracked tree.
- RENAME the ugly canonical names with `git mv` and update assemble.py in the SAME commit:
  * calibration_vanderlinden_v3 -> calibration_vanderlinden
  * validation_trusted -> validation   (only after the legacy `validation` dir is archived)
  * leave elasticity_tuned / decompose_tuned / control_tuned and the sweep_* names as-is (already clean).
  Do NOT rename asset files inside assets/ — only the source output dirs — to avoid churn in the .qmd.
- Also fix the CLI tag doubling for future runs: adjust so `<command>_<experiment>` no longer repeats
  the command (rename the experiment overlays so e.g. decompose + "sectors" -> decompose_sectors, not
  decompose_decomposition_sectors); this only affects NEW runs, so it is safe.
- Commit: "outputs: archive superseded/quick/diagnostic runs (_archive, gitignored); rename canonical runs; de-double CLI tags".

PHASE 5 - rewire the report + VERIFY (the critical gate)
- Update assemble.py RUNS to point only at the KEEP/renamed dirs; DELETE the dead entries (the ones
  its own comments mark SUPERSEDED / STALE / RETIRED / legacy, now archived). Update any direct
  asset-path references in report.qmd / supplementary.qmd if a rename requires it.
- Re-run assemble.py; `quarto render` BOTH report.qmd and supplementary.qmd; confirm they build with
  NO unresolved crossrefs and NO missing figures/tables.
- DIFF against the Phase-0 baseline PDFs (page count + spot-check the figures/tables that come from
  renamed dirs). The rendered content must be UNCHANGED. If anything is missing or changed, FIX the
  path wiring (do not alter scientific content). Report the diff result.
- Commit: "report: rewire assemble.py to curated output dirs; prune dead entries; re-render verified".

PHASE 6 - prompts + docs
- Create prompts/archive/. Move EXECUTED historical prompts into it with `git mv`, keeping only
  active/pending prompts at the top of prompts/ (at minimum keep any not-yet-run prompt, e.g. the
  control figure/table consistency fix, and this restructure prompt). Use prompts/README.md's index /
  the git log to decide what has run; when unsure, archive (they stay indexed).
- Update prompts/README.md so the index reflects prompts/archive/ vs active.
- Update the top-level README.md "structure" section to the new layout (configs/, scripts/,
  outputs/_archive, DLTKcat/ + refs/ documented as external/gitignored).
- Commit: "prompts: archive executed prompts; docs: update README structure + prompts index".

PHASE 7 - finish
- Print a SUMMARY: what moved/renamed/archived/removed, the config path changes, the KEEP set derived
  from the trace, the CLI command used to verify config resolution, and the PDF-diff result
  (baseline vs rebuilt = identical). Leave everything on the `repo-restructure` branch for the user
  to review and merge; do NOT merge to main or push.

VERIFY (report all)
1. Baseline built before changes; final report + supplementary rebuild identically (PDF diff clean).
2. Config resolves from the new configs/ locations (CLI command ran clean); superseded config.yaml
   removed only after confirming nothing loads it.
3. tpc_pipeline/ and outputs/toy_run removed; nothing referenced them.
4. Stale/quick/diagnostic outputs archived under _archive/ (on disk, untracked); canonical dirs
   renamed; assemble.py points only at live dirs with no dead entries.
5. Executed prompts under prompts/archive/; READMEs updated; all on the repo-restructure branch, not
   merged/pushed.

CONSTRAINTS
- Branch only; no merge, no push. Preserve history with git mv. Commit per phase.
- NO change to model behaviour, analysis numbers, or report scientific content — paths/files only.
- The report MUST render identically to baseline; that diff is the gate. If a phase can't verify,
  STOP and report rather than pressing on.
- Raw data (model/, media/, thermal/, proteomics/, dltkcat/ inputs) and calibration_vanderlinden_v3
  contents are untouched (v3 is only renamed, not modified). DLTKcat/ kept in place.
```
