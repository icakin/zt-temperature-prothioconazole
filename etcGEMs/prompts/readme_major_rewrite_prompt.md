# Claude Code prompt — major rewrite of the top-level README: how the repo works (practical + architecture), with a high-level science overview that points to the report (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). DOCS-ONLY: rewrite README.md. No code, model,
analysis, or report-content changes. The README's job is to explain HOW THE REPOSITORY WORKS — the
code, the CLI, the config model, how the modules link, how to run and reproduce things — with only a
brief high-level science overview that points to the report for the detail.

NOTE TO USER: launch in an auto-approving mode.

CONTEXT: the repo was just restructured (configs/ consolidated, outputs curated, legacy removed,
prompts archived). The README must reflect the CURRENT post-restructure layout, and must be ACCURATE
to the code (read it — do not guess command names, module roles, or paths).

---

```
Work AUTONOMOUSLY; commit at the end; print a summary. READ FIRST so every claim is accurate:
- pyproject.toml (package name, the `etcgem` console entry point, deps), requirements.txt.
- src/etcgem/cli.py (ALL subcommands, their args, what each writes) and config.py (the config
  resolution: defaults.yaml + strain.yaml + experiment overlay -> resolved_config.yaml).
- Every module in src/etcgem/*.py — open each and note its ONE-LINE role and what it imports/feeds,
  so the module map and data-flow are real: mmrt, unfolding, dltkcat, providers, enzyme_cost,
  sectors, proteome_alloc, tpc, sensitivity, decomposition, control, calibration (+ calibration_multi),
  validation, thermal_sampling, dissect, plotting, config, cli.
- configs/defaults.yaml + one configs/experiments/*.yaml + strains/eciML1515/strain.yaml (the 3-layer
  config model).
- reports/etcgem/assemble.py + report.qmd (how outputs -> assets -> rendered report).
- The current README.md and docs/RUNBOOK.md (salvage anything still correct; supersede the rest).
- Confirm the live directory layout (post-restructure) before drawing the tree.

WRITE README.md with these sections (prose + tables + one diagram; keep the science SHORT):

1. TITLE + one-paragraph WHAT THIS IS: an enzyme- and temperature-constrained genome-scale metabolic
   model (etcGEM) of E. coli K-12 MG1655 (eciML1515) that predicts growth thermal performance curves
   (TPCs). One or two sentences on the aim.

2. SCIENTIFIC OVERVIEW (HIGH LEVEL, ~2-3 short paragraphs, NOT the detail): what an etcGEM is in one
   line; the two-stage logic (an "emergent" a-priori model built from independent data -> Bayesian
   calibration to a measured TPC to see what the data demand); and one sentence on the headline
   result. Then explicitly POINT to reports/etcgem/report.qmd and supplementary.qmd for the full
   science, equations, and results. Do NOT reproduce equations or detailed findings here.

3. REPOSITORY LAYOUT: an accurate directory tree (top level + key subdirs) with a one-line purpose
   for each: configs/ (defaults, examples, experiments), src/etcgem/ (the package), strains/
   eciML1515/ (data: model, media, thermal, proteomics, dltkcat; + outputs/ + _archive/), reports/
   etcgem/ (qmd + assemble.py + assets), prompts/ (+ archive/), docs/, and the external/gitignored
   DLTKcat/ and refs/.

4. HOW IT WORKS (the core section — practical + architecture):
   a. THE CONFIG MODEL: explain the 3 layers — configs/defaults.yaml (method defaults) +
      strains/NAME/strain.yaml (organism) + configs/experiments/EXP.yaml (method overlay) — merged by
      the CLI into the resolved_config.yaml dumped in each output dir. State precedence/how to add a
      new strain or experiment.
   b. THE CLI: a TABLE of every `etcgem` subcommand (from cli.py) — command, what it does, key args
      (--strain/--experiment/--fits/etc.), and the output dir it writes. Cover the strain-only
      commands and the strain+experiment commands. Give copy-pasteable example invocations.
   c. MODULE MAP + DATA FLOW: a TABLE of each src/etcgem module with its one-line role, THEN a short
      data-flow description (and a Mermaid diagram) showing how they link, e.g.: providers (load
      GECKO model + set medium) -> thermal parameterisation [mmrt kcat(T) + unfolding f_N(T)/Tm +
      dltkcat] -> enzyme_cost (proteome pool) + sectors (sector partition + growth law) + proteome_
      alloc -> tpc (TPC engine + descriptors) -> {sensitivity, decomposition, control} and calibration
      (emcee) with validation against measured curves; plotting + config + cli as support. Make the
      arrows reflect the ACTUAL imports.

5. INSTALL & QUICKSTART: `pip install -e .` (src layout), the console command `etcgem`, a minimal
   worked example (build or a quick sweep) and where its outputs land. Note the OPTIONAL Gurobi
   solver for speed (academic licence; big speedup on the calibration) and the GLPK fallback.

6. REPRODUCING THE REPORT: the outputs -> assemble.py -> `quarto render` flow; which output dirs feed
   the report; how to rebuild report.pdf + the supplement. Mention that the canonical current outputs
   are the ones assemble.py points at (emergent + tuned analyses, the Van Derlinden calibration).

7. DEVELOPMENT WORKFLOW (brief): the project is built by writing prompts in prompts/ that are executed
   autonomously; executed prompts live in prompts/archive/ (see prompts/README.md). One or two lines.

8. DATA & EXTERNAL DEPENDENCIES (brief): where the raw data lives (strains/eciML1515/{model,media,
   thermal,proteomics,dltkcat}); DLTKcat/ and refs/ are external/gitignored; point to the report for
   data provenance (kcat(T) via DLTKcat, Topt via Li-Engqvist, Tm via the Leuenberger meltome, etc.).

STYLE: practical and precise; the reader should understand how the pieces connect and how to run
them. Keep the science section short and defer to the report. Use tables for the CLI and module map.
Use fenced code blocks for commands. Keep paths post-restructure-correct. Do not invent commands,
modules, or options — everything traceable to the code.

VERIFY (report all)
1. Every CLI subcommand and module described actually exists (cross-checked against cli.py and
   src/etcgem/*.py); no invented commands/paths.
2. The directory tree matches the current post-restructure layout.
3. The config-model, CLI table, module map + data-flow diagram, install/quickstart, and reproduce-
   the-report sections are all present and accurate.
4. Science overview is high-level and points to reports/etcgem/report.qmd + supplementary.qmd.
5. Any Mermaid diagram renders (valid syntax).

CONSTRAINTS
- README.md only (may lightly update docs/RUNBOOK.md to cross-reference if it now overlaps; otherwise
  leave it). No code/model/report-content changes.
- Accurate to the code; post-restructure paths; no invented details.
- Autonomous; single commit: "docs: major README rewrite — repo architecture, CLI, config model, module map + data flow; high-level science pointing to the report".
```
