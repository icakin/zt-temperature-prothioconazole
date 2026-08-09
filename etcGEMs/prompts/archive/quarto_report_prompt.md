# Claude Code prompt — set up Quarto reporting (PDF/docx) for the etcGEM project

Self-contained. Run from the project root (`.../MICROADAPT/etcGEMs`). Sets up a
Quarto report that embeds the pipeline's figures and tables and exports to PDF
(and docx) for sharing with collaborators. Mirror the writing conventions already
used in the sibling project at
`.../MICROADAPT/mulltistress_evolution/multistress-evolution/manuscript/natcomms/`
(read it if accessible; the essentials are reproduced below). Do NOT change any
pipeline/scientific code — this only adds a `reports/` tree.

---

```
Set up a Quarto reporting system for this project so results can be written up as a
report, with embedded figures and tables, and rendered to PDF (and docx) to share
with collaborators. Follow the conventions of the sibling project
`../mulltistress_evolution/multistress-evolution/manuscript/natcomms/` (read its
_quarto.yml, manuscript.qmd, supplementary.qmd, header.tex if you can access that
path). Key conventions to reuse:
- a per-target folder holding _quarto.yml + report.qmd (+ supplementary.qmd) +
  header.tex + references.bib + a CSL style; output-dir "_output"; formats pdf +
  docx; crossref with fig-prefix "Fig.".
- FIGURES are pre-generated static PNGs referenced by relative path (NOT recomputed
  in the document). TABLES are rendered by executable chunks that read result CSVs.
- captions are bold-led: "![**Bold title.** caption ...](path){#fig-x fig-pos=htbp}".

DIFFERENCE FROM THE SIBLING: that project is R; ours is Python. Use the Quarto
Jupyter (python3) engine for the table chunks, reading the pipeline's CSV outputs
with pandas. The document must NOT re-run the heavy pipeline — it only embeds
already-generated figures and reads CSVs. Rendering must be fast.

CREATE  reports/etcgem/  with:
- _quarto.yml:
    project: {type: default, output-dir: _output, render: [report.qmd, supplementary.qmd]}
    format:
      pdf: {documentclass: article, colorlinks: true, geometry: margin=2.5cm,
            include-in-header: {file: header.tex}}
      docx: default
    bibliography: references.bib
    csl: nature-communications.csl
    crossref: {fig-prefix: "Fig.", tbl-prefix: "Table"}
    execute: {echo: false, warning: false}
    jupyter: python3
- header.tex: mirror the sibling's — \usepackage{lineno}\linenumbers (keep, useful
  for collaborator review), \usepackage{float}, the float-placement \renewcommands,
  and the caption styling:
    \usepackage{caption}
    \DeclareCaptionLabelSeparator{bar}{\textbf{\enspace\textbar\enspace}}
    \captionsetup[figure]{labelsep=bar, labelfont=bf}
  and a simple title block (title + "Gabriel Yvon-Durocher, University of Exeter"
  + date). No need to copy the sibling's full custom \@maketitle unless trivial.
- references.bib: seed with the method references so citations resolve. Include
  (fetch correct BibTeX/DOIs): Li 2021 Bayesian etc-GEM (Nat Commun) [Li2021];
  Pettersen & Almaas 2023 parameter inference (Sci Rep) [Pettersen2023]; Sanchez
  2017 GECKO (Mol Syst Biol) [Sanchez2017]; Basan 2015 (Nature) [Basan2015]; Scott
  2010 (Science) [Scott2010]; Heckmann 2020 (PNAS) [Heckmann2020]; Davidi 2016
  (PNAS) [Davidi2016]; Hobbs 2013 MMRT (ACS Chem Biol) [Hobbs2013]; Qiu 2024 DLTKcat
  (Brief Bioinform) [Qiu2024]; Machado 2021 CarveMe/polarization [Machado2021].
- nature-communications.csl: copy it from the sibling manuscript folder if
  accessible, else download the Nature Communications CSL.
- assets/  (populated by assemble.py below): figures/ and tables/.
- report.qmd and supplementary.qmd: skeletons described below.

ASSEMBLE SCRIPT  reports/etcgem/assemble.py :
  Copies the chosen run's figures and CSVs into reports/etcgem/assets/ under stable
  names, so the .qmd references stable filenames regardless of experiment/run-dir
  names (mirrors the sibling's curated results/figures/manuscript/ subset). Make the
  source run dirs configurable at the top of the script; defaults for eciML1515:
    sweep     = strains/eciML1515/outputs/default
    dltkcat   = strains/eciML1515/outputs/dltkcat_ext
    decompose = strains/eciML1515/outputs/decompose_decomposition_quick
    control   = strains/eciML1515/outputs/control_control_quick
    calibrated= strains/eciML1515/outputs/sweep_calibrated      # M1.2 (see PREREQS)
    sectors   = strains/eciML1515/outputs/sweep_sectors         # optional (see PREREQS)
  Copy to assets/figures/ (stable names):
    tpc_ensemble.png            <- sweep/tpc_ensemble.png
    sensitivity_heatmap.png     <- sweep/sensitivity_heatmap.png
    descriptor_distributions.png<- sweep/descriptor_distributions.png
    decomp_achievable.png       <- decompose/achievable_ranges.png
    decomp_variance.png         <- decompose/variance_partition.png
    decomp_shapley.png          <- decompose/shapley_effects.png
    control_thermal.png         <- control/thermal_control_bar.png
    control_bottleneck.png      <- control/bottleneck_vs_temperature.png
    control_identifiability.png <- control/identifiability_hist.png
    dltkcat_ensemble.png        <- dltkcat/tpc_ensemble.png
    calibrated_ensemble.png     <- calibrated/tpc_ensemble.png   # calibrated-uncertainty fan (M1.2)
    sectors_sensitivity.png     <- sectors/sensitivity_heatmap.png
  Copy to assets/tables/ (stable names): sweep/descriptors.csv,
    sweep/sensitivity_spearman.csv, sweep/summary.json,
    decompose/decomposition_table.csv, control/thermal_control.csv,
    control/identifiability.csv, calibrated/descriptors.csv,
    calibrated/summary.json, sectors/samples.csv, sectors/descriptors.csv.
  Skip missing sources with a warning; print what was copied.

  DERIVED FIGURES (the pipeline does NOT emit these; generate them in assemble.py
  with matplotlib from the CSVs above, save to assets/figures/):
    calibrated_vs_default.png : for the key descriptors (Topt_C, rmax, CTmax_C,
      niche_width_C) overlay the DEFAULT run's descriptor distribution (hand-set
      LHS ranges) against the CALIBRATED run's (correlated / DLTKcat-posterior
      per-enzyme sampling) — the headline M1.2 comparison of nominal vs calibrated
      uncertainty. Read assets/tables/descriptors.csv and calibrated/descriptors.csv.
    sector_tradeoff.png : from the sectors run, scatter rmax (and/or CTmax_C) vs
      f_metab coloured by f_maint (join samples.csv columns f_metab/f_maint with
      descriptors.csv by row) to show the interior growth optimum (metabolic vs
      biosynthesis co-limitation) that is the whole point of the sector model.
  If a source run is missing, skip its derived figure with a warning (do not fail).
  (An `etcgem report assemble` CLI subcommand wrapping this is a nice-to-have but
  optional.)

report.qmd  (front matter + skeleton; figures via assets, tables via python chunks):
  ---
  title: "Enzyme- and temperature-constrained modelling of thermal performance
    curves in *E. coli* (eciML1515): sensitivity and identifiability"
  abstract: |
    (one paragraph placeholder)
  ---
  Sections (leave prose as short placeholders for the author to expand; wire every
  figure and table to real assets so it renders now):
  1. Objective & model summary (MMRT kcat(T), proteome pool, peak-normalisation).
  2. Nominal TPC on eciML1515 — figure tpc_ensemble.png (@fig-tpc); a table of the
     nominal descriptors read from assets/tables/summary.json (Topt, rmax, CTmin,
     CTmax, breadth, Ea).
  3. What generates TPC variation (global sensitivity) — figures sensitivity_heatmap
     .png (@fig-sens) and descriptor_distributions.png; a table of Spearman indices
     from sensitivity_spearman.csv (@tbl-sens). Interpret vs H1.1 (envelope) / H1.2
     (allocation).
  4. Allocation vs envelope decomposition (H1.3) — figures decomp_achievable.png,
     decomp_variance.png, decomp_shapley.png; a table from decomposition_table.csv
     (@tbl-decomp) giving allocation/envelope/interaction and Shapley fractions.
  5. Per-enzyme control & identifiability — figures control_thermal.png,
     control_bottleneck.png, control_identifiability.png; a top-N table from
     thermal_control.csv and an identifiability summary from identifiability.csv.
     Interpret vs H1.1 (few thermal determinants) and the omics need (H1.2).
  6. DLTKcat-calibrated envelope — figure dltkcat_ensemble.png; brief compare to the
     default-parameter run.
  7. Calibrated thermal-envelope uncertainty (M1.2) — figures calibrated_ensemble.png
     (per-enzyme correlated / DLTKcat-posterior uncertainty fan) and the derived
     calibrated_vs_default.png (nominal hand-set ranges vs calibrated uncertainty for
     Topt/rmax/CTmax/breadth); a table of descriptor median + IQR from the calibrated
     run (@tbl-calibrated, read calibrated/descriptors.csv or summary.json). State the
     shared_fraction (rho) and that mode=posterior draws per-enzyme mean+sd from the
     DLTKcat fits. This is the move from a nominal scan to calibrated uncertainty.
  8. Proteome-sector allocation trade-off (optional; include only if the sectors run
     exists) — figures sector_tradeoff.png (interior growth optimum: metabolic vs
     biosynthesis co-limitation) and sectors_sensitivity.png (Spearman indices over
     f_metab/f_maint/dTopt/dCp_scale). Interpret as the mechanistic allocation axis
     that speaks WP2's maintenance->biosynthesis language.
  9. Interpretation for H1.1-H1.3 and caveats (structural/in-silico; input-range
     dependence of variance fractions; DLTKcat R^2~0.6; one-factor rho is a simple
     stand-in for the true thermostability covariance; sector couplings are
     order-of-magnitude / auto-calibrated).
  Use crossrefs (@fig-..., @tbl-...) and citations ([@Li2021] etc.). Example table
  chunk pattern:
    ```{python}
    #| label: tbl-decomp
    #| tbl-cap: "Variance of each TPC descriptor partitioned into allocation,
    #|   envelope and interaction (Shapley split), from the H1.3 decomposition."
    #| output: asis
    import pandas as pd
    df = pd.read_csv("assets/tables/decomposition_table.csv")
    print(df.round(3).to_markdown(index=False))
    ```
  (Requires `tabulate` for to_markdown; a great_tables version is an optional upgrade.)

supplementary.qmd: the full sensitivity_spearman table, the full thermal_control and
  identifiability tables, and the run's resolved_config.yaml echoed for provenance.
  Same front matter minus the abstract; can reuse header.tex.

PREREQUISITES & RENDER (document in reports/etcgem/README.md and the project RUNBOOK)
- Install: Quarto CLI; a LaTeX toolchain via `quarto install tinytex`; and in the
  project venv `pip install jupyter pandas tabulate` (great_tables optional).
- Workflow: run the pipeline (build/sweep/decompose/control/dltkcat) so the outputs
  exist -> `python reports/etcgem/assemble.py` -> `cd reports/etcgem && quarto
  render` -> _output/report.pdf and report.docx to share.
- IMPORTANT run notes for the new experiments (the report needs their PNG figures):
  * The `*_quick` experiments run with --no-plots and produce NO figures. For the
    report, run the FULL experiments WITH plots so tpc_ensemble.png etc. exist, e.g.
    `etcgem sweep --strain eciML1515 --experiment calibrated` (produces
    strains/eciML1515/outputs/sweep_calibrated with figures). Update assemble.py's
    `calibrated`/`sectors` source dirs to whatever run dirs you actually produce.
  * The `sectors` experiment requires the STRAIN to enable sectors first: set
    `proteome_sectors: {enabled: true, ...}` in strains/eciML1515/strain.yaml (or a
    sector-enabled strain variant) before `etcgem sweep --strain eciML1515
    --experiment sectors`. If you do not run sectors, the report simply omits
    section 8 and its two figures (assemble.py warns and skips).
- Add reports/etcgem/_output/ to .gitignore (keep the .qmd, header.tex,
  references.bib, csl, and assemble.py tracked).

VERIFY (do these; report)
1. `python reports/etcgem/assemble.py` copies the figures + tables that exist and
   warns about any missing source.
2. `cd reports/etcgem && quarto render report.qmd --to pdf` produces
   _output/report.pdf with figures embedded and the tables rendered; then render
   the full project (report + supplementary) to pdf and docx.
3. Cross-references resolve (Fig. 1..n, Table 1..n) and citations render from
   references.bib with the CSL. Fix any unresolved refs.
4. Print the path to the rendered PDF.

CONSTRAINTS
- Report EMBEDS pre-generated pipeline outputs; do not re-run the pipeline in the
  document and do not modify pipeline/scientific code.
- Keep prose as short placeholders (the author will expand); the value here is a
  working, reproducible report that renders current results to a shareable PDF.
- Commit as "add Quarto reporting (PDF/docx) for etcGEM results".
```
