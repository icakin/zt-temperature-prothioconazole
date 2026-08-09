# Temperature x prothioconazole in *Zymoseptoria tritici*

Physiology, transcriptomics and genome-scale metabolic modelling showing that
supra-optimal temperature antagonises azole fungicide action through a shared
thermal-chemical stress response (project 12129, isolate IPO323).

## Layout

    manuscript/   Quarto project: manuscript.qmd + supplementary.qmd + references.bib,
                  and the rendered outputs (Manuscript.docx/pdf, Supplementary_Materials.docx/pdf).
                  Re-render with:  quarto render manuscript/manuscript.qmd
                  (needs Quarto >= 1.4; for PDF: quarto install tinytex; figures load from ../figures)
    data/         Inputs: oxygen/ (PreSens time series), rnaseq/ (counts, metadata,
                  salmon_quant/), reference/ (genome + annotation; large indexes
                  gitignored), qc/ (sequencing + quant QC; gitignored)
    scripts/      The complete numbered pipeline (below). Run everything from the
                  PROJECT ROOT.
    tables/       All derived tables: physiology/, rnaseq/, revision/, gem/
    models/       Fitted objects: physiology/ (brms .rds), gem/ (SBML models,
                  scaffolds, CarveFungi drafts)
    figures/      Final figures (Figure1-6, FigureS1-S11; png + pdf) and
                  diagnostics/ intermediates
    etcGEMs/      Enzyme- and temperature-constrained GEM pipeline (has its own
                  git history; add as a submodule or strip its .git before pushing)
    docs/         Runbook, GEM model report, session logs
    _to_delete/   Superseded files parked during reorganisation - review and delete

## Pipeline (run from the project root)

    01-09   Physiology (R): oxygen -> rates -> Bayesian TPC/Hill/CUE models
            (Rscript scripts/run_all.R runs 01-08; 03 is an interactive trim
             selector; 09/09b/09c/10 are sensitivity + posterior checks)
    11-13   RNA-seq preprocessing (shell): reference download, Salmon, eggNOG
    14-17   RNA-seq statistics (R/Python): DESeq2, GSEA, modules, sensitivity
    18-23   GEM (Python): build ztGEM -> GECKO/thermal strain -> condition
            constraints -> Monte-Carlo -> CYP51 safety margin
            (thermal grounding + calibration run via the etcGEMs CLI:
             see docs/RUNBOOK_mac.md)
    24-36   All figures:  bash scripts/run_all.sh  regenerates every figure
            in figures/ (6 main + 11 supplementary)

## Manuscript

The manuscript is authored in Quarto with references managed in references.bib
(54 entries, author-date citations via @keys). Editing workflow: edit
manuscript.qmd, then `quarto render manuscript/manuscript.qmd` regenerates
Manuscript.docx and Manuscript.pdf with formatted citations. Sequencing reads
are NOT stored in this repository (deposit at ENA; accession in the paper).

## Requirements

R (tidyverse, brms, ggridges, patchwork, minpack.lm, DESeq2, tximport, ashr),
Python 3 (numpy, pandas, scipy, matplotlib, cobra, gseapy, adjustText, emcee),
Salmon, eggNOG-mapper. The etcGEMs package installs with `pip install -e etcGEMs`.
