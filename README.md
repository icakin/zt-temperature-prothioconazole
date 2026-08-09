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
            in figures/ (6 main + 14 supplementary)
    37-46   Robustness and identifiability: full-surface prediction, medium
            and anchoring sensitivity, leave-one-temperature-out, AOX
            identifiability, enzyme/proteome cost, shared-control test
            (44_verify_all.py re-derives every modelling number in the paper)
    47-48   Reference-free test of the temperature x dose interaction
            (GAM separable vs tensor-interaction, exact LOO) -> Fig. S14
    49-50   Climate exposure analysis (NOT used in the manuscript; retained
            for reference): ERA5 hourly 2 m air temperature over the European
            fungicide-spray window, binned into the assayed thermal bands.
            The historical change is small relative to the total window, so
            the analysis was kept out of the paper. Requires network access
            on first run; responses cache under data/climate/cache/

## Manuscript

The manuscript is authored in Quarto with references managed in references.bib
(56 entries, author-date citations via @keys). Editing workflow: edit
manuscript.qmd, then `quarto render manuscript/manuscript.qmd` regenerates
Manuscript.docx and Manuscript.pdf with formatted citations. Sequencing reads
are NOT stored in this repository (deposit at ENA; accession in the paper).

## Requirements

R (tidyverse, brms, ggridges, patchwork, minpack.lm, DESeq2, tximport, ashr),
Python 3 (numpy, pandas, scipy, matplotlib, cobra, gseapy, adjustText, emcee),
R (mgcv, for the reference-free interaction test in script 47),
Salmon, eggNOG-mapper. The etcGEMs package installs with `pip install -e etcGEMs`.
