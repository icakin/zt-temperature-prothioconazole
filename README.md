# Temperature x prothioconazole in *Zymoseptoria tritici*

Physiology, transcriptomics and a targeted respiratory-inhibitor experiment,
showing that supra-optimal temperature antagonises azole fungicide action
through a shared thermal-chemical stress response (project 12129, isolate
IPO323).

## Layout

    manuscript/   Quarto project: manuscript.qmd + supplementary.qmd + references.bib,
                  and the rendered outputs (Manuscript.docx/pdf,
                  Supplementary_Materials.docx/pdf).
                  Re-render with:  quarto render manuscript/manuscript.qmd
                  (needs Quarto >= 1.4; for PDF: quarto install tinytex;
                   figures load from ../figures)
    data/         Inputs: oxygen/ (PreSens time series for the dose-response
                  experiment), aox/ (SHAM and nPG inhibitor traces),
                  rnaseq/ (counts, metadata, salmon_quant/),
                  reference/ (genome + annotation; large indexes gitignored),
                  qc/ (sequencing + quant QC; gitignored)
    scripts/      The numbered pipeline (below). Run everything from the PROJECT ROOT.
    tables/       Derived tables: physiology/, rnaseq/, revision/, aox/
    models/       Fitted objects: physiology/ (brms .rds)
    figures/      Final figures: Figure1-6 and FigureS1-S19 (png + pdf)
    docs/         Runbook, pre-registrations, session logs

## Pipeline (run from the project root)

    00-16   Physiology (R): oxygen traces -> growth/respiration rates ->
            Bayesian TPC, Hill dose-response and CUE models, plus sensitivity
            and posterior checks.
            Rscript scripts/run_all.R runs 01-08.
            03 is an interactive Shiny trim selector: it is run by hand and
            writes the fitting intervals the model fits depend on.
            14 is the reference-free test of the temperature x dose
            interaction (GAM, separable vs tensor interaction, exact LOO).

    20-32   RNA-seq (shell / R / Python): reference download, Salmon
            quantification, eggNOG annotation, DESeq2, KEGG GSEA tables,
            module scores, EC50 target audit, CYP51 transcript check.
            25 writes the pathway tables the figure scripts read; it needs the
            eggNOG database, which is gitignored. Its output,
            data/reference/gene_annotation.csv, IS committed, so the figures
            rebuild without re-running annotation.

    40-49   AOX inhibitor experiments.
            41 fits the dissolved-oxygen model to the SHAM and nPG traces on
            hand-selected intervals and writes the culture-level rate tables:
              Rscript scripts/41_sham_rates.R . sham
              Rscript scripts/41_sham_rates.R . npg
            It requires tables/aox/{sham,npg}_fit_windows.csv, exported from
            the trim selector. Without them it stops rather than guessing.
            42 is the pre-specified, hand-free reanalysis of both compounds
            (rule-based windows, swept reliability filter) and writes
            tables/aox/uniform_window_{fits,analysis}.csv.
            43 extracts initial O2-consumption rates from the azide / SHAM
            respirometry plates (Figure 6c-d) and writes
            tables/aox/azide_well_rates.csv and azide_plate_means.csv:
              Rscript scripts/43_azide_rates.R
            aox_common.R holds the SensorDish reader, the oxygen-model fit
            and the rate fit shared by 39-49, 56, 73 and 74.
            39 turns the SHAM and nPG PreSens exports into the wide oxygen
            tables the rest of the AOX pipeline reads, and records the plate
            layout (dose by column, replicate by row) in data/aox/*_layout.csv.
            40 builds the trim-selector inputs for the SHAM experiment so both
            inhibitor experiments and the factorial are fitted by one procedure.
            45-49 are the prothioconazole x SHAM factorial (Figure 6d-f). No
            fitting interval is chosen by hand:
              45  PreSens exports -> wide traces + the ~41 h handling-step table
              46  -> the analysis-ready series: step offset removed, 4 h moving
                  average on the 27 C plates only, and the interval set by rule
                  (36 h from each well's own post-equilibration maximum).
                  Writes ptc_sham_prepared_traces.csv.gz and ptc_sham_windows.csv
              47  fits the model to those series; rate, condition-mean and
                  interaction tables (tables/aox/ptc_sham_*.csv)
              47  refits all 132 curves under eleven alternative interval rules
                  and reports the contrast and a symmetry diagnostic per rule
              49  Shiny reviewer: inspect every fitted curve, reject a failed
                  well, or override an interval by hand. Overrides carry a
                  reason and 46 reports the contrast with and without them.
                    Rscript scripts/49_ptc_sham_review.R
            44 tests whether the periodic incubator component on the 27 C
            SHAM plate affects the rates: every curve is refit after a 2 h
            and a 4 h moving average and the Figure 6b contrast recomputed
            (tables/aox/sham_denoise_{sensitivity,summary}.csv).

    50-71   Figures. Every manuscript figure is generated by R.
              Rscript scripts/make_main_figures.R     # Figures 1-7
              Rscript scripts/make_supp_figures.R     # Figures S1-S19
            Both accept figure numbers to rebuild a subset, e.g.
              Rscript scripts/make_main_figures.R 3 6
            Each runs its steps in separate R sessions and prints a status
            table naming any figure that failed or is missing.
            bash scripts/run_all.sh rebuilds all figures;
            bash scripts/run_all.sh --all runs the physiology pipeline first.

## Reproducing Figure 6

Figure 6 holds the two intervention experiments in the paper, and both rebuild
from the raw traces:

    Panels a-b (SHAM dose response, 15 and 27 C)
    data/aox/sham_{15,27}_Oxygen.csv        raw dissolved oxygen, 21 curves each
    tables/aox/sham_fit_windows.csv         hand-selected fitting intervals
      -> scripts/41_sham_rates.R            fits O(t) = O2_0 + (K/r)(1 - exp(rt))
    tables/aox/sham_culture_rates.csv       one row per curve

    Panels c-d (azide / SHAM respirometry, 27 C)
    data/oxygen/azide_rep{1,2,3}_Oxygen.xlsx  PreSens exports, one per plate
    data/oxygen/azide_layout.csv              well -> condition
      -> scripts/43_azide_rates.R             linear initial rate, peak -> +10 h
    tables/aox/azide_well_rates.csv           one row per well (72)
    tables/aox/azide_plate_means.csv          plate x condition (the paired unit)

      -> scripts/56_figure6.R               culture-level slopes, paired test,
                                          panels a-c of the figure
      -> scripts/72_figureS16.R             every SHAM trace with window and model fit
      -> scripts/73_figureS17.R             every azide trace with window and linear fit

## Reproducing Figure 6d-f

    Prothioconazole x SHAM factorial (15 and 27 C, three independent cultures)
    data/oxygen/ptc_sham_rep{1,2,3}_{15,27}_Oxygen.xlsx   PreSens exports
      -> scripts/45_ptc_sham_prepare.R        wide traces + handling-step table
      -> scripts/46_ptc_sham_prepared_traces.R  step offset, denoising and the
                                              rule-based interval, in one place
    tables/aox/ptc_sham_prepared_traces.csv.gz  analysis-ready series
    tables/aox/ptc_sham_windows.csv             the interval of every well
      -> scripts/47_ptc_sham_rates.R          oxygen-model fits, plate means,
                                              interaction contrast I and I27 - I15
      -> scripts/48_ptc_sham_window_rules.R   the same result under eleven rules
      -> scripts/56_figure6.R                 panels d-f of the figure
      -> scripts/74_figureS18.R               every trace with interval and fit
      -> scripts/75_figureS19.R               the contrast under eleven interval rules
      -> scripts/49_ptc_sham_review.R         curve-by-curve review (Shiny)

The statistics in the figure are recomputed from the rate tables rather than
hard-coded. Per-curve fits, their intervals and diagnostics for all 78
SHAM and nPG curves are in tables/aox/TableS2_per_curve_fits.csv.

## Manuscript

Authored in Quarto, references in references.bib (author-date via @keys).
Edit manuscript.qmd, then `quarto render manuscript/manuscript.qmd` regenerates
Manuscript.docx and Manuscript.pdf. Sequencing reads are NOT stored here
(deposited at ENA; accession in the paper).

## Requirements

R: ggplot2, patchwork, dplyr, tidyr, readr, readxl, ggrepel, scales, minpack.lm,
   brms, posterior, ggdist, ggridges, cowplot, mgcv, DESeq2, tximport, ashr
Python 3: numpy, pandas, scipy, gseapy
Salmon, eggNOG-mapper (RNA-seq preprocessing only)

The figure scripts need only the ggplot2 stack and minpack.lm; brms and the
RNA-seq tools are required to regenerate the upstream tables, not the figures.

## Note on history

Earlier versions of this project included an enzyme- and temperature-
constrained genome-scale metabolic model (ztGEM/ecZtGEM) and its
preregistrations. That analysis was removed from the manuscript; its code,
model files, preregistrations, reports and derived tables are not in the
working tree. They are preserved on the branch `archive/gem` and in the git
history before commit 03621fe.
