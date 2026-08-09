#!/usr/bin/env bash
# =============================================================================
# run_all.sh — regenerate ALL manuscript figures (6 main + 11 supplementary).
# Usage:  bash scripts/run_all.sh        (from the project root, or anywhere)
# Requires: R (tidyverse, ggridges, patchwork, brms)
#           Python 3 (numpy, pandas, scipy, matplotlib, adjustText)
# Outputs: figures/Figure{1..6}.{png,pdf} and figures/FigureS{1..11}.{png,pdf}
#
# Full pipeline order (see README.md):
#   01-10  physiology  (R;   scripts/run_all.R runs 01-08)
#   11-17  RNA-seq     (shell + R + Python; needs reads/reference)
#   18-23  GEM         (Python; + etcGEMs CLI steps, see docs/RUNBOOK_mac.md)
#   24-36  figures     (this script)
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
S="$ROOT/scripts"
FIG="$ROOT/figures"
DIAG="$FIG/diagnostics"
mkdir -p "$FIG" "$DIAG"
cd "$ROOT"

echo "== Figure 1 (R): thermal performance + CUE =="
Rscript "$S/24_figure1.R"
echo "== Figure 2 (R): dose-response x temperature + Bliss =="
Rscript "$S/25_figure2.R"
echo "== Figure 3 (Py): drug transcriptional programme =="
python3 "$S/26_figure3.py"
cp "$DIAG/Figure4_drug_mechanism.png" "$FIG/Figure3.png"
cp "$DIAG/Figure4_drug_mechanism.pdf" "$FIG/Figure3.pdf"
echo "== Figure 4 (Py): temperature alone =="
python3 "$S/27_figure4.py"
cp "$DIAG/Figure5_temperature.png" "$FIG/Figure4.png"
cp "$DIAG/Figure5_temperature.pdf" "$FIG/Figure4.pdf"
echo "== Figure 5 + S8 (Py): drug-warming convergence =="
python3 "$S/28_figure5_S8.py"
cp "$DIAG/Figure6_interaction.png" "$FIG/Figure5.png"
cp "$DIAG/Figure6_interaction.pdf" "$FIG/Figure5.pdf"
cp "$DIAG/FigureS3_gene_convergence.png" "$FIG/FigureS8.png"
cp "$DIAG/FigureS3_gene_convergence.pdf" "$FIG/FigureS8.pdf"
echo "== Figure 6 (Py): genome-scale metabolic model =="
python3 "$S/29_figure6.py"
echo "== Figures S1-S3 (R) =="
Rscript "$S/30_figureS1_S3.R"
echo "== Figures S4, S5 (R) =="
Rscript "$S/31_figureS4.R"
Rscript "$S/32_figureS5.R"
echo "== Figures S6, S7 (Py) =="
python3 "$S/33_figureS6_S7.py"
cp "$DIAG/FigureS6_kegg_pathways.png" "$FIG/FigureS6.png"
cp "$DIAG/FigureS6_kegg_pathways.pdf" "$FIG/FigureS6.pdf"
cp "$DIAG/FigureS7_gene_grid.png" "$FIG/FigureS7.png"
cp "$DIAG/FigureS7_gene_grid.pdf" "$FIG/FigureS7.pdf"
echo "== Figures S9, S10 (R) =="
Rscript "$S/34_figureS9.R"
Rscript "$S/35_figureS10.R"
echo "== Figure S11 (Py): GEM diagnostics =="
python3 "$S/36_figureS11.py"

echo "All figures (6 main + 11 supplementary) written to $FIG"
