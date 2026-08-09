## =============================================================================
## 54_cyp51_abundance_check.R — Stage 1, S1.4
## Stress test of the constant-CYP51-abundance assumption stated in the Methods.
##
## INTERPRETATION IS EXPLICITLY ASYMMETRIC (spec section 5): a strong temperature
## effect would CHALLENGE the assumption. Absence of an effect does NOT establish
## constant CYP51 protein abundance — transcript is not protein, three
## temperatures is not seven, and the test has limited power. A null result
## merely fails to contradict. This check carries no decision gate.
##
## RUN FROM the master "X0123 copy" folder:
##     Rscript scripts/54_cyp51_abundance_check.R
## Outputs -> tables/revision/predictive/
## =============================================================================
suppressPackageStartupMessages({library(dplyr); library(readr)})

OUT <- "tables/revision/predictive"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "54_cyp51_abundance_check.log"), split = TRUE)

## CYP51/ERG11 in Z. tritici IPO323. The GEM reaction r_0317 is carried by
## Mycgr3G110231 (from the DLTKcat enzyme mapping).
GENE <- "Mycgr3G110231"

contrasts <- list("21 vs 15 C (control)" = "tables/rnaseq/DE_temp_21vs15_proth0.csv",
                  "27 vs 15 C (control)" = "tables/rnaseq/DE_temp_27vs15_proth0.csv")

rows <- list()
for (nm in names(contrasts)) {
  D <- read_csv(contrasts[[nm]], show_col_types = FALSE)
  idcol <- names(D)[1]
  hit <- D[grepl(GENE, D[[idcol]], fixed = TRUE), ]
  if (nrow(hit) == 0) {
    cat("!! ", GENE, " not found in ", contrasts[[nm]], "\n", sep = "")
    next
  }
  lfc <- hit$log2FoldChange[1]; padj <- hit$padj[1]
  ## where does this effect sit in the genome-wide distribution?
  gw <- D$log2FoldChange[is.finite(D$log2FoldChange)]
  pct <- 100 * mean(abs(gw) < abs(lfc))
  cat("\n=== ", nm, " ===\n", sep = "")
  cat("  gene            :", hit[[idcol]][1], "\n")
  cat("  log2FC          :", round(lfc, 4), "  (fold change ",
      round(2^lfc, 3), ")\n", sep = "")
  cat("  adjusted P      :", format.pval(padj, digits = 3), "\n")
  cat("  |log2FC| exceeds", round(pct, 1), "% of genome-wide temperature effects\n")
  cat("  genome-wide |log2FC| median", round(median(abs(gw)), 3),
      " 90th pct", round(quantile(abs(gw), 0.9), 3), "\n")
  rows[[nm]] <- data.frame(contrast = nm, gene = hit[[idcol]][1],
                           log2FC = lfc, fold_change = 2^lfc, padj = padj,
                           pct_of_genomewide_effects = pct,
                           significant_at_0.05 = isTRUE(padj < 0.05))
}
R <- bind_rows(rows)
write_csv(R, file.path(OUT, "cyp51_transcript_temperature.csv"))

cat("\n=== interpretation (asymmetric, pre-registered) ===\n")
if (any(R$significant_at_0.05 & abs(R$log2FC) >= 1, na.rm = TRUE)) {
  cat("A strong temperature effect is present. This CHALLENGES the\n",
      "constant-CYP51-abundance assumption and requires a Methods caveat.\n", sep = "")
} else {
  cat("No strong temperature effect on CYP51 transcript is detected.\n",
      "This FAILS TO CONTRADICT the constant-abundance assumption. It does NOT\n",
      "establish it: transcript is not protein, three temperatures is not seven,\n",
      "and the test has limited power. No positive claim follows.\n", sep = "")
}
sink()
