#!/usr/bin/env Rscript
###############################################################################
## 06_sensitivity_outliers.R
##   Sensitivity analysis: re-run the prothioconazole (drug-vs-control) DE
##   contrasts with the lowest-quality libraries removed, and compare against
##   the full-data results to confirm the conclusions are robust.
##
##   Scenarios:
##     full         all 45 samples (reference)
##     drop_S19     - 12129_15_2_R3_S19      (worst on every metric, ~1M reads)
##     drop_S19_36_49 + 12129_21_2_R5_S36, 12129_27_2_R3_S49 (the three weakest)
##
##   For each scenario it refits DESeq2 (~ group) and extracts the 2-vs-0
##   contrast at 15, 21, 27 C, then reports: number of DE genes, AOX induction,
##   and the per-gene log2FC correlation with the full-data result.
##
##   RUN (from project root):  Rscript scripts/06_sensitivity_outliers.R
##   Outputs -> tables/rnaseq/sensitivity/
###############################################################################
suppressPackageStartupMessages(library(DESeq2))
DE <- "tables/rnaseq"
OUT <- file.path(DE, "sensitivity"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dds0 <- readRDS(file.path(DE, "dds.rds"))          # full, tximport-based, filtered
if (!"group" %in% names(SummarizedExperiment::colData(dds0)))
  stop("dds.rds has no 'group' column in colData")

## AOX gene (for a single-gene readout)
anno <- read.csv("data/reference/gene_annotation.csv", stringsAsFactors = FALSE)
aox <- anno$gene_id[grepl("alternative oxidase", tolower(anno$description))][1]

contrasts <- list(proth_at_15 = c("15_2","15_0"),
                  proth_at_21 = c("21_2","21_0"),
                  proth_at_27 = c("27_2","27_0"))
ALPHA <- 0.05; LFC <- 1

run_scenario <- function(drop) {
  d <- dds0[, !(colnames(dds0) %in% drop)]
  d$group <- droplevels(d$group)
  DESeq2::design(d) <- ~ group
  d <- DESeq2::DESeq(d, quiet = TRUE)
  lapply(contrasts, function(ab) {
    res <- DESeq2::results(d, contrast = c("group", ab[1], ab[2]), alpha = ALPHA)
    shr <- tryCatch(DESeq2::lfcShrink(d, contrast = c("group", ab[1], ab[2]),
                                      res = res, type = "ashr", quiet = TRUE),
                    error = function(e) res)
    data.frame(gene_id = rownames(res),
               log2FoldChange = shr$log2FoldChange, padj = res$padj,
               row.names = rownames(res), stringsAsFactors = FALSE)
  })
}

scenarios <- list(
  full          = character(0),
  drop_S19      = c("12129_15_2_R3_S19"),
  drop_S19_36_49 = c("12129_15_2_R3_S19","12129_21_2_R5_S36","12129_27_2_R3_S49"))

message("Refitting DESeq2 for ", length(scenarios), " scenarios ...")
fits <- lapply(scenarios, run_scenario)

## ---- build comparison table ----------------------------------------------
n_sig <- function(df) sum(!is.na(df$padj) & df$padj < ALPHA & abs(df$log2FoldChange) >= LFC)
n_dir <- function(df, s) sum(!is.na(df$padj) & df$padj < ALPHA & df$log2FoldChange*s >= LFC)
rows <- list()
for (sc in names(fits)) for (cn in names(contrasts)) {
  df <- fits[[sc]][[cn]]; ref <- fits[["full"]][[cn]]
  g <- intersect(rownames(df), rownames(ref))
  r <- if (sc == "full") 1 else suppressWarnings(cor(df[g,"log2FoldChange"], ref[g,"log2FoldChange"], use="complete.obs"))
  rows[[length(rows)+1]] <- data.frame(
    scenario = sc, contrast = cn,
    n_DE = n_sig(df), n_up = n_dir(df, 1), n_down = n_dir(df, -1),
    AOX_log2FC = round(df[aox, "log2FoldChange"], 2),
    cor_log2FC_vs_full = round(r, 3))
}
summ <- do.call(rbind, rows)
write.csv(summ, file.path(OUT, "sensitivity_summary.csv"), row.names = FALSE)
message("\n=== Sensitivity summary (DE genes, AOX induction, log2FC concordance with full) ===")
print(summ, row.names = FALSE)

## ---- scatter: full vs drop_S19_36_49 log2FC per contrast ------------------
pdf(file.path(OUT, "sensitivity_log2FC_scatter.pdf"), width = 10, height = 3.6)
par(mfrow = c(1, 3), mar = c(4,4,3,1))
for (cn in names(contrasts)) {
  ref <- fits[["full"]][[cn]]; alt <- fits[["drop_S19_36_49"]][[cn]]
  g <- intersect(rownames(ref), rownames(alt))
  x <- ref[g,"log2FoldChange"]; y <- alt[g,"log2FoldChange"]
  plot(x, y, pch = 16, cex = 0.3, col = "#3333aa55",
       xlab = "log2FC (all 45)", ylab = "log2FC (drop 3 outliers)",
       main = sprintf("%s  (r = %.3f)", cn, cor(x, y, use="complete.obs")))
  abline(0, 1, col = "red", lty = 2)
}
dev.off()

message("\nSaved -> ", OUT, "/sensitivity_summary.csv and sensitivity_log2FC_scatter.pdf")
message("Interpretation: if n_DE and AOX_log2FC stay similar and cor_log2FC_vs_full is high")
message("(> ~0.95), the conclusions are robust to the low-quality libraries.")
