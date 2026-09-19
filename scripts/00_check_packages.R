#!/usr/bin/env Rscript
# =============================================================================
# 00_check_packages.R -- report which R packages the pipeline needs and which
# are missing, so a failure is a list rather than a stack trace half way through
# a figure build.
#
#   Rscript scripts/00_check_packages.R           # report only
#   Rscript scripts/00_check_packages.R --install # install what is missing
#
# The figures and the AOX analyses need only the CRAN set. The Bayesian
# physiology stage (06-14, Figures 1, 2 and S1-S3, S9-S11) additionally needs
# brms and its Stan backend, and the RNA-seq stage (23-27, Figures 3-5 and
# S4-S8) needs the Bioconductor set.
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE)
GROUPS <- list(
  `core (all figures)` = c("ggplot2", "dplyr", "tidyr", "patchwork", "scales", "ggrepel",
                           "readr", "readxl", "stringr", "grid", "cowplot", "zoo",
                           "minpack.lm", "ggridges"),
  `Bayesian physiology` = c("brms", "posterior", "bayesplot", "tidybayes", "ggdist", "mgcv"),
  `RNA-seq (Bioconductor)` = c("DESeq2", "clusterProfiler", "EnhancedVolcano", "pheatmap", "UpSetR"),
  `interactive review only` = c("shiny", "rstudioapi"))
missing <- character(0)
for (g in names(GROUPS)) {
  have <- vapply(GROUPS[[g]], requireNamespace, logical(1), quietly = TRUE)
  cat(sprintf("\n%s\n", g))
  for (i in seq_along(GROUPS[[g]]))
    cat(sprintf("  %-18s %s\n", GROUPS[[g]][i], if (have[i]) "ok" else "MISSING"))
  missing <- c(missing, GROUPS[[g]][!have])
}
cat(sprintf("\n%d of %d packages missing\n", length(missing), length(unlist(GROUPS))))
if (length(missing)) {
  cat("  ", paste(missing, collapse = ", "), "\n")
  if ("--install" %in% ARGS) {
    bioc <- intersect(missing, GROUPS[["RNA-seq (Bioconductor)"]])
    cran <- setdiff(missing, bioc)
    if (length(cran)) install.packages(cran, repos = "https://cloud.r-project.org")
    if (length(bioc)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(bioc, ask = FALSE)
    }
  } else cat("  rerun with --install to install them\n")
} else cat("  nothing to do\n")
cat("\nNote: brms also needs a working Stan backend (cmdstanr or rstan).\n")
cat("Figures also need a working cairo device for the PDF output; on macOS that means XQuartz.\n")
