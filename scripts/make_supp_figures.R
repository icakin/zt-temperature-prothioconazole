#!/usr/bin/env Rscript
# =============================================================================
# make_supp_figures.R -- build every supplementary figure, in order.
#
#   Rscript scripts/make_supp_figures.R            # run from the repository root
#   Rscript scripts/make_supp_figures.R 12 13 14   # or rebuild only these
#
# Each step runs in its own R session. Every figure ends up in figures/ as
# FigureSN.png / FigureSN.pdf, which is what supplementary.qmd points at.
# =============================================================================
STEPS <- list(
  list(figs = c(1, 2, 3), scripts = "60_figureS1_S3.R"),
  list(figs = 4,          scripts = "61_figureS4.R"),
  list(figs = 5,          scripts = "62_figureS5.R"),
  list(figs = c(6, 7),    scripts = "63_figureS6_S7.R"),
  list(figs = 8,          scripts = "64_figureS8.R"),
  list(figs = 9,          scripts = "65_figureS9.R"),
  list(figs = 10,         scripts = "66_figureS10.R"),
  list(figs = 11,         scripts = "67_figureS11.R"),
  list(figs = 12,         scripts = "68_figureS12.R"),
  list(figs = 13,         scripts = "69_figureS13.R"),
  list(figs = 14,         scripts = "70_figureS14.R"),
  list(figs = 15,         scripts = "71_figureS15.R"),
  list(figs = 16,         scripts = "72_figureS16.R"),
  list(figs = 17,         scripts = "73_figureS17.R"))

args <- commandArgs(trailingOnly = TRUE)
all_figs <- sort(unique(unlist(lapply(STEPS, `[[`, "figs"))))
want <- if (length(args)) as.integer(args) else all_figs

if (!dir.exists("scripts") || !dir.exists("figures"))
  stop("run this from the repository root (the folder containing scripts/ and figures/)")

run_one <- function(s) {
  cat(sprintf("\n--- %s ---\n", s))
  st <- system2("Rscript", c(file.path("scripts", s), "."), stdout = "", stderr = "")
  if (st != 0) stop(sprintf("%s exited with status %d", s, st))
}

ok <- 0; bad <- character(0)
for (S in STEPS) {
  if (!any(S$figs %in% want)) next
  res <- try(for (s in S$scripts) run_one(s), silent = TRUE)
  if (inherits(res, "try-error")) { bad <- c(bad, S$scripts)
    cat(sprintf("!! %s FAILED: %s\n", S$scripts, conditionMessage(attr(res, "condition"))))
  } else ok <- ok + 1
}

cat("\n============ supplementary figures ============\n")
for (n in all_figs) {
  if (!n %in% want) next
  f <- sprintf("figures/FigureS%d.png", n)
  cat(sprintf("  Figure S%-3d %-24s %s\n", n, basename(f),
      if (file.exists(f)) format(file.mtime(f), "%Y-%m-%d %H:%M") else "MISSING"))
}
cat(sprintf("\n%d step(s) built, %d failed\n", ok, length(bad)))
if (length(bad)) quit(status = 1)
