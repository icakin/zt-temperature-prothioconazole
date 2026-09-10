#!/usr/bin/env Rscript
# =============================================================================
# make_main_figures.R -- build every main-text figure, in order.
#
#   Rscript scripts/make_main_figures.R            # run from the repository root
#   Rscript scripts/make_main_figures.R 3 6        # or rebuild only these figures
#
# Each step runs in its own R session, so one failure does not poison the rest
# and package namespaces cannot collide. Every figure ends up in figures/ under
# its canonical name (FigureN.png / FigureN.pdf), which is what manuscript.qmd
# points at.
# =============================================================================
STEPS <- list(
  list(fig = 1, scripts = c("50_fig1_Ea_panel.R", "51_figure1.R"),
       src = "tables/revision/fig1G/Figure1_withG"),
  list(fig = 2, scripts = "52_figure2.R", src = "tables/revision/fig2/Figure2_redesign"),
  list(fig = 3, scripts = "53_figure3.R", src = "figures/Figure3"),
  list(fig = 4, scripts = "54_figure4.R", src = "figures/Figure4"),
  list(fig = 5, scripts = "55_figure5.R", src = "figures/Figure5"),
  list(fig = 6, scripts = "56_figure6.R", src = "figures/Figure6"))

args <- commandArgs(trailingOnly = TRUE)
want <- if (length(args)) as.integer(args) else vapply(STEPS, `[[`, 0, "fig")

if (!dir.exists("scripts") || !dir.exists("figures"))
  stop("run this from the repository root (the folder containing scripts/ and figures/)")
dir.create("figures", showWarnings = FALSE)

run_one <- function(s) {
  cat(sprintf("\n--- %s ---\n", s))
  st <- system2("Rscript", c(file.path("scripts", s), "."), stdout = "", stderr = "")
  if (st != 0) stop(sprintf("%s exited with status %d", s, st))
}
publish <- function(src, fig) {                 # copy to figures/FigureN.{png,pdf}
  n <- 0
  for (ext in c("png", "pdf")) {
    from <- paste0(src, ".", ext); to <- sprintf("figures/Figure%d.%s", fig, ext)
    if (file.exists(from) && normalizePath(from) != normalizePath(to, mustWork = FALSE)) {
      file.copy(from, to, overwrite = TRUE); n <- n + 1
    } else if (file.exists(to)) n <- n + 1
  }
  n
}

ok <- character(0); bad <- character(0)
for (S in STEPS) {
  if (!S$fig %in% want) next
  res <- try({ for (s in S$scripts) run_one(s); publish(S$src, S$fig) }, silent = TRUE)
  if (inherits(res, "try-error")) { bad <- c(bad, sprintf("Figure %d", S$fig))
    cat(sprintf("!! Figure %d FAILED: %s\n", S$fig, conditionMessage(attr(res, "condition"))))
  } else ok <- c(ok, sprintf("Figure %d", S$fig))
}

cat("\n================ main figures ================\n")
for (S in STEPS) {
  if (!S$fig %in% want) next
  f <- sprintf("figures/Figure%d.png", S$fig)
  cat(sprintf("  Figure %d  %-28s %s\n", S$fig, basename(f),
      if (file.exists(f)) format(file.mtime(f), "%Y-%m-%d %H:%M") else "MISSING"))
}
cat(sprintf("\n%d built, %d failed\n", length(ok), length(bad)))
if (length(bad)) quit(status = 1)
