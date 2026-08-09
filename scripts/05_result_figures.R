# =============================================================================
# 05_result_figures.R — Descriptive result figures (BEFORE the Bayesian models)
# -----------------------------------------------------------------------------
# Renders the raw/derived result figures so you can eyeball them before spending
# hours on MCMC in 04. Uses only 03's outputs (results_plot.rds) — no models.
#
# Produces, for growth, respiration, biomass-corrected growth/respiration, CUE,
# and respiration/growth:
#   - TPC scatter (rate vs Temperature, coloured by Dose, with per-dose mean line)
#   - Boxplots by Dose, faceted by Temperature
# Everything is written to figures_dir as PNGs plus a combined multi-page PDF:
#   result_figures_preBayes.pdf
#
# RUN (after 03, before 04):
#   source("05_result_figures.R")
# =============================================================================

# Source shared config (works interactively in RStudio, via source(), or from CLI)
.this_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                 rstudioapi::isAvailable() &&
                 nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}
source(file.path(.this_dir, "00_config.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---- Load 03 outputs --------------------------------------------------------
rp_path <- file.path(models_dir, "results_plot.rds")
if (!file.exists(rp_path))
  stop("results_plot.rds not found in ", models_dir,
       ". Run 04_oxygen_fits.R first.")
results_plot <- readRDS(rp_path)

# Dose levels / colours: reuse 03's saved ones for consistency, else rebuild.
dose_levels <- tryCatch(readRDS(file.path(models_dir, "dose_levels.rds")),
                        error = function(e) NULL)
if (is.null(dose_levels))
  dose_levels <- {
    dl <- dose_levels_from_data(as.character(results_plot$Dose))
    c(intersect("Control", dl), setdiff(dl, "Control"))
  }
dose_cols <- tryCatch(readRDS(file.path(models_dir, "dose_cols.rds")),
                      error = function(e) make_dose_colors(dose_levels))

results_plot <- results_plot %>%
  mutate(Dose = factor(as.character(Dose), levels = dose_levels))

# ---- Plot helpers -----------------------------------------------------------
# TPC: rate vs Temperature, coloured by Dose, with a per-dose mean line.
tpc_plot <- function(df, ycol, ylab, title) {
  d <- df %>% filter(is.finite(.data[[ycol]]), .data[[ycol]] > 0)
  ggplot(d, aes(T, .data[[ycol]], colour = Dose, group = Dose)) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.7, alpha = 0.85) +
    geom_point(size = 2, alpha = 0.85) +
    scale_colour_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = title, x = "Temperature (°C)", y = ylab, colour = "Dose") +
    theme_classic(12)
}

# Boxplot by Dose, faceted by Temperature. Each point is labelled with its
# replicate (R1/R2/R3); a fixed-seed jitter keeps label and point aligned.
box_plot <- function(df, ycol, ylab, title) {
  d <- df %>% filter(is.finite(.data[[ycol]]), .data[[ycol]] > 0)
  pj <- position_jitter(width = 0.15, height = 0, seed = 42)
  ggplot(d, aes(Dose, .data[[ycol]], fill = Dose)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_point(position = pj, size = 1.4, alpha = 0.8) +
    geom_text(aes(label = Replicate), position = pj, size = 2.4,
              vjust = -0.7, colour = "grey20", show.legend = FALSE) +
    facet_wrap(~T, scales = "free_y") +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = title, x = NULL, y = ylab) +
    theme_classic(11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# metric column, file tag, axis label, human title
metrics <- list(
  list("growth_fgC_h",          "growth",         "Growth (fg C h^-1)",            "Growth"),
  list("respiration_fgC_h",     "respiration",    "Respiration (fg C h^-1)",       "Respiration"),
  list("growth_C_per_C_h",      "growth_biomass", "Growth (C per C h^-1)",         "Biomass-corrected growth"),
  list("respiration_C_per_C_h", "resp_biomass",   "Respiration (C per C h^-1)",    "Biomass-corrected respiration"),
  list("CUE",                   "cue",            "CUE",                           "Carbon Use Efficiency (CUE)"),
  list("resp_over_growth",      "resp_over_growth","Respiration / Growth",         "Respiration / Growth")
)

# ---- Build, save PNGs, and collect for a combined PDF -----------------------
all_plots <- list()
for (m in metrics) {
  ycol <- m[[1]]; tag <- m[[2]]; ylab <- m[[3]]; ttl <- m[[4]]
  if (!ycol %in% names(results_plot)) {
    message("Skipping ", ycol, " (column not present in results_plot).")
    next
  }
  p_tpc <- tpc_plot(results_plot, ycol, ylab, paste0(ttl, " vs Temperature (raw data)"))
  p_box <- box_plot(results_plot, ycol, ylab, paste0(ttl, " by Dose (raw data)"))

  ggsave(file.path(figures_dir, paste0("preBayes_tpc_", tag, ".png")),
         p_tpc, width = 7.2, height = 4.6, dpi = 150)
  ggsave(file.path(figures_dir, paste0("preBayes_box_", tag, ".png")),
         p_box, width = 8.5, height = 6.0, dpi = 150)

  all_plots[[length(all_plots) + 1L]] <- p_tpc
  all_plots[[length(all_plots) + 1L]] <- p_box
}

# Combined multi-page PDF for quick flip-through
pdf_out <- file.path(figures_dir, "result_figures_preBayes.pdf")
pdf(pdf_out, width = 9, height = 6.2)
for (p in all_plots) print(p)
dev.off()

message("\nPre-Bayesian result figures written to: ", figures_dir)
message("  - preBayes_tpc_*.png / preBayes_box_*.png")
message("  - ", basename(pdf_out), " (all figures, one per page)")
message("Review these before running 06_bayesian_models.R.")
