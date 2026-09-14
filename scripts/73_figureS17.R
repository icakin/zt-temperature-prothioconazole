#!/usr/bin/env Rscript
# =============================================================================
# 73_figureS17.R -- Figure S17: every dissolved-oxygen trace behind Figure 6c-d.
# Three independent plates (rows; biological replicates) x three conditions
# (columns), all 24 wells per plate. Light line: full raw trace. Solid line: the
# initial-rate window (post-equilibration maximum -> +10 h). Dashed line: the
# linear fit over that window, whose slope is the rate in Figure 6d.
# Inputs : data/oxygen/azide_rep{1,2,3}_Oxygen.xlsx, data/oxygen/azide_layout.csv,
#          tables/aox/azide_well_rates.csv (43_azide_rates.R)
# Outputs: figures/FigureS17.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))
INK <- "#1B2420"; INK2 <- "#555555"

tr <- load_azide_plates(ROOT)
W  <- read.csv(file.path(ROOT, "tables/aox/azide_well_rates.csv"), stringsAsFactors = FALSE)
tr <- tr %>% left_join(W %>% select(plate, well, window_start_h, window_end_h, rate_mgL_h), by = c("plate", "well")) %>%
  mutate(inwin = t_h >= window_start_h & t_h < window_start_h + WIN_H, g = paste(plate, well))

# per-well linear fit over the window (same fit as 43_azide_rates.R)
fitl <- tr %>% filter(inwin, is.finite(o2)) %>% group_by(plate, well, condition, g) %>%
  group_modify(~{ f <- lm(o2 ~ t_h, .x); data.frame(t_h = range(.x$t_h), o2 = predict(f, data.frame(t_h = range(.x$t_h)))) }) %>%
  ungroup()

lab <- W %>% group_by(plate, condition) %>%
  summarise(n = n(), m = mean(rate_mgL_h), s = sd(rate_mgL_h), r2 = mean(r2), .groups = "drop") %>%
  mutate(txt = sprintf("n = %d wells\nrate %.3f ± %.3f mg/L/h\nmean R² = %.2f", n, m, s, r2),
         condition = factor(condition, levels = AZ_ORDER))

fac <- function(df) { df$condition <- factor(df$condition, levels = AZ_ORDER, labels = AZ_LAB[AZ_ORDER])
                      df$plate <- factor(df$plate, labels = paste("Replicate", 1:3)); df }
tr <- fac(tr); fitl <- fac(fitl); lab <- fac(lab)

th <- theme_classic(base_size = 8) +
  theme(axis.line = element_line(linewidth = 0.28, colour = INK),
        axis.ticks = element_line(linewidth = 0.28, colour = INK),
        axis.text = element_text(size = 6.5, colour = INK),
        axis.title = element_text(size = 8, colour = INK),
        strip.background = element_blank(),
        strip.text = element_text(size = 7.5, face = "bold", colour = INK),
        strip.text.y = element_text(angle = 0),
        panel.spacing = unit(2.5, "mm"),
        legend.position = "top", legend.text = element_text(size = 7),
        legend.key.width = unit(6, "mm"), legend.margin = margin(0, 0, 0, 0))

p <- ggplot() +
  geom_line(data = tr, aes(t_h, o2, group = g, colour = condition, linetype = "raw"), linewidth = 0.3, alpha = 0.28) +
  geom_line(data = subset(tr, inwin), aes(t_h, o2, group = g, colour = condition, linetype = "win"), linewidth = 0.45, alpha = 0.9) +
  geom_line(data = fitl, aes(t_h, o2, group = g, linetype = "fit"), colour = INK, linewidth = 0.35) +
  geom_text(data = lab, aes(x = 25.5, y = 6.5, label = txt), hjust = 1, vjust = 0, size = 2.0, colour = INK2, lineheight = 1.0) +
  facet_grid(plate ~ condition) +
  scale_colour_manual(values = setNames(AZ_COL[AZ_ORDER], AZ_LAB[AZ_ORDER]), guide = "none") +
  scale_linetype_manual(name = NULL, values = c(raw = "solid", win = "solid", fit = "22"),
                        breaks = c("raw", "win", "fit"),
                        labels = c("raw trace (all 24 wells per plate)", "initial-rate window (peak → +10 h)", "linear fit")) +
  guides(linetype = guide_legend(override.aes = list(
    colour = c("#888888", "#888888", INK), alpha = c(0.35, 1, 1), linewidth = c(0.5, 0.8, 0.6)))) +
  scale_x_continuous(limits = c(0, 26), breaks = seq(0, 24, 6)) +
  scale_y_continuous(limits = c(6.3, 12.4), breaks = 7:12) +
  labs(x = "Time (h)", y = expression("Dissolved O"[2]*" (mg L"^-1*")")) + th

OUT <- file.path(ROOT, "figures"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "FigureS17.png"), p, width = 7.2, height = 5.0, dpi = 400, bg = "white")
ggsave(file.path(OUT, "FigureS17.pdf"), p, width = 7.2, height = 5.0, device = cairo_pdf, bg = "white")
cat("Figure S17:", nrow(lab), "panels,", length(unique(tr$g)), "wells\n")
