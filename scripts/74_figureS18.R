#!/usr/bin/env Rscript
# =============================================================================
# 74_figureS18.R -- Figure S18: every dissolved-oxygen trace behind Figure 7.
# Six plates (rows: temperature x independent culture) x six conditions
# (columns), all 132 fitted wells. Light line: full raw trace (offset for the
# ~41 h handling step removed). Solid line: the fitting interval, drawn on the
# denoised trace when a moving average was chosen in the trim selector (what
# was fitted). Dashed line: the fitted oxygen model
# O(t) = O2_0 + (K/r)(1 - exp(rt)), refitted here with the procedure of
# 46_ptc_sham_rates.R, so the r values printed in each panel are those in
# tables/aox/ptc_sham_well_rates.csv.
# Inputs : data/aox/ptc_sham_{15,27}_Oxygen.csv, tables/aox/ptc_sham_well_rates.csv,
#          tables/aox/ptc_sham_steps.csv
# Outputs: figures/FigureS18.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))
C15 <- "#2166AC"; C27 <- "#B2182B"; INK <- "#1B2420"; INK2 <- "#555555"

W <- read.csv(file.path(ROOT, "tables/aox/ptc_sham_well_rates.csv"), stringsAsFactors = FALSE)
steps <- read.csv(file.path(ROOT, "tables/aox/ptc_sham_steps.csv"), stringsAsFactors = FALSE)
raw <- list(); fitl <- list(); lab <- list()
for (temp in c(15, 27)) {
  d <- read.csv(file.path(ROOT, sprintf("data/aox/ptc_sham_%d_Oxygen.csv", temp)), check.names = FALSE)
  for (i in which(W$temp == temp)) {
    w <- W[i, ]; cv <- sprintf("%s_R%d%s", w$condition, w$replicate, w$well)
    tt <- d$Time
    st <- steps[steps$T == temp & steps$replicate == w$replicate & steps$well == w$well, ]
    yr <- step_correct(tt, d[[cv]], st$step_time_h, st$step_mgL)
    y  <- smooth_ma(tt, yr, w$smooth_h)                       # what 46 fitted
    m <- tt >= w$fit_start & tt <= w$fit_end & is.finite(y)
    f <- fit_o2_model(tt[m], y[m])
    raw[[length(raw) + 1]]  <- data.frame(temp, replicate = w$replicate, condition = w$condition, well = w$well,
                                          t_h = tt/60, o2 = yr, o2s = y, inwin = m)
    fitl[[length(fitl) + 1]] <- data.frame(temp, replicate = w$replicate, condition = w$condition, well = w$well,
                                           t_h = tt[m]/60, o2 = f$pred)
  }
}
raw <- bind_rows(raw); fitl <- bind_rows(fitl)
lab <- W %>% group_by(temp, replicate, condition) %>%
  summarise(txt = { v <- sprintf("%.3f", r_per_h)
                    if (length(v) > 3) paste(paste(v[1:3], collapse = " "), paste(v[-(1:3)], collapse = " "), sep = "\n") else paste(v, collapse = " ") },
            .groups = "drop")
fac <- function(df) {
  df$condition <- factor(df$condition, levels = PS_ORDER, labels = PS_LAB[PS_ORDER])
  df$plate <- factor(sprintf("%d °C, culture %d", df$temp, df$replicate),
                     levels = sprintf("%d °C, culture %d", rep(c(15, 27), each = 3), 1:3))
  df$temp <- factor(df$temp); df
}
raw <- fac(raw); fitl <- fac(fitl); lab <- fac(lab)
raw$g <- paste(raw$plate, raw$condition, raw$well); fitl$g <- paste(fitl$plate, fitl$condition, fitl$well)

th <- theme_classic(base_size = 8) +
  theme(axis.line = element_line(linewidth = 0.28, colour = INK),
        axis.ticks = element_line(linewidth = 0.28, colour = INK),
        axis.text = element_text(size = 6.2, colour = INK),
        axis.title = element_text(size = 8, colour = INK),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.8, face = "bold", colour = INK),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(1.2, "mm"), panel.spacing.y = unit(2, "mm"),
        legend.position = "top", legend.text = element_text(size = 7),
        legend.key.width = unit(6, "mm"), legend.margin = margin(0, 0, 0, 0))

sm_tab <- aggregate(smooth_h ~ temp, W, function(x) paste(unique(x), collapse = "/"))
sm_txt <- paste(sprintf("%s h at %d °C", sm_tab$smooth_h, sm_tab$temp), collapse = ", ")
win_lab <- if (any(W$smooth_h > 0)) sprintf("fitting interval (moving average: %s)", sm_txt) else "fitting interval"
p <- ggplot() +
  geom_line(data = raw, aes(t_h, o2, group = g, colour = temp, linetype = "raw"), linewidth = 0.28, alpha = 0.28) +
  geom_line(data = subset(raw, inwin), aes(t_h, o2s, group = g, colour = temp, linetype = "win"), linewidth = 0.4, alpha = 0.85) +
  geom_line(data = fitl, aes(t_h, o2, group = g, linetype = "fit"), colour = INK, linewidth = 0.32) +
  geom_text(data = lab, aes(x = 64, y = 16.0, label = txt), hjust = 1, vjust = 1, size = 1.7, colour = INK2, lineheight = 0.9) +
  facet_grid(plate ~ condition) +
  scale_colour_manual(values = c(`15` = C15, `27` = C27), guide = "none") +
  scale_linetype_manual(name = NULL, values = c(raw = "solid", win = "solid", fit = "22"),
                        breaks = c("raw", "win", "fit"),
                        labels = c("raw trace (step offset removed)", win_lab, expression("fitted model   "*O(t)==O[0]+(K/r)*(1-e^{rt})*"   (numbers: "*italic(r)*", "*h^-1*", per well)"))) +
  guides(linetype = guide_legend(override.aes = list(colour = c("#888888", "#888888", INK), alpha = c(0.35, 1, 1), linewidth = c(0.5, 0.8, 0.6)))) +
  scale_x_continuous(limits = c(0, 65), breaks = c(0, 24, 48)) +
  scale_y_continuous(limits = c(0, 16.2), breaks = c(0, 5, 10, 15)) +
  labs(x = "Time (h)", y = expression("Dissolved O"[2]*" (mg L"^-1*")")) + th

OUT <- file.path(ROOT, "figures"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "FigureS18.png"), p, width = 7.2, height = 7.4, dpi = 400, bg = "white")
ggsave(file.path(OUT, "FigureS18.pdf"), p, width = 7.2, height = 7.4, device = cairo_pdf, bg = "white")
cat("Figure S18:", nrow(lab), "panels,", nrow(W), "wells\n")
