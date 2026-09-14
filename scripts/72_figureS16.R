#!/usr/bin/env Rscript
# =============================================================================
# 72_figureS16.R -- Figure S16: every dissolved-oxygen trace behind Figure 6a-b.
# Two temperatures (rows) x six SHAM concentrations (columns) x three cultures.
# Light line: full raw trace. Solid line: the hand-selected fitting interval.
# Dashed line: the fitted model O(t) = O2_0 + (K/r)(1 - exp(rt)), refitted here
# with exactly the procedure of 41_sham_rates.R (bounded nlsLM), so the r values
# printed in each panel are the ones in tables/aox/sham_culture_rates.csv.
# 0.35 mM is omitted, as in the main analysis.
# Inputs : data/aox/sham_{15,27}_Oxygen.csv, tables/aox/sham_fit_windows.csv
# Outputs: figures/FigureS16.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(minpack.lm)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
C15 <- "#2166AC"; C27 <- "#B2182B"; INK <- "#1B2420"; INK2 <- "#555555"
DOSES <- c("Control", "0.0125", "0.025", "0.05", "0.1", "0.2")
DLAB  <- c(Control = "Solvent control", `0.0125` = "0.0125 mM", `0.025` = "0.025 mM",
           `0.05` = "0.05 mM", `0.1` = "0.1 mM", `0.2` = "0.2 mM")

resp_model <- function(t, r, K, O2_0) O2_0 + (K/r) * (1 - exp(r * t))
fit_one <- function(tt, yy) {                       # identical to 41_sham_rates.R
  t0 <- tt - min(tt); sl <- median(diff(yy) / diff(t0))
  st <- list(r = 1e-3, K = min(max(abs(sl), 1e-6), 1), O2_0 = yy[1])
  f <- nlsLM(yy ~ resp_model(t0, r, K, O2_0), start = st,
             lower = c(r = 1e-6, K = 1e-10, O2_0 = min(yy) - 1),
             upper = c(r = 0.15,  K = 1,     O2_0 = max(yy) + 1),
             control = nls.lm.control(maxiter = 200, ftol = 1e-12, ptol = 1e-12))
  p <- coef(f); list(r = p[["r"]], pred = resp_model(t0, p[["r"]], p[["K"]], p[["O2_0"]]))
}

win <- read.csv(file.path(ROOT, "tables/aox/sham_fit_windows.csv"), stringsAsFactors = FALSE)
raw <- list(); fitl <- list(); lab <- list()
for (temp in c(15, 27)) {
  d <- read.csv(file.path(ROOT, sprintf("data/aox/sham_%d_Oxygen.csv", temp)), check.names = FALSE)
  for (dose in DOSES) {
    rs <- numeric(0)
    for (rep in 1:3) {
      cv <- sprintf("%s_R%d", dose, rep); y <- d[[cv]]
      w <- win[win$T == temp & win$Dose == dose & win$Replicate == paste0("R", rep), ][1, ]
      m <- d$Time >= w$fit_start & d$Time <= w$fit_end & is.finite(y)
      f <- fit_one(d$Time[m], y[m]); rs <- c(rs, f$r * 60)
      raw[[length(raw) + 1]]   <- data.frame(temp, dose, rep, t_h = d$Time/60, o2 = y, inwin = m)
      fitl[[length(fitl) + 1]] <- data.frame(temp, dose, rep, t_h = d$Time[m]/60, o2 = f$pred)
    }
    lab[[length(lab) + 1]] <- data.frame(temp, dose, vals = paste(sprintf("%.3f", rs), collapse = "  "))
  }
}
fac <- function(df) { df$dose <- factor(df$dose, levels = DOSES, labels = DLAB[DOSES])
                      df$temp <- factor(df$temp, levels = c(15, 27), labels = c("15 °C", "27 °C")); df }
raw <- fac(bind_rows(raw)); fitl <- fac(bind_rows(fitl)); lab <- fac(bind_rows(lab))

th <- theme_classic(base_size = 8) +
  theme(axis.line = element_line(linewidth = 0.28, colour = INK),
        axis.ticks = element_line(linewidth = 0.28, colour = INK),
        axis.text = element_text(size = 6.5, colour = INK),
        axis.title = element_text(size = 8, colour = INK),
        strip.background = element_blank(),
        strip.text = element_text(size = 7.5, face = "bold", colour = INK),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(1.5, "mm"), panel.spacing.y = unit(3, "mm"),
        legend.position = "top", legend.text = element_text(size = 7),
        legend.key.width = unit(6, "mm"), legend.margin = margin(0, 0, 0, 0))

p <- ggplot() +
  geom_line(data = raw, aes(t_h, o2, group = rep, colour = temp, linetype = "raw"), linewidth = 0.3, alpha = 0.25) +
  geom_line(data = subset(raw, inwin), aes(t_h, o2, group = rep, colour = temp, linetype = "win"), linewidth = 0.45, alpha = 0.8) +
  geom_line(data = fitl, aes(t_h, o2, group = rep, linetype = "fit"), colour = INK, linewidth = 0.4) +
  geom_text(data = lab, aes(x = 68, y = 17.6), label = "italic(r)~(h^-1)", parse = TRUE,
            hjust = 1, vjust = 1, size = 2.0, colour = INK2) +
  geom_text(data = lab, aes(x = 68, y = 15.9, label = vals), hjust = 1, vjust = 1, size = 2.0, colour = INK2) +
  facet_grid(temp ~ dose) +
  scale_colour_manual(values = c(`15 °C` = C15, `27 °C` = C27), guide = "none") +
  scale_linetype_manual(name = NULL, values = c(raw = "solid", win = "solid", fit = "22"),
                        breaks = c("raw", "win", "fit"),
                        labels = c("raw trace", "fit window",
                                   expression("fitted model   "*O(t)==O[0]+(K/r)*(1-e^{rt})))) +
  guides(linetype = guide_legend(override.aes = list(
    colour = c("#888888", "#888888", INK), alpha = c(0.35, 1, 1), linewidth = c(0.5, 0.8, 0.6)))) +
  scale_x_continuous(limits = c(0, 69), breaks = c(0, 24, 48)) +
  scale_y_continuous(limits = c(0, 18), breaks = c(0, 5, 10, 15)) +
  labs(x = "Time (h)", y = expression("Dissolved O"[2]*" (mg L"^-1*")")) + th

OUT <- file.path(ROOT, "figures"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "FigureS16.png"), p, width = 7.2, height = 3.6, dpi = 400, bg = "white")
ggsave(file.path(OUT, "FigureS16.pdf"), p, width = 7.2, height = 3.6, device = cairo_pdf, bg = "white")
cat("Figure S16:", nrow(lab), "panels, 36 curves\n")
