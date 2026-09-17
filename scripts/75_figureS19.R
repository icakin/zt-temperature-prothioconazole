#!/usr/bin/env Rscript
# =============================================================================
# 75_figureS19.R -- Figure S19: the prothioconazole x SHAM interaction under
# every fitting rule.
#
# The size of the interaction depends on the interval over which each oxygen
# trace is fitted, because the vehicle wells reach oxygen depletion sooner than
# the slow drug + SHAM wells. This figure refits all 132 curves under a set of
# rules that make no reference to the result, and reports the contrast under
# each alongside the hand-selected intervals used in Figure 6.
#
#   left   I27 - I15 at 2 mg/L under each rule, with its 95% CI
#   right  the symmetry diagnostic: the median fitted interval of the 27 C
#          drug + SHAM wells divided by that of the vehicle wells. A ratio above
#          one means the slow wells were fitted over a longer span than the
#          wells they are compared with, which inflates the contrast.
#
# Inputs : tables/aox/ptc_sham_window_rules.csv (48_ptc_sham_window_rules.R)
# Outputs: figures/FigureS19.png / .pdf
# Run from the repository root:  Rscript scripts/75_figureS19.R
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
INK <- "#1B2420"; INK2 <- "#555555"

f <- file.path(ROOT, "tables/aox/ptc_sham_window_rules.csv")
if (!file.exists(f)) stop("run scripts/48_ptc_sham_window_rules.R first")
S <- read.csv(f, stringsAsFactors = FALSE)
names(S) <- tolower(names(S))
dcol <- grep("^i27_i15_2$|^d_2$|^diff", names(S), value = TRUE)[1]
pcol <- grep("^p_2$|^p2$|^p$", names(S), value = TRUE)[1]
rcol <- grep("len_ratio|ratio", names(S), value = TRUE)[1]
ncol_ <- grep("^rule$|name", names(S), value = TRUE)[1]
S$rule <- factor(S[[ncol_]], levels = S[[ncol_]][order(S[[dcol]])])

th <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.28, colour = INK),
        axis.ticks = element_line(linewidth = 0.28, colour = INK),
        axis.text = element_text(size = 7.4, colour = INK),
        axis.title = element_text(size = 8.6, colour = INK),
        plot.margin = margin(4, 4, 2, 4))

pA <- ggplot(S, aes(S[[dcol]], rule)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "#B0B4B8", linewidth = 0.3) +
  geom_point(aes(colour = S[[pcol]] < 0.05), size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#8899A6"),
                      labels = c(`TRUE` = "P < 0.05", `FALSE` = "P >= 0.05"), name = NULL) +
  labs(x = expression(italic(I)[27]-italic(I)[15]~"at 2 mg L"^-1), y = NULL) + th +
  theme(legend.position = "top", legend.text = element_text(size = 7))

pB <- ggplot(S, aes(S[[rcol]], rule)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "#B0B4B8", linewidth = 0.3) +
  geom_point(size = 2.2, colour = INK2) +
  labs(x = "interval length, drug + SHAM / vehicle (27 °C)", y = NULL) + th +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

fig <- (pA | pB) + plot_layout(widths = c(1, 0.72)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", colour = INK))

OUT <- file.path(ROOT, "figures"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "FigureS19.png"), fig, width = 7.2, height = 3.6, dpi = 400, bg = "white")
ggsave(file.path(OUT, "FigureS19.pdf"), fig, width = 7.2, height = 3.6, device = cairo_pdf, bg = "white")
cat(sprintf("Figure S19: %d rules; %d of them negative at 2 mg/L\n", nrow(S), sum(S[[dcol]] < 0)))
