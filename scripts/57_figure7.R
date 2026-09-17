#!/usr/bin/env Rscript
# =============================================================================
# 57_figure7.R -- Figure 7: blocking the alternative oxidase potentiates
# prothioconazole at 27 C but not at 15 C.
#   (a) 15 C and (b) 27 C: growth rate relative to the vehicle control for SHAM
#       alone, prothioconazole alone (2 and 4 mg/L) and the combinations; per
#       plate (three independent cultures) with mean +- s.d. For each
#       combination the open symbol is the multiplicative prediction from the
#       two single agents on the same plate (independence); the arrow is the
#       departure from it.
#   (c) Interaction contrast I = log r_PS - log r_P - log r_S + log r_V per
#       replicate at each temperature, paired within replicate, with the mean
#       and 95% CI; 2 mg/L (primary) and 4 mg/L.
#   (d) Respiration parameter K against growth rate r, both relative to vehicle,
#       plate x condition: at 27 C the combinations keep or raise K while r
#       falls, i.e. cells respire but do not multiply.
# Inputs : tables/aox/ptc_sham_condition_means.csv, ptc_sham_interaction.csv,
#          ptc_sham_interaction_diff.csv, ptc_sham_interaction_tests.csv (46)
# Outputs: figures/Figure7.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))
C15 <- "#2166AC"; C27 <- "#B2182B"; INK <- "#1B2420"; INK2 <- "#555555"

th <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.28, colour = INK),
        axis.ticks = element_line(linewidth = 0.28, colour = INK),
        axis.text = element_text(size = 7.6, colour = INK),
        axis.title = element_text(size = 8.8, colour = INK),
        legend.text = element_text(size = 7.2, colour = INK),
        legend.key.size = unit(3.5, "mm"), legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8.4, face = "bold", colour = INK),
        plot.margin = margin(2, 2, 2, 2))

M <- read.csv(file.path(ROOT, "tables/aox/ptc_sham_condition_means.csv"), stringsAsFactors = FALSE)
I <- read.csv(file.path(ROOT, "tables/aox/ptc_sham_interaction.csv"), stringsAsFactors = FALSE)
D <- read.csv(file.path(ROOT, "tables/aox/ptc_sham_interaction_diff.csv"), stringsAsFactors = FALSE)
S <- read.csv(file.path(ROOT, "tables/aox/ptc_sham_interaction_tests.csv"), stringsAsFactors = FALSE)

## ---- a, b: relative growth with independence predictions -------------------
XPOS <- c(S = 1, P2 = 2.2, P2S = 3.0, P4 = 4.2, P4S = 5.0)
XLAB <- c(S = "SHAM", P2 = "2", P2S = "2 + SHAM", P4 = "4", P4S = "4 + SHAM")
rel <- M %>% filter(condition != "V") %>% mutate(x = XPOS[condition])
pred <- I %>% select(temp, replicate, P2S = pred_rel_P2S, P4S = pred_rel_P4S) %>%
  tidyr::pivot_longer(c(P2S, P4S), names_to = "condition", values_to = "pred") %>%
  mutate(x = XPOS[condition])
obs  <- rel %>% filter(condition %in% c("P2S", "P4S")) %>% select(temp, replicate, condition, obs = r_rel)
arr  <- merge(pred, obs) %>% group_by(temp, condition, x) %>%
  summarise(pred = mean(pred), obs = mean(obs), .groups = "drop")
ms <- rel %>% group_by(temp, condition, x) %>% summarise(m = mean(r_rel), s = sd(r_rel), .groups = "drop")
SHP <- c(`1` = 21, `2` = 22, `3` = 24)

panel_rel <- function(tp, tag_col, show_y = TRUE) {
  ggplot() +
    geom_hline(yintercept = 1, linetype = "22", colour = "#B0B4B8", linewidth = 0.3) +
    annotate("rect", xmin = 1.7, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "#F6F6F6") +
    annotate("rect", xmin = 3.7, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = "#F6F6F6") +
    annotate("text", x = c(2.6, 4.6), y = 1.34, label = c("prothioconazole 2 mg/L", "prothioconazole 4 mg/L"),
             size = 2.35, colour = INK2) +
    geom_segment(data = subset(arr, temp == tp), aes(x = x + 0.22, xend = x + 0.22, y = pred, yend = obs, colour = condition),
                 linewidth = 0.5, arrow = arrow(length = unit(1.4, "mm"), type = "closed"), show.legend = FALSE) +
    geom_point(data = subset(arr, temp == tp), aes(x + 0.22, pred), shape = 23, colour = INK, fill = "white", size = 2.3, stroke = 0.5) +
    geom_point(data = subset(rel, temp == tp), aes(x, r_rel, fill = condition, shape = factor(replicate)),
               colour = INK, size = 1.5, stroke = 0.3, alpha = 0.85) +
    geom_errorbar(data = subset(ms, temp == tp), aes(x, ymin = m - s, ymax = m + s, colour = condition), width = 0.16, linewidth = 0.45) +
    geom_point(data = subset(ms, temp == tp), aes(x, m, fill = condition), shape = 21, colour = INK, size = 2.6, stroke = 0.5) +
    annotate("text", x = 0.62, y = 1.34, label = sprintf("%d °C", tp), colour = tag_col, size = 3.2, fontface = "bold", hjust = 0) +
    scale_fill_manual(values = PS_COL, guide = "none") + scale_colour_manual(values = PS_COL, guide = "none") +
    scale_shape_manual(values = SHP, name = "Culture") +
    scale_x_continuous(breaks = XPOS, labels = XLAB, limits = c(0.5, 5.55), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.42), breaks = seq(0, 1.25, 0.25)) +
    labs(x = NULL, y = if (show_y) expression("Growth rate relative to vehicle,"~italic(r)/italic(r)[V]) else NULL) + th +
    theme(legend.position = "none", axis.text.x = element_text(size = 7.0))
}
pA <- panel_rel(15, C15); pB <- panel_rel(27, C27, show_y = FALSE) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

## ---- c: interaction contrast, paired within replicate -----------------------
Il <- I %>% select(temp, replicate, `2` = I_2, `4` = I_4) %>%
  tidyr::pivot_longer(c(`2`, `4`), names_to = "dose", values_to = "I") %>%
  mutate(dose = factor(paste0(dose, " mg/L"), levels = c("2 mg/L", "4 mg/L")),
         xt = ifelse(temp == 15, 1, 2), temp = factor(temp))
Is <- Il %>% group_by(dose, temp, xt) %>%
  summarise(m = mean(I), lo = m + qt(0.025, 2) * sd(I) / sqrt(3), hi = m - qt(0.025, 2) * sd(I) / sqrt(3), .groups = "drop")
lab <- data.frame(dose = factor(c("2 mg/L", "4 mg/L"), levels = c("2 mg/L", "4 mg/L")),
                  txt = c(sprintf("Delta == %.2f", S$mean[1]), sprintf("Delta == %.2f", S$mean[2])),
                  txt2 = c(sprintf("'paired '*italic(P)*' = %.2f'", S$P[1]), sprintf("'paired '*italic(P)*' = %.2f'", S$P[2])))
pC <- ggplot() +
  geom_hline(yintercept = 0, linetype = "22", colour = "#B0B4B8", linewidth = 0.3) +
  geom_line(data = Il, aes(xt, I, group = replicate), colour = "#C4C4C4", linewidth = 0.4) +
  geom_point(data = Il, aes(xt, I, colour = temp, shape = factor(replicate)), fill = "white", size = 1.6, stroke = 0.45) +
  geom_errorbar(data = Is, aes(xt + 0.3, ymin = lo, ymax = hi, colour = temp), width = 0.12, linewidth = 0.5) +
  geom_point(data = Is, aes(xt + 0.3, m, colour = temp), size = 2.4) +
  geom_text(data = lab, aes(x = 0.62, y = -1.62, label = txt), parse = TRUE, hjust = 0, size = 2.6, colour = INK2) +
  geom_text(data = lab, aes(x = 0.62, y = -1.85, label = txt2), parse = TRUE, hjust = 0, size = 2.6, colour = INK2) +
  facet_wrap(~dose, nrow = 1) +
  scale_colour_manual(values = c(`15` = C15, `27` = C27), guide = "none") +
  scale_shape_manual(values = SHP, name = "Culture") +
  scale_x_continuous(breaks = c(1.15, 2.15), labels = c("15 °C", "27 °C"), limits = c(0.55, 2.6)) +
  scale_y_continuous(limits = c(-1.95, 0.45), breaks = seq(-1.5, 0.5, 0.5)) +
  labs(x = NULL, y = expression(atop("Interaction contrast "*italic(I), "(0 = independent; < 0 = SHAM potentiates)"))) +
  th + theme(legend.position = "none", panel.spacing.x = unit(3, "mm"))

## ---- d: respiration parameter vs growth rate ---------------------------------
Md <- M %>% filter(condition != "V") %>% mutate(temp = factor(temp))
pD <- ggplot(Md, aes(r_rel, K_rel)) +
  geom_abline(slope = 0, intercept = 1, linetype = "22", colour = "#B0B4B8", linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = "22", colour = "#B0B4B8", linewidth = 0.3) +
  geom_point(aes(fill = condition, shape = temp), colour = INK, size = 2.0, stroke = 0.35, alpha = 0.9) +
  annotate("text", x = 0.02, y = 1.86, label = "respiring, not multiplying", size = 2.4, colour = INK2, hjust = 0) +
  scale_fill_manual(values = PS_COL, labels = PS_LAB, name = NULL) +
  scale_shape_manual(values = c(`15` = 21, `27` = 24), labels = c("15 °C", "27 °C"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.35), breaks = seq(0, 1.25, 0.25)) +
  scale_y_continuous(limits = c(0.3, 1.92)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2.2), order = 1),
         shape = guide_legend(override.aes = list(fill = "#BBBBBB"), order = 2)) +
  labs(x = expression("Growth rate, "*italic(r)/italic(r)[V]), y = expression("Respiration parameter, "*italic(K)/italic(K)[V])) +
  th + theme(legend.position = c(1.0, 1.0), legend.justification = c(1, 1),
             legend.key.height = unit(2.8, "mm"), legend.spacing.y = unit(0, "mm"),
             legend.text = element_text(size = 6.6))

## ---- layout ---------------------------------------------------------------------
row1 <- pA + pB + plot_layout(widths = c(1, 0.92))
row2 <- pC + pD + plot_layout(widths = c(1.15, 1))
fig <- (row1 / row2) + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", colour = INK), plot.tag.position = "topleft",
        plot.margin = margin(5, 3, 2, 3))
OUT <- file.path(ROOT, "figures"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "Figure7.png"), fig, width = 7.2, height = 5.6, dpi = 600, bg = "white")
ggsave(file.path(OUT, "Figure7.pdf"), fig, width = 7.2, height = 5.6, device = cairo_pdf, bg = "white")
cat(sprintf("Figure 7: I27-I15 at 2 mg/L = %+.3f (95%% CI %+.2f to %+.2f, paired P = %.3f), %d/3 negative\n",
            S$mean[1], S$lo[1], S$hi[1], S$P[1], S$n_same_sign[1]))
