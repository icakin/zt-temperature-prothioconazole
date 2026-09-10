#!/usr/bin/env Rscript
# =============================================================================
# 71_figureS15.R -- Figure S15: n-propyl gallate, the second AOX inhibitor.
# Same design, same cultures, same model and the same culture-level analysis as
# Figure 6, over the dose range nPG solubility allows (0.0125-0.1 mM).
# Inputs : tables/aox/npg_culture_rates.csv   (written by 41_sham_rates.R . npg)
#          tables/aox/sham_culture_rates.csv  (for the SHAM effect size marker)
# Outputs: figures/FigureS15.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
C15 <- "#2166AC"; C27 <- "#B2182B"; INK <- "#1B2420"

slopes_of <- function(d, drop = NULL) {
  s <- d[d$dose_mM > 0 & !(d$dose_mM %in% drop), ]
  do.call(rbind, lapply(split(s, list(s$temp, s$culture), drop = TRUE), function(g)
    data.frame(temp = g$temp[1], culture = g$culture[1],
               slope = unname(coef(lm(r_per_h ~ log10(dose_mM), g))[2]),
               icpt  = unname(coef(lm(r_per_h ~ log10(dose_mM), g))[1]))))
}
d  <- read.csv(file.path(ROOT, "tables/aox/npg_culture_rates.csv"))
sh <- read.csv(file.path(ROOT, "tables/aox/sham_culture_rates.csv"))
d$temp <- factor(d$temp, levels = c(15, 27))
ctrl <- subset(d, dose_mM == 0); dose <- subset(d, dose_mM > 0)
S  <- slopes_of(d); S$temp <- factor(S$temp, levels = c(15, 27))
SH <- slopes_of(sh, drop = 0.35)
SHAM_DIFF <- mean(SH$slope[SH$temp == 27]) - mean(SH$slope[SH$temp == 15])
tt <- t.test(S$slope[S$temp == 27], S$slope[S$temp == 15], var.equal = FALSE)
DIFF <- mean(S$slope[S$temp == 27]) - mean(S$slope[S$temp == 15])

th <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.28, colour = INK),
        axis.ticks = element_line(linewidth = 0.28, colour = INK),
        axis.text = element_text(size = 7.6, colour = INK),
        axis.title = element_text(size = 8.8, colour = INK))
COLS <- c(`15` = C15, `27` = C27); SHP <- c(`15` = 16, `27` = 17)

cs <- ctrl %>% group_by(temp) %>% summarise(m = mean(r_per_h), s = sd(r_per_h), .groups = "drop")
ds <- dose %>% group_by(temp, dose_mM) %>% summarise(m = mean(r_per_h), s = sd(r_per_h), .groups = "drop")
xs <- 10^seq(log10(min(dose$dose_mM) * 0.9), log10(max(dose$dose_mM) * 1.1), length.out = 120)
pc <- do.call(rbind, lapply(seq_len(nrow(S)), function(i)
  data.frame(temp = S$temp[i], x = xs, y = S$icpt[i] + S$slope[i]*log10(xs),
             g = paste0(S$temp[i], "_", S$culture[i]))))

pA <- ggplot() +
  geom_hline(data = cs, aes(yintercept = m, colour = temp), linetype = "22", linewidth = 0.25, alpha = 0.30) +
  geom_line(data = pc, aes(x, y, colour = temp, group = g), linewidth = 0.45, alpha = 0.75) +
  geom_point(data = dose, aes(dose_mM, r_per_h, colour = temp, shape = temp,
                              group = interaction(temp, culture)), size = 1.5, alpha = 0.85) +
  annotate("text", x = min(dose$dose_mM)*0.92, y = c(0.0125, 0.0645), label = c("15 °C", "27 °C"),
           colour = c(C15, C27), size = 3.2, fontface = "bold", hjust = 0) +
  scale_colour_manual(values = COLS, guide = "none") + scale_shape_manual(values = SHP, guide = "none") +
  scale_x_log10(breaks = sort(unique(dose$dose_mM)),
                labels = c("0.0125", "0.025", "0.05", "0.075", "0.1")) +
  coord_cartesian(ylim = c(0, 0.066)) +
  labs(x = "n-propyl gallate (mM)", y = expression(italic(r)~"("*h^-1*")"),
       subtitle = "one line per culture; no aggregate trend is drawn") + th +
  theme(plot.subtitle = element_text(size = 6.6, colour = "#666666"))

sm <- S %>% group_by(temp) %>%
  summarise(m = mean(slope), lo = m + qt(0.025, 2)*sd(slope)/sqrt(3),
            hi = m - qt(0.025, 2)*sd(slope)/sqrt(3), .groups = "drop")
pB <- ggplot() +
  geom_hline(yintercept = 0, linetype = "33", colour = "#B0B4B8", linewidth = 0.26) +
  geom_hline(yintercept = SHAM_DIFF, linetype = "53", colour = C27, linewidth = 0.4) +
  annotate("text", x = 1.5, y = SHAM_DIFF, vjust = -0.55, hjust = 0.5, size = 2.5, colour = C27,
           label = "SHAM contrast") +
  geom_point(data = S, aes(temp, slope, colour = temp), size = 1.2, alpha = 0.42) +
  geom_errorbar(data = sm, aes(temp, ymin = lo, ymax = hi, colour = temp), width = 0.10, linewidth = 0.45) +
  geom_point(data = sm, aes(temp, m, colour = temp, shape = temp), size = 2.0) +
  annotate("text", x = 0.42, y = 0.062, hjust = 0, vjust = 1, size = 2.7, colour = "#333333",
           label = sprintf("Δ slope = %+.4f\n95%% CI [%.3f, %.3f]\nWelch P = %.2f",
                           DIFF, tt$conf.int[1], tt$conf.int[2], tt$p.value)) +
  annotate("text", x = 0.42, y = -0.072, hjust = 0, vjust = 0, size = 2.5, colour = "#777777",
           label = "cultures disagree in sign at\nboth temperatures (2 of 3 negative)") +
  scale_colour_manual(values = COLS, guide = "none") + scale_shape_manual(values = SHP, guide = "none") +
  scale_x_discrete(labels = c("15 °C", "27 °C"), expand = expansion(add = 0.7)) +
  coord_cartesian(ylim = c(-0.075, 0.064)) +
  labs(x = NULL, y = expression(atop(Delta~"slope", "("*h^-1*" per tenfold)"))) + th

fig <- (pA | pB) + plot_layout(widths = c(1, 0.75)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1.02))
dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(ROOT, "figures/FigureS15.png"), fig, width = 6.3, height = 2.9, dpi = 600, bg = "white")
ggsave(file.path(ROOT, "figures/FigureS15.pdf"), fig, width = 6.3, height = 2.9, device = cairo_pdf, bg = "white")
cat(sprintf("  drop most extreme culture -> 15C %+.5f, 27C %+.5f\n",
  mean(S$slope[S$temp==15][-which.max(abs(S$slope[S$temp==15]))]),
  mean(S$slope[S$temp==27][-which.max(abs(S$slope[S$temp==27]))])))
cat(sprintf("FigS15 (nPG): 15C %+.5f (P=%.3f) | 27C %+.5f (P=%.3f) | diff %+.4f [%.4f, %.4f] P=%.3f\n",
  mean(S$slope[S$temp==15]), t.test(S$slope[S$temp==15])$p.value,
  mean(S$slope[S$temp==27]), t.test(S$slope[S$temp==27])$p.value,
  DIFF, tt$conf.int[1], tt$conf.int[2], tt$p.value))
