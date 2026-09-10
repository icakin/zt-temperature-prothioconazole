# =============================================================================
# 69_figureS13.R -- Figure S13: AOX expression tracks the loss of carbon-use
# efficiency across the six transcriptomics conditions.
# Inputs : tables/revision/integration/{module_scores_by_condition_extended,
#          physiology_condition_values}.csv ; tables/revision/fig5/fig5_loo_sensitivity.csv
# Outputs: figures/FigureS13.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/fig_style.R"))

M <- read.csv(file.path(ROOT, "tables/revision/integration/module_scores_by_condition_extended.csv"),
              check.names = FALSE)
P <- read.csv(file.path(ROOT, "tables/revision/integration/physiology_condition_values.csv"))
L <- read.csv(file.path(ROOT, "tables/revision/fig5/fig5_loo_sensitivity.csv"))
D <- merge(P, M[, c("cond", "AOX")], by = "cond")
D$temp <- factor(D$T, levels = c(15, 21, 27))
D$dose <- factor(D$dose, levels = c(0, 2))
loo <- subset(L, module == "AOX" & trait == "CUE")

fit <- lm(CUE ~ AOX, D); rr <- cor(D$AOX, D$CUE)
gx <- seq(min(D$AOX), max(D$AOX), length.out = 50)
ln <- data.frame(AOX = gx, CUE = predict(fit, data.frame(AOX = gx)))

pA <- ggplot(D, aes(AOX, CUE)) +
  geom_line(data = ln, colour = "#BBBBBB", linewidth = 0.4) +
  geom_line(aes(group = temp, colour = temp), linewidth = 0.3, alpha = 0.55) +
  geom_point(aes(colour = temp, shape = dose), size = 2.2) +
  geom_text(aes(label = sprintf("%s °C, %s", T, dose)), size = 2.1,
            colour = INK2, hjust = -0.16, vjust = -0.5) +
  annotate("text", x = min(D$AOX), y = min(D$CUE),
           label = sprintf("r = %.2f (descriptive, n = 6)\nleave-one-out %.2f to %.2f",
                           loo$r_full[1], loo$r_loo_min[1], loo$r_loo_max[1]),
           hjust = 0, vjust = 0, size = 2.3, colour = INK2) +
  scale_colour_manual(values = TEMP3, name = "Temperature") +
  scale_shape_manual(values = c(`0` = 16, `2` = 17),
                     name = expression("Prothioconazole (mg L"^-1*")")) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.30))) +
  labs(x = "AOX module score (VST z)", y = "Carbon-use efficiency") +
  theme_pub() + theme(legend.position = "right", legend.key.height = unit(3.2, "mm"))

pB <- ggplot(D, aes(temp, AOX, colour = temp, shape = dose, group = dose)) +
  geom_line(colour = "#CCCCCC", linewidth = 0.3) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = TEMP3, guide = "none") +
  scale_shape_manual(values = c(`0` = 16, `2` = 17), guide = "none") +
  labs(x = "Temperature (°C)", y = "AOX module score (VST z)") + theme_pub()

fig <- (pA | pB) + plot_layout(widths = c(1, 0.52)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1.02))
dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(ROOT, "figures/FigureS13.png"), fig, width = 6.3, height = 2.9, dpi = 600, bg = "white")
ggsave(file.path(ROOT, "figures/FigureS13.pdf"), fig, width = 6.3, height = 2.9, device = cairo_pdf, bg = "white")
cat(sprintf("FigS13: AOX-CUE r = %.3f (table %.3f), LOO %.3f to %.3f\n",
            rr, loo$r_full[1], loo$r_loo_min[1], loo$r_loo_max[1]))
