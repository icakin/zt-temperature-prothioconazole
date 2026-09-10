# =============================================================================
# 68_figureS12.R -- Figure S12: SHAM temperature contrast vs fitting-interval rule
# Inputs : tables/aox/window_robustness.csv (rule, drop035, s15, s27, diff, p)
# Outputs: figures/FigureS12.png, figures/FigureS12.pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
C15 <- "#2166AC"; C27 <- "#B2182B"; INK <- "#1B2420"
PRIMARY <- -0.0238                    # hand-selected windows, Fig. 6

d <- read.csv(file.path(ROOT, "tables/aox/window_robustness.csv"), stringsAsFactors = FALSE)
d <- d[as.logical(d$drop035), ]       # 0.35 mM excluded, as in the main analysis
d$lab <- ifelse(grepl("^fixed", d$rule), paste0(d$rule, " min"), d$rule)
d$lab <- factor(d$lab, levels = rev(d$lab))

th <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.28, colour = INK),
        axis.ticks = element_line(linewidth = 0.28, colour = INK),
        axis.text = element_text(size = 7.4, colour = INK),
        axis.title = element_text(size = 8.4, colour = INK),
        plot.margin = margin(2, 3, 2, 2))

lng <- rbind(data.frame(lab = d$lab, v = d$s15, temp = "15"),
             data.frame(lab = d$lab, v = d$s27, temp = "27"))
pA <- ggplot() +
  geom_vline(xintercept = 0, linetype = "33", colour = "#B0B4B8", linewidth = 0.26) +
  geom_segment(data = d, aes(x = s15, xend = s27, y = lab, yend = lab),
               colour = "#C7CDD1", linewidth = 0.3) +
  geom_point(data = lng, aes(v, lab, colour = temp, shape = temp), size = 1.9) +
  scale_colour_manual(values = c(`15` = C15, `27` = C27), labels = c("15 °C", "27 °C"), name = NULL) +
  scale_shape_manual(values = c(`15` = 16, `27` = 17), labels = c("15 °C", "27 °C"), name = NULL) +
  scale_x_continuous(limits = c(-0.0215, 0.0122)) +
  labs(x = expression("slope of "*italic(r)*" on "*log[10]*"[SHAM]  ("*h^-1*" per tenfold)"), y = NULL) +
  th + theme(legend.position = c(0.14, 0.12), legend.text = element_text(size = 7.4),
             legend.key.height = unit(3.2, "mm"), legend.background = element_blank())

pB <- ggplot(d) +
  geom_vline(xintercept = 0, linetype = "33", colour = "#B0B4B8", linewidth = 0.26) +
  geom_vline(xintercept = PRIMARY, linetype = "53", colour = C27, linewidth = 0.38) +
  geom_point(aes(diff, lab), shape = 18, size = 2.4, colour = "#3D4A52") +
  annotate("text", x = PRIMARY, y = 0.62, label = "  primary (Fig. 6)",
           colour = C27, size = 2.6, hjust = 0, vjust = 0) +
  scale_x_continuous(limits = c(-0.030, 0.002)) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.4))) +
  labs(x = "27 °C − 15 °C difference", y = NULL) +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.line.y = element_blank())
pA <- pA + scale_y_discrete(expand = expansion(add = c(0.9, 0.4)))

fig <- (pA | pB) + plot_layout(widths = c(1, 0.85)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1.02))

dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(ROOT, "figures/FigureS12.png"), fig, width = 6.3, height = 2.7, dpi = 600, bg = "white")
ggsave(file.path(ROOT, "figures/FigureS12.pdf"), fig, width = 6.3, height = 2.7, device = cairo_pdf, bg = "white")
cat(sprintf("FigS12: %d rules | 27C steeper in %d/%d | diff range %.4f to %.4f\n",
    nrow(d), sum(d$s27 < d$s15), nrow(d), min(d$diff), max(d$diff)))
