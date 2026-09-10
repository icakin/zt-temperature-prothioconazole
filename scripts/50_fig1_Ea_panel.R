# 50_fig1_Ea_panel.R — Figure 1 panel G: dose-specific activation energies (E_a)
# for growth and respiration (Sharpe-Schoolfield fits). Right-side bracket: the
# respiration-growth gap equals the independently fitted CUE thermal slope.
suppressPackageStartupMessages({library(ggplot2)})
source("scripts/fig_style.R")

dr <- readRDS("tables/revision/v3build/Ea_draws.rds")
key <- readRDS("models/physiology/dose_key_tbl.rds")
doses <- key$Dose

summ <- function(df, trait) {
  do.call(rbind, lapply(seq_along(df), function(i) {
    x <- df[[i]]
    data.frame(Dose = doses[i], trait = trait, med = median(x),
               lo95 = quantile(x, .025), hi95 = quantile(x, .975),
               lo66 = quantile(x, .17),  hi66 = quantile(x, .83))
  }))
}
S <- rbind(summ(dr$gE, "Growth"), summ(dr$rE, "Respiration"))
S$Dose <- factor(S$Dose, levels = DOSES)
S$x <- as.numeric(S$Dose) + ifelse(S$trait == "Growth", -0.17, 0.17)

cueE  <- abs(median(dr$cE[[1]]))
gmean <- mean(S$med[S$trait == "Growth"]); rmean <- mean(S$med[S$trait == "Respiration"])
gap   <- rmean - gmean
bx <- 9.05  # bracket x position

p <- ggplot(S, aes(x = x, colour = Dose)) +
  geom_segment(aes(x = 0.5, xend = 8.6, y = gmean, yend = gmean),
               colour = INK2, linetype = "13", linewidth = 0.3, inherit.aes = FALSE) +
  geom_segment(aes(x = 0.5, xend = 8.6, y = rmean, yend = rmean),
               colour = INK2, linetype = "13", linewidth = 0.3, inherit.aes = FALSE) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.35) +
  geom_linerange(aes(ymin = lo66, ymax = hi66), linewidth = 0.95) +
  geom_point(aes(y = med, shape = trait, fill = Dose), size = 1.9,
             colour = INK, stroke = 0.35) +
  annotate("segment", x = bx, xend = bx, y = gmean, yend = rmean,
           colour = "#D64B21", linewidth = 0.5) +
  annotate("segment", x = bx - 0.12, xend = bx, y = gmean, yend = gmean, colour = "#D64B21", linewidth = 0.5) +
  annotate("segment", x = bx - 0.12, xend = bx, y = rmean, yend = rmean, colour = "#D64B21", linewidth = 0.5) +
  annotate("text", x = bx + 0.18, y = (gmean + rmean)/2,
           label = sprintf("Delta*italic(E)[a] == %.2f~eV", gap), parse = TRUE,
           angle = 90, size = 2.3, colour = "#D64B21", family = FONT) +
  annotate("text", x = bx + 0.62, y = (gmean + rmean)/2,
           label = sprintf("group('|',italic(E)[CUE],'|') == %.2f~eV", cueE), parse = TRUE,
           angle = 90, size = 2.3, colour = INK2, family = FONT) +
  scale_shape_manual(values = c(Growth = 21, Respiration = 24), name = NULL) +
  scale_colour_manual(values = DOSE_RAMP, limits = DOSES, guide = "none") +
  scale_fill_manual(values = DOSE_RAMP, limits = DOSES, guide = "none") +
  scale_x_continuous(breaks = 1:8, labels = DOSES, limits = c(0.4, 9.8)) +
  labs(x = expression("Prothioconazole (mg "*L^-1*")"),
       y = expression("Activation energy "*italic(E)[a]*" (eV)")) +
  theme_pub() +
  theme(legend.position = c(0.16, 0.95), legend.direction = "horizontal",
        legend.key.size = unit(2.8, "mm"))

ggsave("tables/revision/fig1G/Fig1_panelG_Ea.png", p, width = 120, height = 62,
       units = "mm", dpi = 400, bg = "white")
ggsave("tables/revision/fig1G/Fig1_panelG_Ea.pdf", p, width = 120, height = 62,
       units = "mm", device = cairo_pdf)
cat(sprintf("gap = %.3f eV vs |E_CUE| = %.3f eV\n", gap, cueE))
