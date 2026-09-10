# =============================================================================
# 56_figure6.R  -- Figure 6: SHAM (AOX inhibitor) dose response at 15 and 27 C
# Inputs : tables/aox/sham_culture_rates.csv  (temp, dose_mM, culture, r_per_h)
# Outputs: figures/Figure6.png, figures/Figure6.pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork); library(grid)})

ARGS <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(ARGS)) ARGS[1] else "."
IN   <- file.path(ROOT, "tables/aox/sham_culture_rates.csv")
OUT  <- file.path(ROOT, "figures")

C15 <- "#2166AC"; C27 <- "#B2182B"
INK <- "#1B2420"; MUT <- "#7A8288"; LINE <- "#A6AEB3"

d <- read.csv(IN, stringsAsFactors = FALSE)
d$temp <- factor(d$temp, levels = c(15, 27))
ctrl <- subset(d, dose_mM == 0)
# 0.35 mM is excluded on ESTIMABILITY grounds: at 27 C two of the three wells
# consume too little oxygen for the model to return an identifiable rate.
# The exclusion is conservative -- those wells lie in the direction of the
# reported effect. See Methods and Fig. S12.
DROP_DOSE <- 0.35
dose <- subset(d, dose_mM > 0 & dose_mM != DROP_DOSE)
dose$ld <- log10(dose$dose_mM)

## ---- culture-level slopes (the unit of replication) -------------------------
slopes <- do.call(rbind, lapply(split(dose, list(dose$temp, dose$culture), drop = TRUE),
  function(g) data.frame(temp = g$temp[1], culture = g$culture[1],
                         slope = unname(coef(lm(r_per_h ~ ld, g))[2]),
                         icpt  = unname(coef(lm(r_per_h ~ ld, g))[1]))))
s15 <- slopes$slope[slopes$temp == 15]; s27 <- slopes$slope[slopes$temp == 27]
tt  <- t.test(s27, s15, var.equal = FALSE)          # Welch: variances differ ~10x
DIFF <- mean(s27) - mean(s15); PVAL <- tt$p.value

## ---- panel A: pathway schematic --------------------------------------------
box <- function(cx, cy, w, h, lab, fill, col, size, tcol = INK, face = "bold") {
  list(annotate("rect", xmin = cx - w/2, xmax = cx + w/2, ymin = cy - h/2, ymax = cy + h/2,
                fill = fill, colour = col, linewidth = 0.38),
       annotate("text", x = cx, y = cy, label = lab, size = size, fontface = face, colour = tcol))
}
LX <- 22; BWL <- 36; RX <- 73; BWR <- 52; BY <- 60; BH <- 15; UQY <- 86
ah <- arrow(length = unit(1.5, "mm"), type = "closed")
pA <- ggplot() + xlim(0, 100) + ylim(0, 100) + theme_void() +
  box(50, UQY, 36, 12, "UQ pool", "#EAEFEC", "#93A399", 3.7) +
  annotate("segment", x = 50, xend = 50, y = UQY - 6, yend = 78, colour = LINE, linewidth = 0.32) +
  annotate("segment", x = LX, xend = RX, y = 78, yend = 78, colour = LINE, linewidth = 0.32) +
  annotate("segment", x = c(LX, RX), xend = c(LX, RX), y = 78, yend = BY + BH/2 + 1.2,
           colour = LINE, linewidth = 0.32, arrow = ah) +
  box(LX, BY, BWL, BH, "AOX", "#FBEAE7", C27, 4.05, tcol = C27) +
  box(RX, BY, BWR, BH, "complex III → IV", "#E9F0F8", "#7FA3C8", 3.4) +
  annotate("segment", x = c(LX, RX), xend = c(LX, RX), y = BY - BH/2 - 1.2, yend = 46,
           colour = LINE, linewidth = 0.32, arrow = ah) +
  annotate("text", x = c(LX, RX), y = 43, label = "O[2]%->%H[2]*O", parse = TRUE,
           size = 3.05, colour = MUT, vjust = 1) +
  annotate("text", x = LX, y = 36, label = "no proton\npumping", size = 3.05,
           colour = MUT, fontface = "italic", vjust = 1, lineheight = 0.95) +
  annotate("text", x = RX, y = 36, label = "proton pumping\n→ ATP", size = 3.05,
           colour = MUT, fontface = "italic", vjust = 1, lineheight = 0.95) +
  annotate("text", x = 4, y = 86, label = "SHAM", hjust = 0, size = 4.05,
           fontface = "bold", colour = C27) +
  annotate("segment", x = 12, xend = 12, y = 81, yend = BY + BH/2 + 2.4, colour = C27, linewidth = 0.6) +
  annotate("segment", x = 5, xend = 19, y = BY + BH/2 + 2.4, yend = BY + BH/2 + 2.4,
           colour = C27, linewidth = 1.0) +
  annotate("text", x = 25, y = 25, label = "hypothesized greater\ncontribution at 27 °C",
           size = 3.05, fontface = "bold", colour = C27, vjust = 1, lineheight = 1.0) +
  coord_cartesian(expand = FALSE, clip = "off")

## ---- shared theme for the data panels --------------------------------------
th <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.28, colour = INK),
        axis.ticks = element_line(linewidth = 0.28, colour = INK),
        axis.text = element_text(size = 7.6, colour = INK),
        axis.title = element_text(size = 8.8, colour = INK),
        plot.margin = margin(2, 2, 2, 2))
COLS <- c(`15` = C15, `27` = C27); SHP <- c(`15` = 16, `27` = 17)

## ---- panel B: control strip + dose response --------------------------------
cs <- ctrl %>% group_by(temp) %>% summarise(m = mean(r_per_h), s = sd(r_per_h), .groups = "drop")
pB1 <- ggplot() +
  geom_point(data = ctrl, aes(1, r_per_h, colour = temp), size = 1.0, alpha = 0.40) +
  geom_errorbar(data = cs, aes(1, ymin = m - s, ymax = m + s, colour = temp), width = 0.10, linewidth = 0.42) +
  geom_point(data = cs, aes(1, m, colour = temp, shape = temp), size = 1.7) +
  scale_colour_manual(values = COLS, guide = "none") + scale_shape_manual(values = SHP, guide = "none") +
  scale_x_continuous(breaks = 1, labels = "ctrl", limits = c(0.72, 1.28)) +
  scale_y_continuous(limits = c(0, 0.062), breaks = seq(0, 0.06, 0.02)) +
  labs(x = NULL, y = expression(italic(r)~"("*h^-1*")")) + th

ds <- dose %>% group_by(temp, dose_mM) %>% summarise(m = mean(r_per_h), s = sd(r_per_h), .groups = "drop")
xs <- 10^seq(log10(0.0115), log10(0.225), length.out = 120)
percult <- do.call(rbind, lapply(seq_len(nrow(slopes)), function(i)
  data.frame(temp = slopes$temp[i], x = xs, y = slopes$icpt[i] + slopes$slope[i]*log10(xs),
             g = paste0(slopes$temp[i], "_", slopes$culture[i]))))
mn <- slopes %>% group_by(temp) %>% summarise(sl = mean(slope), .groups = "drop")
mn <- merge(mn, dose %>% group_by(temp) %>% summarise(my = mean(r_per_h), mx = mean(ld), .groups = "drop"))
meanline <- do.call(rbind, lapply(seq_len(nrow(mn)), function(i)
  data.frame(temp = mn$temp[i], x = xs, y = mn$my[i] + mn$sl[i]*(log10(xs) - mn$mx[i]))))

pB2 <- ggplot() +
  geom_hline(data = cs, aes(yintercept = m, colour = temp), linetype = "22", linewidth = 0.25, alpha = 0.30) +
  geom_point(data = dose, aes(dose_mM*ifelse(temp == "15", 0.965, 1.035), r_per_h, colour = temp),
             size = 0.9, alpha = 0.32) +
  geom_line(data = percult, aes(x, y, colour = temp, group = g), linewidth = 0.25, alpha = 0.30) +
  geom_line(data = meanline, aes(x, y, colour = temp), linewidth = 0.62) +
  geom_errorbar(data = ds, aes(dose_mM, ymin = m - s, ymax = m + s, colour = temp), width = 0, linewidth = 0.42) +
  geom_point(data = ds, aes(dose_mM, m, colour = temp, shape = temp), size = 1.7) +
  annotate("text", x = 0.228, y = c(0.0535, 0.0205), label = c("15 °C", "27 °C"),
           colour = c(C15, C27), size = 3.2, fontface = "bold", hjust = 1) +
  scale_colour_manual(values = COLS, guide = "none") + scale_shape_manual(values = SHP, guide = "none") +
  scale_x_log10(breaks = sort(unique(dose$dose_mM)),
                labels = format(sort(unique(dose$dose_mM)), trim = TRUE, scientific = FALSE),
                limits = range(dose$dose_mM) * c(0.88, 1.175)) +
  scale_y_continuous(limits = c(0, 0.062)) +
  labs(x = "SHAM (mM)", y = NULL) + th +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
pB <- pB1 + pB2 + plot_layout(widths = c(0.14, 1))

## ---- panel C: culture-level slopes -----------------------------------------
sm <- slopes %>% group_by(temp) %>%
  summarise(m = mean(slope), lo = m + qt(0.025, 2)*sd(slope)/sqrt(3),
            hi = m - qt(0.025, 2)*sd(slope)/sqrt(3), .groups = "drop")
lab <- sprintf("Δ slope = −%.4f\nWelch %s = %.3f", abs(DIFF), "italic(P)", PVAL)
pC <- ggplot() +
  geom_hline(yintercept = 0, linetype = "33", colour = "#B0B4B8", linewidth = 0.26) +
  geom_point(data = slopes, aes(temp, slope, colour = temp), size = 1.2, alpha = 0.42) +
  geom_errorbar(data = sm, aes(temp, ymin = lo, ymax = hi, colour = temp), width = 0.10, linewidth = 0.45) +
  geom_point(data = sm, aes(temp, m, colour = temp, shape = temp), size = 2.0) +
  annotate("text", x = 0.42, y = -0.0405,
           label = sprintf("Delta*' slope = −%.4f'", abs(DIFF)), parse = TRUE,
           hjust = 0, vjust = 0, size = 2.85, colour = "#333333") +
  annotate("text", x = 0.42, y = -0.0455, label = sprintf("'Welch '*italic(P)*' = %.3f'", PVAL),
           parse = TRUE, hjust = 0, vjust = 0, size = 2.85, colour = "#333333") +
  scale_colour_manual(values = COLS, guide = "none") + scale_shape_manual(values = SHP, guide = "none") +
  scale_x_discrete(labels = c("15 °C", "27 °C"), expand = expansion(add = 0.7)) +
  scale_y_continuous(limits = c(-0.048, 0.012)) +
  labs(x = NULL, y = expression(atop(Delta~"slope", "("*h^-1*" per tenfold)"))) + th

## ---- assemble ---------------------------------------------------------------
right <- pB / pC + plot_layout(heights = c(1, 0.90))
fig <- pA | right
fig <- fig + plot_layout(widths = c(0.475, 0.525)) +
  plot_annotation(tag_levels = list(c("A", "B", "", "C"))) &
  theme(plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(-0.02, 1.02))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "Figure6.png"), fig, width = 6.3, height = 4.44, dpi = 600, bg = "white")
ggsave(file.path(OUT, "Figure6.pdf"), fig, width = 6.3, height = 4.44, device = cairo_pdf, bg = "white")
cat(sprintf("Figure6: 15C slope %+.5f (P=%.3f) | 27C slope %+.5f (P=%.3f) | diff %+.4f Welch P=%.3f\n",
    mean(s15), t.test(s15)$p.value, mean(s27), t.test(s27)$p.value, DIFF, PVAL))
