#!/usr/bin/env Rscript
# =============================================================================
# 56_figure6.R -- Figure 6: the alternative oxidase mediates the temperature
# dependence of prothioconazole efficacy.
#
#   a  growth rate against SHAM concentration at 15 and 27 C, solvent controls
#      in the narrow strip at left. Faint points are the three cultures, thin
#      lines their within-culture regressions, thick lines the mean slope. Open
#      symbols mark rates whose oxygen trace was indistinguishable from linear,
#      so r is not identifiable there (r2_gain < 0.005; 41_sham_rates.R).
#   b  per-culture dose slopes at each temperature, mean and 95% CI (Welch t).
#   c  initial O2 consumption with cytochrome oxidase blocked by 0.5 mM azide,
#      with and without 1 mM SHAM; cell-free wells as the background band.
#   d  growth relative to vehicle at 15 C for SHAM, prothioconazole 2 mg/L and
#      the combination. Open diamond: the independence prediction from the two
#      single agents on the same plate; arrow runs prediction -> observed.
#   e  the same at 27 C.
#   f  interaction contrast I = log r_PS - log r_P - log r_S + log r_V, paired
#      within culture across the two temperatures (0 = independent action).
#
# Panels a-c come from the SHAM and azide inhibitor experiments, d-f from the
# prothioconazole x SHAM factorial. These are separate culture batches, so every
# quantity is referenced to the control or vehicle on its own plate; absolute
# rates are not comparable between the two halves.
#
# Inputs : tables/aox/sham_culture_rates.csv   (41_sham_rates.R)
#          tables/aox/azide_well_rates.csv     (43_azide_rates.R)
#          data/oxygen/azide_rep{1,2,3}_Oxygen.xlsx + azide_layout.csv
#          tables/aox/ptc_sham_condition_means.csv, ptc_sham_interaction.csv,
#          ptc_sham_interaction_diff.csv, ptc_sham_interaction_tests.csv
#                                              (47_ptc_sham_rates.R)
# Outputs: figures/Figure6.png / .pdf
# Run from the repository root:  Rscript scripts/56_figure6.R
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork); library(grid)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))

inh <- local({
  C15 <- "#2166AC"; C27 <- "#B2182B"; INK <- "#1B2420"; INK2 <- "#555555"

  ## ---- shared theme -----------------------------------------------------------
  th <- theme_classic(base_size = 9) +
    theme(axis.line = element_line(linewidth = 0.28, colour = INK),
          axis.ticks = element_line(linewidth = 0.28, colour = INK),
          axis.text = element_text(size = 7.6, colour = INK),
          axis.title = element_text(size = 8.8, colour = INK),
          legend.text = element_text(size = 7.2, colour = INK),
          legend.key.size = unit(3.5, "mm"), legend.background = element_blank(),
          plot.margin = margin(2, 2, 2, 2))

  ## =========================================================================== a, b
  d <- read.csv(file.path(ROOT, "tables/aox/sham_culture_rates.csv"), stringsAsFactors = FALSE)
  d$temp <- factor(d$temp, levels = c(15, 27))
  ctrl <- subset(d, dose_mM == 0)
  DROP_DOSE <- 0.35    # not estimable at 27 C (see Methods, Fig. S12)
  # Identifiability. r is the curvature of the oxygen decline, so a trace that is
  # essentially linear carries no growth signal whatever its R^2. r2_gain
  # (41_sham_rates.R) is how much the curved model improves on a straight line on
  # the same points; below GAIN_MIN the growth rate is not determined by the data.
  # The criterion is applied identically at both temperatures; flagged points are
  # drawn open in panel a and the slopes in b are unchanged.
  GAIN_MIN <- 0.005
  d$weak <- is.finite(d$r2_gain) & d$r2_gain < GAIN_MIN
  dose <- subset(d, dose_mM > 0 & dose_mM != DROP_DOSE); dose$ld <- log10(dose$dose_mM)

  slopes <- do.call(rbind, lapply(split(dose, list(dose$temp, dose$culture), drop = TRUE),
    function(g) data.frame(temp = g$temp[1], culture = g$culture[1],
                           slope = unname(coef(lm(r_per_h ~ ld, g))[2]),
                           icpt  = unname(coef(lm(r_per_h ~ ld, g))[1]))))
  s15 <- slopes$slope[slopes$temp == 15]; s27 <- slopes$slope[slopes$temp == 27]
  tt   <- t.test(s27, s15, var.equal = FALSE)
  DIFF <- mean(s27) - mean(s15); PVAL <- tt$p.value
  COLS <- c(`15` = C15, `27` = C27); SHP <- c(`15` = 16, `27` = 17)

  cs <- ctrl %>% group_by(temp) %>% summarise(m = mean(r_per_h), s = sd(r_per_h), .groups = "drop")
  pA1 <- ggplot() +
    geom_point(data = ctrl, aes(1, r_per_h, colour = temp), size = 1.0, alpha = 0.40) +
    geom_errorbar(data = cs, aes(1, ymin = m - s, ymax = m + s, colour = temp), width = 0.10, linewidth = 0.42) +
    geom_point(data = cs, aes(1, m, colour = temp, fill = temp, shape = temp), size = 1.7, stroke = 0.5) +
    scale_colour_manual(values = COLS, guide = "none") + scale_fill_manual(values = COLS, guide = "none") +
    scale_shape_manual(values = c(`15` = 21, `27` = 24), guide = "none") +
    scale_x_continuous(breaks = 1, labels = "ctrl", limits = c(0.72, 1.28)) +
    scale_y_continuous(limits = c(0, 0.062), breaks = seq(0, 0.06, 0.02)) +
    labs(x = NULL, y = expression(italic(r)~"("*h^-1*")")) + th

  ds <- dose %>% group_by(temp, dose_mM) %>%
    summarise(m = mean(r_per_h), s = sd(r_per_h), weak = any(weak), .groups = "drop")
  xs <- 10^seq(log10(0.0115), log10(0.225), length.out = 120)
  percult <- do.call(rbind, lapply(seq_len(nrow(slopes)), function(i)
    data.frame(temp = slopes$temp[i], x = xs, y = slopes$icpt[i] + slopes$slope[i]*log10(xs),
               g = paste0(slopes$temp[i], "_", slopes$culture[i]))))
  mn <- slopes %>% group_by(temp) %>% summarise(sl = mean(slope), .groups = "drop")
  mn <- merge(mn, dose %>% group_by(temp) %>% summarise(my = mean(r_per_h), mx = mean(ld), .groups = "drop"))
  meanline <- do.call(rbind, lapply(seq_len(nrow(mn)), function(i)
    data.frame(temp = mn$temp[i], x = xs, y = mn$my[i] + mn$sl[i]*(log10(xs) - mn$mx[i]))))

  pA2 <- ggplot() +
    geom_hline(data = cs, aes(yintercept = m, colour = temp), linetype = "22", linewidth = 0.25, alpha = 0.30) +
    geom_point(data = subset(dose, !weak), aes(dose_mM*ifelse(temp == "15", 0.965, 1.035), r_per_h, colour = temp),
               size = 0.9, alpha = 0.32) +
    geom_point(data = subset(dose, weak), aes(dose_mM*ifelse(temp == "15", 0.965, 1.035), r_per_h, colour = temp),
               size = 0.9, alpha = 0.55, shape = 1, stroke = 0.4) +
    geom_line(data = percult, aes(x, y, colour = temp, group = g), linewidth = 0.25, alpha = 0.30) +
    geom_line(data = meanline, aes(x, y, colour = temp), linewidth = 0.62) +
    geom_errorbar(data = ds, aes(dose_mM, ymin = m - s, ymax = m + s, colour = temp), width = 0, linewidth = 0.42) +
    geom_point(data = subset(ds, !weak), aes(dose_mM, m, colour = temp, fill = temp, shape = temp), size = 1.7, stroke = 0.5) +
    geom_point(data = subset(ds, weak), aes(dose_mM, m, colour = temp, shape = temp),
               size = 1.7, fill = "white", stroke = 0.7) +
    annotate("text", x = 0.228, y = c(0.0535, 0.0205), label = c("15 °C", "27 °C"),
             colour = c(C15, C27), size = 3.2, fontface = "bold", hjust = 1) +
    scale_colour_manual(values = COLS, guide = "none") + scale_fill_manual(values = COLS, guide = "none") +
    scale_shape_manual(values = c(`15` = 21, `27` = 24), guide = "none") +
    scale_x_log10(breaks = sort(unique(dose$dose_mM)), labels = c("0.0125", "", "0.05", "0.1", "0.2"),
                  limits = range(dose$dose_mM) * c(0.88, 1.175)) +
    scale_y_continuous(limits = c(0, 0.062)) +
    labs(x = "SHAM (mM)", y = NULL) + th +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

  sm <- slopes %>% group_by(temp) %>%
    summarise(m = mean(slope), lo = m + qt(0.025, 2)*sd(slope)/sqrt(3),
              hi = m - qt(0.025, 2)*sd(slope)/sqrt(3), .groups = "drop")
  pB <- ggplot() +
    geom_hline(yintercept = 0, linetype = "33", colour = "#B0B4B8", linewidth = 0.26) +
    geom_point(data = slopes, aes(temp, slope, colour = temp), size = 1.2, alpha = 0.42) +
    geom_errorbar(data = sm, aes(temp, ymin = lo, ymax = hi, colour = temp), width = 0.10, linewidth = 0.45) +
    geom_point(data = sm, aes(temp, m, colour = temp, shape = temp), size = 2.0) +
    annotate("text", x = 0.44, y = -0.0408, label = sprintf("Delta*' = −%.4f'", abs(DIFF)), parse = TRUE,
             hjust = 0, vjust = 0, size = 2.85, colour = "#333333") +
    annotate("text", x = 0.44, y = -0.0462, label = sprintf("'Welch '*italic(P)*' = %.3f'", PVAL),
             parse = TRUE, hjust = 0, vjust = 0, size = 2.85, colour = "#333333") +
    scale_colour_manual(values = COLS, guide = "none") + scale_shape_manual(values = SHP, guide = "none") +
    scale_x_discrete(labels = c("15 °C", "27 °C"), expand = expansion(add = 0.7)) +
    scale_y_continuous(limits = c(-0.050, 0.012)) +
    labs(x = NULL, y = expression(atop("Dose slope", "("*h^-1*" per tenfold)"))) + th

  ## =========================================================================== c, d
  tr <- load_azide_plates(ROOT)
  W  <- read.csv(file.path(ROOT, "tables/aox/azide_well_rates.csv"), stringsAsFactors = FALSE)
  W$condition <- factor(W$condition, levels = AZ_ORDER)
  M  <- W %>% group_by(plate, condition) %>% summarise(rate = mean(rate_mgL_h), .groups = "drop")
  Mw <- reshape(as.data.frame(M), idvar = "plate", timevar = "condition", direction = "wide")
  names(Mw) <- sub("rate\\.", "", names(Mw))
  P_AZ <- t.test(Mw$azide, Mw$`azide+SHAM`, paired = TRUE)$p.value
  TPK  <- mean(W$window_start_h)

  # (c) per-plate condition means of the raw traces (no smoothing)
  cm <- tr %>% filter(is.finite(o2)) %>%
    mutate(tb = round(t_h * 60)) %>%                       # one value per minute
    group_by(plate, condition, tb) %>% summarise(t_h = mean(t_h), o2 = mean(o2), .groups = "drop") %>%
    mutate(g = paste(plate, condition))
  pC <- ggplot() +
    annotate("rect", xmin = TPK, xmax = TPK + WIN_H, ymin = -Inf, ymax = Inf, fill = "#F3F3F3") +
    annotate("text", x = TPK + WIN_H/2, y = 12.55, label = "initial-rate window", size = 2.2, colour = INK2, vjust = 1) +
    geom_line(data = cm, aes(t_h, o2, colour = condition, group = g), linewidth = 0.45, alpha = 0.9) +
    annotate("text", x = 25.3, y = 12.6, label = "27 °C", colour = C27, size = 2.8, fontface = "bold", hjust = 1, vjust = 1) +
    scale_colour_manual(values = AZ_COL, labels = AZ_LAB, name = NULL) +
    scale_x_continuous(limits = c(0, 25.5), breaks = seq(0, 24, 6), expand = c(0, 0)) +
    scale_y_continuous(limits = c(7.6, 12.7), breaks = 8:12) +
    labs(x = "Time (h)", y = expression("Dissolved O"[2]*" (mg L"^-1*")")) + th +
    theme(legend.position = c(0.99, 0.80), legend.justification = c(1, 1),
          legend.key.height = unit(3, "mm"), legend.spacing.y = unit(0.5, "mm")) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.2)))

  # (d) initial rates
  set.seed(1)
  cellW <- W %>% filter(condition != "blank") %>%
    mutate(x = as.integer(condition) - 2 + (plate - 2) * 0.20 + runif(n(), -0.05, 0.05))
  cellM <- M %>% filter(condition != "blank") %>% mutate(x = as.integer(condition) - 2 + (plate - 2) * 0.20)
  bm <- M$rate[M$condition == "blank"]
  MK <- c(`1` = 21, `2` = 22, `3` = 24)
  pD <- ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = mean(bm) - sd(bm), ymax = mean(bm) + sd(bm), fill = "#E8E8E8") +
    geom_hline(yintercept = mean(bm), colour = "#A8A8A8", linewidth = 0.3, linetype = "32") +
    annotate("text", x = -0.55, y = mean(bm) + sd(bm) + 0.002, label = "cell-free\nbackground",
             hjust = 0, vjust = 0, size = 2.1, colour = INK2, lineheight = 0.9) +
    geom_point(data = cellW, aes(x, rate_mgL_h, fill = condition, colour = condition, shape = factor(plate)),
               size = 1.2, alpha = 0.5, stroke = 0.3) +
    geom_line(data = cellM, aes(x, rate, group = plate), colour = "#BDBDBD", linewidth = 0.35) +
    geom_point(data = cellM, aes(x, rate, fill = condition, shape = factor(plate)),
               size = 2.4, colour = INK, stroke = 0.4) +
    annotate("segment", x = 0, xend = 0, y = 0.134, yend = 0.138, colour = INK2, linewidth = 0.3) +
    annotate("segment", x = 1, xend = 1, y = 0.134, yend = 0.138, colour = INK2, linewidth = 0.3) +
    annotate("segment", x = 0, xend = 1, y = 0.138, yend = 0.138, colour = INK2, linewidth = 0.3) +
    annotate("text", x = 0.5, y = 0.140, label = sprintf("italic(P)*' = %.3f'", P_AZ), parse = TRUE,
             size = 2.5, vjust = 0, colour = INK) +
    annotate("text", x = 0.5, y = 0.150, label = "paired t, n = 3", size = 2.3, vjust = 0, colour = INK) +
    scale_colour_manual(values = AZ_COL, guide = "none") +
    scale_fill_manual(values = AZ_COL, guide = "none") +
    scale_shape_manual(values = MK, labels = paste("Replicate", 1:3), name = NULL) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Azide", "Azide\n+ SHAM"), limits = c(-0.6, 1.6)) +
    scale_y_continuous(limits = c(0, 0.158), breaks = seq(0, 0.14, 0.02)) +
    labs(x = NULL, y = expression(atop("Initial O"[2]*" consumption", "(mg L"^-1*" h"^-1*")"))) + th +
    theme(legend.position = c(1.0, 0.66), legend.justification = c(1, 0.5),
          legend.key.height = unit(3, "mm")) +
    guides(shape = guide_legend(override.aes = list(fill = "white", size = 2.0, alpha = 1)))

  list(pA1 = pA1, pA2 = pA2, pB = pB, pC = pC, pD = pD, s15 = s15, s27 = s27, DIFF = DIFF, PVAL = PVAL, P_AZ = P_AZ, azide = Mw$azide, azideSHAM = Mw$`azide+SHAM`)
})

fac <- local({
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
  XPOS <- c(S = 1, P2 = 2.4, P2S = 3.6)
  XLAB <- c(S = "SHAM", P2 = "2 mg/L", P2S = "2 mg/L + SHAM")
  rel <- M %>% filter(condition %in% c("S", "P2", "P2S")) %>% mutate(x = XPOS[condition])
  pred <- I %>% select(temp, replicate, P2S = pred_rel_P2S) %>%
    tidyr::pivot_longer(c(P2S), names_to = "condition", values_to = "pred") %>%
    mutate(x = XPOS[condition])
  obs  <- rel %>% filter(condition %in% c("P2S")) %>% select(temp, replicate, condition, obs = r_rel)
  arr  <- merge(pred, obs) %>% group_by(temp, condition, x) %>%
    summarise(pred = mean(pred), obs = mean(obs), .groups = "drop")
  ms <- rel %>% group_by(temp, condition, x) %>% summarise(m = mean(r_rel), s = sd(r_rel), .groups = "drop")
  SHP <- c(`1` = 21, `2` = 22, `3` = 24)

  panel_rel <- function(tp, tag_col, show_y = TRUE) {
    ggplot() +
      geom_hline(yintercept = 1, linetype = "22", colour = "#B0B4B8", linewidth = 0.3) +
      annotate("rect", xmin = 1.8, xmax = 4.3, ymin = -Inf, ymax = Inf, fill = "#F6F6F6") +
      annotate("text", x = 3.05, y = 1.34, label = "prothioconazole 2 mg/L",
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
      scale_x_continuous(breaks = XPOS, labels = XLAB, limits = c(0.5, 4.35), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1.42), breaks = seq(0, 1.25, 0.25)) +
      labs(x = NULL, y = if (show_y) expression("Growth rate relative to vehicle,"~italic(r)/italic(r)[V]) else NULL) + th +
      theme(legend.position = "none", axis.text.x = element_text(size = 7.0))
  }
  pA <- panel_rel(15, C15); pB <- panel_rel(27, C27, show_y = FALSE) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

  ## ---- c: interaction contrast, paired within replicate -----------------------
  Il <- I %>% select(temp, replicate, `2` = I_2) %>%
    tidyr::pivot_longer(c(`2`), names_to = "dose", values_to = "I") %>%
    mutate(dose = factor(paste0(dose, " mg/L"), levels = c("2 mg/L")),
           xt = ifelse(temp == 15, 1, 2), temp = factor(temp))
  Is <- Il %>% group_by(dose, temp, xt) %>%
    summarise(m = mean(I), lo = m + qt(0.025, 2) * sd(I) / sqrt(3), hi = m - qt(0.025, 2) * sd(I) / sqrt(3), .groups = "drop")
  lab <- data.frame(dose = factor("2 mg/L", levels = "2 mg/L"),
                    txt = sprintf("Delta == %.2f", S$mean[1]),
                    txt2 = sprintf("'paired '*italic(P)*' = %.3f'", S$P[1]))
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
    scale_y_continuous(breaks = seq(-1.5, 0.5, 0.5)) +
    coord_cartesian(ylim = c(-1.95, 0.5)) +
    labs(x = NULL, y = expression(atop("Interaction contrast "*italic(I), "(0 = independent; < 0 = SHAM potentiates)"))) +
    th + theme(legend.position = "none", panel.spacing.x = unit(3, "mm"))

  ## ---- d: respiration parameter vs growth rate ---------------------------------
  Md <- M %>% filter(condition %in% c("S", "P2", "P2S")) %>% mutate(temp = factor(temp))
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

  list(pA = pA, pB = pB, pC = pC, pD = pD, S = S)
})

## ---- assemble ---------------------------------------------------------------
row1 <- inh$pA1 + inh$pA2 + inh$pB + inh$pD + plot_layout(widths = c(0.13, 1.00, 0.52, 0.62))
row2 <- fac$pA + fac$pB + fac$pC + plot_layout(widths = c(1.00, 0.92, 1.05))
fig <- (row1 / row2) + plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = list(c("a", "", "b", "c", "d", "e", "f"))) &
  theme(plot.tag = element_text(size = 11, face = "bold"), plot.tag.position = "topleft",
        plot.margin = margin(5, 3, 2, 3))

OUT <- file.path(ROOT, "figures"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "Figure6.png"), fig, width = 10.6, height = 6.4, dpi = 500, bg = "white")
ggsave(file.path(OUT, "Figure6.pdf"), fig, width = 10.6, height = 6.4, device = cairo_pdf, bg = "white")

cat(sprintf("Figure 6\n  (b) 15 C slopes %s | 27 C slopes %s | diff %+.5f | Welch P = %.4f\n",
    paste(sprintf("%+.5f", inh$s15), collapse = " "), paste(sprintf("%+.5f", inh$s27), collapse = " "),
    inh$DIFF, inh$PVAL))
cat(sprintf("  (c) azide %s | azide+SHAM %s | paired P = %.4f\n",
    paste(sprintf("%.4f", inh$azide), collapse = " "), paste(sprintf("%.4f", inh$azideSHAM), collapse = " "), inh$P_AZ))
cat(sprintf("  (f) I27-I15 at 2 mg/L = %+.3f (95%% CI %+.2f to %+.2f, paired P = %.3f), %d/3 cultures negative\n",
    fac$S$mean[1], fac$S$lo[1], fac$S$hi[1], fac$S$P[1], fac$S$n_same_sign[1]))
