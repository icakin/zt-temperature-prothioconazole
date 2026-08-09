# =============================================================================
# figure1_final.R — NEW Figure 1: whole-organism thermal performance.
#   A growth TPC by dose (Sharpe–Schoolfield posterior)   B lnB0 ridges
#   C respiration TPC by dose                             D lnB0 ridges
#   E CUE Boltzmann–Arrhenius (common slope)              F alpha ridges
# Curves recomputed from the SAVED posterior draws with the pipeline's own
# equations (Tref = 293.15 K, 1 eV/k = 11604.51812 K, kB = 8.617e-5 eV/K).
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(ggridges); library(patchwork)
})
source("scripts/fig_style.R")   # run from project root

DAT  <- "tables/physiology"
UPL  <- "tables/physiology"
OUT  <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

KEV   <- 11604.51812; TREF <- 293.15; KB <- 8.617e-5
dose_levels <- c("Control","0.06","0.12","0.25","0.5","1","2","4")
ridge_levels <- rev(dose_levels)           # Control on top

norm_dose <- function(x) {
  x <- as.character(x)
  x[x %in% c("0","Control","control")] <- "Control"
  sub("^0\\.50$","0.5", sub("^1\\.00$","1", x))
}

# ---- data -------------------------------------------------------------------
raw <- read_csv(file.path(DAT, "derived_N0_R_results_with_carbon.csv"), show_col_types = FALSE) %>%
  mutate(Dose = factor(norm_dose(Dose), levels = dose_levels))
dg <- read_csv(file.path(DAT, "posterior_sharpe_schoolfield_tpc_growth_C_per_C_h_by_dose.csv"),
               show_col_types = FALSE) %>% mutate(Dose = factor(norm_dose(Dose), levels = dose_levels))
dr <- read_csv(file.path(UPL, "posterior_sharpe_schoolfield_tpc_respiration_C_per_C_h_by_dose.csv"),
               show_col_types = FALSE) %>% mutate(Dose = factor(norm_dose(Dose), levels = dose_levels))
dc <- read_csv(file.path(UPL, "posterior_arrhenius_log_CUE_common_slope.csv"),
               show_col_types = FALSE) %>% mutate(Dose = factor(norm_dose(Dose), levels = dose_levels))

# ---- SS curve quantiles from draws ------------------------------------------
Tg <- seq(14.8, 28.6, length.out = 140)
ss_band <- function(d) {
  bind_rows(lapply(split(d, d$Dose), function(dd) {
    TK <- Tg + 273.15
    ln <- outer(seq_len(nrow(dd)), seq_along(TK), function(i, j)
      dd$lnB0[i] - dd$E[i] * KEV * (1/TK[j] - 1/TREF) -
        log1p(exp(dd$Eh[i] * KEV * (1/dd$Th[i] - 1/TK[j]))))
    q <- apply(exp(ln), 2, quantile, c(.025, .5, .975), na.rm = TRUE)
    tibble(Dose = dd$Dose[1], T = Tg, lo = q[1,], med = q[2,], hi = q[3,])
  }))
}
set.seed(1)
sub_d <- function(d, n = 2000) d %>% group_by(Dose) %>% slice_sample(n = n) %>% ungroup()
bg <- ss_band(sub_d(dg)); br <- ss_band(sub_d(dr))

dose_scale_c <- scale_colour_manual(values = DOSE_RAMP, limits = dose_levels,
                                    name = expression("Prothioconazole (mg L"^-1*")"))
dose_scale_f <- scale_fill_manual(values = DOSE_RAMP, limits = dose_levels,
                                  name = expression("Prothioconazole (mg L"^-1*")"))

tpc_panel <- function(band, pts, ycol, ylab) {
  ggplot() +
    geom_ribbon(data = band, aes(T, ymin = lo, ymax = hi, fill = Dose), alpha = 0.11) +
    geom_line(data = band, aes(T, med, colour = Dose), linewidth = 0.55) +
    geom_point(data = pts, aes(T, .data[[ycol]], colour = Dose),
               size = 0.75, alpha = 0.75, stroke = 0) +
    scale_y_log10() + dose_scale_c + dose_scale_f +
    guides(fill = "none", colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linewidth = 1.4, size = 0))) +
    labs(x = "Temperature (\u00b0C)", y = ylab) + theme_pub()
}
pA <- tpc_panel(bg, raw, "growth_C_per_C_h",
                expression("Growth rate (C C"^-1*" h"^-1*")"))
pC <- tpc_panel(br, raw, "respiration_C_per_C_h",
                expression("Respiration rate (C C"^-1*" h"^-1*")"))

# ---- lnB0 / alpha ridge panels ----------------------------------------------
ridge_panel <- function(d, vcol, xlab) {
  ctrl_med <- median(d[[vcol]][d$Dose == "Control"])
  ggplot(d, aes(.data[[vcol]], factor(Dose, levels = ridge_levels),
                fill = Dose, colour = Dose)) +
    geom_density_ridges(scale = 1.25, alpha = 0.55, linewidth = 0.35,
                        quantile_lines = TRUE, quantiles = 2,
                        rel_min_height = 0.005) +
    geom_vline(xintercept = ctrl_med, linetype = "22", colour = INK2, linewidth = 0.35) +
    dose_scale_c + dose_scale_f + guides(colour = "none", fill = "none") +
    labs(x = xlab, y = NULL) + theme_pub() + theme(legend.position = "none")
}
pB <- ridge_panel(dg, "lnB0", expression("ln(B"[0]*"), growth"))
pD <- ridge_panel(dr, "lnB0", expression("ln(B"[0]*"), respiration"))

# ---- CUE Boltzmann (common slope) -------------------------------------------
xg <- 1/(KB*TREF) - 1/(KB*(Tg + 273.15))
cue_band <- bind_rows(lapply(split(dc, dc$Dose), function(dd) {
  ln <- outer(dd$alpha, xg, function(a, x) a) + outer(dd$E, xg, function(e, x) e * x)
  q <- apply(exp(ln), 2, quantile, c(.025, .5, .975), na.rm = TRUE)
  tibble(Dose = dd$Dose[1], x = xg, lo = q[1,], med = q[2,], hi = q[3,])
}))
raw_cue <- raw %>% mutate(x = 1/(KB*TREF) - 1/(KB*(T + 273.15)))
pE <- ggplot() +
  geom_ribbon(data = cue_band, aes(x, ymin = lo, ymax = hi, fill = Dose), alpha = 0.11) +
  geom_line(data = cue_band, aes(x, med, colour = Dose), linewidth = 0.55) +
  geom_point(data = raw_cue, aes(x, CUE, colour = Dose), size = 0.75, alpha = 0.75, stroke = 0) +
  scale_y_log10() + dose_scale_c + dose_scale_f +
  guides(fill = "none", colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linewidth = 1.4, size = 0))) +
  labs(x = expression(frac(1, k*T[ref]) - frac(1, k*T)), y = "CUE") + theme_pub()
pF <- ridge_panel(dc, "alpha", expression(alpha*"  (equiv. ln(CUE) at T"[ref]*")"))

# ---- assemble ---------------------------------------------------------------
fig <- (pA + pB) / (pC + pD) / (pE + pF) +
  plot_layout(guides = "collect", widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.title = element_text(size = BASE),
        plot.tag = element_text(size = LABELPT, face = "bold", family = FONT))
ggsave(file.path(OUT, "Figure1.png"), fig, width = W2COL_MM, height = 205,
       units = "mm", dpi = 600, bg = "white")
ggsave(file.path(OUT, "Figure1.pdf"), fig, width = W2COL_MM, height = 205,
       units = "mm", bg = "white", device = cairo_pdf)
cat("Figure 1 done\n")
