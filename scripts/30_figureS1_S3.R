# =============================================================================
# figureS1_S3.R — Supplementary Figures S1–S3 (run from project root).
#   S1  Sharpe–Schoolfield fit diagnostics: observed replicates vs posterior-
#       median TPC per dose, growth (top) and respiration (bottom).
#   S2  Posterior-predictive checks: observed vs replicated densities for the
#       Hill growth model and the log-linear respiration model, + the
#       pipeline-wide PPC coverage summary.
#   S3  Growth dose–response by temperature (Hill posterior bands, facetted).
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork)
  library(brms)
})
source("scripts/fig_style.R")
OUT <- "figures"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
KEV <- 11604.51812; TREF <- 293.15
dose_levels <- c("Control","0.06","0.12","0.25","0.5","1","2","4")
norm_dose <- function(x){x<-as.character(x); x[x %in% c("0","Control")]<-"Control"; x}

raw <- read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv", show_col_types = FALSE) %>%
  mutate(Dose = factor(norm_dose(Dose), levels = dose_levels))

ss_med_curve <- function(post_csv) {
  d <- read_csv(post_csv, show_col_types = FALSE) %>%
    mutate(Dose = factor(norm_dose(Dose), levels = dose_levels)) %>%
    group_by(Dose) %>%
    summarise(lnB0 = median(lnB0), E = median(E), Eh = median(Eh), Th = median(Th), .groups="drop") %>%
    crossing(T = seq(14.8, 28.6, length.out = 120)) %>%
    mutate(rate = exp(lnB0 - E*KEV*(1/(T+273.15) - 1/TREF) -
                        log1p(exp(Eh*KEV*(1/Th - 1/(T+273.15))))))
  d
}

# ---------------- S1: fit diagnostics ----------------------------------------
diag_panel <- function(curve, ycol, ylab) {
  ggplot() +
    geom_point(data = raw, aes(T, .data[[ycol]]), size = 0.6, alpha = 0.7, colour = INK2) +
    geom_line(data = curve, aes(T, rate), colour = "#27519E", linewidth = 0.5) +
    facet_wrap(~Dose, nrow = 2) + scale_y_log10() +
    labs(x = "Temperature (\u00b0C)", y = ylab) + theme_pub()
}
s1 <- diag_panel(ss_med_curve("tables/physiology/posterior_sharpe_schoolfield_tpc_growth_C_per_C_h_by_dose.csv"),
                 "growth_C_per_C_h", expression("Growth rate (C C"^-1*" h"^-1*")")) /
      diag_panel(ss_med_curve("tables/physiology/posterior_sharpe_schoolfield_tpc_respiration_C_per_C_h_by_dose.csv"),
                 "respiration_C_per_C_h", expression("Respiration rate (C C"^-1*" h"^-1*")")) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = LABELPT, face = "bold", family = FONT))
ggsave(file.path(OUT, "FigureS1.png"), s1, width = W2COL_MM, height = 150, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(OUT, "FigureS1.pdf"), s1, width = W2COL_MM, height = 150, units = "mm", bg = "white", device = cairo_pdf)
cat("S1 done\n")

# ---------------- S2: posterior-predictive checks ----------------------------
set.seed(7)
ppc_dens <- function(fit, obs, xlab, n_rep = 60, log10_obs = FALSE) {
  yrep <- posterior_predict(fit, ndraws = n_rep)
  dd <- as_tibble(t(yrep)) %>% mutate(.obs = obs) %>%
    pivot_longer(-.obs, names_to = "rep", values_to = "y")
  ggplot() +
    geom_density(data = dd, aes(y, group = rep), colour = "#A9C4E4", linewidth = 0.25, alpha = 0.6) +
    geom_density(aes(obs), colour = INK, linewidth = 0.8) +
    labs(x = xlab, y = "Density") + theme_pub()
}
fit_hill <- readRDS("models/physiology/growth/fit_hill.rds")
fit_resp <- readRDS("models/physiology/respiration/fit_loglinear_respiration.rds")
pS2a <- ppc_dens(fit_hill, fit_hill$data$rate, expression("Growth rate (C C"^-1*" h"^-1*")"))
pS2b <- ppc_dens(fit_resp, fit_resp$data$log_rate, "log(respiration rate)")
cov <- read_csv("tables/revision/ppc/ppc_summary.csv", show_col_types = FALSE)
pS2c <- ggplot(cov, aes(model, coverage95)) +
  geom_hline(yintercept = 0.95, linetype = "22", colour = "#C8C8C8", linewidth = 0.4) +
  geom_col(fill = "#457BBE", width = 0.5) +
  geom_text(aes(label = sprintf("%.2f", coverage95)), vjust = -0.4, size = 2.4, colour = INK) +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  labs(x = NULL, y = "Empirical 95% coverage") + theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
s2 <- (pS2a + pS2b + pS2c) + plot_annotation(tag_levels = "A",
        caption = "Light lines: 60 posterior-predictive replicates; dark line: observed. C: coverage of the 95% predictive interval per model.") &
  theme(plot.tag = element_text(size = LABELPT, face = "bold", family = FONT),
        plot.caption = element_text(size = BASE - 1, colour = INK2))
ggsave(file.path(OUT, "FigureS2.png"), s2, width = W2COL_MM, height = 70, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(OUT, "FigureS2.pdf"), s2, width = W2COL_MM, height = 70, units = "mm", bg = "white", device = cairo_pdf)
cat("S2 done\n")

# ---------------- S3: dose–response by temperature ---------------------------
temps <- c(15,18,21,24,26,27,28)
nd <- expand_grid(temperature = factor(temps), conc = seq(0, 4, length.out = 121))
ep <- fitted(fit_hill, newdata = nd, re_formula = NA, probs = c(.025, .975))
band <- bind_cols(nd, as_tibble(ep)) %>%
  rename(med = Estimate, lo = `Q2.5`, hi = `Q97.5`) %>%
  mutate(temperature = as.character(temperature))
pts <- as_tibble(fit_hill$data) %>% mutate(temperature = as.character(temperature))
s3 <- ggplot() +
  geom_ribbon(data = band, aes(conc, ymin = lo, ymax = hi, fill = temperature), alpha = 0.20) +
  geom_line(data = band, aes(conc, med, colour = temperature), linewidth = 0.55) +
  geom_point(data = pts, aes(conc, rate, colour = temperature), size = 0.7, alpha = 0.8, stroke = 0) +
  facet_wrap(~factor(temperature, levels = temps), nrow = 2) +
  scale_colour_manual(values = TEMP_RAMP, limits = as.character(temps), guide = "none") +
  scale_fill_manual(values = TEMP_RAMP, limits = as.character(temps), guide = "none") +
  labs(x = expression("Prothioconazole (mg L"^-1*")"),
       y = expression("Growth rate (C C"^-1*" h"^-1*")")) + theme_pub()
ggsave(file.path(OUT, "FigureS3.png"), s3, width = W2COL_MM, height = 100, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(OUT, "FigureS3.pdf"), s3, width = W2COL_MM, height = 100, units = "mm", bg = "white", device = cairo_pdf)
cat("S3 done\n")
