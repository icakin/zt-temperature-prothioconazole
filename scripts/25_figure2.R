# =============================================================================
# figure2_final.R — NEW Figure 2: drug sensitivity and its temperature dependence.
#   A growth dose–response (Hill, posterior epred bands)   B pEC50(T)
#   C respiration vs dose (log-linear posterior)           D dose slope(T)
#   E Bliss deviation vs temperature (by dose)             F Bliss heatmap
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork)
  library(brms)
})
source("scripts/fig_style.R")   # run from project root
# (paths now absolute from project root)
OUT <- "figures"

temps <- c(15,18,21,24,26,27,28)
temp_scale_c <- scale_colour_manual(values = TEMP_RAMP, limits = as.character(temps),
                                    name = "Temperature (\u00b0C)")
temp_scale_f <- scale_fill_manual(values = TEMP_RAMP, limits = as.character(temps),
                                  name = "Temperature (\u00b0C)")

# ---- A: Hill epred bands ----------------------------------------------------
fit <- readRDS(file.path("models/physiology/growth/fit_hill.rds"))
nd <- expand_grid(temperature = factor(temps, levels = temps),
                  conc = c(seq(0, 4, length.out = 161)))
ep <- fitted(fit, newdata = nd, re_formula = NA, probs = c(.025, .975))
band <- bind_cols(nd, as_tibble(ep)) %>%
  rename(med = Estimate, lo = `Q2.5`, hi = `Q97.5`) %>%
  mutate(temperature = as.character(temperature))
pts <- as_tibble(fit$data) %>% mutate(temperature = as.character(temperature))

pA <- ggplot() +
  geom_ribbon(data = band, aes(conc, ymin = lo, ymax = hi, fill = temperature), alpha = 0.13) +
  geom_line(data = band, aes(conc, med, colour = temperature), linewidth = 0.55) +
  geom_point(data = pts, aes(conc, rate, colour = temperature), size = 0.8, alpha = 0.7, stroke = 0) +
  temp_scale_c + temp_scale_f +
  guides(fill = "none", colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linewidth = 1.4, size = 0))) +
  labs(x = expression("Prothioconazole (mg L"^-1*")"), y = expression("Growth rate (C C"^-1*" h"^-1*")")) +
  theme_pub()

# ---- B: pEC50 by temperature ------------------------------------------------
pec <- read_csv(file.path("tables/physiology/growth/04_pec50_by_temperature.csv"), show_col_types = FALSE) %>%
  mutate(temperature = as.character(temperature))
glob <- read_csv(file.path("tables/physiology/growth/02_hill_model_parameters.csv"), show_col_types = FALSE)
gp <- if ("pEC50_median" %in% names(glob)) glob else NULL
pB <- ggplot(pec, aes(factor(temperature, levels = temps), pEC50_median, colour = temperature)) +
  geom_hline(yintercept = median(pec$pEC50_median), linetype = "22", colour = "#CFCFCF", linewidth = 0.35) +
  geom_linerange(aes(ymin = pEC50_lower, ymax = pEC50_upper), linewidth = 0.9) +
  geom_point(size = 1.9) +
  annotate("text", x = 7.35, y = max(pec$pEC50_upper), label = '"More sensitive" * symbol("\\255")',
           hjust = 1, vjust = 0, size = 2.3, colour = INK2, fontface = "italic", parse = TRUE) +
  temp_scale_c + guides(colour = "none") +
  labs(x = "Temperature (\u00b0C)", y = expression("pEC"[50]*" (-log"[10]*" M)")) +
  theme_pub()

# ---- C: respiration vs dose (log-linear model, posterior lines) --------------
rl <- read_csv(file.path("tables/physiology/respiration/03_loglinear_params_by_temperature.csv"),
               show_col_types = FALSE) %>% mutate(temperature = as.character(temperature))
dr_pts <- read_csv(file.path("tables/physiology/fig1_extra_dose_response_by_temp.csv"), show_col_types = FALSE) %>%
  filter(grepl("Respiration", Trait)) %>%
  transmute(temperature = as.character(Temperature_C), conc = Prothioconazole_mg_L,
            rate = Rate_raw)
cg <- seq(0, 4, length.out = 100)
rline <- rl %>% crossing(conc = cg) %>%
  mutate(med = exp(log(a) + slope * conc),
         lo  = exp(log(a_lower) + slope_lower * conc),
         hi  = exp(log(a_upper) + slope_upper * conc))
pC <- ggplot() +
  geom_line(data = rline, aes(conc, med, colour = temperature), linewidth = 0.55) +
  geom_point(data = dr_pts, aes(conc, rate, colour = temperature), size = 0.8, alpha = 0.7, stroke = 0) +
  scale_y_log10() + temp_scale_c +
  guides(colour = "none") +
  labs(x = expression("Prothioconazole (mg L"^-1*")"), y = expression("Respiration rate (C C"^-1*" h"^-1*")")) +
  theme_pub()

# ---- D: respiration dose slope by temperature -------------------------------
pD <- ggplot(rl, aes(factor(temperature, levels = temps), slope, colour = temperature)) +
  geom_hline(yintercept = 0, linetype = "13", colour = "#BEBEBE", linewidth = 0.35) +
  geom_linerange(aes(ymin = slope_lower, ymax = slope_upper), linewidth = 0.9) +
  geom_point(size = 1.9) +
  temp_scale_c + guides(colour = "none") +
  labs(x = "Temperature (\u00b0C)", y = expression("Respiration dose slope (per mg L"^-1*")")) +
  theme_pub()

# ---- E: Bliss deviation vs temperature --------------------------------------
bl <- read_csv(file.path("tables/physiology/growth/08_bliss_deviation_pointwise.csv"), show_col_types = FALSE) %>%
  mutate(Dose = as.character(conc))
dose_nz <- c("0.06","0.12","0.25","0.5","1","2","4")
bl_scale_c <- scale_colour_manual(values = DOSE_RAMP[dose_nz], limits = dose_nz,
                                  name = expression("Prothioconazole (mg L"^-1*")"))
bl_scale_f <- scale_fill_manual(values = DOSE_RAMP[dose_nz], limits = dose_nz,
                                name = expression("Prothioconazole (mg L"^-1*")"))
pE <- ggplot(bl, aes(temperature, delta_median, colour = Dose, fill = Dose)) +
  geom_hline(yintercept = 0, linetype = "22", colour = INK2, linewidth = 0.35) +
  geom_ribbon(aes(ymin = delta_lower, ymax = delta_upper), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.1) +
  bl_scale_c + bl_scale_f +
  guides(fill = "none", colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linewidth = 1.4))) +
  annotate("text", x = 15.1, y = max(bl$delta_upper)*0.97, label = '"Antagonism" * symbol("\\255")',
           hjust = 0, size = 2.3, colour = INK2, fontface = "italic", parse = TRUE) +
  annotate("text", x = 15.1, y = min(bl$delta_lower)*0.97, label = '"Synergy" * symbol("\\257")',
           hjust = 0, size = 2.3, colour = INK2, fontface = "italic", parse = TRUE) +
  labs(x = "Temperature (\u00b0C)", y = expression(Delta*" (Bliss deviation)")) +
  theme_pub()

# ---- F: Bliss heatmap -------------------------------------------------------
hm <- bl %>% mutate(Tf = factor(temperature, levels = temps),
                    Df = factor(Dose, levels = dose_nz))
pF <- ggplot(hm, aes(Df, Tf, fill = delta_mean)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", delta_mean),
                colour = abs(delta_mean) > 0.17), size = 1.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0, guide = "none") +
  labs(x = expression("Prothioconazole (mg L"^-1*")"), y = "Temperature (\u00b0C)") +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank())

# ---- assemble ---------------------------------------------------------------
fig <- (pA + pB) / (pC + pD) / (pE + pF) +
  plot_layout(guides = "collect", widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(1, 1, 1, 1),
        plot.tag = element_text(size = LABELPT, face = "bold", family = FONT))
ggsave(file.path(OUT, "Figure2.png"), fig, width = W2COL_MM, height = 200,
       units = "mm", dpi = 600, bg = "white")
ggsave(file.path(OUT, "Figure2.pdf"), fig, width = W2COL_MM, height = 200,
       units = "mm", bg = "white", device = cairo_pdf)
cat("Figure 2 done\n")
