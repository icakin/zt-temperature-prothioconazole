# =============================================================================
# 52_figure2.R — Figure 2 redesign (5 panels) + supplementary
# respiration small-multiples.
#   A growth dose-response (Hill epred bands)      [unchanged]
#   B pEC50(T)                                     [unchanged]
#   C NEW: drug effect on respiration as % change per mg/L, by temperature
#          (100*[exp(beta_T)-1], posterior medians + 95% CrI, zero line;
#          replaces old C spaghetti + old D slopes)
#   D Bliss deviation vs temperature (old E)
#   E Bliss heatmap, full width (old F)
# Supp: respiration vs dose small multiples (7 facets, shared log y,
#       points + posterior median + 95% ribbon) — model-adequacy diagnostic.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork)
})
source("scripts/fig_style.R")
OUT <- "tables/revision/fig2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
temps <- c(15,18,21,24,26,27,28)
temp_scale_c <- scale_colour_manual(values = TEMP_RAMP, limits = as.character(temps),
                                    name = "Temperature (\u00b0C)")
temp_scale_f <- scale_fill_manual(values = TEMP_RAMP, limits = as.character(temps),
                                  name = "Temperature (\u00b0C)")

# ---- A: Hill epred bands (cached CSV preferred; brms fallback) ----
if (file.exists(file.path(OUT, "fig2A_hill_band.csv"))) {
  band <- read_csv(file.path(OUT, "fig2A_hill_band.csv"), show_col_types = FALSE) %>%
    mutate(temperature = as.character(as.integer(temperature)))
  pts  <- read_csv(file.path(OUT, "fig2A_points.csv"), show_col_types = FALSE) %>%
    mutate(temperature = as.character(as.integer(temperature)))
} else {
  library(brms)
  fit <- readRDS("models/physiology/growth/fit_hill.rds")
  nd <- expand_grid(temperature = factor(temps, levels = temps),
                    conc = seq(0, 4, length.out = 161))
  ep <- fitted(fit, newdata = nd, re_formula = NA, probs = c(.025, .975))
  band <- bind_cols(nd, as_tibble(ep)) %>%
    rename(med = Estimate, lo = `Q2.5`, hi = `Q97.5`) %>%
    mutate(temperature = as.character(temperature))
  pts <- as_tibble(fit$data) %>% mutate(temperature = as.character(temperature))
}
pA <- ggplot() +
  geom_ribbon(data = band, aes(conc, ymin = lo, ymax = hi, fill = temperature), alpha = 0.13) +
  geom_line(data = band, aes(conc, med, colour = temperature), linewidth = 0.55) +
  geom_point(data = pts, aes(conc, rate, colour = temperature), size = 0.8, alpha = 0.7, stroke = 0) +
  temp_scale_c + temp_scale_f +
  guides(fill = "none", colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linewidth = 1.4, size = 0))) +
  labs(x = expression("Prothioconazole (mg L"^-1*")"),
       y = expression("Growth rate (C C"^-1*" h"^-1*")"), tag = "A") + theme_pub()

# ---- B: pEC50 (verbatim) ----
pec <- read_csv("tables/physiology/growth/04_pec50_by_temperature.csv", show_col_types = FALSE) %>%
  mutate(temperature = as.character(temperature))
pB <- ggplot(pec, aes(factor(temperature, levels = temps), pEC50_median, colour = temperature)) +
  geom_hline(yintercept = median(pec$pEC50_median), linetype = "22", colour = "#CFCFCF", linewidth = 0.35) +
  geom_linerange(aes(ymin = pEC50_lower, ymax = pEC50_upper), linewidth = 0.9) +
  geom_point(size = 1.9) +
  annotate("text", x = 7.35, y = max(pec$pEC50_upper), label = '"More sensitive" * symbol("\\255")',
           hjust = 1, vjust = 0, size = 2.3, colour = INK2, fontface = "italic", parse = TRUE) +
  temp_scale_c + guides(colour = "none") +
  labs(x = "Temperature (\u00b0C)", y = expression("pEC"[50]*" (-log"[10]*" M)"), tag = "B") +
  theme_pub()

# ---- C NEW: respiration dose effect as % per mg/L ----
fr <- readRDS("models/physiology/respiration/fit_loglinear_respiration.rds")
.sim <- attr(fr$fit, "sim"); .pn <- .sim$fnames_oi; .wu <- .sim$warmup2[1]
.gd <- function(nm){ i <- match(nm, .pn)
  unlist(lapply(.sim$samples, function(ch) ch[[i]][(.wu+1):length(ch[[i]])])) }
dw <- data.frame(b_conc = .gd("b_conc"), check.names = FALSE)
for (tt in c("18","21","24","26","27","28")) dw[[paste0("b_conc:temperature",tt)]] <- .gd(paste0("b_conc:temperature",tt))
for (tt in c("18","21","24","26","27","28")) dw[[paste0("b_temperature",tt)]] <- .gd(paste0("b_temperature",tt))
dw[["b_Intercept"]] <- .gd("b_Intercept")
beta <- tibble(`15` = dw$b_conc,
               `18` = dw$b_conc + dw$`b_conc:temperature18`,
               `21` = dw$b_conc + dw$`b_conc:temperature21`,
               `24` = dw$b_conc + dw$`b_conc:temperature24`,
               `26` = dw$b_conc + dw$`b_conc:temperature26`,
               `27` = dw$b_conc + dw$`b_conc:temperature27`,
               `28` = dw$b_conc + dw$`b_conc:temperature28`)
pct <- beta %>% mutate(across(everything(), ~ 100*(exp(.x)-1)))
Csum <- tibble(temperature = as.character(temps),
               med  = sapply(pct, median),
               lo95 = sapply(pct, quantile, .025),
               hi95 = sapply(pct, quantile, .975),
               lo66 = sapply(pct, quantile, .17),
               hi66 = sapply(pct, quantile, .83),
               p_pos = sapply(beta, function(x) mean(x > 0)))
p_max24 <- mean(apply(as.matrix(beta), 1, which.max) == 4)
cat("== new panel C stats ==\n"); print(Csum, n = 8)
cat(sprintf("P(beta_24 = max of 7) = %.3f\n", p_max24))
write_csv(Csum, file.path(OUT, "fig2C_effect_summary.csv"))
pC <- ggplot(Csum, aes(factor(temperature, levels = temps))) +
  geom_hline(yintercept = 0, linetype = "13", colour = "#BEBEBE", linewidth = 0.35) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), colour = INK2, linewidth = 0.35, alpha = 0.8) +
  geom_linerange(aes(ymin = lo66, ymax = hi66), colour = INK2, linewidth = 0.95, alpha = 0.9) +
  geom_point(aes(y = med, fill = temperature), shape = 21, size = 2.1, colour = INK, stroke = 0.35) +
  scale_fill_manual(values = TEMP_RAMP, limits = as.character(temps), guide = "none") +
  labs(x = "Temperature (\u00b0C)",
       y = expression("Respiration change (% per mg L"^-1*")"), tag = "C") +
  theme_pub()

# ---- D: Bliss deviation vs temperature (old E, verbatim) ----
bl <- read_csv("tables/physiology/growth/08_bliss_deviation_pointwise.csv", show_col_types = FALSE) %>%
  mutate(Dose = as.character(conc))
dose_nz <- c("0.06","0.12","0.25","0.5","1","2","4")
bl_scale_c <- scale_colour_manual(values = DOSE_RAMP[dose_nz], limits = dose_nz,
                                  name = expression("Prothioconazole (mg L"^-1*")"))
bl_scale_f <- scale_fill_manual(values = DOSE_RAMP[dose_nz], limits = dose_nz,
                                name = expression("Prothioconazole (mg L"^-1*")"))
pD <- ggplot(bl, aes(temperature, delta_median, colour = Dose, fill = Dose)) +
  geom_hline(yintercept = 0, linetype = "22", colour = INK2, linewidth = 0.35) +
  geom_ribbon(aes(ymin = delta_lower, ymax = delta_upper), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.1) +
  bl_scale_c + bl_scale_f +
  guides(fill = "none", colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linewidth = 1.4))) +
  annotate("text", x = 15.1, y = max(bl$delta_upper)*0.97, label = '"Antagonism" * symbol("\\255")',
           hjust = 0, size = 2.3, colour = INK2, fontface = "italic", parse = TRUE) +
  annotate("text", x = 15.1, y = min(bl$delta_lower)*0.97, label = '"Synergy" * symbol("\\257")',
           hjust = 0, size = 2.3, colour = INK2, fontface = "italic", parse = TRUE) +
  labs(x = "Temperature (\u00b0C)", y = expression(Delta*" (Bliss deviation)"), tag = "D") +
  theme_pub()

# ---- E: Bliss heatmap, full width (old F) ----
hm <- bl %>% mutate(Tf = factor(temperature, levels = temps),
                    Df = factor(Dose, levels = dose_nz))
pE <- ggplot(hm, aes(Df, Tf, fill = delta_mean)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", delta_mean),
                colour = abs(delta_mean) > 0.17), size = 2.0, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0, guide = "none") +
  coord_fixed(ratio = 1) +
  labs(x = expression("Prothioconazole (mg L"^-1*")"), y = "Temperature (\u00b0C)", tag = "E") +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank())

fig <- (pA + pB) / (pC + pD) / pE +
  plot_layout(guides = "collect", heights = c(1, 1, 0.95)) &
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(1, 1, 1, 1),
        plot.tag = element_text(size = LABELPT, face = "bold", family = FONT))
ggsave(file.path(OUT, "Figure2_redesign.png"), fig, width = 180, height = 205,
       units = "mm", dpi = 450, bg = "white")
ok <- tryCatch({ ggsave(file.path(OUT, "Figure2_redesign.pdf"), fig, width = 180, height = 205,
         units = "mm", bg = "white", device = cairo_pdf); TRUE }, error = function(e) FALSE)
if (!ok) message("PDF export failed in this R build - PNG written; PDF available in repo.")
# ---- Supplementary: respiration small multiples ----
dr_pts <- read_csv("tables/physiology/fig1_extra_dose_response_by_temp.csv", show_col_types = FALSE) %>%
  filter(grepl("Respiration", Trait)) %>%
  transmute(temperature = as.character(Temperature_C), conc = Prothioconazole_mg_L, rate = Rate_raw)
cg <- seq(0, 4, length.out = 80)
tint <- list(`15`=0, `18`=dw$b_temperature18, `21`=dw$b_temperature21, `24`=dw$b_temperature24,
             `26`=dw$b_temperature26, `27`=dw$b_temperature27, `28`=dw$b_temperature28)
bandS <- bind_rows(lapply(as.character(temps), function(tt){
  mu0 <- dw$b_Intercept + tint[[tt]]; b <- beta[[tt]]
  ln <- outer(mu0, rep(1, length(cg))) + outer(b, cg)
  q <- apply(exp(ln), 2, quantile, c(.025, .5, .975))
  tibble(temperature = tt, conc = cg, lo = q[1,], med = q[2,], hi = q[3,])}))
labfun <- function(x) paste0(x, " \u00b0C")
pS <- ggplot() +
  geom_ribbon(data = bandS, aes(conc, ymin = lo, ymax = hi, fill = temperature), alpha = 0.18) +
  geom_line(data = bandS, aes(conc, med, colour = temperature), linewidth = 0.6) +
  geom_point(data = dr_pts, aes(conc, rate, colour = temperature), size = 0.9, alpha = 0.75, stroke = 0) +
  facet_wrap(~factor(temperature, levels = temps), nrow = 2, labeller = as_labeller(labfun)) +
  scale_y_log10() + temp_scale_c + temp_scale_f + guides(colour = "none", fill = "none") +
  labs(x = expression("Prothioconazole (mg L"^-1*")"),
       y = expression("Respiration rate (C C"^-1*" h"^-1*")")) +
  theme_pub()
ggsave(file.path(OUT, "FigureS_resp_dose_multiples.png"), pS, width = 180, height = 85,
       units = "mm", dpi = 450, bg = "white")
ok <- tryCatch({ ggsave(file.path(OUT, "FigureS_resp_dose_multiples.pdf"), pS, width = 180, height = 85,
         units = "mm", bg = "white", device = cairo_pdf); TRUE }, error = function(e) FALSE)
if (!ok) message("PDF export failed in this R build - PNG written; PDF available in repo.")
cat("Figure 2 redesign done\n")
