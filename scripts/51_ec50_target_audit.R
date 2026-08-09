## =============================================================================
## 51_ec50_target_audit.R — Stage 1, S1.1
## EC50 target-integrity and dose-range audit, and export of JOINT posterior
## draws for the resampling analysis in script 53.
##
## Pre-registered in Stage1_predictive_spec.md revision 3. Diagnostic only.
## Touches no manuscript output and no existing model.
##
## RUN FROM the master "X0123 copy" folder:
##     Rscript scripts/51_ec50_target_audit.R
## Outputs -> tables/revision/predictive/
## =============================================================================
suppressPackageStartupMessages({
  library(brms); library(posterior); library(dplyr); library(readr); library(tidyr)
})

OUT <- "tables/revision/predictive"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "51_ec50_target_audit.log"), split = TRUE)

TEMPS <- c(15, 18, 21, 24, 26, 27, 28)
MAX_DOSE <- 4.0            # maximum tested prothioconazole concentration, mg/L

## ---- 1. joint posterior draws ----------------------------------------------
fit <- readRDS("models/physiology/growth/fit_hill.rds")
D <- as_draws_df(fit)
cat("posterior draws available:", ndraws(D), "\n\n")

## temperature enters as intercept + contrasts, so all seven temperatures share
## the intercept. Draws MUST be taken jointly; independent per-temperature
## resampling would destroy this induced correlation (spec section 4.1).
grab <- function(par) {
  int <- D[[paste0("b_", par, "_Intercept")]]
  out <- lapply(TEMPS, function(t) {
    if (t == 15) int else int + D[[paste0("b_", par, "_temperature", t)]]
  })
  names(out) <- as.character(TEMPS)
  as.data.frame(out, check.names = FALSE)
}
ec50_d <- grab("ec50")
r0_d   <- grab("r0")
n_d    <- grab("n")

draws <- bind_rows(lapply(seq_along(TEMPS), function(i) {
  data.frame(.draw = seq_len(nrow(ec50_d)),
             temperature = TEMPS[i],
             ec50 = ec50_d[[i]], r0 = r0_d[[i]], n = n_d[[i]])
}))
write_csv(draws, file.path(OUT, "hill_joint_draws.csv"))
cat("wrote joint draws:", nrow(draws), "rows\n\n")

## correlation of log10 EC50 across temperatures, to document that the draws
## are NOT independent
lg <- as.data.frame(lapply(ec50_d, function(x) log10(pmax(x, 1e-9))))
names(lg) <- as.character(TEMPS)
cat("=== correlation of log10 EC50 across temperatures (joint posterior) ===\n")
print(round(cor(lg), 3))
cat("\nmean off-diagonal correlation:",
    round(mean(cor(lg)[upper.tri(cor(lg))]), 3), "\n\n")

## ---- 2. per-temperature target audit ---------------------------------------
hill <- read_csv("tables/physiology/growth/03_hill_parameters_by_temperature.csv",
                 show_col_types = FALSE)
pec  <- read_csv("tables/physiology/growth/04_pec50_by_temperature.csv",
                 show_col_types = FALSE)

## interpolation flag: was >=50% inhibition relative to the dose-zero baseline
## actually OBSERVED at any tested dose at this temperature?
der <- read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv",
                show_col_types = FALSE) %>%
  mutate(Dose = ifelse(Dose == "Control", "0", Dose),
         conc = as.numeric(Dose), Temp = as.numeric(T),
         g = growth_C_per_C_h) %>%
  filter(is.finite(g), g > 0)

grp <- der %>% group_by(Temp, conc) %>%
  summarise(g_mean = mean(g), n_rep = n(), .groups = "drop")
ctrl <- grp %>% filter(conc == 0) %>% select(Temp, g_ctrl = g_mean)
grp <- grp %>% left_join(ctrl, by = "Temp") %>% mutate(rel = g_mean / g_ctrl)

interp <- grp %>% filter(conc > 0) %>% group_by(Temp) %>%
  summarise(min_rel_observed = min(rel),
            dose_at_min_rel  = conc[which.min(rel)],
            reached_50pc     = min(rel) <= 0.5,
            .groups = "drop")

audit <- hill %>%
  left_join(pec %>% select(temperature, pEC50_sd), by = "temperature") %>%
  left_join(interp, by = c("temperature" = "Temp")) %>%
  transmute(
    temperature,
    EC50 = EC50, EC50_lower, EC50_upper,
    log10_EC50 = log10(EC50),
    sd_log10_EC50 = pEC50_sd,
    ci_halfwidth_log10 = (log10(EC50_upper) - log10(EC50_lower)) / 2,
    hill_n = n, hill_n_lower = n_lower, hill_n_upper = n_upper,
    median_above_max_dose = EC50 > MAX_DOSE,
    upper_above_max_dose  = EC50_upper > MAX_DOSE,
    min_rel_observed = round(min_rel_observed, 3),
    dose_at_min_rel,
    interpolation_flag = reached_50pc          # TRUE = EC50 is interpolated
  )
write_csv(audit, file.path(OUT, "ec50_target_audit.csv"))

cat("=== per-temperature EC50 target audit ===\n")
print(as.data.frame(audit %>% select(temperature, EC50, EC50_lower, EC50_upper,
                                     sd_log10_EC50, median_above_max_dose,
                                     upper_above_max_dose, min_rel_observed,
                                     interpolation_flag)), row.names = FALSE)

n_interp <- sum(audit$interpolation_flag)
cat("\ntemperatures satisfying the >=50%-inhibition interpolation flag:",
    n_interp, "of 7\n")
cat("temperatures with median EC50 beyond the max tested dose (",
    MAX_DOSE, "mg/L):", sum(audit$median_above_max_dose), "\n")
cat("temperatures with upper 95% bound beyond the max tested dose:",
    sum(audit$upper_above_max_dose), "\n")

## ---- 3. target-uncertainty scale (TUS) -------------------------------------
## NOTE: this is a SCALE, not a floor on achievable prediction error. LOTO error
## is measured against posterior medians, which a model could in principle match
## arbitrarily closely regardless of posterior width. See spec section 2.
TUS <- sqrt(mean(audit$sd_log10_EC50^2))
signal_span <- max(audit$log10_EC50) - min(audit$log10_EC50)
cat("\n=== target-uncertainty scale ===\n")
cat("TUS (RMS posterior SD of log10 EC50) =", round(TUS, 4), "log10 units\n")
cat("signal span (max - min log10 EC50)   =", round(signal_span, 4), "log10 units\n")
cat("ratio TUS / signal span              =", round(TUS / signal_span, 3), "\n")
cat("[scale for judging error differences; not a floor]\n")

write_csv(data.frame(quantity = c("TUS_log10", "signal_span_log10",
                                  "ratio", "n_interpolated", "n_temperatures"),
                     value = c(TUS, signal_span, TUS / signal_span,
                               n_interp, nrow(audit))),
          file.path(OUT, "target_uncertainty_scale.csv"))
sink()
