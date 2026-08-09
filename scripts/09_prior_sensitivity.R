# =============================================================================
# 09_prior_sensitivity.R  — OPTIONAL, RUN-ONCE supplementary analysis
# -----------------------------------------------------------------------------
# Refits the GROWTH Sharpe-Schoolfield TPC under alternative priors for E (the
# activation energy) and Th (the high-temperature deactivation reference), and
# tabulates how E, lnB0 and the derived thermal optimum (Topt) shift.
#
# Purpose (for the supplementary): demonstrate that the DATA-driven parameters
# (E, lnB0) are stable to the prior, while the weakly-observed decline
# (Th / Eh / Topt) tracks its prior. This is the honest robustness check a
# reviewer will ask for, given the temperatures barely reach the optimum.
#
# NOTE: this is SLOW — it runs several MCMC fits. Run it deliberately, once:
#   source("09_prior_sensitivity.R")
# It does NOT affect the main pipeline. Output: tables/supp_prior_sensitivity.csv
# =============================================================================

.this_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                 rstudioapi::isAvailable() &&
                 nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}
source(file.path(.this_dir, "00_config.R"))
suppressPackageStartupMessages({
  library(brms); library(posterior)
})

# Growth biomass TPC data (written by 04_oxygen_fits.R)
growth_biomass_dat <- readRDS(file.path(models_dir, "growth_biomass_dat.rds"))
kB     <- 8.617e-5
kB_inv <- 11604.51812   # 1 / kB, matches 06_bayesian_models.R

# Fit the SS TPC with overridable E and Th priors (everything else as in 04).
fit_ss_prior <- function(dat, e_mean, e_sd, th_offset_c, th_sd) {
  max_TK  <- max(dat$TK, na.rm = TRUE)
  lnB0_pm <- log(stats::median(dat$y_raw, na.rm = TRUE))
  th_pm   <- max_TK + th_offset_c
  bf_ss <- brms::bf(
    y ~ lnB0 -
      (E * 11604.51812) * ((1 / TK) - (1 / 293.15)) -
      log(1 + exp((Eh * 11604.51812) * ((1 / Th) - (1 / TK)))),
    lnB0 ~ 0 + Dose_key, E ~ 0 + Dose_key, Eh ~ 0 + Dose_key, Th ~ 0 + Dose_key,
    nl = TRUE
  )
  pri <- c(
    brms::prior_string(sprintf("normal(%.6f, 2)", lnB0_pm), nlpar = "lnB0"),
    brms::prior_string(sprintf("normal(%g, %g)", e_mean, e_sd),
                       nlpar = "E", lb = 0, ub = 3),
    brms::prior_string(sprintf("normal(%g, %g)", SS_PRIOR_EH_MEAN, SS_PRIOR_EH_SD),
                       nlpar = "Eh", lb = 0, ub = 12),
    brms::prior_string(sprintf("normal(%.6f, %g)", th_pm, th_sd),
                       nlpar = "Th",
                       lb = max_TK + SS_PRIOR_TH_LB_C, ub = max_TK + SS_PRIOR_TH_UB_C),
    brms::prior_string(sprintf("exponential(%g)", SS_PRIOR_SIGMA_RATE), class = "sigma")
  )
  brms::brm(bf_ss, data = dat, family = bayes_family(), prior = pri,
            iter = BAYES_ITER, warmup = BAYES_WARMUP, chains = BAYES_CHAINS,
            seed = BAYES_SEED,
            control = list(adapt_delta = BAYES_ADAPT, max_treedepth = BAYES_MAX_TD),
            backend = "rstan", refresh = 0, init = 0)
}

# Prior settings to compare (baseline = the values in 00_config.R).
settings <- data.frame(
  label   = c("baseline", "E_wide", "E_neutral", "Th_offset_0", "Th_offset_3", "Th_wide"),
  e_mean  = c(0.65,        0.65,     0.00,        0.65,          0.65,          0.65),
  e_sd    = c(0.30,        0.60,     1.00,        0.30,          0.30,          0.30),
  th_off  = c(1.0,         1.0,      1.0,         0.0,           3.0,           1.0),
  th_sd   = c(1.5,         1.5,      1.5,         1.5,           1.5,           3.0),
  stringsAsFactors = FALSE
)

rows <- list()
for (i in seq_len(nrow(settings))) {
  st <- settings[i, ]
  message(sprintf("[%d/%d] fitting prior setting '%s' ...",
                  i, nrow(settings), st$label))
  fit <- fit_ss_prior(growth_biomass_dat, st$e_mean, st$e_sd, st$th_off, st$th_sd)
  ps  <- brms::posterior_summary(fit)
  gm  <- function(pref) mean(ps[grepl(pref, rownames(ps)), "Estimate"])
  E_m <- gm("^b_E_Dose_key"); lnB0_m <- gm("^b_lnB0_Dose_key")
  Eh_m <- gm("^b_Eh_Dose_key"); Th_m <- gm("^b_Th_Dose_key")
  Topt_m <- if (Eh_m > E_m)
    (Eh_m * Th_m) / (Eh_m + kB * Th_m * log(Eh_m / E_m - 1)) - 273.15 else NA_real_
  rows[[i]] <- data.frame(
    setting = st$label, E_prior = sprintf("N(%.2f,%.2f)", st$e_mean, st$e_sd),
    Th_prior = sprintf("N(maxT+%.0f,%.1f)", st$th_off, st$th_sd),
    E_mean = round(E_m, 3), lnB0_mean = round(lnB0_m, 3),
    Th_mean_C = round(Th_m - 273.15, 2), Topt_C = round(Topt_m, 2)
  )
}
out <- do.call(rbind, rows)
readr::write_csv(out, file.path(tables_dir, "supp_prior_sensitivity.csv"))
message("\nSaved: ", file.path(tables_dir, "supp_prior_sensitivity.csv"))
print(out)
cat("\nInterpretation: E and lnB0 stable across settings => data-driven.",
    "\n Topt / Th shifting with the Th prior => optimum is prior-informed.\n")
