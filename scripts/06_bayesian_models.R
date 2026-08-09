# ==============================================================================
# 06_bayesian_models.R
# Stage 2: Bayesian Thermal Model Fitting
#
# This script:
#  - Fits Bayesian Arrhenius models (separate and hierarchical)
#  - Fits Bayesian Sharpe-Schoolfield thermal performance curves (TPCs)
#  - Extracts posterior summaries and prediction grids
#  - Saves all models and posterior data for downstream plotting
#
# NO plotting output—only model fitting and data extraction.
# ==============================================================================

# Source shared config (works interactively in RStudio, via source(), or from CLI)
.this_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                 rstudioapi::isAvailable() &&
                 nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}
source(file.path(.this_dir, "00_config.R"))

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(grid)
})

options(mc.cores = max(1, parallel::detectCores() - 1))

# Lift R's vector-memory cap on macOS (large models + workspace save).
if (Sys.info()[["sysname"]] == "Darwin") {
  try(suppressWarnings(mem.maxVSize(vsize = Inf)), silent = TRUE)
}

# ===== Reload Stage 1 outputs ===================================================
results_plot     <- readRDS(file.path(models_dir, "results_plot.rds"))
ratio_dat        <- readRDS(file.path(models_dir, "ratio_dat.rds"))
dose_levels      <- readRDS(file.path(models_dir, "dose_levels.rds"))
dose_cols        <- readRDS(file.path(models_dir, "dose_cols.rds"))
dose_key_tbl     <- readRDS(file.path(models_dir, "dose_key_tbl.rds"))
growth_dat       <- readRDS(file.path(models_dir, "growth_dat.rds"))
resp_dat         <- readRDS(file.path(models_dir, "resp_dat.rds"))
ratio_dat_arr    <- readRDS(file.path(models_dir, "ratio_dat_arr.rds"))
cue_dat          <- readRDS(file.path(models_dir, "cue_dat.rds"))
growth_biomass_dat <- readRDS(file.path(models_dir, "growth_biomass_dat.rds"))
resp_biomass_dat <- readRDS(file.path(models_dir, "resp_biomass_dat.rds"))

# ===== Helper function: posterior summary =====================================
posterior_summary_df <- function(fit, pars = NULL, probs = c(0.025, 0.5, 0.975)) {
  as.data.frame(
    posterior::summarise_draws(
      posterior::as_draws_df(fit, variable = pars),
      mean,
      sd,
      ~posterior::quantile2(.x, probs = probs)
    )
  )
}

# ===== Helper function: TPC check plot ========================================
make_tpc_check_plot <- function(dat, grid, title_txt, ylab_txt, log_scale = FALSE) {
  p <- ggplot(dat, aes(T, y_raw, color = Dose, fill = Dose)) +
    geom_point(size = 2.0, alpha = 0.9) +
    {
      if (nrow(grid) > 0) {
        geom_ribbon(
          data = grid,
          aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose),
          inherit.aes = FALSE,
          alpha = 0.18,
          colour = NA
        )
      }
    } +
    {
      if (nrow(grid) > 0) {
        geom_line(
          data = grid,
          aes(x = T, y = raw_q50, color = Dose),
          inherit.aes = FALSE,
          linewidth = 1
        )
      }
    } +
    facet_wrap(~Dose, scales = "free_y") +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(
      title = title_txt,
      x = "Temperature (°C)",
      y = ylab_txt,
      color = "Dose",
      fill = "Dose"
    ) +
    theme_classic(12) +
    theme(legend.position = "none")

  if (isTRUE(log_scale)) p <- p + scale_y_log10()
  p
}

# ===== Descriptive CUE plot on Boltzmann x-axis ===============================
p_cue_boltz_scatter <- ggplot(cue_dat, aes(boltz_shift, y_raw, color = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_y_log10() +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(
    title = "CUE vs Boltzmann temperature scale",
    x = expression(frac(1, k[T[ref]]) - frac(1, k*T)),
    y = "CUE (log scale)",
    color = "Dose"
  ) +
  theme_classic(12)

# ===== Bayesian Arrhenius model ==============================================
fit_bayes_arr_by_dose <- function(dat) {
  if (nrow(dat) < 8) stop("Too few rows for Bayesian Arrhenius fit.")
  if (dplyr::n_distinct(dat$Dose_key) < 1) stop("No dose levels found.")
  if (dplyr::n_distinct(dat$TK) < 3) stop("Too few distinct temperatures for Arrhenius fit.")

  bf_arr <- brms::bf(
    y ~ alpha + E * boltz_shift,
    alpha ~ 0 + Dose_key,
    E     ~ 0 + Dose_key,
    nl = TRUE
  )

  pri <- c(
    brms::prior("normal(0, 5)", nlpar = "alpha"),
    brms::prior("normal(0.65, 0.35)", nlpar = "E", lb = 0, ub = 3),
    brms::prior("exponential(1)", class = "sigma")
  )

  brms::brm(
    formula = bf_arr,
    data = dat,
    family = bayes_family(),   # gaussian or robust student-t (see BAYES_LIKELIHOOD)
    prior = pri,
    iter = BAYES_ITER,
    warmup = BAYES_WARMUP,
    chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = BAYES_ADAPT, max_treedepth = BAYES_MAX_TD),
    backend = "rstan",
    refresh = 0,
    init = 0
  )
}

fit_bayes_arr_cue_common_slope <- function(dat) {
  if (nrow(dat) < 8) stop("Too few rows for Bayesian CUE Arrhenius fit.")
  if (dplyr::n_distinct(dat$Dose_key) < 1) stop("No dose levels found.")
  if (dplyr::n_distinct(dat$TK) < 3) stop("Too few distinct temperatures for CUE Arrhenius fit.")

  bf_arr_cue <- brms::bf(
    y ~ alpha + E * boltz_shift,
    alpha ~ 0 + Dose_key,
    E ~ 1,
    nl = TRUE
  )

  pri <- c(
    brms::prior("normal(0, 5)", nlpar = "alpha"),
    brms::prior("normal(0, 1)", nlpar = "E"),
    brms::prior("exponential(1)", class = "sigma")
  )

  brms::brm(
    formula = bf_arr_cue,
    data = dat,
    family = bayes_family(),   # gaussian or robust student-t (see BAYES_LIKELIHOOD)
    prior = pri,
    iter = BAYES_ITER,
    warmup = BAYES_WARMUP,
    chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = BAYES_ADAPT, max_treedepth = BAYES_MAX_TD),
    backend = "rstan",
    refresh = 0,
    init = 0
  )
}

# ===== Separate hierarchical Arrhenius sensitivity analysis ===================
fit_bayes_arr_by_dose_hier <- function(dat) {
  if (nrow(dat) < 8) stop("Too few rows for hierarchical Bayesian Arrhenius fit.")
  if (dplyr::n_distinct(dat$Dose_key) < 1) stop("No dose levels found.")
  if (dplyr::n_distinct(dat$TK) < 3) stop("Too few distinct temperatures for Arrhenius fit.")
  if (!"cond_id" %in% names(dat)) stop("cond_id column not found.")

  bf_arr_hier <- brms::bf(
    y ~ alpha + E * boltz_shift,
    alpha ~ 0 + Dose_key + (1 | cond_id),
    E     ~ 0 + Dose_key,
    nl = TRUE
  )

  pri <- c(
    brms::prior("normal(0, 5)", nlpar = "alpha"),
    brms::prior("normal(0.65, 0.35)", nlpar = "E", lb = 0, ub = 3),
    brms::prior("student_t(3, 0, 2.5)", class = "sd", nlpar = "alpha", group = "cond_id"),
    brms::prior("exponential(1)", class = "sigma")
  )

  brms::brm(
    formula = bf_arr_hier,
    data = dat,
    family = bayes_family(),   # gaussian or robust student-t (see BAYES_LIKELIHOOD)
    prior = pri,
    iter = BAYES_ITER,
    warmup = BAYES_WARMUP,
    chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = 0.99, max_treedepth = BAYES_MAX_TD),
    backend = "rstan",
    refresh = 0,
    init = 0
  )
}

# ===== Bayesian full TPC: Sharpe–Schoolfield =================================
fit_bayes_ss_tpc_by_dose <- function(dat) {
  if (nrow(dat) < 8) stop("Too few rows for Bayesian Sharpe–Schoolfield TPC fit.")
  if (dplyr::n_distinct(dat$Dose_key) < 1) stop("No dose levels found.")
  if (dplyr::n_distinct(dat$T) < 4) stop("Too few distinct temperatures for TPC fit.")

  max_TK_obs <- max(dat$TK, na.rm = TRUE)

  lnB0_prior_mean <- log(stats::median(dat$y_raw, na.rm = TRUE))
  Th_prior_mean   <- max_TK_obs + SS_PRIOR_TH_OFFSET_C   # offset in K (== deg C step)
  Th_lb           <- max_TK_obs + SS_PRIOR_TH_LB_C
  Th_ub           <- max_TK_obs + SS_PRIOR_TH_UB_C

  prior_lnB0 <- sprintf("normal(%0.8f, 2)", lnB0_prior_mean)
  prior_Th   <- sprintf("normal(%0.8f, %0.4f)", Th_prior_mean, SS_PRIOR_TH_SD)

  bf_ss <- brms::bf(
    y ~ lnB0 -
      (E * 11604.51812) * ((1 / TK) - (1 / 293.15)) -
      log(1 + exp((Eh * 11604.51812) * ((1 / Th) - (1 / TK)))),
    lnB0 ~ 0 + Dose_key,
    E    ~ 0 + Dose_key,
    Eh   ~ 0 + Dose_key,
    Th   ~ 0 + Dose_key,
    nl = TRUE
  )

  prior_E     <- sprintf("normal(%g, %g)", SS_PRIOR_E_MEAN,  SS_PRIOR_E_SD)
  prior_Eh    <- sprintf("normal(%g, %g)", SS_PRIOR_EH_MEAN, SS_PRIOR_EH_SD)
  prior_sigma <- sprintf("exponential(%g)", SS_PRIOR_SIGMA_RATE)

  pri <- c(
    do.call(brms::prior, list(prior = prior_lnB0, nlpar = "lnB0")),
    do.call(brms::prior, list(prior = prior_E,  nlpar = "E",  lb = SS_PRIOR_E_LB,  ub = SS_PRIOR_E_UB)),
    do.call(brms::prior, list(prior = prior_Eh, nlpar = "Eh", lb = SS_PRIOR_EH_LB, ub = SS_PRIOR_EH_UB)),
    do.call(brms::prior, list(prior = prior_Th, nlpar = "Th", lb = Th_lb, ub = Th_ub)),
    do.call(brms::prior, list(prior = prior_sigma, class = "sigma"))
  )

  brms::brm(
    formula = bf_ss,
    data = dat,
    family = bayes_family(),   # gaussian or robust student-t (see BAYES_LIKELIHOOD)
    prior = pri,
    iter = BAYES_ITER,
    warmup = BAYES_WARMUP,
    chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = 0.99, max_treedepth = BAYES_MAX_TD),
    backend = "rstan",
    refresh = 0
  )
}

# ===== Fit Bayesian thermal models ===========================================
fit_growth_ok <- nrow(growth_dat) >= 8 && dplyr::n_distinct(growth_dat$TK) >= 3 && dplyr::n_distinct(growth_dat$Dose_key) >= 1
fit_resp_ok   <- nrow(resp_dat)   >= 8 && dplyr::n_distinct(resp_dat$TK)   >= 3 && dplyr::n_distinct(resp_dat$Dose_key)   >= 1
fit_ratio_ok  <- nrow(ratio_dat_arr) >= 8 && dplyr::n_distinct(ratio_dat_arr$TK) >= 3 && dplyr::n_distinct(ratio_dat_arr$Dose_key) >= 1
fit_cue_ok    <- nrow(cue_dat) >= 8 && dplyr::n_distinct(cue_dat$TK) >= 3 && dplyr::n_distinct(cue_dat$Dose_key) >= 1

fit_growth_biomass_ok <- nrow(growth_biomass_dat) >= 8 && dplyr::n_distinct(growth_biomass_dat$TK) >= 3 && dplyr::n_distinct(growth_biomass_dat$Dose_key) >= 1
fit_resp_biomass_ok   <- nrow(resp_biomass_dat)   >= 8 && dplyr::n_distinct(resp_biomass_dat$TK)   >= 3 && dplyr::n_distinct(resp_biomass_dat$Dose_key)   >= 1

# ESSENTIAL models (needed by the kept figures) — always fit:
#   SS TPC biomass growth/respiration, and the common-slope CUE model.
fit_cue_brm    <- if (fit_cue_ok)    fit_bayes_arr_cue_common_slope(cue_dat) else NULL
fit_growth_tpc_biomass_brm <- if (fit_growth_biomass_ok) fit_bayes_ss_tpc_by_dose(growth_biomass_dat) else NULL
fit_resp_tpc_biomass_brm   <- if (fit_resp_biomass_ok)   fit_bayes_ss_tpc_by_dose(resp_biomass_dat) else NULL

# NON-ESSENTIAL models — only fit when FIT_ALL_MODELS is TRUE (see 00_config.R).
# These fed figures that were removed; all downstream use is NULL-guarded, so
# skipping them just no-ops those plots and saves the bulk of the MCMC time.
fit_growth_brm <- if (FIT_ALL_MODELS && fit_growth_ok) fit_bayes_arr_by_dose(growth_dat) else NULL
fit_resp_brm   <- if (FIT_ALL_MODELS && fit_resp_ok)   fit_bayes_arr_by_dose(resp_dat) else NULL
fit_ratio_brm  <- if (FIT_ALL_MODELS && fit_ratio_ok)  fit_bayes_arr_by_dose(ratio_dat_arr) else NULL

fit_growth_biomass_brm <- if (FIT_ALL_MODELS && fit_growth_biomass_ok) fit_bayes_arr_by_dose(growth_biomass_dat) else NULL
fit_resp_biomass_brm   <- if (FIT_ALL_MODELS && fit_resp_biomass_ok)   fit_bayes_arr_by_dose(resp_biomass_dat) else NULL

fit_growth_tpc_brm <- if (FIT_ALL_MODELS && fit_growth_ok) fit_bayes_ss_tpc_by_dose(growth_dat) else NULL
fit_resp_tpc_brm   <- if (FIT_ALL_MODELS && fit_resp_ok)   fit_bayes_ss_tpc_by_dose(resp_dat) else NULL

fit_growth_hier_brm <- if (FIT_ALL_MODELS && fit_growth_ok) fit_bayes_arr_by_dose_hier(growth_dat) else NULL
fit_resp_hier_brm   <- if (FIT_ALL_MODELS && fit_resp_ok)   fit_bayes_arr_by_dose_hier(resp_dat) else NULL
fit_ratio_hier_brm  <- if (FIT_ALL_MODELS && fit_ratio_ok)  fit_bayes_arr_by_dose_hier(ratio_dat_arr) else NULL

fit_growth_biomass_hier_brm <- if (FIT_ALL_MODELS && fit_growth_biomass_ok) fit_bayes_arr_by_dose_hier(growth_biomass_dat) else NULL
fit_resp_biomass_hier_brm   <- if (FIT_ALL_MODELS && fit_resp_biomass_ok)   fit_bayes_arr_by_dose_hier(resp_biomass_dat) else NULL

# ---- SUPPLEMENTARY: MCMC convergence diagnostics (Rhat / ESS / divergences) --
# Written for the models the paper uses. A reviewer will want to see these.
.diag_one <- function(fit, name) {
  if (is.null(fit)) return(NULL)
  tryCatch({
    dr <- posterior::as_draws(fit)
    s  <- posterior::summarise_draws(dr, "rhat", "ess_bulk", "ess_tail")
    np <- brms::nuts_params(fit)
    ndiv <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
    data.frame(model = name, variable = s$variable, rhat = s$rhat,
               ess_bulk = s$ess_bulk, ess_tail = s$ess_tail, n_divergent = ndiv)
  }, error = function(e) NULL)
}
.conv <- dplyr::bind_rows(
  .diag_one(fit_growth_tpc_biomass_brm, "SS_TPC_growth"),
  .diag_one(fit_resp_tpc_biomass_brm,   "SS_TPC_respiration"),
  .diag_one(fit_cue_brm,                "CUE_common_slope")
)
if (nrow(.conv) > 0) {
  readr::write_csv(.conv, file.path(tables_dir, "supp_bayes_convergence.csv"))
  message("Saved MCMC convergence diagnostics: ",
          file.path(tables_dir, "supp_bayes_convergence.csv"),
          "  (max Rhat = ", sprintf("%.3f", max(.conv$rhat, na.rm = TRUE)),
          ", total divergences = ", sum(unique(.conv[, c("model","n_divergent")])$n_divergent), ")")
}

if (!is.null(fit_growth_brm)) saveRDS(fit_growth_brm, brms_growth_rds)
if (!is.null(fit_resp_brm)) saveRDS(fit_resp_brm, brms_resp_rds)
if (!is.null(fit_ratio_brm)) saveRDS(fit_ratio_brm, brms_ratio_rds)
if (!is.null(fit_cue_brm)) saveRDS(fit_cue_brm, brms_cue_rds)
if (!is.null(fit_growth_biomass_brm)) saveRDS(fit_growth_biomass_brm, brms_growth_biomass_rds)
if (!is.null(fit_resp_biomass_brm)) saveRDS(fit_resp_biomass_brm, brms_resp_biomass_rds)

if (!is.null(fit_growth_tpc_brm)) saveRDS(fit_growth_tpc_brm, brms_tpc_growth_rds)
if (!is.null(fit_resp_tpc_brm)) saveRDS(fit_resp_tpc_brm, brms_tpc_resp_rds)
if (!is.null(fit_growth_tpc_biomass_brm)) saveRDS(fit_growth_tpc_biomass_brm, brms_tpc_growth_biomass_rds)
if (!is.null(fit_resp_tpc_biomass_brm)) saveRDS(fit_resp_tpc_biomass_brm, brms_tpc_resp_biomass_rds)

if (!is.null(fit_growth_hier_brm)) saveRDS(fit_growth_hier_brm, brms_growth_hier_rds)
if (!is.null(fit_resp_hier_brm)) saveRDS(fit_resp_hier_brm, brms_resp_hier_rds)
if (!is.null(fit_ratio_hier_brm)) saveRDS(fit_ratio_hier_brm, brms_ratio_hier_rds)
if (!is.null(fit_growth_biomass_hier_brm)) saveRDS(fit_growth_biomass_hier_brm, brms_growth_biomass_hier_rds)
if (!is.null(fit_resp_biomass_hier_brm)) saveRDS(fit_resp_biomass_hier_brm, brms_resp_biomass_hier_rds)

# ===== Posterior summaries ====================================================
get_key_par_names_arr <- function(keys) {
  c(
    paste0("b_alpha_Dose_key", keys),
    paste0("b_E_Dose_key", keys)
  )
}

get_key_par_names_ss <- function(keys) {
  c(
    paste0("b_lnB0_Dose_key", keys),
    paste0("b_E_Dose_key", keys),
    paste0("b_Eh_Dose_key", keys),
    paste0("b_Th_Dose_key", keys)
  )
}

extract_arr_posterior <- function(fit, dose_key_tbl) {
  draws <- posterior::as_draws_df(fit) %>% as_tibble()

  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]
    tibble(
      Dose = d,
      Dose_key = dk,
      .draw = seq_len(nrow(draws)),
      alpha = draws[[paste0("b_alpha_Dose_key", dk)]],
      E     = draws[[paste0("b_E_Dose_key", dk)]],
      A     = exp(alpha)
    )
  })
}

extract_arr_cue_common_posterior <- function(fit, dose_key_tbl) {
  draws <- posterior::as_draws_df(fit) %>% as_tibble()

  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]
    tibble(
      Dose = d,
      Dose_key = dk,
      .draw = seq_len(nrow(draws)),
      alpha = draws[[paste0("b_alpha_Dose_key", dk)]],
      E     = draws[["b_E_Intercept"]],
      A     = exp(alpha)
    )
  })
}

extract_ss_posterior <- function(fit, dose_key_tbl) {
  draws <- posterior::as_draws_df(fit) %>% as_tibble()

  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]

    tibble(
      Dose = d,
      Dose_key = dk,
      .draw = seq_len(nrow(draws)),
      lnB0 = draws[[paste0("b_lnB0_Dose_key", dk)]],
      E    = draws[[paste0("b_E_Dose_key", dk)]],
      Eh   = draws[[paste0("b_Eh_Dose_key", dk)]],
      Th   = draws[[paste0("b_Th_Dose_key", dk)]],
      Th_C = draws[[paste0("b_Th_Dose_key", dk)]] - 273.15
    )
  })
}

if (!is.null(fit_growth_brm)) {
  pars_growth <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_growth <- posterior_summary_df(fit_growth_brm, pars = pars_growth)
  readr::write_csv(summ_growth, summary_growth_csv)
  readr::write_csv(extract_arr_posterior(fit_growth_brm, dose_key_tbl), post_growth_csv)
}

if (!is.null(fit_resp_brm)) {
  pars_resp <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_resp <- posterior_summary_df(fit_resp_brm, pars = pars_resp)
  readr::write_csv(summ_resp, summary_resp_csv)
  readr::write_csv(extract_arr_posterior(fit_resp_brm, dose_key_tbl), post_resp_csv)
}

if (!is.null(fit_ratio_brm)) {
  pars_ratio <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_ratio <- posterior_summary_df(fit_ratio_brm, pars = pars_ratio)
  readr::write_csv(summ_ratio, summary_ratio_csv)
  readr::write_csv(extract_arr_posterior(fit_ratio_brm, dose_key_tbl), post_ratio_csv)
}

if (!is.null(fit_cue_brm)) {
  pars_cue <- c(paste0("b_alpha_Dose_key", dose_key_tbl$Dose_key), "b_E_Intercept", "sigma")
  summ_cue <- posterior_summary_df(fit_cue_brm, pars = pars_cue)
  readr::write_csv(summ_cue, summary_cue_csv)
  readr::write_csv(extract_arr_cue_common_posterior(fit_cue_brm, dose_key_tbl), post_cue_csv)
}

if (!is.null(fit_growth_biomass_brm)) {
  pars_growth_biomass <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_growth_biomass <- posterior_summary_df(fit_growth_biomass_brm, pars = pars_growth_biomass)
  readr::write_csv(summ_growth_biomass, summary_growth_biomass_csv)
  readr::write_csv(extract_arr_posterior(fit_growth_biomass_brm, dose_key_tbl), post_growth_biomass_csv)
}

if (!is.null(fit_resp_biomass_brm)) {
  pars_resp_biomass <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_resp_biomass <- posterior_summary_df(fit_resp_biomass_brm, pars = pars_resp_biomass)
  readr::write_csv(summ_resp_biomass, summary_resp_biomass_csv)
  readr::write_csv(extract_arr_posterior(fit_resp_biomass_brm, dose_key_tbl), post_resp_biomass_csv)
}

if (!is.null(fit_growth_tpc_brm)) {
  pars_tpc_g <- c(get_key_par_names_ss(dose_key_tbl$Dose_key), "sigma")
  summ_tpc_g <- posterior_summary_df(fit_growth_tpc_brm, pars = pars_tpc_g)
  readr::write_csv(summ_tpc_g, summary_tpc_growth_csv)
  readr::write_csv(extract_ss_posterior(fit_growth_tpc_brm, dose_key_tbl), post_tpc_growth_csv)
}

if (!is.null(fit_resp_tpc_brm)) {
  pars_tpc_r <- c(get_key_par_names_ss(dose_key_tbl$Dose_key), "sigma")
  summ_tpc_r <- posterior_summary_df(fit_resp_tpc_brm, pars = pars_tpc_r)
  readr::write_csv(summ_tpc_r, summary_tpc_resp_csv)
  readr::write_csv(extract_ss_posterior(fit_resp_tpc_brm, dose_key_tbl), post_tpc_resp_csv)
}

if (!is.null(fit_growth_tpc_biomass_brm)) {
  pars_tpc_g_bio <- c(get_key_par_names_ss(dose_key_tbl$Dose_key), "sigma")
  summ_tpc_g_bio <- posterior_summary_df(fit_growth_tpc_biomass_brm, pars = pars_tpc_g_bio)
  readr::write_csv(summ_tpc_g_bio, summary_tpc_growth_biomass_csv)
  readr::write_csv(extract_ss_posterior(fit_growth_tpc_biomass_brm, dose_key_tbl), post_tpc_growth_biomass_csv)
}

if (!is.null(fit_resp_tpc_biomass_brm)) {
  pars_tpc_r_bio <- c(get_key_par_names_ss(dose_key_tbl$Dose_key), "sigma")
  summ_tpc_r_bio <- posterior_summary_df(fit_resp_tpc_biomass_brm, pars = pars_tpc_r_bio)
  readr::write_csv(summ_tpc_r_bio, summary_tpc_resp_biomass_csv)
  readr::write_csv(extract_ss_posterior(fit_resp_tpc_biomass_brm, dose_key_tbl), post_tpc_resp_biomass_csv)
}

if (!is.null(fit_growth_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_growth_hier_brm, pars = pars), summary_growth_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_growth_hier_brm, dose_key_tbl), post_growth_hier_csv)
}

if (!is.null(fit_resp_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_resp_hier_brm, pars = pars), summary_resp_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_resp_hier_brm, dose_key_tbl), post_resp_hier_csv)
}

if (!is.null(fit_ratio_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_ratio_hier_brm, pars = pars), summary_ratio_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_ratio_hier_brm, dose_key_tbl), post_ratio_hier_csv)
}

if (!is.null(fit_growth_biomass_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_growth_biomass_hier_brm, pars = pars), summary_growth_biomass_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_growth_biomass_hier_brm, dose_key_tbl), post_growth_biomass_hier_csv)
}

if (!is.null(fit_resp_biomass_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_resp_biomass_hier_brm, pars = pars), summary_resp_biomass_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_resp_biomass_hier_brm, dose_key_tbl), post_resp_biomass_hier_csv)
}

# ===== Posterior prediction grids ============================================
make_arr_grid_from_posterior <- function(fit, dat, dose_key_tbl) {
  dr <- posterior::as_draws_df(fit) %>% as_tibble()

  grid_base <- dat %>%
    group_by(Dose, Dose_key) %>%
    summarise(T_min = min(T), T_max = max(T), .groups = "drop") %>%
    right_join(dose_key_tbl, by = c("Dose", "Dose_key")) %>%
    mutate(
      T_min = ifelse(is.na(T_min), min(dat$T, na.rm = TRUE), T_min),
      T_max = ifelse(is.na(T_max), max(dat$T, na.rm = TRUE), T_max)
    )

  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]

    T_seq <- seq(
      grid_base$T_min[grid_base$Dose_key == dk],
      grid_base$T_max[grid_base$Dose_key == dk],
      length.out = 200
    )
    TK_seq <- T_seq + 273.15
    boltz_shift <- (1 / (k_B * T_ref)) - (1 / (k_B * TK_seq))

    alpha <- dr[[paste0("b_alpha_Dose_key", dk)]]
    E     <- dr[[paste0("b_E_Dose_key", dk)]]

    pred_log_mat <- sapply(seq_along(boltz_shift), function(j) alpha + E * boltz_shift[j])
    pred_raw_mat <- exp(pred_log_mat)

    tibble(
      Dose = d,
      Dose_key = dk,
      T = T_seq,
      boltz_shift = boltz_shift,
      log_q2.5  = apply(pred_log_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      log_q50   = apply(pred_log_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      log_q97.5 = apply(pred_log_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      raw_q2.5  = apply(pred_raw_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      raw_q50   = apply(pred_raw_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      raw_q97.5 = apply(pred_raw_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
  })
}

make_arr_grid_from_posterior_common_slope <- function(fit, dat, dose_key_tbl) {
  dr <- posterior::as_draws_df(fit) %>% as_tibble()

  grid_base <- dat %>%
    group_by(Dose, Dose_key) %>%
    summarise(T_min = min(T), T_max = max(T), .groups = "drop") %>%
    right_join(dose_key_tbl, by = c("Dose", "Dose_key")) %>%
    mutate(
      T_min = ifelse(is.na(T_min), min(dat$T, na.rm = TRUE), T_min),
      T_max = ifelse(is.na(T_max), max(dat$T, na.rm = TRUE), T_max)
    )

  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]

    T_seq <- seq(
      grid_base$T_min[grid_base$Dose_key == dk],
      grid_base$T_max[grid_base$Dose_key == dk],
      length.out = 200
    )
    TK_seq <- T_seq + 273.15
    boltz_shift <- (1 / (k_B * T_ref)) - (1 / (k_B * TK_seq))

    alpha <- dr[[paste0("b_alpha_Dose_key", dk)]]
    E     <- dr[["b_E_Intercept"]]

    pred_log_mat <- sapply(seq_along(boltz_shift), function(j) alpha + E * boltz_shift[j])
    pred_raw_mat <- exp(pred_log_mat)

    tibble(
      Dose = d,
      Dose_key = dk,
      T = T_seq,
      boltz_shift = boltz_shift,
      log_q2.5  = apply(pred_log_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      log_q50   = apply(pred_log_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      log_q97.5 = apply(pred_log_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      raw_q2.5  = apply(pred_raw_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      raw_q50   = apply(pred_raw_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      raw_q97.5 = apply(pred_raw_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
  })
}

make_ss_grid_from_posterior <- function(fit, dat, dose_key_tbl) {
  dr <- posterior::as_draws_df(fit) %>% as_tibble()

  grid_base <- dat %>%
    group_by(Dose, Dose_key) %>%
    summarise(T_min = min(T), T_max = max(T), .groups = "drop") %>%
    right_join(dose_key_tbl, by = c("Dose", "Dose_key")) %>%
    mutate(
      T_min = ifelse(is.na(T_min), min(dat$T, na.rm = TRUE), T_min),
      T_max = ifelse(is.na(T_max), max(dat$T, na.rm = TRUE), T_max)
    )

  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]

    T_seq <- seq(
      grid_base$T_min[grid_base$Dose_key == dk],
      grid_base$T_max[grid_base$Dose_key == dk],
      length.out = 200
    )
    TK_seq <- T_seq + 273.15

    lnB0 <- dr[[paste0("b_lnB0_Dose_key", dk)]]
    E    <- dr[[paste0("b_E_Dose_key", dk)]]
    Eh   <- dr[[paste0("b_Eh_Dose_key", dk)]]
    Th   <- dr[[paste0("b_Th_Dose_key", dk)]]

    pred_log_mat <- sapply(seq_along(TK_seq), function(j) {
      lnB0 -
        (E * 11604.51812) * ((1 / TK_seq[j]) - (1 / 293.15)) -
        log(1 + exp((Eh * 11604.51812) * ((1 / Th) - (1 / TK_seq[j]))))
    })
    pred_raw_mat <- exp(pred_log_mat)

    tibble(
      Dose = d,
      Dose_key = dk,
      T = T_seq,
      log_q2.5  = apply(pred_log_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      log_q50   = apply(pred_log_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      log_q97.5 = apply(pred_log_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      raw_q2.5  = apply(pred_raw_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      raw_q50   = apply(pred_raw_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      raw_q97.5 = apply(pred_raw_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
  })
}

arr_grid_growth <- if (!is.null(fit_growth_brm)) make_arr_grid_from_posterior(fit_growth_brm, growth_dat, dose_key_tbl) else tibble()
arr_grid_resp   <- if (!is.null(fit_resp_brm))   make_arr_grid_from_posterior(fit_resp_brm, resp_dat, dose_key_tbl) else tibble()
arr_grid_ratio  <- if (!is.null(fit_ratio_brm))  make_arr_grid_from_posterior(fit_ratio_brm, ratio_dat_arr, dose_key_tbl) else tibble()
arr_grid_cue    <- if (!is.null(fit_cue_brm))    make_arr_grid_from_posterior_common_slope(fit_cue_brm, cue_dat, dose_key_tbl) else tibble()

arr_grid_growth_biomass <- if (!is.null(fit_growth_biomass_brm)) make_arr_grid_from_posterior(fit_growth_biomass_brm, growth_biomass_dat, dose_key_tbl) else tibble()
arr_grid_resp_biomass   <- if (!is.null(fit_resp_biomass_brm))   make_arr_grid_from_posterior(fit_resp_biomass_brm, resp_biomass_dat, dose_key_tbl) else tibble()

tpc_grid_growth <- if (!is.null(fit_growth_tpc_brm)) make_ss_grid_from_posterior(fit_growth_tpc_brm, growth_dat, dose_key_tbl) else tibble()
tpc_grid_resp   <- if (!is.null(fit_resp_tpc_brm))   make_ss_grid_from_posterior(fit_resp_tpc_brm, resp_dat, dose_key_tbl) else tibble()

tpc_grid_growth_biomass <- if (!is.null(fit_growth_tpc_biomass_brm)) make_ss_grid_from_posterior(fit_growth_tpc_biomass_brm, growth_biomass_dat, dose_key_tbl) else tibble()
tpc_grid_resp_biomass   <- if (!is.null(fit_resp_tpc_biomass_brm))   make_ss_grid_from_posterior(fit_resp_tpc_biomass_brm, resp_biomass_dat, dose_key_tbl) else tibble()

arr_grid_growth_hier <- if (!is.null(fit_growth_hier_brm)) make_arr_grid_from_posterior(fit_growth_hier_brm, growth_dat, dose_key_tbl) else tibble()
arr_grid_resp_hier   <- if (!is.null(fit_resp_hier_brm))   make_arr_grid_from_posterior(fit_resp_hier_brm, resp_dat, dose_key_tbl) else tibble()
arr_grid_ratio_hier  <- if (!is.null(fit_ratio_hier_brm))  make_arr_grid_from_posterior(fit_ratio_hier_brm, ratio_dat_arr, dose_key_tbl) else tibble()

arr_grid_growth_biomass_hier <- if (!is.null(fit_growth_biomass_hier_brm)) make_arr_grid_from_posterior(fit_growth_biomass_hier_brm, growth_biomass_dat, dose_key_tbl) else tibble()
arr_grid_resp_biomass_hier   <- if (!is.null(fit_resp_biomass_hier_brm))   make_arr_grid_from_posterior(fit_resp_biomass_hier_brm, resp_biomass_dat, dose_key_tbl) else tibble()

# ===== Save workspace for downstream plotting scripts =========================
# Save an EXPLICIT, curated list of objects rather than save.image(). save.image()
# serialises the entire global environment, so any large objects left over from an
# interactive session (e.g. per-draw prediction objects) get swept into the file
# and can bloat it to several GB. The list below is exactly what 07 consumes:
# datasets, dose metadata, the results plot, the TPC-check plotting helper, the
# Arrhenius/Sharpe-Schoolfield prediction grids, and the fitted brms models
# (which are also written individually as .rds above). intersect(..., ls()) makes
# the save robust to any model that was skipped and therefore never created.
stage2_keep <- c(
  # datasets
  "growth_dat", "resp_dat", "growth_biomass_dat", "resp_biomass_dat",
  "cue_dat", "ratio_dat_arr",
  # dose metadata
  "dose_levels", "dose_cols", "dose_key_tbl",
  # results plot + plotting helper used by 07
  "results_plot", "make_tpc_check_plot",
  # Arrhenius prediction grids
  "arr_grid_growth", "arr_grid_resp", "arr_grid_ratio", "arr_grid_cue",
  "arr_grid_growth_biomass", "arr_grid_resp_biomass",
  "arr_grid_growth_hier", "arr_grid_resp_hier", "arr_grid_ratio_hier",
  "arr_grid_growth_biomass_hier", "arr_grid_resp_biomass_hier",
  # Sharpe-Schoolfield TPC grids
  "tpc_grid_growth", "tpc_grid_resp",
  "tpc_grid_growth_biomass", "tpc_grid_resp_biomass",
  # fitted brms models (also stored individually as .rds)
  "fit_growth_brm", "fit_resp_brm", "fit_ratio_brm", "fit_cue_brm",
  "fit_growth_biomass_brm", "fit_resp_biomass_brm",
  "fit_growth_tpc_brm", "fit_resp_tpc_brm",
  "fit_growth_tpc_biomass_brm", "fit_resp_tpc_biomass_brm",
  "fit_growth_hier_brm", "fit_resp_hier_brm", "fit_ratio_hier_brm",
  "fit_growth_biomass_hier_brm", "fit_resp_biomass_hier_brm"
)
save(list = intersect(stage2_keep, ls()),
     file = file.path(models_dir, "stage2_workspace.RData"))
message("06_bayesian_models.R complete. Workspace saved (curated object list; ",
        "models also stored separately as .rds).")
