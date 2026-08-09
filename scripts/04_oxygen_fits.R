# ===== Stage 1: Oxygen model fits, carbon unit conversions, and descriptive plots
#
# This script:
#  - Fits exponential-decay respiration models to dissolved oxygen time series data
#  - Computes initial cell counts (N0) and carbon-based metrics from oxygen consumption
#  - Derives growth rate, respiration rate, and CUE (Carbon Use Efficiency)
#  - Filters outlier points and generates descriptive plots
#  - Prepares data for downstream Bayesian thermal models
#
# Output files are saved to tables_dir for use by downstream scripts

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
  library(minpack.lm)
  library(brms)
  library(posterior)
  library(grid)
})

options(mc.cores = max(1, parallel::detectCores() - 1))

# ===== Read oxygen data =======================================================
colspec <- list(
  T         = readr::col_double(),
  Dose      = readr::col_character(),
  Replicate = readr::col_character(),
  Time      = readr::col_double(),
  Oxygen    = readr::col_double()
)

o2_file <- readr::read_csv(IN_CSV, col_types = do.call(readr::cols, colspec))

# ── Apply interactive 02b trim windows (so you don't have to re-run 02) ────────
# 03_trim_selector.R writes tables/manual_fit_windows.csv (T, Dose, Replicate,
# fit_start, fit_end). 02 bakes those windows into Oxygen_Data_Filtered.csv only
# when it re-runs; to keep the pipeline linear (01 -> 02 -> 02b -> 03, each once),
# 03 applies any window that differs from what 02 already baked. It re-derives ONLY
# the changed curves by cropping the RAW long series to the window (no smoothing —
# everything downstream uses raw oxygen). Untouched curves keep their filtered data.
manual_start_override <- tibble(T = numeric(0), Dose = character(0),
                                Replicate = character(0), delta_start = numeric(0))
.mfw_path  <- file.path(tables_dir, "manual_fit_windows.csv")
.long_path <- file.path(tables_dir, "Oxygen_All_Long.csv")
if (file.exists(.mfw_path) && file.exists(TRIM_META_CSV) && file.exists(.long_path)) {
  mfw <- suppressWarnings(readr::read_csv(.mfw_path, show_col_types = FALSE))
  if (all(c("T", "Dose", "Replicate", "fit_start", "fit_end") %in% names(mfw)) && nrow(mfw) > 0) {
    metawin <- suppressWarnings(readr::read_csv(TRIM_META_CSV, show_col_types = FALSE)) %>%
      transmute(T = as.numeric(T), Dose = norm_dose(Dose),
                Replicate = toupper(as.character(Replicate)),
                meta_start = as.numeric(peak_time), meta_end = as.numeric(chosen_end_time))
    mfw <- mfw %>%
      transmute(T = as.numeric(T), Dose = norm_dose(Dose),
                Replicate = toupper(as.character(Replicate)),
                fit_start = suppressWarnings(as.numeric(fit_start)),
                fit_end   = suppressWarnings(as.numeric(fit_end))) %>%
      left_join(metawin, by = c("T", "Dose", "Replicate")) %>%
      mutate(eff_start = ifelse(is.finite(fit_start), fit_start, meta_start),
             eff_end   = ifelse(is.finite(fit_end),   fit_end,   meta_end)) %>%
      filter(is.finite(eff_start), is.finite(eff_end), eff_end > eff_start,
             (abs(eff_start - meta_start) > 0.5 | abs(eff_end - meta_end) > 0.5 |
                is.na(meta_start) | is.na(meta_end)))
    if (nrow(mfw) > 0) {
      rawlong <- readr::read_csv(.long_path, show_col_types = FALSE) %>%
        mutate(T = as.numeric(T), Dose = norm_dose(Dose),
               Replicate = toupper(as.character(Replicate)))
      changed_keys <- paste(mfw$T, mfw$Dose, mfw$Replicate)
      o2_file <- o2_file %>%
        mutate(.k = paste(as.numeric(T), norm_dose(Dose), toupper(Replicate))) %>%
        filter(!.k %in% changed_keys) %>% select(-.k)
      rebuilt <- list()
      for (i in seq_len(nrow(mfw))) {
        rw  <- mfw[i, ]
        seg <- rawlong %>%
          filter(T == rw$T, Dose == rw$Dose, Replicate == rw$Replicate,
                 is.finite(Time), is.finite(Oxygen),
                 Time >= rw$eff_start, Time <= rw$eff_end) %>%
          arrange(Time)
        if (nrow(seg) >= 2) {
          rebuilt[[length(rebuilt) + 1L]] <-
            seg %>% transmute(T, Dose, Replicate, Time, Oxygen)
        }
      }
      if (length(rebuilt) > 0) {
        add <- bind_rows(rebuilt)
        for (cc in setdiff(names(o2_file), names(add))) add[[cc]] <- NA
        o2_file <- bind_rows(o2_file, add[names(o2_file)])
      }
      # The N0 anchor (delta_Ninoc_to_N0_min) is the window start; update it for
      # changed curves so N0/R use the new start (matches 02's manual path).
      manual_start_override <- mfw %>% transmute(T, Dose, Replicate, delta_start = eff_start)
      message(sprintf(
        "03: applied %d interactive 02b window(s) from manual_fit_windows.csv.", nrow(mfw)))
    }
  }
}

o2f_raw <- o2_file %>%
  mutate(
    T         = as.numeric(T),
    Dose      = norm_dose(Dose),
    Replicate = toupper(as.character(Replicate)),
    series_id = paste0("T=", T, " | Dose=", Dose, " | Rep=", Replicate),
    # Use the RAW oxygen measurements everywhere — no spline/smoothing is used
    # for fitting or plotting in this pipeline.
    Oxygen_used = Oxygen
  ) %>%
  arrange(T, Dose, Replicate, Time)

stopifnot(all(c("T", "Dose", "Replicate", "Time", "Oxygen_used") %in% names(o2f_raw)))

if (is.null(allowed_doses)) {
  allowed_doses <- dose_levels_from_data(o2f_raw$Dose)
} else {
  allowed_doses <- norm_dose(allowed_doses)
}

o2f <- o2f_raw %>% filter(Dose %in% allowed_doses)
if (nrow(o2f) == 0) stop("No rows remain after filtering by allowed_doses.")

dose_levels <- dose_levels_from_data(o2f$Dose)
dose_levels <- c(intersect("Control", dose_levels), setdiff(dose_levels, "Control"))
dose_cols   <- make_dose_colors(dose_levels)

dose_key_tbl <- tibble(
  Dose = dose_levels,
  Dose_key = paste0("D", seq_along(dose_levels))
)
readr::write_csv(dose_key_tbl, dose_key_csv)

ymin_all <- suppressWarnings(min(o2f$Oxygen_used, na.rm = TRUE))
ymax_all <- suppressWarnings(max(o2f$Oxygen_used, na.rm = TRUE))
pad_all  <- 0.02 * (ymax_all - ymin_all)
if (!is.finite(pad_all)) pad_all <- 0.1
Y_LIMITS_SERIES <- c(ymin_all - pad_all, ymax_all + pad_all)

# ===== Read trimming metadata =================================================
trim_meta <- readr::read_csv(TRIM_META_CSV, show_col_types = FALSE) %>%
  mutate(
    T = as.numeric(T),
    Dose = norm_dose(Dose),
    Replicate = toupper(as.character(Replicate)),
    delta_Ninoc_to_N0_min = as.numeric(delta_Ninoc_to_N0_min)
  ) %>%
  select(T, Dose, Replicate, delta_Ninoc_to_N0_min) %>%
  distinct()

# For curves re-windowed from 02b above, the N0 anchor moves to the new start.
if (nrow(manual_start_override) > 0) {
  trim_meta <- trim_meta %>%
    left_join(manual_start_override, by = c("T", "Dose", "Replicate")) %>%
    mutate(delta_Ninoc_to_N0_min = ifelse(!is.na(delta_start),
                                           delta_start, delta_Ninoc_to_N0_min)) %>%
    select(-delta_start)
}

stopifnot(all(c("T", "Dose", "Replicate", "delta_Ninoc_to_N0_min") %in% names(trim_meta)))

group_lookup <- o2f %>%
  distinct(T, Dose, Replicate) %>%
  left_join(trim_meta, by = c("T", "Dose", "Replicate")) %>%
  mutate(N_inoculation_cells_per_L = as.numeric(N_inoculation_cells_per_L)) %>%
  arrange(T, Dose, Replicate)

readr::write_csv(group_lookup, group_lookup_csv)

# ── Sync intermediate files to the windows 03 ACTUALLY fit ───────────────────
# 02 writes Oxygen_Data_Filtered.csv + Oxygen_Trimmed_Series_Metadata.csv with
# the AUTO window. Here 03 rewrites both to the EFFECTIVE window it used (your
# 02b window where set, else the auto window), so no intermediate file ever shows
# a stale auto window. Also writes effective_fit_windows_used.csv as a plain record.
effective_windows_used <- o2f_raw %>%
  group_by(T, Dose, Replicate) %>%
  summarise(fit_start_used = min(Time, na.rm = TRUE),
            fit_end_used   = max(Time, na.rm = TRUE),
            window_len_min = max(Time, na.rm = TRUE) - min(Time, na.rm = TRUE),
            n_points_fit   = dplyr::n(), .groups = "drop") %>%
  left_join(group_lookup %>% select(T, Dose, Replicate, delta_Ninoc_to_N0_min),
            by = c("T", "Dose", "Replicate")) %>%
  arrange(T, Dose, Replicate)
readr::write_csv(effective_windows_used,
                 file.path(tables_dir, "effective_fit_windows_used.csv"))

# Rewrite Oxygen_Data_Filtered.csv to the effective (02b) windows — raw oxygen.
tryCatch({
  ckey <- readr::read_csv(file.path(tables_dir, "Oxygen_Curve_Code_Key.csv"),
                          show_col_types = FALSE) %>%
    transmute(T = as.numeric(T), Dose = norm_dose(Dose),
              Replicate = toupper(as.character(Replicate)), curve_code)
  o2f_raw %>%
    transmute(T, Dose, Replicate, Time, Oxygen) %>%
    left_join(ckey, by = c("T", "Dose", "Replicate")) %>%
    select(curve_code, T, Dose, Replicate, Time, Oxygen) %>%
    arrange(T, Dose, Replicate, Time) %>%
    readr::write_csv(IN_CSV)
}, error = function(e)
  warning("Could not rewrite Oxygen_Data_Filtered.csv: ", conditionMessage(e)))

# Update the window columns in the metadata to the effective windows.
tryCatch({
  eff <- effective_windows_used %>%
    transmute(T, Dose, Replicate, .es = fit_start_used, .ee = fit_end_used,
              .len = window_len_min)
  mm <- readr::read_csv(TRIM_META_CSV, show_col_types = FALSE) %>%
    mutate(T = as.numeric(T), Dose = norm_dose(Dose),
           Replicate = toupper(as.character(Replicate))) %>%
    left_join(eff, by = c("T", "Dose", "Replicate")) %>%
    mutate(
      peak_time             = ifelse(!is.na(.es), .es, peak_time),
      main_run_start_time   = ifelse(!is.na(.es), .es, main_run_start_time),
      chosen_end_time       = ifelse(!is.na(.ee), .ee, chosen_end_time),
      main_run_end_time     = ifelse(!is.na(.ee), .ee, main_run_end_time),
      fit_duration_min      = ifelse(!is.na(.len), .len, fit_duration_min),
      delta_Ninoc_to_N0_min = ifelse(!is.na(.es), .es, delta_Ninoc_to_N0_min)
    ) %>%
    select(-.es, -.ee, -.len)
  readr::write_csv(mm, TRIM_META_CSV)
}, error = function(e)
  warning("Could not update metadata windows: ", conditionMessage(e)))

message("Synced Oxygen_Data_Filtered.csv + metadata to the effective (02b) windows; ",
        "wrote effective_fit_windows_used.csv.")

if (any(is.na(group_lookup$delta_Ninoc_to_N0_min))) {
  warning("Missing delta_Ninoc_to_N0_min for some groups in trimming metadata. N0 and R will be NA for those groups.")
}

# ===== Fit one oxygen series ==================================================
fit_one <- function(df, y_limits = NULL, rmse_keep_threshold = 0.06) {
  df0 <- df %>%
    arrange(Time) %>%
    mutate(Time0 = Time - min(Time, na.rm = TRUE))

  y <- df0$Oxygen_used

  base_plot <- function(subtitle_txt) {
    p <- ggplot(df0, aes(Time0, Oxygen_used)) +
      geom_point(size = 1.3) +
      labs(
        title = df$series_id[1],
        subtitle = subtitle_txt,
        x = "Time (min, rebased)",
        y = "O₂ (mg/L)"
      ) +
      theme_classic(12)

    if (!is.null(y_limits) && all(is.finite(y_limits))) {
      p <- p + coord_cartesian(ylim = y_limits)
    }
    p
  }

  if (nrow(df0) < 6 || any(!is.finite(y))) {
    return(list(
      coefs = tibble(
        parameter = c("O2_0", "r", "K"),
        Estimate  = c(NA_real_, NA_real_, NA_real_),
        SE        = NA_real_,
        p_value   = NA_real_
      ),
      metrics = tibble(
        T = df$T[1], Dose = df$Dose[1], Replicate = df$Replicate[1],
        n = nrow(df0), r2 = NA_real_, rmse = NA_real_, rss = NA_real_,
        aic = NA_real_, aicc = NA_real_, T_end_min = NA_real_, keep = FALSE
      ),
      keep = FALSE,
      plot = base_plot("Too few/invalid points")
    ))
  }

  # Starting values and bounds are kept IDENTICAL to the guide-curve fit in
  # 03_trim_selector.R (fit_resp / winfit) so the blue guide line you see in 02b
  # IS exactly the model 03 fits. Same model (resp_model), same time rebasing
  # (Time0), same optimiser (nlsLM). If you change these, change them in 02b too.
  t0     <- df0$Time0
  slope0 <- suppressWarnings(median(diff(y) / diff(t0), na.rm = TRUE))
  K_start    <- if (is.finite(slope0)) max(abs(slope0), 1e-6) else 1e-3
  r_start    <- 1e-3
  O2_0_start <- y[1]

  r_lower  <- 1e-6;       r_upper  <- 0.15
  K_lower  <- 1e-10;      K_upper  <- 1
  O2_lower <- min(y) - 1; O2_upper <- max(y) + 1

  fit <- try(
    nlsLM(
      Oxygen_used ~ resp_model(r, K, Time0, O2_0),
      data = df0,
      start = list(r = r_start, K = K_start, O2_0 = O2_0_start),
      lower = c(r = r_lower, K = K_lower, O2_0 = O2_lower),
      upper = c(r = r_upper, K = K_upper, O2_0 = O2_upper),
      control = nls.lm.control(maxiter = 1500, ftol = 1e-12, ptol = 1e-12)
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(list(
      coefs = tibble(
        parameter = c("O2_0", "r", "K"),
        Estimate  = c(NA_real_, NA_real_, NA_real_),
        SE        = NA_real_,
        p_value   = NA_real_
      ),
      metrics = tibble(
        T = df$T[1], Dose = df$Dose[1], Replicate = df$Replicate[1],
        n = nrow(df0), r2 = NA_real_, rmse = NA_real_, rss = NA_real_,
        aic = NA_real_, aicc = NA_real_, T_end_min = NA_real_, keep = FALSE
      ),
      keep = FALSE,
      plot = base_plot("Fit failed")
    ))
  }

  preds <- as.numeric(predict(fit, df0))
  n <- nrow(df0)

  rss <- sum((df0$Oxygen_used - preds)^2, na.rm = TRUE)
  rmse <- sqrt(rss / n)

  ss_tot <- sum((df0$Oxygen_used - mean(df0$Oxygen_used, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- if (is.finite(ss_tot) && ss_tot > 1e-12) 1 - rss / ss_tot else NA_real_

  k_param <- 3L
  aic <- n * log(rss / n) + 2 * k_param
  aicc <- if (n > k_param + 1) aic + (2 * k_param * (k_param + 1)) / (n - k_param - 1) else NA_real_

  T_end_min <- suppressWarnings(max(df0$Time0, na.rm = TRUE))

  co <- coef(fit)
  low_vec <- c(r = r_lower, K = K_lower, O2_0 = O2_lower)
  up_vec  <- c(r = r_upper, K = K_upper, O2_0 = O2_upper)
  on_boundary <- any(abs(co - low_vec[names(co)]) < 1e-10 | abs(co - up_vec[names(co)]) < 1e-10, na.rm = TRUE)

  # Include ALL successfully-fitted curves — no quality filtering (RMSE, min
  # points, or parameter-boundary). Diagnostics (rmse, r2, on_boundary) are still
  # recorded in fit_metrics.csv so you can flag questionable curves yourself.
  quality_ok  <- (rmse < rmse_keep_threshold) && (n >= 6) && !on_boundary
  keep <- TRUE

  co_sum <- as.data.frame(summary(fit)$parameters) %>%
    tibble::rownames_to_column("parameter") %>%
    as_tibble()

  nm <- names(co_sum)
  if ("Std. Error" %in% nm) nm[nm == "Std. Error"] <- "SE"
  if ("Pr(>|t|)"  %in% nm) nm[nm == "Pr(>|t|)"]    <- "p_value"
  names(co_sum) <- nm

  if (!"SE" %in% names(co_sum)) co_sum$SE <- NA_real_
  if (!"p_value" %in% names(co_sum)) co_sum$p_value <- NA_real_

  metrics <- tibble(
    T = df$T[1], Dose = df$Dose[1], Replicate = df$Replicate[1],
    n = n, r2 = r2, rmse = rmse, rss = rss,
    aic = aic, aicc = aicc, T_end_min = T_end_min,
    on_boundary = on_boundary, quality_ok = quality_ok, keep = keep
  )

  subtitle_txt <- paste0(
    "R²=", ifelse(is.na(r2), "NA", sprintf("%.3f", r2)),
    " | RMSE=", sprintf("%.3g", rmse),
    " | AIC=", sprintf("%.1f", aic),
    if (!quality_ok) " — low quality (kept)" else ""
  )

  p <- ggplot(df0, aes(Time0, Oxygen_used)) +
    geom_point(size = 1.3) +
    geom_line(aes(y = preds), linewidth = 0.9, color = "red") +
    labs(
      title = df$series_id[1],
      subtitle = subtitle_txt,
      x = "Time (min, rebased)",
      y = "O₂ (mg/L)"
    ) +
    theme_classic(12)

  if (!is.null(y_limits) && all(is.finite(y_limits))) {
    p <- p + coord_cartesian(ylim = y_limits)
  }

  list(coefs = co_sum, metrics = metrics, keep = keep, plot = p)
}

# ===== Run oxygen fits ========================================================
groups <- o2f %>%
  group_by(T, Dose, Replicate) %>%
  group_split()

all_coef_rows <- list()
all_metrics   <- list()

pdf(pdf_path, width = 6.8, height = 4.6)
for (g in groups) {
  res <- fit_one(g, y_limits = Y_LIMITS_SERIES, rmse_keep_threshold = RMSE_KEEP_THRESHOLD)
  all_metrics[[length(all_metrics) + 1]] <- res$metrics

  if (isTRUE(res$keep)) {
    print(res$plot)
    all_coef_rows[[length(all_coef_rows) + 1]] <-
      tibble(T = g$T[1], Dose = g$Dose[1], Replicate = g$Replicate[1]) %>%
      bind_cols(res$coefs) %>%
      mutate(T_end_min = res$metrics$T_end_min[1])
  }
}
dev.off()

coef_out <- bind_rows(all_coef_rows)
fit_metrics_out <- bind_rows(all_metrics)

readr::write_csv(coef_out, coef_csv)
readr::write_csv(fit_metrics_out, fit_metrics_csv)

# ===== Wide coefficients ======================================================
coef_wide <- coef_out %>%
  select(T, Dose, Replicate, parameter, Estimate, T_end_min) %>%
  tidyr::pivot_wider(names_from = parameter, values_from = Estimate) %>%
  group_by(T, Dose, Replicate) %>%
  summarise(
    r = first(r),
    K = first(K),
    O2_0 = first(O2_0),
    T_end_min = first(T_end_min),
    .groups = "drop"
  ) %>%
  arrange(T, Dose, Replicate)

readr::write_csv(coef_wide, coef_wide_csv)

# ===== Compute sample-specific N0, respiration, carbon units, CUE ============
results <- coef_wide %>%
  left_join(group_lookup, by = c("T", "Dose", "Replicate")) %>%
  mutate(
    N0_cells_per_L = dplyr::if_else(
      is.finite(N_inoculation_cells_per_L) & N_inoculation_cells_per_L > 0 &
        is.finite(delta_Ninoc_to_N0_min) & delta_Ninoc_to_N0_min >= 0 &
        is.finite(r) & r > 0,
      N_inoculation_cells_per_L * exp(r * delta_Ninoc_to_N0_min),
      NA_real_
    ),
    C_tot_O2_mg_per_L = dplyr::if_else(
      is.finite(K) & is.finite(r) & r > 0 &
        is.finite(T_end_min) & T_end_min > 0,
      (K / r) * (exp(r * T_end_min) - 1),
      NA_real_
    ),
    biomass_integral_cells_min_per_L = dplyr::if_else(
      is.finite(N0_cells_per_L) & N0_cells_per_L > 0 &
        is.finite(r) & r > 0 &
        is.finite(T_end_min) & T_end_min > 0,
      N0_cells_per_L * (exp(r * T_end_min) - 1) / r,
      NA_real_
    ),
    R_O2_mg_cell_min = dplyr::if_else(
      is.finite(C_tot_O2_mg_per_L) & C_tot_O2_mg_per_L > 0 &
        is.finite(biomass_integral_cells_min_per_L) & biomass_integral_cells_min_per_L > 0,
      C_tot_O2_mg_per_L / biomass_integral_cells_min_per_L,
      NA_real_
    ),
    cell_volume_um3 = CELL_VOLUME_UM3,
    cell_carbon_fg = CELL_CARBON_FG_PER_CELL,
    growth_fgC_h = dplyr::if_else(
      is.finite(r) & r > 0,
      r * cell_carbon_fg * MIN_TO_H,
      NA_real_
    ),
    respiration_fgC_h = dplyr::if_else(
      is.finite(R_O2_mg_cell_min) & R_O2_mg_cell_min > 0,
      R_O2_mg_cell_min * MG_TO_FG * RESPIRATORY_QUOTIENT * MIN_TO_H,
      NA_real_
    ),
    growth_C_per_C_h = dplyr::if_else(
      is.finite(growth_fgC_h) & growth_fgC_h > 0 &
        is.finite(cell_carbon_fg) & cell_carbon_fg > 0,
      growth_fgC_h / cell_carbon_fg,
      NA_real_
    ),
    respiration_C_per_C_h = dplyr::if_else(
      is.finite(respiration_fgC_h) & respiration_fgC_h > 0 &
        is.finite(cell_carbon_fg) & cell_carbon_fg > 0,
      respiration_fgC_h / cell_carbon_fg,
      NA_real_
    ),
    CUE = dplyr::if_else(
      is.finite(growth_fgC_h) & growth_fgC_h > 0 &
        is.finite(respiration_fgC_h) & respiration_fgC_h > 0,
      growth_fgC_h / (growth_fgC_h + respiration_fgC_h),
      NA_real_
    ),
    resp_over_growth = dplyr::if_else(
      is.finite(respiration_fgC_h) & respiration_fgC_h > 0 &
        is.finite(growth_fgC_h) & growth_fgC_h > 0,
      respiration_fgC_h / growth_fgC_h,
      NA_real_
    )
  ) %>%
  arrange(T, Dose, Replicate)

# ===== Exclude specific outlier points =======================================
# Excluded samples are chosen INTERACTIVELY in 03_trim_selector.R, which writes
# tables/plot_exclude_points.csv (columns: T, Dose, Replicate). Select/deselect
# there — do NOT hardcode the list here. If the file is missing or empty, no
# samples are excluded (clean slate).
exclude_points_csv <- file.path(tables_dir, "plot_exclude_points.csv")
EXCLUDE_POINTS <- if (file.exists(exclude_points_csv)) {
  ep <- readr::read_csv(exclude_points_csv, show_col_types = FALSE)
  if (all(c("T", "Dose", "Replicate") %in% names(ep)) && nrow(ep) > 0) {
    tibble(T         = as.numeric(ep$T),
           Dose      = as.character(ep$Dose),
           Replicate = toupper(as.character(ep$Replicate)))
  } else {
    tibble(T = numeric(0), Dose = character(0), Replicate = character(0))
  }
} else {
  tibble(T = numeric(0), Dose = character(0), Replicate = character(0))
}

results <- results %>%
  anti_join(EXCLUDE_POINTS, by = c("T", "Dose", "Replicate"))

message(sprintf(
  "Exclusion filter applied: %d point(s) removed (from %s). %d rows retained.",
  nrow(EXCLUDE_POINTS), basename(exclude_points_csv), nrow(results)
))
# =============================================================================

readr::write_csv(results, derived_csv)

# ===== Descriptive plots ======================================================
results_plot <- results %>%
  filter(Dose %in% allowed_doses) %>%
  mutate(Dose = factor(Dose, levels = dose_levels))

# ===== Replication summary ====================================================
replication_summary <- results_plot %>%
  count(T, Dose, name = "n_replicates") %>%
  arrange(Dose, T)

readr::write_csv(replication_summary, replication_summary_csv)

p_box_growth <- ggplot(results_plot %>% filter(is.finite(growth_fgC_h)), aes(Dose, growth_fgC_h, fill = Dose)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~T, scales = "free_y") +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Growth in carbon units", x = NULL, y = "Growth (fg C h^-1)") +
  theme_classic(12) +
  theme(legend.position = "none")

p_box_resp <- ggplot(results_plot %>% filter(is.finite(respiration_fgC_h)), aes(Dose, respiration_fgC_h, fill = Dose)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~T, scales = "free_y") +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Respiration in carbon units", x = NULL, y = "Respiration (fg C h^-1)") +
  theme_classic(12) +
  theme(legend.position = "none")

ratio_dat <- results_plot %>%
  filter(
    is.finite(growth_fgC_h), growth_fgC_h > 0,
    is.finite(respiration_fgC_h), respiration_fgC_h > 0
  ) %>%
  mutate(
    log_resp_over_growth = log(resp_over_growth),
    TK = T + 273.15
  )

p_rg_t <- ggplot(ratio_dat, aes(T, resp_over_growth, color = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Respiration / Growth vs Temperature", x = "Temperature (°C)", y = "Respiration / Growth", color = "Dose") +
  theme_classic(12)

p_growth_vs_T <- ggplot(results_plot %>% filter(is.finite(growth_fgC_h), growth_fgC_h > 0), aes(T, growth_fgC_h, color = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Growth vs Temperature", x = "Temperature (°C)", y = "Growth (fg C h^-1)", color = "Dose") +
  theme_classic(12)

p_resp_vs_T <- ggplot(results_plot %>% filter(is.finite(respiration_fgC_h), respiration_fgC_h > 0), aes(T, respiration_fgC_h, color = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Respiration vs Temperature", x = "Temperature (°C)", y = "Respiration (fg C h^-1)", color = "Dose") +
  theme_classic(12)

p_box_growth_biomass <- ggplot(
  results_plot %>% filter(is.finite(growth_C_per_C_h), growth_C_per_C_h > 0),
  aes(Dose, growth_C_per_C_h, fill = Dose)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~T, scales = "free_y") +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Biomass-corrected growth", x = NULL, y = "Growth (C per C per h)") +
  theme_classic(12) +
  theme(legend.position = "none")

p_box_resp_biomass <- ggplot(
  results_plot %>% filter(is.finite(respiration_C_per_C_h), respiration_C_per_C_h > 0),
  aes(Dose, respiration_C_per_C_h, fill = Dose)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~T, scales = "free_y") +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Biomass-corrected respiration", x = NULL, y = "Respiration (C per C per h)") +
  theme_classic(12) +
  theme(legend.position = "none")

p_growth_biomass_vs_T <- ggplot(
  results_plot %>% filter(is.finite(growth_C_per_C_h), growth_C_per_C_h > 0),
  aes(T, growth_C_per_C_h, color = Dose)
) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Biomass-corrected growth vs Temperature", x = "Temperature (°C)", y = "Growth (C per C per h)", color = "Dose") +
  theme_classic(12)

p_resp_biomass_vs_T <- ggplot(
  results_plot %>% filter(is.finite(respiration_C_per_C_h), respiration_C_per_C_h > 0),
  aes(T, respiration_C_per_C_h, color = Dose)
) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Biomass-corrected respiration vs Temperature", x = "Temperature (°C)", y = "Respiration (C per C per h)", color = "Dose") +
  theme_classic(12)

# ===== Prepare data for Bayesian thermal models ===============================
growth_dat <- results_plot %>%
  filter(is.finite(growth_fgC_h), growth_fgC_h > 0) %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  transmute(
    T = as.numeric(T),
    TK = T + 273.15,
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * (T + 273.15))),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(growth_fgC_h)),
    y_raw = as.numeric(growth_fgC_h)
  )

resp_dat <- results_plot %>%
  filter(is.finite(respiration_fgC_h), respiration_fgC_h > 0) %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  transmute(
    T = as.numeric(T),
    TK = T + 273.15,
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * (T + 273.15))),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(respiration_fgC_h)),
    y_raw = as.numeric(respiration_fgC_h)
  )

growth_biomass_dat <- results_plot %>%
  filter(is.finite(growth_C_per_C_h), growth_C_per_C_h > 0) %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  transmute(
    T = as.numeric(T),
    TK = T + 273.15,
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * (T + 273.15))),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(growth_C_per_C_h)),
    y_raw = as.numeric(growth_C_per_C_h)
  )

resp_biomass_dat <- results_plot %>%
  filter(is.finite(respiration_C_per_C_h), respiration_C_per_C_h > 0) %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  transmute(
    T = as.numeric(T),
    TK = T + 273.15,
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * (T + 273.15))),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(respiration_C_per_C_h)),
    y_raw = as.numeric(respiration_C_per_C_h)
  )

ratio_dat_arr <- ratio_dat %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  filter(is.finite(log_resp_over_growth)) %>%
  transmute(
    T = as.numeric(T),
    TK = as.numeric(TK),
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * TK)),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = as.numeric(log_resp_over_growth),
    y_raw = as.numeric(resp_over_growth)
  )

cue_dat <- ratio_dat %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  filter(is.finite(CUE), CUE > 0) %>%
  transmute(
    T = as.numeric(T),
    TK = as.numeric(TK),
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * TK)),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(CUE)),
    y_raw = as.numeric(CUE)
  )

# Save intermediates for downstream scripts
saveRDS(results_plot, file.path(models_dir, "results_plot.rds"))
saveRDS(ratio_dat,    file.path(models_dir, "ratio_dat.rds"))
saveRDS(dose_levels,  file.path(models_dir, "dose_levels.rds"))
saveRDS(dose_cols,    file.path(models_dir, "dose_cols.rds"))
saveRDS(dose_key_tbl, file.path(models_dir, "dose_key_tbl.rds"))

# Save data prepared for Bayesian models
saveRDS(growth_dat,          file.path(models_dir, "growth_dat.rds"))
saveRDS(resp_dat,            file.path(models_dir, "resp_dat.rds"))
saveRDS(ratio_dat_arr,       file.path(models_dir, "ratio_dat_arr.rds"))
saveRDS(cue_dat,             file.path(models_dir, "cue_dat.rds"))
saveRDS(growth_biomass_dat,  file.path(models_dir, "growth_biomass_dat.rds"))
saveRDS(resp_biomass_dat,    file.path(models_dir, "resp_biomass_dat.rds"))

message("04_oxygen_fits.R complete.")
