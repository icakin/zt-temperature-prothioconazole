# =============================================================================
# aox_common.R -- shared helpers for the azide/SHAM respirometry (Figure 6c-d,
# Figure S17). Sourced by 43_azide_rates.R, 56_figure6.R and 73_figureS17.R.
#
#   read_sdr_xlsx(path)    PreSens SensorDish export (sheet "Oxygen Orginal")
#                          -> long data.frame: plate, well, t_h, o2
#   initial_rate(t_h, o2)  linear O2-consumption rate over the initial window:
#                          from the post-equilibration maximum (15-point moving
#                          average, first third of the trace) to +WIN_H hours.
#                          Returns rate (mg L^-1 h^-1, positive = consumption),
#                          R^2, window start and end (h).
# =============================================================================
suppressPackageStartupMessages(library(readxl))

WIN_H <- 10          # initial-rate window length (h) after the peak
AZ_COND <- c(B = "blank", A = "azide", AS = "azide+SHAM")
AZ_LAB  <- c(blank = "No cells", azide = "Azide 0.5 mM", `azide+SHAM` = "Azide 0.5 mM + SHAM 1 mM")
AZ_COL  <- c(blank = "#9A9A9A", azide = "#E08214", `azide+SHAM` = "#433D84")
AZ_ORDER <- c("blank", "azide", "azide+SHAM")

read_sdr_xlsx <- function(path, plate = NA) {
  raw <- suppressMessages(read_excel(path, sheet = "Oxygen Orginal", skip = 12,
                                     col_types = "text", .name_repair = "minimal"))
  nm <- names(raw)
  tcol <- which(nm == "Time/Min.")[1]
  wells <- nm[grepl("^[A-D][1-6]$", nm)]
  t_h <- suppressWarnings(as.numeric(raw[[tcol]])) / 60
  keep <- is.finite(t_h)
  out <- lapply(wells, function(w) data.frame(
    plate = plate, well = w, t_h = t_h[keep],
    o2 = suppressWarnings(as.numeric(raw[[w]][keep])), stringsAsFactors = FALSE))
  do.call(rbind, out)
}

peak_index <- function(o2, k = 15) {
  ok <- is.finite(o2); sm <- o2
  sm[ok] <- stats::filter(o2[ok], rep(1 / k, k), sides = 2)
  sm[is.na(sm)] <- -Inf
  n <- length(sm); which.max(sm[seq_len(max(3, n %/% 3))])
}

initial_rate <- function(t_h, o2, win = WIN_H) {
  ok <- is.finite(o2) & is.finite(t_h); t_h <- t_h[ok]; o2 <- o2[ok]
  pk <- peak_index(o2); t0 <- t_h[pk]
  m <- t_h >= t0 & t_h < t0 + win & o2 > 0.5
  f <- lm(o2 ~ t_h, data = data.frame(t_h = t_h[m], o2 = o2[m]))
  list(rate = -unname(coef(f)[2]), r2 = summary(f)$r.squared,
       t_start = t0, t_end = max(t_h[m]), intercept = unname(coef(f)[1]))
}

load_azide_plates <- function(root = ".") {
  lay <- read.csv(file.path(root, "data/oxygen/azide_layout.csv"), stringsAsFactors = FALSE)
  files <- c(`1` = "azide_rep1_Oxygen.xlsx", `2` = "azide_rep2_Oxygen.xlsx", `3` = "azide_rep3_Oxygen.xlsx")
  tr <- do.call(rbind, lapply(names(files), function(p)
    read_sdr_xlsx(file.path(root, "data/oxygen", files[[p]]), plate = as.integer(p))))
  tr$condition <- lay$condition[match(tr$well, lay$well)]
  tr$condition <- factor(tr$condition, levels = AZ_ORDER)
  tr
}

# --- the manuscript oxygen model, fitted exactly as in 04_oxygen_fits.R and
# 41_sham_rates.R (bounded nlsLM); shared by 46_ptc_sham_rates.R and the
# figure scripts. tt in minutes; returns r in h^-1.
resp_model <- function(t, r, K, O2_0) O2_0 + (K/r) * (1 - exp(r * t))
fit_o2_model <- function(tt, yy) {
  if (!requireNamespace("minpack.lm", quietly = TRUE)) stop("minpack.lm is required")
  t0 <- tt - min(tt); sl <- median(diff(yy) / diff(t0))
  st <- list(r = 1e-3, K = min(max(abs(sl), 1e-6), 1), O2_0 = yy[1])
  f <- try(minpack.lm::nlsLM(yy ~ resp_model(t0, r, K, O2_0), start = st,
                 lower = c(r = 1e-6, K = 1e-10, O2_0 = min(yy) - 1),
                 upper = c(r = 0.15,  K = 1,     O2_0 = max(yy) + 1),
                 control = minpack.lm::nls.lm.control(maxiter = 200, ftol = 1e-12, ptol = 1e-12)),
           silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  p <- coef(f); pred <- resp_model(t0, p[["r"]], p[["K"]], p[["O2_0"]])
  rss <- sum((yy - pred)^2); sst <- sum((yy - mean(yy))^2)
  list(r_per_h = unname(p[["r"]]) * 60, K = unname(p[["K"]]), O2_0 = unname(p[["O2_0"]]),
       R2 = if (sst > 1e-12) 1 - rss/sst else NA_real_, RMSE = sqrt(rss/length(yy)),
       n_pts = length(yy), pred = pred)
}

# prothioconazole x SHAM factorial (45/46/57/74)
PS_ORDER <- c("V", "S", "P2", "P2S", "P4", "P4S")
PS_LAB   <- c(V = "Vehicle", S = "SHAM", P2 = "PTC 2 mg/L", P2S = "PTC 2 + SHAM",
              P4 = "PTC 4 mg/L", P4S = "PTC 4 + SHAM")
PS_COL   <- c(V = "#8C8C8C", S = "#433D84", P2 = "#E5A24A", P2S = "#7B5EA7",
              P4 = "#B2182B", P4S = "#2A1A5E")

# Trace preparation for the prothioconazole x SHAM plates, applied identically
# before trimming (45 -> selector), fitting (46) and plotting (74):
#   1. handling-step offset: readings after step_time_h are shifted by -step_mgL
#      (the plates were briefly removed for an interim export at ~41 h);
#   2. denoising: centred moving average of width PS_SMOOTH_H hours. The 27 C
#      plates carry a 3-4 h re-aeration sawtooth (cell wells only) after the
#      plates were moved to a second incubator; a 4 h average removes a 4 h
#      cycle exactly and most of a 3 h one, while leaving the multi-hour growth
#      trend untouched (cf. 44_sham_denoise_sensitivity.R). Set PS_SMOOTH_H <- 0
#      to fit the raw traces.
PS_SMOOTH_H <- 4
prep_trace <- function(t_min, o2, step_time_h = NA, step_mgL = 0, smooth_h = PS_SMOOTH_H) {
  y <- o2
  if (is.finite(step_time_h) && is.finite(step_mgL) && step_mgL != 0)
    y[t_min > step_time_h * 60] <- y[t_min > step_time_h * 60] - step_mgL
  if (smooth_h > 0) {
    dt <- median(diff(t_min), na.rm = TRUE); k <- max(3L, as.integer(round(smooth_h * 60 / dt)))
    if (k %% 2 == 0) k <- k + 1L
    ok <- is.finite(y); ys <- rep(NA_real_, length(y))
    ys[ok] <- as.numeric(stats::filter(y[ok], rep(1 / k, k), sides = 2))
    y <- ys
  }
  y
}
