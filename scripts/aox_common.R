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
# 41_sham_rates.R (bounded nlsLM); shared by 47_ptc_sham_rates.R and the
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

# --- trace helpers shared by 03_trim_selector.R, 45, 46 and 74 -----------------
# Offset correction for a handling step (plates briefly removed and re-inserted):
# readings after step_time_h are shifted by -step_mgL. Offset only.
step_correct <- function(t_min, y, step_time_h = NA, step_mgL = 0) {
  if (is.finite(step_time_h) && is.finite(step_mgL) && step_mgL != 0)
    y[t_min > step_time_h * 60] <- y[t_min > step_time_h * 60] - step_mgL
  y
}
# Denoising: centred moving average of width smooth_h hours (0 = none). For the
# model O(t) = O2_0 + (K/r)(1 - exp(rt)) a centred moving average leaves r
# unchanged and multiplies K by sinh(rW/2)/(rW/2) (1.002 for r = 0.06 h^-1,
# W = 4 h), so the growth rate is unbiased by construction. Returns NA at the
# edges (half a window at each end).
smooth_ma <- function(t_min, y, smooth_h = 0) {
  if (!is.finite(smooth_h) || smooth_h <= 0) return(y)
  dt <- median(diff(t_min), na.rm = TRUE); k <- max(3L, as.integer(round(smooth_h * 60 / dt)))
  if (k %% 2 == 0) k <- k + 1L
  ok <- is.finite(y); ys <- rep(NA_real_, length(y))
  ys[ok] <- as.numeric(stats::filter(y[ok], rep(1 / k, k), sides = 2))
  ys
}

# =============================================================================
# Rule-based fitting intervals for the prothioconazole x SHAM factorial (45-47,
# 57, 74). No interval is chosen by hand: ps_window() derives one from the trace
# by a named rule, so the same code produces the primary intervals in
# 47_ptc_sham_rates.R and every alternative in 48_ptc_sham_window_rules.R.
#
# Bounds common to every rule, applied first:
#   start >= PS_MIN_START   the first ~3 h is the vial equilibrating (O2 falls
#                           from ~15 to ~11 mg/L and then flattens); that fall
#                           is decelerating while the model is accelerating, so
#                           including it biases r down and K up.
#   end    <  the time O2 first drops below PS_O2_FLOOR, i.e. before the trace
#                           becomes a plateau that carries no rate information.
#
# Rules:
#   stable_r      each curve's own exponential phase: fit r for a series of end
#                 points and keep the longest run over which r stays within
#                 PS_R_TOL of that run's median. Adapts to lag and to early
#                 stationary phase, and stops where the curve leaves exponential
#                 growth, without any reference to the result.
#   fixed(a,b)    the same clock interval for every curve (minutes).
#   draw(lo,hi)   between lo and hi of that curve's own total oxygen drawdown.
#   peak(dur)     from the post-equilibration maximum for dur minutes.
# =============================================================================
PS_MIN_START <- 180      # minutes
PS_O2_FLOOR  <- 2.5      # mg/L
PS_R_TOL     <- 0.05     # stable_r: fractional tolerance on r
PS_MIN_PTS   <- 120      # minimum points in a fitted interval

ps_bounds <- function(t_min, o2) {
  ok <- is.finite(o2) & is.finite(t_min)
  t_min <- t_min[ok]; o2 <- o2[ok]
  lo <- max(PS_MIN_START, min(t_min))
  below <- which(o2 < PS_O2_FLOOR & t_min > lo)
  hi <- if (length(below)) t_min[below[1]] else max(t_min)
  c(lo, hi)
}
# post-equilibration maximum, searched inside the bounds
ps_peak <- function(t_min, o2, bd) {
  cand <- which(t_min >= bd[1] & t_min <= bd[1] + (bd[2] - bd[1]) / 3 & is.finite(o2))
  if (!length(cand)) return(bd[1])
  t_min[cand[which.max(o2[cand])]]
}
ps_window <- function(t_min, o2, rule = list(type = "stable_r")) {
  ok <- is.finite(o2) & is.finite(t_min); t_min <- t_min[ok]; o2 <- o2[ok]
  if (length(t_min) < PS_MIN_PTS) return(c(NA_real_, NA_real_, NA))
  bd <- ps_bounds(t_min, o2)
  if (!is.finite(bd[2]) || bd[2] - bd[1] < 60) return(c(NA_real_, NA_real_, NA))
  fallback <- FALSE
  w <- switch(rule$type,
    fixed = c(max(bd[1], rule$a), min(bd[2], rule$b)),
    peak  = { p <- ps_peak(t_min, o2, bd); c(p, min(bd[2], p + rule$dur)) },
    draw  = { m <- t_min >= bd[1] & t_min <= bd[2]
              top <- max(o2[m]); bot <- min(o2[m]); rng <- top - bot
              if (rng <= 0) c(NA, NA) else {
                s <- t_min[m][which(o2[m] <= top - rule$a * rng)[1]]
                e <- t_min[m][which(o2[m] <= top - rule$b * rng)[1]]
                c(if (is.na(s)) bd[1] else s, if (is.na(e)) bd[2] else e) } },
    stable_r = {
      s <- ps_peak(t_min, o2, bd)
      idx <- which(t_min > s & t_min <= bd[2])
      res <- c(NA, NA)
      if (length(idx) >= PS_MIN_PTS) {
        cand <- idx[seq(PS_MIN_PTS, length(idx), length.out = min(25, length(idx) %/% 20))]
        rr <- sapply(cand, function(j) {
          m <- t_min >= s & t_min <= t_min[j]
          f <- fit_o2_model(t_min[m], o2[m]); if (is.null(f)) NA_real_ else f$r_per_h })
        good <- is.finite(rr)
        best <- -Inf
        for (a in seq_along(cand)) if (good[a]) for (b in length(cand):a) {
          if (!all(good[a:b])) next
          seg <- rr[a:b]; md <- median(seg)
          if (all(abs(seg - md) <= PS_R_TOL * abs(md))) {
            len <- t_min[cand[b]] - s
            if (len > best) { best <- len; res <- c(s, t_min[cand[b]]) }
            break
          }
        }
      }
      if (!is.finite(res[1])) { fallback <- TRUE; res <- bd }   # no stable stretch: plate bounds
      res
    },
    stop("unknown rule type: ", rule$type))
  if (!all(is.finite(w)) || w[2] - w[1] < 60) { fallback <- TRUE; w <- bd }
  c(w[1], w[2], fallback)
}
