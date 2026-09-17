#!/usr/bin/env Rscript
# =============================================================================
# 44_sham_denoise_sensitivity.R -- does the periodic incubator component on the
# 27 C SHAM plate affect the growth rates or the temperature contrast?
#
# The 27 C inhibitor plates carry a periodic component (period ~1.5-4 h, up to
# ~2 mg/L peak-to-peak late in the run) that is absent from the 15 C plate and
# from the seven dose-response plates; it is common to every well and is read as
# incubator temperature cycling picked up by the oxygen sensor spots. This script
# removes it and refits every SHAM curve with exactly the procedure of
# 41_sham_rates.R, on the same hand-selected intervals.
#
# Denoising: centred moving average whose width equals one cycle. A moving
# average of width W cancels a periodic component of period W exactly and
# attenuates nearby periods strongly, while leaving the multi-hour growth trend
# untouched. Two widths are used (1.5 h and 2 h) to bracket the observed period.
#
# Inputs : data/aox/sham_{15,27}_Oxygen.csv, tables/aox/sham_fit_windows.csv
# Outputs: tables/aox/sham_denoise_sensitivity.csv  per curve: raw and denoised r
#          tables/aox/sham_denoise_summary.csv      per plate: residual noise,
#                                                  control CV, slopes, contrast, P
# Run from the repository root:  Rscript scripts/44_sham_denoise_sensitivity.R
# =============================================================================
suppressPackageStartupMessages(library(minpack.lm))
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
WIDTHS_H <- c(2, 4, 12)        # moving-average widths (h): bracket the cycle and
                               # the width used in the primary analysis
DROP_DOSE <- 0.35              # as in 56_figure6.R

resp_model <- function(t, r, K, O2_0) O2_0 + (K/r) * (1 - exp(r * t))
fit_one <- function(tt, yy) {                       # identical to 41_sham_rates.R
  t0 <- tt - min(tt); sl <- median(diff(yy) / diff(t0))
  st <- list(r = 1e-3, K = min(max(abs(sl), 1e-6), 1), O2_0 = yy[1])
  f <- try(nlsLM(yy ~ resp_model(t0, r, K, O2_0), start = st,
                 lower = c(r = 1e-6, K = 1e-10, O2_0 = min(yy) - 1),
                 upper = c(r = 0.15,  K = 1,     O2_0 = max(yy) + 1),
                 control = nls.lm.control(maxiter = 200, ftol = 1e-12, ptol = 1e-12)), silent = TRUE)
  if (inherits(f, "try-error")) return(NA_real_)
  unname(coef(f)[["r"]]) * 60
}
movavg <- function(y, k) {                          # centred, NA at the edges
  k <- max(3L, as.integer(k)); if (k %% 2 == 0) k <- k + 1L
  as.numeric(stats::filter(y, rep(1/k, k), sides = 2))
}
# residual noise after removing the multi-hour trend (6 h running median), and
# the dominant period of that residual between 0.5 and 12 h
noise_stats <- function(t_h, y) {
  ok <- is.finite(y); y <- y[ok]; t_h <- t_h[ok]
  dt <- median(diff(t_h)); w <- as.integer(6/dt); if (w %% 2 == 0) w <- w + 1L
  tr <- stats::runmed(y, w, endrule = "median"); r <- y - tr
  F <- abs(fft(r - mean(r)))[seq_len(length(r) %/% 2)]
  fr <- seq_along(F) / (length(r) * dt)
  m <- fr > 1/12 & fr < 1/0.5
  c(resid_sd = sd(r), period_h = 1/fr[m][which.max(F[m])],
    late_p2p = diff(range(r[seq(floor(2*length(r)/3), length(r))])))
}

# Intervals: the hand-selected windows used by 41_sham_rates.R, so this test is
# run on exactly the curves and intervals behind Figure 6a,b.
MANW <- file.path(ROOT, "tables/aox/sham/manual_fit_windows.csv")
win  <- read.csv(if (file.exists(MANW)) MANW else
                 file.path(ROOT, "tables/aox/sham_fit_windows.csv"), stringsAsFactors = FALSE)
rows <- list(); noise <- list()
for (temp in c(15, 27)) {
  d <- read.csv(file.path(ROOT, sprintf("data/aox/sham_%d_Oxygen.csv", temp)), check.names = FALSE)
  dt_min <- median(diff(d$Time))
  for (cv in setdiff(names(d), c("Time", "T"))) {
    dose <- sub("_R[0-9]+$", "", cv); rep <- as.integer(sub("^.*_R", "", cv))
    w <- win[win$T == temp & win$Dose == dose & win$Replicate == paste0("R", rep), ][1, ]
    if (!nrow(w) || is.na(w$fit_start)) next
    y <- d[[cv]]
    ns <- noise_stats(d$Time/60, y)
    noise[[length(noise) + 1]] <- data.frame(temp, dose, culture = rep, t(ns))
    m <- d$Time >= w$fit_start & d$Time <= w$fit_end
    r_raw <- fit_one(d$Time[m & is.finite(y)], y[m & is.finite(y)])
    r_dn <- sapply(WIDTHS_H, function(h) {
      ys <- movavg(y, h*60/dt_min); mm <- m & is.finite(ys)
      fit_one(d$Time[mm], ys[mm])
    })
    rows[[length(rows) + 1]] <- data.frame(
      temp, dose_mM = if (dose == "Control") 0 else as.numeric(dose), culture = rep,
      r_raw = r_raw, r_dn_2h = r_dn[1], r_dn_4h = r_dn[2], r_dn_12h = r_dn[3],
      pct_change_2h = 100*(r_dn[1]/r_raw - 1), pct_change_4h = 100*(r_dn[2]/r_raw - 1),
      pct_change_12h = 100*(r_dn[3]/r_raw - 1))
  }
}
res <- do.call(rbind, rows); res <- res[order(res$temp, res$dose_mM, res$culture), ]
noise <- do.call(rbind, noise)
write.csv(res, file.path(ROOT, "tables/aox/sham_denoise_sensitivity.csv"), row.names = FALSE)

# --- the Figure 6b statistic, raw vs denoised -------------------------------
contrast <- function(col) {
  dd <- res[res$dose_mM > 0 & res$dose_mM != DROP_DOSE, ]
  dd$ld <- log10(dd$dose_mM)
  sl <- sapply(split(dd, list(dd$temp, dd$culture), drop = TRUE),
               function(g) c(temp = g$temp[1], slope = unname(coef(lm(g[[col]] ~ g$ld))[2])))
  s15 <- sl["slope", sl["temp", ] == 15]; s27 <- sl["slope", sl["temp", ] == 27]
  tt <- t.test(s27, s15, var.equal = FALSE)
  cvc <- sapply(c(15, 27), function(tp) { x <- res[[col]][res$temp == tp & res$dose_mM == 0]; 100*sd(x)/mean(x) })
  data.frame(rates = col, slope15 = mean(s15), slope27 = mean(s27), diff = mean(s27) - mean(s15),
             welch_P = tt$p.value, n27_negative = sum(s27 < 0), ctrl_CV15 = cvc[1], ctrl_CV27 = cvc[2])
}
summ <- do.call(rbind, lapply(c("r_raw", "r_dn_2h", "r_dn_4h", "r_dn_12h"), contrast))
nz <- aggregate(cbind(resid_sd, period_h, late_p2p) ~ temp, noise, median)
write.csv(summ, file.path(ROOT, "tables/aox/sham_denoise_summary.csv"), row.names = FALSE)

cat("Residual noise after 6 h trend removal (plate medians):\n"); print(nz, row.names = FALSE)
ok <- res$dose_mM != DROP_DOSE      # 0.35 mM: non-identifiable at 27 C, excluded in the paper
cat("\nPer-curve change in r after denoising (curves entering Figure 6, 12 h average):\n")
print(summary(res$pct_change_12h[ok]))
cat(sprintf("  |change| > 1%%: %d of %d curves;  > 5%%: %d\n",
            sum(abs(res$pct_change_12h[ok]) > 1), sum(ok), sum(abs(res$pct_change_12h[ok]) > 5)))
cat("\nFigure 6b statistic:\n"); print(summ, row.names = FALSE, digits = 4)
