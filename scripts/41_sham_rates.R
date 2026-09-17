#!/usr/bin/env Rscript
# =============================================================================
# 41_sham_rates.R -- growth rates for the AOX-inhibitor experiments.
#   Rscript scripts/41_sham_rates.R . sham     (default)
#   Rscript scripts/41_sham_rates.R . npg
#
# This is the step that turns raw dissolved-oxygen traces into the table Figure 6
# is drawn from. It uses the SAME model, bounds and fitting procedure as the rest
# of the paper (cf. 04_oxygen_fits.R):
#
#       O(t) = O2_0 + (K/r) * (1 - exp(r*t))
#
# fitted by bounded nonlinear least squares over an interval chosen by hand in
# the trim selector (03_trim_selector.R), one interval per curve. r is returned
# per minute by the fit and converted to h^-1 here.
#
# Inputs
#   data/aox/sham_{15,27}_Oxygen.csv   wide traces: Time (min), T, then one
#                                      column per curve named <Dose>_R<n>,
#                                      with Dose = "Control" or the mM value
#   tables/aox/sham_fit_windows.csv    T, Dose, Replicate, fit_start, fit_end
#                                      (minutes) -- exported from 03_trim_selector
# Output
#   tables/aox/sham_culture_rates.csv  temp, dose_mM, culture, r_per_h, K, R2,
#                                      RMSE, n_pts, fit_start, fit_end
#
# Run from the repository root:  Rscript scripts/41_sham_rates.R
# =============================================================================
suppressPackageStartupMessages(library(minpack.lm))
smooth_ma <- function(t_min, y, smooth_h = 0) {
  if (!is.finite(smooth_h) || smooth_h <= 0) return(y)
  dt <- median(diff(t_min), na.rm = TRUE); k <- max(3L, as.integer(round(smooth_h * 60 / dt)))
  if (k %% 2 == 0) k <- k + 1L
  ok <- is.finite(y); ys <- rep(NA_real_, length(y))
  ys[ok] <- as.numeric(stats::filter(y[ok], rep(1 / k, k), sides = 2))
  ys
}
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
INH  <- if (length(ARGS) > 1) ARGS[2] else "sham"      # "sham" or "npg"
# Interval source. If tables/aox/<inh>/manual_fit_windows.csv exists (written by
# 03_trim_selector.R after 40_aox_trim_inputs.R), it is used in preference to
# the older <inh>_fit_windows.csv, together with the per-curve denoising chosen
# in the selector. This is what puts this experiment through the same procedure
# as the prothioconazole x SHAM factorial (46_ptc_sham_prepared_traces.R).
MANW <- file.path(ROOT, sprintf("tables/aox/%s/manual_fit_windows.csv", INH))
WIN  <- if (file.exists(MANW)) MANW else file.path(ROOT, sprintf("tables/aox/%s_fit_windows.csv", INH))
if (!file.exists(WIN))
  stop("missing ", WIN, "\n  Set the intervals in 03_trim_selector.R for the SHAM run\n",
       "  and save its window table here before running this script.")

resp_model <- function(t, r, K, O2_0) O2_0 + (K/r) * (1 - exp(r * t))

fit_one <- function(tt, yy) {
  t0 <- tt - min(tt)
  sl <- median(diff(yy) / diff(t0))
  st <- list(r = 1e-3, K = min(max(abs(sl), 1e-6), 1), O2_0 = yy[1])
  f <- try(nlsLM(yy ~ resp_model(t0, r, K, O2_0), start = st,
                 lower = c(r = 1e-6, K = 1e-10, O2_0 = min(yy) - 1),
                 upper = c(r = 0.15,  K = 1,     O2_0 = max(yy) + 1),
                 control = nls.lm.control(maxiter = 200, ftol = 1e-12, ptol = 1e-12)),
           silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  p <- coef(f); pred <- resp_model(t0, p[["r"]], p[["K"]], p[["O2_0"]])
  rss <- sum((yy - pred)^2); sst <- sum((yy - mean(yy))^2)
  se <- tryCatch(unname(summary(f)$coefficients["r", "Std. Error"]), error = function(e) NA_real_)
  # identifiability diagnostic: how much the curved model buys over a straight
  # line on the same points. r is the CURVATURE of the oxygen decline, so a
  # trace that is essentially linear carries no growth signal and its r is not
  # determined by the data, however well the model happens to fit.
  rss_lin <- sum(residuals(lm(yy ~ t0))^2)
  list(r_per_h = unname(p[["r"]]) * 60, r_se_per_h = se * 60, K = unname(p[["K"]]),
       R2 = if (sst > 1e-12) 1 - rss/sst else NA_real_,
       r2_gain = if (sst > 1e-12) (rss_lin - rss)/sst else NA_real_,
       RMSE = sqrt(rss/length(yy)), n_pts = length(yy))
}

win <- read.csv(WIN, stringsAsFactors = FALSE)
out <- list()
for (temp in c(15, 27)) {
  f <- file.path(ROOT, sprintf("data/aox/%s_%d_Oxygen.csv", INH, temp))
  d <- read.csv(f, check.names = FALSE)
  curves <- setdiff(names(d), c("Time", "T"))
  for (cv in curves) {
    dose <- sub("_R[0-9]+$", "", cv); rep <- sub("^.*_R", "", cv)
    w <- win[win$T == temp & win$Dose == dose & sub("^R", "", win$Replicate) %in% c(rep, LETTERS[as.integer(rep)]), ]
    if (!nrow(w)) w <- win[win$T == temp & win$Dose == dose & win$Replicate == LETTERS[as.integer(rep)], ]
    if (!nrow(w)) { message("no window for ", temp, " ", cv, " - skipped"); next }
    sm <- if ("smooth_h" %in% names(w) && is.finite(w$smooth_h[1])) w$smooth_h[1] else 0
    yy <- smooth_ma(d$Time, d[[cv]], sm)
    m <- d$Time >= w$fit_start[1] & d$Time <= w$fit_end[1] & is.finite(yy)
    fit <- fit_one(d$Time[m], yy[m])
    if (is.null(fit)) { message("fit failed: ", temp, " ", cv); next }
    out[[length(out) + 1]] <- data.frame(
      temp = temp, dose_mM = if (dose == "Control") 0 else as.numeric(dose),
      culture = as.integer(rep), r_per_h = fit$r_per_h, r_se_per_h = fit$r_se_per_h,
      K = fit$K, R2 = fit$R2, r2_gain = fit$r2_gain,
      RMSE = fit$RMSE, n_pts = fit$n_pts, fit_start = w$fit_start[1], fit_end = w$fit_end[1],
      smooth_h = sm)
  }
}
res <- do.call(rbind, out)
res <- res[order(res$temp, res$dose_mM, res$culture), ]
dir.create(file.path(ROOT, "tables/aox"), showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(ROOT, sprintf("tables/aox/%s_culture_rates.csv", INH)), row.names = FALSE)

cat(sprintf("intervals from %s\n", basename(WIN)))
cat(sprintf("fitted %d curves | median R2 = %.4f | worst R2 = %.4f\n",
            nrow(res), median(res$R2, na.rm = TRUE), min(res$R2, na.rm = TRUE)))
cat(sprintf("curvature over a straight line (r2_gain): median %.4f | %d of %d curves below %.3f (r not identified)\n",
            median(res$r2_gain, na.rm = TRUE), sum(res$r2_gain < 0.005, na.rm = TRUE), nrow(res), 0.005))
for (tp in c(15, 27)) {
  s <- res[res$temp == tp & res$dose_mM > 0 & res$dose_mM < 0.35, ]   # 0.35 mM: not estimable at 27 C
  sl <- sapply(split(s, s$culture), function(g) coef(lm(r_per_h ~ log10(dose_mM), g))[2])
  cat(sprintf("%d C: per-culture slopes %s | mean %+.5f\n", tp,
              paste(sprintf("%+.5f", sl), collapse = " "), mean(sl)))
}
