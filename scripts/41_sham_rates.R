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
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
INH  <- if (length(ARGS) > 1) ARGS[2] else "sham"      # "sham" or "npg"
WIN <- file.path(ROOT, sprintf("tables/aox/%s_fit_windows.csv", INH))
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
  list(r_per_h = unname(p[["r"]]) * 60, K = unname(p[["K"]]),
       R2 = if (sst > 1e-12) 1 - rss/sst else NA_real_,
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
    m <- d$Time >= w$fit_start[1] & d$Time <= w$fit_end[1] & is.finite(d[[cv]])
    fit <- fit_one(d$Time[m], d[[cv]][m])
    if (is.null(fit)) { message("fit failed: ", temp, " ", cv); next }
    out[[length(out) + 1]] <- data.frame(
      temp = temp, dose_mM = if (dose == "Control") 0 else as.numeric(dose),
      culture = as.integer(rep), r_per_h = fit$r_per_h, K = fit$K, R2 = fit$R2,
      RMSE = fit$RMSE, n_pts = fit$n_pts, fit_start = w$fit_start[1], fit_end = w$fit_end[1])
  }
}
res <- do.call(rbind, out)
res <- res[order(res$temp, res$dose_mM, res$culture), ]
dir.create(file.path(ROOT, "tables/aox"), showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(ROOT, sprintf("tables/aox/%s_culture_rates.csv", INH)), row.names = FALSE)

cat(sprintf("fitted %d curves | median R2 = %.4f | worst R2 = %.4f\n",
            nrow(res), median(res$R2, na.rm = TRUE), min(res$R2, na.rm = TRUE)))
for (tp in c(15, 27)) {
  s <- res[res$temp == tp & res$dose_mM > 0 & res$dose_mM < 0.35, ]   # 0.35 mM: not estimable at 27 C
  sl <- sapply(split(s, s$culture), function(g) coef(lm(r_per_h ~ log10(dose_mM), g))[2])
  cat(sprintf("%d C: per-culture slopes %s | mean %+.5f\n", tp,
              paste(sprintf("%+.5f", sl), collapse = " "), mean(sl)))
}
