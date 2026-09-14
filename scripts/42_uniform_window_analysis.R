#!/usr/bin/env Rscript
# =============================================================================
# 42_uniform_window_analysis.R
#
# Pre-specified, hand-free reanalysis of BOTH AOX-inhibitor experiments.
#
# Motivation: in 41_sham_rates.R each fit window is chosen by hand in the trim
# selector. That is defensible but not auditable, and for nPG the fitted r is
# strongly associated with window length (27 C) or with total oxygen drawdown
# (15 C). This script removes analyst choice entirely: every window is set by a
# rule applied identically to every curve of both compounds, and every
# reliability filter is a property of the trace, never of the result.
#
# Rules (all automated):
#   fixed_<a>_<b>      absolute time window, minutes
#   draw_<lo>_<hi>     between lo% and hi% of that curve's total drawdown
#
# Reliability filter: minimum oxygen consumed across the window. Curves with
# little drawdown are close to linear, so r and K are not separately
# identifiable. Threshold is swept, not chosen.
#
# Output: tables/aox/uniform_window_analysis.csv
# =============================================================================
suppressPackageStartupMessages(library(minpack.lm))
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."

resp_model <- function(t, r, K, O2_0) O2_0 + (K/r) * (1 - exp(r * t))

fit_one <- function(tt, yy) {
  if (length(yy) < 30) return(NULL)
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
       used = max(yy) - min(yy), n_pts = length(yy),
       t0 = min(tt), t1 = max(tt))
}

RULES <- list(
  list(id = "fixed 400-3200",  type = "fixed", a = 400,  b = 3200),
  list(id = "fixed 600-3200",  type = "fixed", a = 600,  b = 3200),
  list(id = "fixed 300-2500",  type = "fixed", a = 300,  b = 2500),
  list(id = "fixed 800-3600",  type = "fixed", a = 800,  b = 3600),
  list(id = "drawdown 5-85%",  type = "draw",  a = 0.05, b = 0.85),
  list(id = "drawdown 10-90%", type = "draw",  a = 0.10, b = 0.90),
  list(id = "drawdown 5-95%",  type = "draw",  a = 0.05, b = 0.95)
)

fit_all <- function(inh, rule) {
  out <- list()
  for (temp in c(15, 27)) {
    f <- file.path(ROOT, sprintf("data/aox/%s_%d_Oxygen.csv", inh, temp))
    if (!file.exists(f)) next
    d <- read.csv(f, check.names = FALSE)
    for (cv in setdiff(names(d), c("Time", "T"))) {
      y <- d[[cv]]; tt <- d$Time; ok <- is.finite(y) & is.finite(tt)
      y <- y[ok]; tt <- tt[ok]
      if (rule$type == "fixed") {
        m <- tt >= rule$a & tt <= rule$b
      } else {
        pk <- which.max(y); y2 <- y[pk:length(y)]; t2 <- tt[pk:length(tt)]
        span <- max(y2) - min(y2)
        if (span <= 0) next
        frac <- (max(y2) - y2) / span
        m <- rep(FALSE, length(y)); m[pk:length(y)] <- frac >= rule$a & frac <= rule$b
      }
      if (sum(m) < 30) next
      fit <- fit_one(tt[m], y[m]); if (is.null(fit)) next
      dose <- sub("_R[0-9]+$", "", cv); rep <- as.integer(sub("^.*_R", "", cv))
      out[[length(out) + 1]] <- data.frame(
        inhibitor = inh, rule = rule$id, temp = temp,
        dose_mM = if (dose == "Control") 0 else as.numeric(dose),
        culture = rep, r_per_h = fit$r_per_h, K = fit$K, R2 = fit$R2,
        O2_used = fit$used, n_pts = fit$n_pts, win_start = fit$t0, win_end = fit$t1)
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

all_fits <- do.call(rbind, lapply(c("sham", "npg"), function(i)
  do.call(rbind, lapply(RULES, function(r) fit_all(i, r)))))
write.csv(all_fits, file.path(ROOT, "tables/aox/uniform_window_fits.csv"), row.names = FALSE)

# ---- culture slopes and the 27-15 contrast, swept over drawdown thresholds ----
summarise <- function(df) {
  res <- list()
  for (thr in c(0, 3, 4, 5, 6)) {
    for (inh in unique(df$inhibitor)) for (rl in unique(df$rule)) {
      s <- df[df$inhibitor == inh & df$rule == rl & df$dose_mM > 0 &
              df$dose_mM != 0.35 & df$O2_used >= thr, ]
      if (!nrow(s)) next
      sl <- lapply(split(s, list(s$temp, s$culture), drop = TRUE), function(g)
        if (length(unique(g$dose_mM)) >= 3)
          data.frame(temp = g$temp[1], culture = g$culture[1],
                     slope = unname(coef(lm(r_per_h ~ log10(dose_mM), g))[2])) else NULL)
      sl <- do.call(rbind, sl); if (is.null(sl)) next
      s15 <- sl$slope[sl$temp == 15]; s27 <- sl$slope[sl$temp == 27]
      if (length(s15) < 2 || length(s27) < 2) next
      tt <- t.test(s27, s15, var.equal = FALSE)
      res[[length(res) + 1]] <- data.frame(
        inhibitor = inh, rule = rl, min_O2 = thr,
        n15 = length(s15), n27 = length(s27),
        n_curves = nrow(s),
        s15 = mean(s15), s27 = mean(s27),
        diff = mean(s27) - mean(s15), p = tt$p.value)
    }
  }
  do.call(rbind, res)
}
S <- summarise(all_fits)
S <- S[order(S$inhibitor, S$min_O2, S$rule), ]
write.csv(S, file.path(ROOT, "tables/aox/uniform_window_analysis.csv"), row.names = FALSE)

cat("\n=== 27C minus 15C slope difference, no hand-chosen windows ===\n")
for (inh in c("sham", "npg")) {
  cat(sprintf("\n--- %s ---\n", toupper(inh)))
  x <- S[S$inhibitor == inh, ]
  for (thr in unique(x$min_O2)) {
    y <- x[x$min_O2 == thr, ]
    cat(sprintf(" min O2 %d mg/L : %d/%d rules negative | diff %+.5f to %+.5f | P %.4f to %.4f | curves %d\n",
                thr, sum(y$diff < 0), nrow(y), min(y$diff), max(y$diff),
                min(y$p), max(y$p), round(mean(y$n_curves))))
  }
}
cat("\nwrote tables/aox/uniform_window_analysis.csv and uniform_window_fits.csv\n")
