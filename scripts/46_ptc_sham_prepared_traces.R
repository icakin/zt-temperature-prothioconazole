#!/usr/bin/env Rscript
# =============================================================================
# 46_ptc_sham_prepared_traces.R -- turn the raw oxygen readings of the
# prothioconazole x SHAM factorial into the analysis-ready series, i.e. exactly
# what the growth model is fitted to in 47_ptc_sham_rates.R. Nothing here looks
# at a growth rate, so no step can be tuned to the result.
#
# Four steps, in order, applied identically to all 132 cell-containing wells:
#
#   1. READ        data/aox/ptc_sham_{15,27}_Oxygen.csv (written by 45 from the
#                  PreSens exports), one column per well, named
#                  <condition>_R<culture><well>; cell-free blanks (bV, bS) are
#                  carried through unfitted for reference.
#   2. STEP        the plates were removed for an interim export at ~41 h and
#                  the reading jumps on re-insertion. Readings after the jump
#                  are shifted by its measured size (tables/aox/ptc_sham_steps.csv).
#                  An offset only: no point is deleted and no slope is altered.
#   3. DENOISE     a 4 h centred moving average on the 27 C plates, none at 15 C.
#                  After those plates were moved to a second incubator their
#                  cell-containing wells (but not the cell-free blanks beside
#                  them) show a 3-4 h re-aeration sawtooth, so the artefact is
#                  specific to those plates. For O(t) = O2_0 + (K/r)(1 - e^rt) a
#                  centred moving average of width W leaves r unchanged and
#                  multiplies K by sinh(rW/2)/(rW/2), 1.002 at r = 0.06 h^-1 and
#                  W = 4 h, so the growth rate is unbiased by construction.
#   4. TRIM        the fitting interval from ps_window() under PRIMARY_RULE:
#                  36 h from each well's own post-equilibration maximum. The
#                  start adapts to a culture that lags; 36 h is the longest span
#                  the fastest wells support, the 27 C vehicle wells approaching
#                  depletion near 40 h; and every condition on a plate gets the
#                  same span, so the vehicle and the slow drug + SHAM wells are
#                  never measured over different phases of growth. Rules without
#                  that symmetry inflate the contrast: see 47 and Fig. S19.
#
# Any interval set by hand in the review app (49_ptc_sham_review.R) is applied
# on top of step 4, flagged in the output and reported separately by 46.
#
# Inputs : data/aox/ptc_sham_{15,27}_Oxygen.csv, tables/aox/ptc_sham_steps.csv,
#          tables/aox/ptc_sham/manual_overrides.csv (optional)
# Outputs: tables/aox/ptc_sham_windows.csv          one row per well: the interval
#          tables/aox/ptc_sham_prepared_traces.csv.gz  long, analysis-ready:
#              temp, replicate, condition, well, drug_mgL, sham_mM, time_min,
#              o2_raw, o2_step, o2_prepared, in_window
#          The rows with in_window = TRUE and the o2_prepared column are what
#          47_ptc_sham_rates.R fits; o2_step is the same series without
#          denoising, used there as the sensitivity check.
# Run from the repository root:  Rscript scripts/46_ptc_sham_prepared_traces.R
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))

PRIMARY_RULE   <- list(type = "peak", dur = 36 * 60)   # minutes
SMOOTH_BY_TEMP <- c(`15` = 0, `27` = 4)                # hours
# WINDOW_SOURCE picks where the fitting intervals come from.
#   "manual"  the intervals set by hand in 03_trim_selector.R
#             (tables/aox/ptc_sham/manual_fit_windows.csv), with the per-curve
#             denoising chosen in the same tool. This matches how the fit
#             windows were set for the dose-response experiment (Figs 1-2).
#   "rule"    PRIMARY_RULE, chosen without seeing any rate.
# The rule interval is recorded as rule_start_min / rule_end_min either way, so
# 48_ptc_sham_window_rules.R can always compare the two.
WINDOW_SOURCE  <- "manual"
DRUG <- c(V = 0, S = 0, P2 = 2, P2S = 2, P4 = 4, P4S = 4, bV = 0, bS = 0)
SHAM <- c(V = 0, S = 0.2, P2 = 0, P2S = 0.2, P4 = 0, P4S = 0.2, bV = 0, bS = 0.2)

steps <- read.csv(file.path(ROOT, "tables/aox/ptc_sham_steps.csv"), stringsAsFactors = FALSE)
steps$Replicate <- sprintf("R%d%s", steps$replicate, steps$well)

# Manual interval overrides recorded in the review app (49_ptc_sham_review.R).
# Each carries a reason, is flagged in the output, and 47_ptc_sham_rates.R
# reports the contrast both with and without them, so an override can never
# change the result unnoticed.
OVR <- file.path(ROOT, "tables/aox/ptc_sham/manual_overrides.csv")
ovr <- if (file.exists(OVR)) read.csv(OVR, stringsAsFactors = FALSE) else NULL
if (!is.null(ovr) && nrow(ovr)) {
  ovr$key <- sprintf("%s_%s_%s", ovr$T, ovr$Dose, toupper(ovr$Replicate))
  if (!"reason" %in% names(ovr)) ovr$reason <- ""
} else ovr <- NULL

MANF <- file.path(ROOT, "tables/aox/ptc_sham/manual_fit_windows.csv")
man  <- if (file.exists(MANF)) read.csv(MANF, stringsAsFactors = FALSE) else NULL
if (!is.null(man)) man$key <- sprintf("%s_%s_%s", man$T, man$Dose, toupper(man$Replicate))
if (WINDOW_SOURCE == "manual" && is.null(man))
  stop("WINDOW_SOURCE is 'manual' but ", MANF, " is missing; run 03_trim_selector.R first")
if (WINDOW_SOURCE != "manual") man <- NULL

long <- list(); wins <- list()
for (temp in c(15, 27)) {
  d  <- read.csv(file.path(ROOT, sprintf("data/aox/ptc_sham_%d_Oxygen.csv", temp)), check.names = FALSE)
  SM <- unname(SMOOTH_BY_TEMP[as.character(temp)])
  tt <- d$Time
  for (cv in setdiff(names(d), c("Time", "T"))) {
    cond <- sub("_R.*$", "", cv); repw <- sub("^.*_", "", cv)
    rep  <- as.integer(substr(repw, 2, 2)); well <- substr(repw, 3, 4)
    st   <- steps[steps$T == temp & steps$Replicate == repw, ][1, ]

    mrow    <- if (!is.null(man)) man[man$key == sprintf("%d_%s_%s", temp, cond, repw), ] else NULL
    has_man <- !is.null(mrow) && nrow(mrow) == 1 &&
               is.finite(mrow$fit_start[1]) && is.finite(mrow$fit_end[1])
    SMc <- if (has_man && "smooth_h" %in% names(mrow) && is.finite(mrow$smooth_h[1])) mrow$smooth_h[1] else SM

    o2_raw  <- d[[cv]]
    o2_step <- step_correct(tt, o2_raw, st$step_time_h, st$step_mgL)   # 2
    o2_prep <- smooth_ma(tt, o2_step, SMc)                             # 3

    if (substr(cond, 1, 1) == "b") {                                   # blanks: not fitted
      inwin <- rep(FALSE, length(tt)); w <- c(NA_real_, NA_real_, NA)
    } else {
      w <- ps_window(tt, o2_prep, PRIMARY_RULE)                        # 4
      rule_start <- w[1]; rule_end <- w[2]; overridden <- FALSE; why <- ""
      if (has_man) { w[1] <- mrow$fit_start[1]; w[2] <- mrow$fit_end[1]; w[3] <- 0 }
      if (!is.null(ovr)) {
        o <- ovr[ovr$key == sprintf("%d_%s_%s", temp, cond, repw), ][1, ]
        if (nrow(o) && is.finite(o$fit_start) && is.finite(o$fit_end)) {
          w[1] <- o$fit_start; w[2] <- o$fit_end; overridden <- TRUE; why <- o$reason
        }
      }
      inwin <- is.finite(o2_prep) & tt >= w[1] & tt <= w[2]
      wins[[length(wins) + 1]] <- data.frame(
        temp = temp, replicate = rep, condition = cond, well = well,
        drug_mgL = DRUG[[cond]], sham_mM = SHAM[[cond]],
        rule = if (has_man) "manual (03_trim_selector.R)" else
                 sprintf("%s %g h", PRIMARY_RULE$type, PRIMARY_RULE$dur / 60),
        smooth_h = SMc, step_time_h = st$step_time_h, step_mgL = st$step_mgL,
        fit_start_min = w[1], fit_end_min = w[2], window_h = (w[2] - w[1]) / 60,
        rule_start_min = rule_start, rule_end_min = rule_end,
        overridden = overridden, override_reason = why,
        rule_fallback = as.logical(w[3]), n_points = sum(inwin),
        o2_at_start = o2_prep[which(inwin)[1]],
        o2_consumed = if (any(inwin)) max(o2_prep[inwin]) - min(o2_prep[inwin]) else NA_real_,
        stringsAsFactors = FALSE)
    }
    long[[length(long) + 1]] <- data.frame(
      temp = temp, replicate = rep, condition = cond, well = well,
      drug_mgL = DRUG[[cond]], sham_mM = SHAM[[cond]], time_min = tt,
      o2_raw = o2_raw, o2_step = o2_step, o2_prepared = o2_prep, in_window = inwin,
      stringsAsFactors = FALSE)
  }
}
L <- do.call(rbind, long); Wn <- do.call(rbind, wins)
Wn <- Wn[order(Wn$temp, Wn$replicate, match(Wn$condition, PS_ORDER), Wn$well), ]

dir.create(file.path(ROOT, "tables/aox"), showWarnings = FALSE, recursive = TRUE)
write.csv(Wn, file.path(ROOT, "tables/aox/ptc_sham_windows.csv"), row.names = FALSE)
gz <- gzfile(file.path(ROOT, "tables/aox/ptc_sham_prepared_traces.csv.gz"), "w")
write.csv(L, gz, row.names = FALSE); close(gz)

cat(sprintf("prepared %d wells (%d fitted, %d cell-free blanks) from %d readings\n",
            length(unique(paste(L$temp, L$replicate, L$well))), nrow(Wn),
            length(unique(paste(L$temp, L$replicate, L$well))) - nrow(Wn), nrow(L)))
cat(sprintf("  step offset applied to %d of %d fitted wells (median |shift| %.2f mg/L)\n",
            sum(abs(Wn$step_mgL) > 0.01), nrow(Wn), median(abs(Wn$step_mgL))))
cat(sprintf("  window source: %s\n", WINDOW_SOURCE))
cat(sprintf("  denoising: %s\n", if (WINDOW_SOURCE == "manual")
      paste(sprintf("%g h on %d wells", sort(unique(Wn$smooth_h)),
                    as.integer(table(Wn$smooth_h)[order(unique(sort(Wn$smooth_h)))])), collapse = ", ")
    else paste(sprintf("%g h at %s C", SMOOTH_BY_TEMP, names(SMOOTH_BY_TEMP)), collapse = ", ")))
cat(sprintf("  intervals from %s: median window %.1f h (rule '%s %g h' recorded alongside)\n",
            if (WINDOW_SOURCE == "manual") "03_trim_selector.R" else "the rule",
            median(Wn$window_h), PRIMARY_RULE$type, PRIMARY_RULE$dur / 60))
if (any(Wn$overridden)) {
  cat(sprintf("  MANUAL OVERRIDES on %d of %d wells (46 reports the contrast with and without them):\n", sum(Wn$overridden), nrow(Wn)))
  o <- Wn[Wn$overridden, ]
  for (i in seq_len(nrow(o))) cat(sprintf("    %d C %-4s culture %d %s: %.0f-%.0f min (rule %.0f-%.0f) - %s\n",
      o$temp[i], o$condition[i], o$replicate[i], o$well[i], o$fit_start_min[i], o$fit_end_min[i],
      o$rule_start_min[i], o$rule_end_min[i], o$override_reason[i]))
} else cat("  manual overrides: none\n")
cat("\nper plate and condition: start (h), O2 consumed over the interval (mg/L), points\n")
s <- aggregate(cbind(fit_start_min / 60, o2_consumed, n_points) ~ temp + replicate + condition, Wn, mean)
names(s)[4:6] <- c("start_h", "o2_used", "n_pts")
s <- s[order(s$temp, s$replicate, match(s$condition, PS_ORDER)), ]
print(s, row.names = FALSE, digits = 3)
