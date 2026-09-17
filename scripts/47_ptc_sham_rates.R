#!/usr/bin/env Rscript
# =============================================================================
# 47_ptc_sham_rates.R -- growth rates and the interaction contrast for the
# prothioconazole x SHAM factorial (Figure 7, Figure S18).
#
# Every curve is fitted with the manuscript oxygen model
#       O(t) = O2_0 + (K/r) * (1 - exp(r*t))
# by bounded nonlinear least squares, as in 04_oxygen_fits.R and 41_sham_rates.R.
#
# Fitting intervals are NOT chosen by hand. Each is derived from the trace by
# ps_window() in aox_common.R under the rule named in PRIMARY_RULE below, so the
# analyst never sees a rate while choosing an interval. The primary rule is
# peak + 36 h: each curve is fitted for 36 h from its own post-equilibration
# maximum. The start adapts to a culture that lags; the duration is the longest
# span the fastest wells can support, since the 27 C vehicle wells approach
# oxygen depletion around 40 h. Above all it is symmetric: every condition on a
# plate is fitted over the same span, so the vehicle and the slow drug + SHAM
# wells are not measured over different phases of growth. 48_ptc_sham_window_rules.R
# shows that rules which are not symmetric inflate the contrast (Fig. S19).
# 48_ptc_sham_window_rules.R refits everything under several alternative rules
# and reports the contrast under each (Fig. S19).
#
# Handling step: the plates were taken out for an interim export at ~41 h and
# the reading jumps when they go back (tables/aox/ptc_sham_steps.csv). The
# post-step readings are shifted by the measured jump (an offset only).
#
# Statistics (prespecified in experiments/ptc_sham_mediation_design.pdf):
#   per plate, condition mean r over wells; relative growth r/r_V;
#   interaction contrast at dose d
#       I_d = log r_PS,d - log r_P,d - log r_S + log r_V
#   (0 = the two agents act independently on log growth rate; < 0 = SHAM
#   potentiates the drug). The temperature question is I_27 - I_15 within each
#   replicate. Primary: 2 mg/L; secondary: 4 mg/L and the two doses pooled.
#
# Inputs : tables/aox/ptc_sham_prepared_traces.csv.gz and ptc_sham_windows.csv
#          (46_ptc_sham_prepared_traces.R); ptc_sham/plot_exclude_points.csv (optional)
# Outputs: tables/aox/ptc_sham_well_rates.csv        one row per fitted well
#          tables/aox/ptc_sham_condition_means.csv   plate x condition
#          tables/aox/ptc_sham_interaction.csv       I per plate and dose
#          tables/aox/ptc_sham_interaction_diff.csv  I27 - I15 per replicate
#          tables/aox/ptc_sham_interaction_tests.csv tests
# Run from the repository root:  Rscript scripts/47_ptc_sham_rates.R
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))

# The analysis-ready series and the fitting intervals come from
# 46_ptc_sham_prepared_traces.R, so the step correction, the denoising and the
# interval rule are defined in exactly one place and this script only fits.
PREP <- file.path(ROOT, "tables/aox/ptc_sham_prepared_traces.csv.gz")
WINS <- file.path(ROOT, "tables/aox/ptc_sham_windows.csv")
if (!file.exists(PREP) || !file.exists(WINS))
  stop("run scripts/46_ptc_sham_prepared_traces.R first")
L  <- read.csv(PREP, stringsAsFactors = FALSE)
Wn <- read.csv(WINS, stringsAsFactors = FALSE)

EXC <- file.path(ROOT, "tables/aox/ptc_sham/plot_exclude_points.csv")
excl <- if (file.exists(EXC)) { e <- read.csv(EXC, stringsAsFactors = FALSE); paste(e$T, e$Dose, toupper(e$Replicate)) } else character(0)

rows <- list()
for (i in seq_len(nrow(Wn))) {
  w <- Wn[i, ]
  g <- L[L$temp == w$temp & L$replicate == w$replicate & L$well == w$well & L$in_window, ]
  k <- paste(w$temp, w$condition, sprintf("R%d%s", w$replicate, w$well))
  f <- fit_o2_model(g$time_min, g$o2_prepared)                 # the analysis-ready series
  if (is.null(f)) { message("fit failed: ", k); next }
  fr <- if (w$smooth_h > 0) fit_o2_model(g$time_min, g$o2_step) else f   # same interval, undenoised
  rows[[length(rows) + 1]] <- data.frame(
    temp = w$temp, replicate = w$replicate, well = w$well, condition = w$condition,
    drug_mgL = w$drug_mgL, sham_mM = w$sham_mM,
    r_per_h = f$r_per_h, K = f$K, R2 = f$R2, RMSE = f$RMSE, n_pts = f$n_pts,
    r_per_h_raw = if (is.null(fr)) NA_real_ else fr$r_per_h, smooth_h = w$smooth_h,
    fit_start = w$fit_start_min, fit_end = w$fit_end_min, o2_used = w$o2_consumed,
    rule = w$rule, rule_fallback = w$rule_fallback, overridden = w$overridden,
    step_shift_mgL = w$step_mgL,
    excluded = k %in% excl, stringsAsFactors = FALSE)
}
W <- do.call(rbind, rows)
W <- W[order(W$temp, W$replicate, match(W$condition, PS_ORDER), W$well), ]
dir.create(file.path(ROOT, "tables/aox"), showWarnings = FALSE, recursive = TRUE)
write.csv(W, file.path(ROOT, "tables/aox/ptc_sham_well_rates.csv"), row.names = FALSE)

# --- plate x condition means -------------------------------------------------
U <- W[!W$excluded, ]
M <- aggregate(cbind(r_per_h, K) ~ temp + replicate + condition, U, mean)
n <- aggregate(r_per_h ~ temp + replicate + condition, U, length); names(n)[4] <- "n_wells"
M <- merge(M, n)
V <- M[M$condition == "V", c("temp", "replicate", "r_per_h", "K")]; names(V)[3:4] <- c("rV", "KV")
M <- merge(M, V); M$r_rel <- M$r_per_h / M$rV; M$K_rel <- M$K / M$KV
M <- M[order(M$temp, M$replicate, match(M$condition, PS_ORDER)), c("temp", "replicate", "condition", "n_wells", "r_per_h", "K", "r_rel", "K_rel")]
write.csv(M, file.path(ROOT, "tables/aox/ptc_sham_condition_means.csv"), row.names = FALSE)

# --- interaction contrasts ---------------------------------------------------
get <- function(tp, rp, cd) M$r_per_h[M$temp == tp & M$replicate == rp & M$condition == cd]
I <- do.call(rbind, lapply(1:3, function(rp) do.call(rbind, lapply(c(15, 27), function(tp) {
  lv <- log(get(tp, rp, "V")); ls <- log(get(tp, rp, "S"))
  data.frame(temp = tp, replicate = rp,
             I_2 = log(get(tp, rp, "P2S")) - log(get(tp, rp, "P2")) - ls + lv,
             I_4 = log(get(tp, rp, "P4S")) - log(get(tp, rp, "P4")) - ls + lv,
             pred_rel_P2S = get(tp, rp, "P2") * get(tp, rp, "S") / get(tp, rp, "V")^2,
             obs_rel_P2S  = get(tp, rp, "P2S") / get(tp, rp, "V"),
             pred_rel_P4S = get(tp, rp, "P4") * get(tp, rp, "S") / get(tp, rp, "V")^2,
             obs_rel_P4S  = get(tp, rp, "P4S") / get(tp, rp, "V"))
}))))
I$I_pooled <- (I$I_2 + I$I_4) / 2
D <- data.frame(replicate = 1:3,
  d_2 = I$I_2[I$temp == 27] - I$I_2[I$temp == 15],
  d_4 = I$I_4[I$temp == 27] - I$I_4[I$temp == 15])
D$d_pooled <- (D$d_2 + D$d_4) / 2
tst <- function(x) { t <- t.test(x); c(mean = mean(x), lo = t$conf.int[1], hi = t$conf.int[2], P = t$p.value) }
S <- rbind(
  data.frame(contrast = "I27 - I15, 2 mg/L (primary)", t(tst(D$d_2))),
  data.frame(contrast = "I27 - I15, 4 mg/L",            t(tst(D$d_4))),
  data.frame(contrast = "I27 - I15, pooled",            t(tst(D$d_pooled))),
  data.frame(contrast = "I15, 2 mg/L vs 0",             t(tst(I$I_2[I$temp == 15]))),
  data.frame(contrast = "I27, 2 mg/L vs 0",             t(tst(I$I_2[I$temp == 27]))),
  data.frame(contrast = "I15, 4 mg/L vs 0",             t(tst(I$I_4[I$temp == 15]))),
  data.frame(contrast = "I27, 4 mg/L vs 0",             t(tst(I$I_4[I$temp == 27]))))
S$n_same_sign <- c(sum(D$d_2 < 0), sum(D$d_4 < 0), sum(D$d_pooled < 0),
                   NA, sum(I$I_2[I$temp == 27] < 0), NA, sum(I$I_4[I$temp == 27] < 0))
out <- list(per_plate = I, per_replicate_difference = D, tests = S)
write.csv(I, file.path(ROOT, "tables/aox/ptc_sham_interaction.csv"), row.names = FALSE)
write.csv(D, file.path(ROOT, "tables/aox/ptc_sham_interaction_diff.csv"), row.names = FALSE)
write.csv(S, file.path(ROOT, "tables/aox/ptc_sham_interaction_tests.csv"), row.names = FALSE)

cat(sprintf("ptc_sham: %d wells fitted (%d excluded), intervals: '%s' from 48 (no hand trimming), mean R2 = %.4f, min R2 = %.4f\n",
            nrow(W), sum(W$excluded), W$rule[1], mean(W$R2), min(W$R2)))
cat(sprintf("  rule fell back to the plate bounds on %d of %d curves\n", sum(W$rule_fallback), nrow(W)))
cat(sprintf("  window length (h): 15 C median %.1f, 27 C median %.1f; O2 consumed median %.2f mg/L\n",
            median((W$fit_end - W$fit_start)[W$temp == 15])/60, median((W$fit_end - W$fit_start)[W$temp == 27])/60,
            median(W$o2_used)))
cat(sprintf("  step correction applied to %d wells (median |shift| %.2f mg/L)\n",
            sum(W$step_shift_mgL != 0), median(abs(W$step_shift_mgL[W$step_shift_mgL != 0]))))
sm_tab <- aggregate(smooth_h ~ temp, W, function(x) paste(unique(x), collapse = "/"))
cat("  denoising (h) by temperature:", paste(sprintf("%d C: %s", sm_tab$temp, sm_tab$smooth_h), collapse = ", "), "\n")
if (any(W$smooth_h > 0)) cat(sprintf("  raw-vs-denoised r on denoised curves: median change %.1f%%, max %.1f%%\n",
            100 * median(abs(W$r_per_h / W$r_per_h_raw - 1)[W$smooth_h > 0], na.rm = TRUE),
            100 * max(abs(W$r_per_h / W$r_per_h_raw - 1)[W$smooth_h > 0], na.rm = TRUE)))
cat("\nrelative growth rate r/r_V (plate means):\n")
rel <- reshape(M[, c("temp", "replicate", "condition", "r_rel")], idvar = c("temp", "replicate"),
               timevar = "condition", direction = "wide")
names(rel) <- sub("r_rel\\.", "", names(rel)); print(round(rel[, c("temp", "replicate", PS_ORDER)], 3), row.names = FALSE)
cat("\ninteraction contrast I (0 = independent, < 0 = SHAM potentiates the drug):\n")
print(round(I[, c("temp", "replicate", "I_2", "I_4", "obs_rel_P2S", "pred_rel_P2S", "obs_rel_P4S", "pred_rel_P4S")], 3), row.names = FALSE)
cat("\nI27 - I15 per replicate:\n"); print(round(D, 3), row.names = FALSE)
# --- with and without any manual interval overrides --------------------------
if (any(W$overridden)) {
  ru <- W; ru$r_per_h <- NA_real_
  for (i in which(W$overridden)) {                       # refit those wells on the rule interval
    w <- Wn[Wn$temp == W$temp[i] & Wn$replicate == W$replicate[i] & Wn$well == W$well[i], ][1, ]
    g <- L[L$temp == w$temp & L$replicate == w$replicate & L$well == w$well &
           L$time_min >= w$rule_start_min & L$time_min <= w$rule_end_min & is.finite(L$o2_prepared), ]
    f <- fit_o2_model(g$time_min, g$o2_prepared); if (!is.null(f)) ru$r_per_h[i] <- f$r_per_h
  }
  ru$r_per_h[!W$overridden] <- W$r_per_h[!W$overridden]
  Mr <- aggregate(r_per_h ~ temp + replicate + condition, ru[!ru$excluded, ], mean)
  gg <- function(tp, rp, cd) Mr$r_per_h[Mr$temp == tp & Mr$replicate == rp & Mr$condition == cd]
  II <- function(tp, rp, dd) log(gg(tp,rp,paste0(dd,"S"))) - log(gg(tp,rp,dd)) - log(gg(tp,rp,"S")) + log(gg(tp,rp,"V"))
  dr <- sapply(1:3, function(rp) II(27,rp,"P2") - II(15,rp,"P2"))
  cat(sprintf("\nmanual overrides on %d wells. I27 - I15 at 2 mg/L: %+.3f with them, %+.3f on the rule intervals alone\n",
              sum(W$overridden), mean(D$d_2), mean(dr)))
  for (i in which(W$overridden)) cat(sprintf("    %d C %-4s culture %d %s: %s\n",
      W$temp[i], W$condition[i], W$replicate[i], W$well[i],
      Wn$override_reason[Wn$temp == W$temp[i] & Wn$replicate == W$replicate[i] & Wn$well == W$well[i]][1]))
} else cat("\nmanual interval overrides: none; every interval came from the rule\n")

cat("\ntests (n = 3 replicates; paired within replicate):\n"); print(S, row.names = FALSE, digits = 3)
