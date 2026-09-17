#!/usr/bin/env Rscript
# =============================================================================
# 46_ptc_sham_rates.R -- growth rates and the interaction contrast for the
# prothioconazole x SHAM factorial (Figure 7, Figure S18).
#
# Every curve is fitted with the manuscript oxygen model
#       O(t) = O2_0 + (K/r) * (1 - exp(r*t))
# by bounded nonlinear least squares over its fitting interval, exactly as in
# 04_oxygen_fits.R and 41_sham_rates.R. Intervals come from the trim selector
# (tables/aox/ptc_sham/manual_fit_windows.csv, written by
# 03_trim_selector.R tables/aox/ptc_sham); if that file is absent the
# rule-based intervals from 45 are used and the output says so.
#
# Handling step: the plates were taken out for an interim export at ~41 h and
# the reading jumps when they go back (tables/aox/ptc_sham_steps.csv). The
# post-step readings are shifted by the measured jump (an offset only; recorded
# per well where the interval spans it).
# Denoising: the moving-average width chosen in the trim selector and saved per
# curve as smooth_h in manual_fit_windows.csv is applied to that curve before fitting
# (smooth_ma in aox_common.R; r is unbiased by a centred moving average under
# this model). The fit of the same window on the undenoised trace is kept as
# r_per_h_raw.
#
# Statistics (prespecified in experiments/ptc_sham_mediation_design.pdf):
#   per plate, condition mean r over wells; relative growth r/r_V;
#   interaction contrast at dose d
#       I_d = log r_PS,d - log r_P,d - log r_S + log r_V
#   (0 = the two agents act independently on log growth rate; < 0 = SHAM
#   potentiates the drug). The temperature question is I_27 - I_15 within each
#   replicate. Primary: 2 mg/L; secondary: 4 mg/L and the two doses pooled.
#
# Inputs : data/aox/ptc_sham_{15,27}_Oxygen.csv, data/aox/ptc_sham_layout.csv,
#          tables/aox/ptc_sham/manual_fit_windows.csv (or ..._fit_windows_auto.csv),
#          tables/aox/ptc_sham/plot_exclude_points.csv (optional),
#          tables/aox/ptc_sham_steps.csv
# Outputs: tables/aox/ptc_sham_well_rates.csv        one row per fitted well
#          tables/aox/ptc_sham_condition_means.csv   plate x condition
#          tables/aox/ptc_sham_interaction.csv       I per plate and dose,
#                                                    I27 - I15 per replicate, tests
# Run from the repository root:  Rscript scripts/46_ptc_sham_rates.R
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))

MAN <- file.path(ROOT, "tables/aox/ptc_sham/manual_fit_windows.csv")
AUTO <- file.path(ROOT, "tables/aox/ptc_sham_fit_windows_auto.csv")
win_src <- if (file.exists(MAN)) "manual (trim selector)" else "rule-based (45_ptc_sham_prepare.R)"
win <- read.csv(if (file.exists(MAN)) MAN else AUTO, stringsAsFactors = FALSE)
win$Replicate <- toupper(win$Replicate)
# denoising width chosen in the trim selector, stored per curve (0 = raw);
# applied to each curve after the step offset
if (!"smooth_h" %in% names(win)) win$smooth_h <- 0
win$smooth_h[!is.finite(win$smooth_h)] <- 0
if (file.exists(MAN)) {                      # fall back per curve to the auto window where none was set
  auto <- read.csv(AUTO, stringsAsFactors = FALSE)
  key <- function(d) paste(d$T, d$Dose, toupper(d$Replicate))
  miss <- auto[!key(auto) %in% key(win), ]
  auto$smooth_h <- 0
  if (nrow(miss)) { win <- rbind(win[, names(auto)], miss[, names(auto)]); message(nrow(miss), " curves without a manual window: auto window used (raw)") }
}
EXC <- file.path(ROOT, "tables/aox/ptc_sham/plot_exclude_points.csv")
excl <- if (file.exists(EXC)) { e <- read.csv(EXC, stringsAsFactors = FALSE); paste(e$T, e$Dose, toupper(e$Replicate)) } else character(0)
steps <- read.csv(file.path(ROOT, "tables/aox/ptc_sham_steps.csv"), stringsAsFactors = FALSE)
steps$Replicate <- sprintf("R%d%s", steps$replicate, steps$well)

rows <- list()
for (temp in c(15, 27)) {
  d <- read.csv(file.path(ROOT, sprintf("data/aox/ptc_sham_%d_Oxygen.csv", temp)), check.names = FALSE)
  for (cv in setdiff(names(d), c("Time", "T"))) {
    cond <- sub("_R.*$", "", cv); repw <- sub("^.*_", "", cv)
    if (substr(cond, 1, 1) == "b") next
    k <- paste(temp, cond, repw)
    w <- win[win$T == temp & win$Dose == cond & win$Replicate == repw, ][1, ]
    if (!nrow(w) || is.na(w$fit_start)) { message("no window: ", k); next }
    tt <- d$Time
    st <- steps[steps$T == temp & steps$Replicate == repw, ][1, ]
    shifted <- if (isTRUE(w$fit_end > st$step_time_h * 60)) st$step_mgL else 0   # window spans the handling step
    yr <- step_correct(tt, d[[cv]], st$step_time_h, st$step_mgL)   # offset only
    SM <- w$smooth_h
    y  <- smooth_ma(tt, yr, SM)                                    # + denoising chosen in the selector for this curve
    m <- tt >= w$fit_start & tt <= w$fit_end & is.finite(y)
    f <- fit_o2_model(tt[m], y[m])
    if (is.null(f)) { message("fit failed: ", k); next }
    mr <- tt >= w$fit_start & tt <= w$fit_end & is.finite(yr)      # sensitivity: same window, raw
    fr <- if (SM > 0) fit_o2_model(tt[mr], yr[mr]) else f
    r_raw <- if (is.null(fr)) NA_real_ else fr$r_per_h
    rows[[length(rows) + 1]] <- data.frame(
      temp = temp, replicate = as.integer(substr(repw, 2, 2)), well = substr(repw, 3, 4), condition = cond,
      drug_mgL = c(V = 0, S = 0, P2 = 2, P2S = 2, P4 = 4, P4S = 4)[[cond]],
      sham_mM = if (cond %in% c("S", "P2S", "P4S")) 0.2 else 0,
      r_per_h = f$r_per_h, K = f$K, R2 = f$R2, RMSE = f$RMSE, n_pts = f$n_pts,
      r_per_h_raw = r_raw, smooth_h = SM,
      fit_start = w$fit_start, fit_end = w$fit_end, step_shift_mgL = shifted,
      excluded = k %in% excl, stringsAsFactors = FALSE)
  }
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

cat(sprintf("ptc_sham: %d wells fitted (%d excluded), windows: %s, mean R2 = %.4f, min R2 = %.4f\n",
            nrow(W), sum(W$excluded), win_src, mean(W$R2), min(W$R2)))
cat(sprintf("  step correction applied to %d wells (median |shift| %.2f mg/L)\n",
            sum(W$step_shift_mgL != 0), median(abs(W$step_shift_mgL[W$step_shift_mgL != 0]))))
sm_tab <- aggregate(smooth_h ~ temp, W, function(x) paste(unique(x), collapse = "/"))
cat("  denoising width (h) by temperature:", paste(sprintf("%d C: %s", sm_tab$temp, sm_tab$smooth_h), collapse = ", "), "\n")
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
cat("\ntests (n = 3 replicates; paired within replicate):\n"); print(S, row.names = FALSE, digits = 3)
