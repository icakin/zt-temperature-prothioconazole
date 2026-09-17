#!/usr/bin/env Rscript
# =============================================================================
# 45_ptc_sham_prepare.R -- prothioconazole x SHAM factorial (Figure 7):
# raw PreSens exports -> pipeline tables.
#
# Design (experiments/ptc_sham_mediation_design.pdf): three independent
# cultures (replicates 1-3), each run as one 15 C and one 27 C SensorDish plate
# from the same cups on the same day. Per plate, six conditions at fixed
# positions:
#   row A  V   vehicle (5 wells) + bV cell-free blank (A6)
#   row B  P2  prothioconazole 2 mg/L (B1-B3), P4 4 mg/L (B4-B6)
#   row C  S   SHAM 0.2 mM (5 wells) + bS cell-free blank with SHAM (C6)
#   row D  P2S 2 mg/L + SHAM (D1-D3), P4S 4 mg/L + SHAM (D4-D6)
# 0.4% DMSO and 0.2% ethanol in every vial; 5.05 mL per vial.
#
# Inputs : data/oxygen/ptc_sham_rep{1,2,3}_{15,27}_Oxygen.xlsx  PreSens exports
# Outputs: data/aox/ptc_sham_{15,27}_Oxygen.csv
#            wide traces in the format of the other oxygen files: Time (min), T,
#            then one column per well named <condition>_R<rep><well>
#            (e.g. P2S_R2D1); blanks as bV_/bS_.
#          tables/aox/ptc_sham/Oxygen_All_Long.csv
#          tables/aox/ptc_sham/Oxygen_Trimmed_Series_Metadata.csv
#            what 03_trim_selector.R needs to show these curves (readings with
#            only the ~41 h handling-step offset removed; denoising is chosen in
#            the selector and saved with the windows):
#              Rscript scripts/03_trim_selector.R tables/aox/ptc_sham
#            It writes tables/aox/ptc_sham/manual_fit_windows.csv, which
#            46_ptc_sham_rates.R uses when present.
#          tables/aox/ptc_sham_fit_windows_auto.csv
#            rule-based windows (post-equilibration peak -> before O2 falls
#            under 2.5 mg/L, and before the 41 h handling step on the 27 C
#            plates); the reference lines in the selector and the fallback
#            for 46 when no manual windows exist.
#          tables/aox/ptc_sham_steps.csv
#            the plates were taken out of the incubator for an interim export
#            at ~41 h; the reading jumps when they go back. Per-well time and
#            size of that step, for the offset correction in 46.
# Run from the repository root:  Rscript scripts/45_ptc_sham_prepare.R
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))

PS_COND  <- c(V = "vehicle", S = "SHAM 0.2 mM", P2 = "prothioconazole 2 mg/L",
              P4 = "prothioconazole 4 mg/L", P2S = "2 mg/L + SHAM", P4S = "4 mg/L + SHAM",
              bV = "blank", bS = "blank + SHAM")
ps_layout <- function() {
  w <- as.vector(outer(LETTERS[1:4], 1:6, paste0))
  cond <- sapply(w, function(x) {
    r <- substr(x, 1, 1); c <- as.integer(substr(x, 2, 2))
    switch(r, A = if (c <= 5) "V" else "bV", C = if (c <= 5) "S" else "bS",
              B = if (c <= 3) "P2" else "P4", D = if (c <= 3) "P2S" else "P4S")
  })
  data.frame(well = w, condition = unname(cond), stringsAsFactors = FALSE)
}
STEP_LO <- 40; STEP_HI <- 43; O2_FLOOR <- 2.5; PRE_STEP_END_27 <- 40.5

lay <- ps_layout()
dir.create(file.path(ROOT, "data/aox"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "tables/aox/ptc_sham"), showWarnings = FALSE, recursive = TRUE)
write.csv(lay, file.path(ROOT, "data/aox/ptc_sham_layout.csv"), row.names = FALSE)

long_all <- list(); meta_all <- list(); win_all <- list(); step_all <- list()
for (temp in c(15, 27)) {
  wide <- NULL
  for (rep in 1:3) {
    f <- file.path(ROOT, sprintf("data/oxygen/ptc_sham_rep%d_%d_Oxygen.xlsx", rep, temp))
    tr <- read_sdr_xlsx(f, plate = rep)
    tr$condition <- lay$condition[match(tr$well, lay$well)]
    tr$Replicate <- sprintf("R%d%s", rep, tr$well)
    tr$Time <- round(tr$t_h * 60, 2)
    # wide block for this plate
    wb <- reshape(tr[, c("Time", "condition", "Replicate", "o2")],
                  idvar = "Time", timevar = "Replicate", direction = "wide",
                  drop = "condition")
    names(wb) <- sub("^o2\\.", "", names(wb))
    cn <- setdiff(names(wb), "Time")
    names(wb)[match(cn, names(wb))] <- paste0(tr$condition[match(cn, tr$Replicate)], "_", cn)
    wide <- if (is.null(wide)) wb else merge(wide, wb, by = "Time", all = TRUE)

    # per-well: long rows, auto window, handling step
    for (w in unique(tr$well)) {
      g <- tr[tr$well == w & is.finite(tr$o2), ]; cond <- g$condition[1]
      t_h <- g$t_h; o2 <- g$o2
      # handling step: largest one-sample rise in 40-43 h, size = median 1 h after - median 1 h before
      m <- t_h > STEP_LO & t_h < STEP_HI
      j <- which(m)[which.max(diff(o2[m]))]
      st_t <- t_h[j]
      st_dy <- median(o2[t_h > st_t & t_h < st_t + 1]) - median(o2[t_h > st_t - 1 & t_h <= st_t])
      step_all[[length(step_all) + 1]] <- data.frame(
        T = temp, replicate = rep, well = w, condition = cond,
        step_time_h = round(st_t, 2), step_mgL = round(st_dy, 3))
      if (substr(cond, 1, 1) == "b") next          # blanks are not fitted
      # selector table: raw readings with only the handling-step offset removed
      # (so a moving average chosen in the selector does not smear the jump)
      o2c <- step_correct(t_h * 60, o2, st_t, st_dy)
      long_all[[length(long_all) + 1]] <- data.frame(
        File = basename(f), Time = round(t_h * 60, 2), T = temp, Dose = cond,
        Replicate = sprintf("R%d%s", rep, w), Oxygen = o2c)
      # start: post-equilibration maximum of a 1 h moving average, searched after
      # the first 2.5 h (the trace starts high and falls while the vial equilibrates)
      sm1 <- as.numeric(stats::filter(o2, rep(1/51, 51), sides = 2)); sm1[is.na(sm1)] <- -Inf
      cand <- which(t_h > 2.5 & t_h < max(t_h)/3)
      pk <- cand[which.max(sm1[cand])]; t0 <- t_h[pk]
      sm <- as.numeric(stats::filter(o2, rep(1/201, 201), sides = 2))
      below <- which(sm < O2_FLOOR & t_h > t0)
      t1 <- if (length(below)) t_h[below[1]] else max(t_h)
      if (temp == 27) t1 <- min(t1, PRE_STEP_END_27)
      meta_all[[length(meta_all) + 1]] <- data.frame(
        T = temp, Dose = cond, Replicate = sprintf("R%d%s", rep, w),
        main_run_start_time = round(t0 * 60, 1), steepest_drop_time = round(t1 * 60, 1))
    }
  }
  wide <- wide[order(wide$Time), ]
  wide <- cbind(Time = wide$Time, T = temp, wide[, setdiff(names(wide), "Time"), drop = FALSE])
  write.csv(wide, file.path(ROOT, sprintf("data/aox/ptc_sham_%d_Oxygen.csv", temp)), row.names = FALSE, na = "")
}
long <- do.call(rbind, long_all); meta <- do.call(rbind, meta_all); steps <- do.call(rbind, step_all)
win_auto <- data.frame(T = meta$T, Dose = meta$Dose, Replicate = meta$Replicate,
                       fit_start = meta$main_run_start_time, fit_end = meta$steepest_drop_time)
write.csv(long,  file.path(ROOT, "tables/aox/ptc_sham/Oxygen_All_Long.csv"), row.names = FALSE)
write.csv(meta,  file.path(ROOT, "tables/aox/ptc_sham/Oxygen_Trimmed_Series_Metadata.csv"), row.names = FALSE)
write.csv(win_auto, file.path(ROOT, "tables/aox/ptc_sham_fit_windows_auto.csv"), row.names = FALSE)
write.csv(steps, file.path(ROOT, "tables/aox/ptc_sham_steps.csv"), row.names = FALSE)

cat(sprintf("ptc_sham: %d curves (%d per temperature) written; blanks: %d wells\n",
            length(unique(paste(long$T, long$Replicate))), length(unique(long$Replicate[long$T == 15])),
            sum(substr(steps$condition, 1, 1) == "b")))
cat("handling step at ~41 h, median size (mg/L) by plate:\n")
print(aggregate(step_mgL ~ T + replicate, steps[substr(steps$condition, 1, 1) != "b", ], median), row.names = FALSE)
cat("Next: Rscript scripts/03_trim_selector.R tables/aox/ptc_sham   (set the windows by hand)\n",
    "      Rscript scripts/46_ptc_sham_rates.R\n")
