#!/usr/bin/env Rscript
# =============================================================================
# 40_aox_trim_inputs.R -- put the SHAM dose-response experiment (Figure 6a,b)
# through the same interactive trim selector used for the prothioconazole x SHAM
# factorial, so both experiments are fitted by one procedure.
#
# Builds the two files 03_trim_selector.R reads, from the wide oxygen exports:
#   tables/aox/sham/Oxygen_All_Long.csv
#     File, Time (min), T, Dose, Replicate, Oxygen  -- one row per reading
#   tables/aox/sham/Oxygen_Trimmed_Series_Metadata.csv
#     T, Dose, Replicate, main_run_start_time, steepest_drop_time
#     seeded from tables/aox/sham_fit_windows.csv, the intervals already chosen
#     for Figure 6, so the selector opens on those rather than on nothing.
#
# These plates were not taken out mid-run, so there is no handling-step offset
# to remove; the readings go through unchanged.
#
# Inputs : data/aox/sham_{15,27}_Oxygen.csv, tables/aox/sham_fit_windows.csv
# Outputs: tables/aox/sham/Oxygen_All_Long.csv
#          tables/aox/sham/Oxygen_Trimmed_Series_Metadata.csv
# Run from the repository root:  Rscript scripts/40_aox_trim_inputs.R
# Then:   Rscript scripts/03_trim_selector.R tables/aox/sham
#         which writes tables/aox/sham/manual_fit_windows.csv
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
OUT <- file.path(ROOT, "tables/aox/sham")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

win <- read.csv(file.path(ROOT, "tables/aox/sham_fit_windows.csv"), stringsAsFactors = FALSE)
long <- list(); meta <- list()
for (temp in c(15, 27)) {
  f <- sprintf("data/aox/sham_%d_Oxygen.csv", temp)
  d <- read.csv(file.path(ROOT, f), check.names = FALSE)
  for (cv in setdiff(names(d), c("Time", "T"))) {
    dose <- sub("_R[0-9]+$", "", cv); rep <- sub("^.*_(R[0-9]+)$", "\\1", cv)
    ok <- is.finite(d[[cv]])
    long[[length(long) + 1]] <- data.frame(
      File = basename(f), Time = d$Time[ok], T = temp, Dose = dose,
      Replicate = rep, Oxygen = d[[cv]][ok], stringsAsFactors = FALSE)
    w <- win[win$T == temp & win$Dose == dose & win$Replicate == rep, ][1, ]
    meta[[length(meta) + 1]] <- data.frame(
      T = temp, Dose = dose, Replicate = rep,
      main_run_start_time = if (nrow(w) && is.finite(w$fit_start)) w$fit_start else min(d$Time[ok]),
      steepest_drop_time  = if (nrow(w) && is.finite(w$fit_end))   w$fit_end   else max(d$Time[ok]),
      stringsAsFactors = FALSE)
  }
}
L <- do.call(rbind, long); M <- do.call(rbind, meta)
write.csv(L, file.path(OUT, "Oxygen_All_Long.csv"), row.names = FALSE)
write.csv(M, file.path(OUT, "Oxygen_Trimmed_Series_Metadata.csv"), row.names = FALSE)
cat(sprintf("SHAM experiment: %d curves (%d per temperature), %d readings\n",
            nrow(M), sum(M$T == 15), nrow(L)))
cat(sprintf("  seeded from sham_fit_windows.csv on %d of %d curves\n",
            sum(!is.na(match(paste(M$T, M$Dose, M$Replicate), paste(win$T, win$Dose, win$Replicate)))), nrow(M)))
cat("Next: Rscript scripts/03_trim_selector.R tables/aox/sham\n")
