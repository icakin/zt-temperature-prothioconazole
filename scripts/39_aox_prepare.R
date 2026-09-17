#!/usr/bin/env Rscript
# =============================================================================
# 39_aox_prepare.R -- AOX-inhibitor dose-response experiments (SHAM, n-propyl
# gallate): raw PreSens exports -> the wide oxygen tables the rest of the
# pipeline reads.
#
# Plate layout. Dose is assigned to plate column and replicate to plate row, so
# the dose axis and the column axis coincide on these plates:
#   SHAM  columns 1-6 = 0.0125, 0.025, 0.05, 0.1, 0.2, 0.35 mM; rows A-C are the
#         three replicate cultures; D1-D3 are the solvent controls.
#   nPG   column 1 = solvent control, columns 2-6 = 0.0125, 0.025, 0.05, 0.075,
#         0.1 mM; rows A-C are the three replicate cultures.
#
# Inputs : data/oxygen/sham_{15,27}_Oxygen.xlsx
#          data/oxygen/npg_{15,27}_Oxygen.xlsx
# Outputs: data/aox/sham_{15,27}_Oxygen.csv
#          data/aox/npg_{15,27}_Oxygen.csv
#            Time (min), T, then one column per curve named <Dose>_R<n>, with
#            Dose = "Control" or the concentration in mM.
#          data/aox/sham_layout.csv, data/aox/npg_layout.csv  well -> dose, replicate
# Run from the repository root:  Rscript scripts/39_aox_prepare.R
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))

layouts <- list(
  sham = list(doses = c(`1` = "0.0125", `2` = "0.025", `3` = "0.05",
                        `4` = "0.1",    `5` = "0.2",   `6` = "0.35"),
              reps  = c(A = 1, B = 2, C = 3),
              ctrl  = c(D1 = 1, D2 = 2, D3 = 3)),
  npg  = list(doses = c(`2` = "0.0125", `3` = "0.025", `4` = "0.05",
                        `5` = "0.075",  `6` = "0.1"),
              reps  = c(A = 1, B = 2, C = 3),
              ctrl  = c(A1 = 1, B1 = 2, C1 = 3)))

dir.create(file.path(ROOT, "data/aox"), showWarnings = FALSE, recursive = TRUE)
for (inh in names(layouts)) {
  L <- layouts[[inh]]
  lay <- do.call(rbind, c(
    lapply(names(L$ctrl), function(w) data.frame(well = w, dose = "Control", replicate = L$ctrl[[w]])),
    lapply(names(L$reps), function(r) do.call(rbind, lapply(names(L$doses), function(cl) {
      w <- paste0(r, cl)
      if (w %in% names(L$ctrl)) NULL
      else data.frame(well = w, dose = L$doses[[cl]], replicate = L$reps[[r]])
    })))))
  lay <- lay[!duplicated(lay$well), ]
  write.csv(lay, file.path(ROOT, sprintf("data/aox/%s_layout.csv", inh)), row.names = FALSE)

  for (temp in c(15, 27)) {
    f <- file.path(ROOT, sprintf("data/oxygen/%s_%d_Oxygen.xlsx", inh, temp))
    tr <- read_sdr_xlsx(f, plate = temp)
    tr <- tr[tr$well %in% lay$well, ]
    tr$curve <- sprintf("%s_R%d", lay$dose[match(tr$well, lay$well)],
                                  lay$replicate[match(tr$well, lay$well)])
    tr$Time <- round(tr$t_h * 60, 2)
    wide <- reshape(tr[, c("Time", "curve", "o2")], idvar = "Time",
                    timevar = "curve", direction = "wide")
    names(wide) <- sub("^o2\\.", "", names(wide))
    ord <- c("Control", L$doses)
    cn <- unlist(lapply(ord, function(d) sprintf("%s_R%d", d, 1:3)))
    wide <- wide[order(wide$Time), c("Time", intersect(cn, names(wide)))]
    wide <- cbind(Time = wide$Time, T = temp, wide[, -1, drop = FALSE])
    write.csv(wide, file.path(ROOT, sprintf("data/aox/%s_%d_Oxygen.csv", inh, temp)),
              row.names = FALSE, na = "")
    cat(sprintf("%s %d C: %d curves, %d readings\n", inh, temp, ncol(wide) - 2, nrow(wide)))
  }
}
