# =============================================================================
# aox_common.R -- shared helpers for the azide/SHAM respirometry (Figure 6c-d,
# Figure S17). Sourced by 43_azide_rates.R, 56_figure6.R and 73_figureS17.R.
#
#   read_sdr_xlsx(path)    PreSens SensorDish export (sheet "Oxygen Orginal")
#                          -> long data.frame: plate, well, t_h, o2
#   initial_rate(t_h, o2)  linear O2-consumption rate over the initial window:
#                          from the post-equilibration maximum (15-point moving
#                          average, first third of the trace) to +WIN_H hours.
#                          Returns rate (mg L^-1 h^-1, positive = consumption),
#                          R^2, window start and end (h).
# =============================================================================
suppressPackageStartupMessages(library(readxl))

WIN_H <- 10          # initial-rate window length (h) after the peak
AZ_COND <- c(B = "blank", A = "azide", AS = "azide+SHAM")
AZ_LAB  <- c(blank = "No cells", azide = "Azide 0.5 mM", `azide+SHAM` = "Azide 0.5 mM + SHAM 1 mM")
AZ_COL  <- c(blank = "#9A9A9A", azide = "#E08214", `azide+SHAM` = "#433D84")
AZ_ORDER <- c("blank", "azide", "azide+SHAM")

read_sdr_xlsx <- function(path, plate = NA) {
  raw <- suppressMessages(read_excel(path, sheet = "Oxygen Orginal", skip = 12,
                                     col_types = "text", .name_repair = "minimal"))
  nm <- names(raw)
  tcol <- which(nm == "Time/Min.")[1]
  wells <- nm[grepl("^[A-D][1-6]$", nm)]
  t_h <- suppressWarnings(as.numeric(raw[[tcol]])) / 60
  keep <- is.finite(t_h)
  out <- lapply(wells, function(w) data.frame(
    plate = plate, well = w, t_h = t_h[keep],
    o2 = suppressWarnings(as.numeric(raw[[w]][keep])), stringsAsFactors = FALSE))
  do.call(rbind, out)
}

peak_index <- function(o2, k = 15) {
  ok <- is.finite(o2); sm <- o2
  sm[ok] <- stats::filter(o2[ok], rep(1 / k, k), sides = 2)
  sm[is.na(sm)] <- -Inf
  n <- length(sm); which.max(sm[seq_len(max(3, n %/% 3))])
}

initial_rate <- function(t_h, o2, win = WIN_H) {
  ok <- is.finite(o2) & is.finite(t_h); t_h <- t_h[ok]; o2 <- o2[ok]
  pk <- peak_index(o2); t0 <- t_h[pk]
  m <- t_h >= t0 & t_h < t0 + win & o2 > 0.5
  f <- lm(o2 ~ t_h, data = data.frame(t_h = t_h[m], o2 = o2[m]))
  list(rate = -unname(coef(f)[2]), r2 = summary(f)$r.squared,
       t_start = t0, t_end = max(t_h[m]), intercept = unname(coef(f)[1]))
}

load_azide_plates <- function(root = ".") {
  lay <- read.csv(file.path(root, "data/oxygen/azide_layout.csv"), stringsAsFactors = FALSE)
  files <- c(`1` = "azide_rep1_Oxygen.xlsx", `2` = "azide_rep2_Oxygen.xlsx", `3` = "azide_rep3_Oxygen.xlsx")
  tr <- do.call(rbind, lapply(names(files), function(p)
    read_sdr_xlsx(file.path(root, "data/oxygen", files[[p]]), plate = as.integer(p))))
  tr$condition <- lay$condition[match(tr$well, lay$well)]
  tr$condition <- factor(tr$condition, levels = AZ_ORDER)
  tr
}
