#!/usr/bin/env Rscript
# =============================================================================
# 48_ptc_sham_window_rules.R -- how much does the prothioconazole x SHAM result
# depend on the rule used to choose the fitting intervals? (Figure S19)
#
# Every one of the 132 curves is refitted under each rule in RULES (ps_window()
# in aox_common.R), the plate x condition means and the interaction contrasts
# are recomputed, and two things are reported per rule:
#
#   the contrast   I_27 - I_15 at each dose, its mean over the three cultures
#                  and the paired t;
#   the symmetry   the ratio of the median fitted interval of the drug + SHAM
#                  wells to that of the vehicle wells at 27 C. A rule that gives
#                  the slow combination wells much longer intervals than the
#                  vehicle is measuring different phases of growth in the two
#                  conditions and will exaggerate the contrast, whatever the
#                  biology. The primary rule should be one whose ratio is near 1.
#
# Outputs: tables/aox/ptc_sham_window_rules.csv
# Run from the repository root:  Rscript scripts/48_ptc_sham_window_rules.R
# =============================================================================
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))
SMOOTH_BY_TEMP <- c(`15` = 0, `27` = 4)

RULES <- list(
  list(id = "stable r (5%)",            rule = list(type = "stable_r")),
  list(id = "peak + 24 h",              rule = list(type = "peak", dur = 24*60)),
  list(id = "peak + 30 h",              rule = list(type = "peak", dur = 30*60)),
  list(id = "peak + 36 h",              rule = list(type = "peak", dur = 36*60)),
  list(id = "fixed 3-36 h",             rule = list(type = "fixed", a = 180, b = 36*60)),
  list(id = "fixed 3-40 h",             rule = list(type = "fixed", a = 180, b = 40*60)),
  list(id = "fixed 3-48 h",             rule = list(type = "fixed", a = 180, b = 48*60)),
  list(id = "fixed 5-40 h",             rule = list(type = "fixed", a = 300, b = 40*60)),
  list(id = "drawdown 5-60%",           rule = list(type = "draw", a = 0.05, b = 0.60)),
  list(id = "drawdown 5-80%",           rule = list(type = "draw", a = 0.05, b = 0.80)),
  list(id = "drawdown 10-90%",          rule = list(type = "draw", a = 0.10, b = 0.90))
)

traces <- list()
steps <- read.csv(file.path(ROOT, "tables/aox/ptc_sham_steps.csv"), stringsAsFactors = FALSE)
steps$Replicate <- sprintf("R%d%s", steps$replicate, steps$well)
for (temp in c(15, 27)) {
  d <- read.csv(file.path(ROOT, sprintf("data/aox/ptc_sham_%d_Oxygen.csv", temp)), check.names = FALSE)
  SM <- unname(SMOOTH_BY_TEMP[as.character(temp)])
  for (cv in setdiff(names(d), c("Time", "T"))) {
    cond <- sub("_R.*$", "", cv); repw <- sub("^.*_", "", cv)
    if (substr(cond, 1, 1) == "b") next
    st <- steps[steps$T == temp & steps$Replicate == repw, ][1, ]
    y <- smooth_ma(d$Time, step_correct(d$Time, d[[cv]], st$step_time_h, st$step_mgL), SM)
    traces[[length(traces) + 1]] <- list(temp = temp, replicate = as.integer(substr(repw, 2, 2)),
      condition = cond, well = substr(repw, 3, 4), t = d$Time, y = y)
  }
}

out <- list()
for (R in RULES) {
  rows <- lapply(traces, function(g) {
    w <- ps_window(g$t, g$y, R$rule)
    if (!all(is.finite(w[1:2]))) return(NULL)
    m <- g$t >= w[1] & g$t <= w[2] & is.finite(g$y)
    f <- fit_o2_model(g$t[m], g$y[m]); if (is.null(f)) return(NULL)
    data.frame(temp = g$temp, replicate = g$replicate, condition = g$condition, well = g$well,
               r_per_h = f$r_per_h, R2 = f$R2, len_h = (w[2] - w[1])/60, fallback = as.logical(w[3]))
  })
  D <- do.call(rbind, rows)
  M <- aggregate(r_per_h ~ temp + replicate + condition, D, mean)
  g <- function(tp, rp, cd) { v <- M$r_per_h[M$temp == tp & M$replicate == rp & M$condition == cd]; if (length(v)) v else NA_real_ }
  I <- function(tp, rp, dd) log(g(tp,rp,paste0(dd,"S"))) - log(g(tp,rp,dd)) - log(g(tp,rp,"S")) + log(g(tp,rp,"V"))
  d2 <- sapply(1:3, function(rp) I(27,rp,"P2") - I(15,rp,"P2"))
  d4 <- sapply(1:3, function(rp) I(27,rp,"P4") - I(15,rp,"P4"))
  lenPS <- median(D$len_h[D$temp == 27 & D$condition %in% c("P2S","P4S")])
  lenV  <- median(D$len_h[D$temp == 27 & D$condition == "V"])
  out[[length(out) + 1]] <- data.frame(
    rule = R$id, len_ratio_PS_over_V = lenPS/lenV, R2_min = min(D$R2), fallback = sum(D$fallback),
    I27_I15_2 = mean(d2), P_2 = t.test(d2)$p.value, n_neg_2 = sum(d2 < 0),
    I27_I15_4 = mean(d4), P_4 = t.test(d4)$p.value, n_neg_4 = sum(d4 < 0),
    d2_1 = d2[1], d2_2 = d2[2], d2_3 = d2[3])
}
S <- do.call(rbind, out)
dir.create(file.path(ROOT, "tables/aox"), showWarnings = FALSE, recursive = TRUE)
write.csv(S, file.path(ROOT, "tables/aox/ptc_sham_window_rules.csv"), row.names = FALSE)
cat("Interval rule sweep (132 curves refitted under each rule)\n")
cat("len ratio = median interval of the 27 C drug+SHAM wells / that of the vehicle wells;\n")
cat("a rule near 1 fits both conditions over comparable spans.\n\n")
print(S[, c("rule","len_ratio_PS_over_V","R2_min","I27_I15_2","P_2","n_neg_2","I27_I15_4","P_4","n_neg_4")],
      row.names = FALSE, digits = 3)
