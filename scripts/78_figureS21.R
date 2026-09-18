#!/usr/bin/env Rscript
# =============================================================================
# 78_figureS21.R -- Figure S21: is the prothioconazole x SHAM interaction an
# artefact of how the 27 C traces are denoised, or of the background oxygen loss?
#
# Two checks, both on the hand-selected fitting intervals, which are held fixed
# throughout so that only the pre-processing varies.
#
#   a  the moving-average width applied to the 27 C traces is varied from 0 to
#      12 h while 15 C is left unsmoothed, and the last row applies 12 h at both
#      temperatures. If the interaction were produced by attenuating the fast
#      vehicle wells more than the slow drug + SHAM wells, it would grow with
#      the width of the average. It does the opposite.
#   b  the cell-free wells lose oxygen slowly over the fitting intervals, by
#      0.021 mg/L/h at 27 C and nothing at 15 C, which is 10% of the decline in
#      the vehicle wells and 16% of that in the drug + SHAM wells. The right
#      panel adds that background back to every well before refitting.
#
# Inputs : data/oxygen/ptc_sham_rep{1,2,3}_{15,27}_Oxygen.xlsx
#          data/aox/ptc_sham_layout.csv
#          tables/aox/ptc_sham/manual_fit_windows.csv
# Outputs: tables/aox/ptc_sham_preprocessing.csv
#          figures/FigureS21.png / .pdf
# Run from the repository root:  Rscript scripts/78_figureS21.R
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/aox_common.R"))
INK <- "#1B2420"; INK2 <- "#555555"; C15 <- "#2166AC"; C27 <- "#B2182B"

lay <- read.csv(file.path(ROOT, "data/aox/ptc_sham_layout.csv"), stringsAsFactors = FALSE)
MW  <- read.csv(file.path(ROOT, "tables/aox/ptc_sham/manual_fit_windows.csv"), stringsAsFactors = FALSE)

TR <- list(); BSL <- list()
for (tp in c(15, 27)) for (rp in 1:3) {
  tr <- read_sdr_xlsx(file.path(ROOT, sprintf("data/oxygen/ptc_sham_rep%d_%d_Oxygen.xlsx", rp, tp)), plate = rp)
  tr$condition <- lay$condition[match(tr$well, lay$well)]
  TR[[paste(tp, rp)]] <- tr
  bw <- unique(tr$well[substr(tr$condition, 1, 1) == "b"])
  BSL[[paste(tp, rp)]] <- mean(sapply(bw, function(w) {       # background loss, mg/L/h
    g <- tr[tr$well == w & is.finite(tr$o2), ]; i0 <- which.max(g$o2)
    -unname(coef(lm(g$o2[i0:nrow(g)] ~ g$t_h[i0:nrow(g)]))[2]) }))
}

refit <- function(sm15, sm27, debackground = FALSE) {
  rows <- list()
  for (tp in c(15, 27)) for (rp in 1:3) {
    tr <- TR[[paste(tp, rp)]]; bsl <- BSL[[paste(tp, rp)]]
    for (w in unique(tr$well)) {
      g <- tr[tr$well == w & is.finite(tr$o2), ]; cond <- g$condition[1]
      if (substr(cond, 1, 1) == "b") next
      mw <- MW[MW$T == tp & MW$Dose == cond & MW$Replicate == sprintf("R%d%s", rp, w), ][1, ]
      tm <- g$t_h * 60; th <- g$t_h; o2 <- g$o2
      m <- th > 40 & th < 43; j <- which(m)[which.max(diff(o2[m]))]; st <- th[j]
      dy <- median(o2[th > st & th < st + 1]) - median(o2[th > st - 1 & th <= st])
      y <- step_correct(tm, o2, st, dy)
      if (debackground) y <- y + bsl * (th - th[1])
      y <- smooth_ma(tm, y, if (tp == 15) sm15 else sm27)
      k <- tm >= mw$fit_start & tm <= mw$fit_end & is.finite(y)
      rows[[length(rows) + 1]] <- data.frame(tp, rp, cond, r = fit_o2_model(tm[k], y[k])$r_per_h)
    } }
  aggregate(r ~ tp + rp + cond, do.call(rbind, rows), mean)
}
I_of <- function(A, t, p, d) { g <- function(c) A$r[A$tp == t & A$rp == p & A$cond == c]
  log(g(paste0("P", d, "S"))) - log(g(paste0("P", d))) - log(g("S")) + log(g("V")) }
summarise_one <- function(A, label) {
  i15 <- sapply(1:3, function(p) I_of(A, 15, p, 2)); i27 <- sapply(1:3, function(p) I_of(A, 27, p, 2))
  d <- i27 - i15; tt <- t.test(d)
  data.frame(label = label, I15 = mean(i15), I27 = mean(i27), P27 = t.test(i27)$p.value,
             diff = mean(d), lo = tt$conf.int[1], hi = tt$conf.int[2], P = tt$p.value)
}
CFG <- list(c(0,0), c(0,2), c(0,4), c(0,8), c(0,12), c(12,12))
S <- do.call(rbind, lapply(CFG, function(cf)
  cbind(summarise_one(refit(cf[1], cf[2]), sprintf("%g h / %g h", cf[1], cf[2])),
        sm15 = cf[1], sm27 = cf[2], background = "as measured")))
S2 <- cbind(summarise_one(refit(0, 12, TRUE), "0 h / 12 h"), sm15 = 0, sm27 = 12,
            background = "added back")
OUT_T <- file.path(ROOT, "tables/aox"); dir.create(OUT_T, showWarnings = FALSE, recursive = TRUE)
write.csv(rbind(S, S2), file.path(OUT_T, "ptc_sham_preprocessing.csv"), row.names = FALSE)

S$label <- factor(S$label, levels = rev(S$label))
th <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.28, colour = INK),
        axis.ticks = element_line(linewidth = 0.28, colour = INK),
        axis.text = element_text(size = 7.4, colour = INK),
        axis.title = element_text(size = 8.6, colour = INK),
        plot.margin = margin(4, 4, 2, 4))

L <- rbind(data.frame(label = S$label, temp = "15", I = S$I15),
           data.frame(label = S$label, temp = "27", I = S$I27))
pA <- ggplot(L, aes(I, label, colour = temp)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "#B0B4B8", linewidth = 0.3) +
  geom_line(aes(group = label), colour = "#C8C8C8", linewidth = 0.4) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(`15` = C15, `27` = C27),
                      labels = c("15 °C", "27 °C"), name = NULL) +
  labs(x = expression("interaction contrast "*italic(I)*" at 2 mg L"^-1),
       y = "moving average, 15 °C / 27 °C") + th +
  theme(legend.position = "top", legend.text = element_text(size = 7.4))

B <- rbind(cbind(S[S$sm27 == 12 & S$sm15 == 0, ], what = "as measured"),
           cbind(S2, what = "background added back"))
pB <- ggplot(B, aes(diff, what)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "#B0B4B8", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.10, linewidth = 0.45, colour = INK) +
  geom_point(size = 2.4, colour = INK) +
  labs(x = expression(italic(I)[27]-italic(I)[15]), y = NULL) + th

fig <- (pA | pB) + plot_layout(widths = c(1, 0.78)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", colour = INK))
OUT <- file.path(ROOT, "figures"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "FigureS21.png"), fig, width = 7.2, height = 3.4, dpi = 400, bg = "white")
ok <- tryCatch({ ggsave(file.path(OUT, "FigureS21.pdf"), fig, width = 7.2, height = 3.4,
                        device = cairo_pdf, bg = "white"); TRUE }, error = function(e) FALSE)
cat("Figure S21\n"); print(rbind(S, S2)[, c("label","background","I15","I27","P27","diff","lo","hi","P")],
                            row.names = FALSE, digits = 3)
cat(sprintf("\nbackground loss over the fitting intervals (mg/L/h): 15 C %.4f, 27 C %.4f\n",
            mean(unlist(BSL[grep("^15", names(BSL))])), mean(unlist(BSL[grep("^27", names(BSL))]))))
if (!ok) cat("note: PDF not written (cairo unavailable)\n")
