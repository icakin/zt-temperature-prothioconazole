## R_figureS5_physiology.R — Supplementary Figure S5 (EC50 extrapolation,
## growth-rate test, Bliss reference scan), ggplot2 port
## RUN FROM the master "X0123 copy" folder:  Rscript revision/scripts_R/R_figureS5_physiology.R
source("scripts/fig_style.R")
suppressMessages({library(ggplot2); library(patchwork); library(dplyr); library(readr); library(tidyr)})

TAB <- "tables/revision"; OUT <- "figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
th <- theme_pub(9) +
  theme(plot.title = element_text(size = 10, face = "plain", hjust = 0))

H  <- read_csv(file.path(TAB, "physiology/growthrate_vs_pec50.csv"), show_col_types = FALSE)
BL <- read_csv(file.path(TAB, "physiology/bliss_reference_scan_summary.csv"), show_col_types = FALSE)
gr <- read_csv(file.path(TAB, "physiology/growthrate_tolerance_stats.csv"), show_col_types = FALSE)
r_p <- gr$value[gr$stat == "pearson_r"]; p_p <- gr$value[gr$stat == "pearson_p"]

## ---- A: pEC50 vs r0 ---------------------------------------------------------
H <- H |> mutate(side = factor(ifelse(supra, "supra", "sub"), levels = c("sub","supra")))
pA <- ggplot(H, aes(r0, pEC50_median)) +
  geom_errorbar(aes(ymin = pEC50_median - 1.96*pEC50_sd, ymax = pEC50_median + 1.96*pEC50_sd),
                colour = "grey80", width = 0) +
  geom_point(aes(colour = side), size = 2.6) +
  geom_text(aes(label = paste0(temperature, "\u00b0")), hjust = -0.35, vjust = -0.35, size = 2.8) +
  scale_colour_manual(values = c(sub = "#0072B2", supra = "#D55E00"), name = NULL,
                      labels = c("\u2264 Topt", "> Topt")) +
  labs(title = sprintf("A  Sensitivity vs baseline growth (r=%.2f, p=%.2f)", r_p, p_p),
       x = expression("untreated growth rate r"[0]*" (h"^-1*")"),
       y = expression("pEC"[50]*" (higher = more sensitive)")) +
  th + theme(legend.position = c(0.85, 0.2))

## ---- B: EC50 vs temperature with tested range ------------------------------
pB <- ggplot(H, aes(temperature, EC50)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 4, fill = "grey93") +
  geom_hline(yintercept = 4, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(ymin = EC50_lower, ymax = EC50_upper), width = 0.35, linewidth = 0.35) +
  geom_point(size = 2) +
  annotate("text", x = 16.6, y = 4.5, hjust = 0, size = 2.7, colour = "grey40",
           label = "highest tested dose (4 mg/L)") +
  labs(title = expression("B  27"*degree*"C  EC"[50]*" exceeds tested range"),
       x = "temperature (\u00b0C)", y = expression("EC"[50]*" (mg L"^-1*")")) + th

## ---- C: Bliss deviation sign vs reference choice ---------------------------
hc <- BL |> select(reference_T, median_delta_15C, median_delta_27C) |>
  pivot_longer(-reference_T, names_to = "at", values_to = "delta") |>
  mutate(at = ifelse(grepl("15", at), "d15", "d27"))
pC <- ggplot(hc, aes(factor(reference_T), delta, fill = at)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c(d15 = "#2166AC", d27 = "#B2182B"), name = NULL,
                    labels = c("\u0394 at 15 \u00b0C", "\u0394 at 27 \u00b0C")) +
  labs(title = "C  Interaction sign vs reference choice",
       x = "reference temperature (\u00b0C)", y = "median Bliss deviation \u0394") +
  th + theme(legend.position = c(0.85, 0.9))

layout <- pA | pB | pC
ggsave(file.path(OUT, "FigureS9.png"), layout, width = 12.5, height = 3.9, dpi = 200)
ggsave(file.path(OUT, "FigureS9.pdf"), layout, width = 12.5, height = 3.9, device = cairo_pdf)
cat("saved", file.path(OUT, "FigureS9.png"), "\n")
