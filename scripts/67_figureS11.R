# =============================================================================
# 67_figureS11.R -- Figure S11: reference-free test of the temperature x dose
# interaction (port of 48_figureS14.py; no Bliss framework, no reference temp).
# Inputs : tables/revision/interaction/{interaction_residual_surface,
#          interaction_modelfree_stats}.csv   (script 47)
# Outputs: figures/FigureS11.png, figures/FigureS11.pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork); library(scales)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/fig_style.R"))
TAB <- file.path(ROOT, "tables/revision/interaction")

G <- read.csv(file.path(TAB, "interaction_residual_surface.csv"))
S0 <- read.csv(file.path(TAB, "interaction_modelfree_stats.csv"))
S  <- setNames(S0$value, S0$stat)
G$dose <- factor(G$conc, levels = sort(unique(G$conc)))
RAMP <- DOSE_RAMP[-1]                                   # drop the control grey
names(RAMP) <- levels(G$dose)

## ---- A: residual against temperature, one line per dose --------------------
ymax <- max(G$delta_med)
pA <- ggplot(G, aes(Temp, delta_med, colour = dose)) +
  geom_hline(yintercept = 0, linetype = "32", colour = "#BBBBBB", linewidth = 0.3) +
  geom_vline(xintercept = 24, colour = "#DDDDDD", linewidth = 0.28) +
  geom_line(linewidth = 0.45) +
  annotate("text", x = 24.2, y = ymax * 0.93, label = "T[opt]", parse = TRUE,
           size = 2.1, colour = INK2, hjust = 0) +
  scale_colour_manual(values = RAMP, name = expression(mg~L^-1)) +
  labs(x = "Temperature (°C)",
       y = expression(atop("Interaction residual", italic(F)*"("*T*","*c*") − [ "*italic(f)*"("*T*") + "*italic(h)*"("*c*") ]  (log scale)"))) +
  theme_pub() + theme(legend.position = c(0.86, 0.16), legend.key.height = unit(2.6, "mm"),
                      legend.text = element_text(size = 5.6), legend.title = element_text(size = 5.8)) +
  guides(colour = guide_legend(ncol = 2))

## ---- B: the same residual as a surface, comparable with Fig. 2F ------------
v <- max(abs(G$delta_med))
pB <- ggplot(G, aes(Temp, dose, fill = delta_med)) +
  geom_raster(interpolate = FALSE) +
  geom_vline(xintercept = 24, colour = "#333333", linewidth = 0.28, linetype = "32") +
  scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0,
                       limits = c(-v, v), name = "Interaction residual\n(log scale)") +
  labs(x = "Temperature (°C)", y = expression("Prothioconazole (mg L"^-1*")"),
       title = "blue, growth below separable (synergy)\nred, growth above separable (antagonism)") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_pub() + theme(plot.title = element_text(size = 5.8, face = "plain", colour = INK2),
                      legend.key.width = unit(2.4, "mm"), legend.key.height = unit(6, "mm"),
                      legend.text = element_text(size = 5.6), legend.title = element_text(size = 5.8))

## ---- C: exact leave-one-out evidence ---------------------------------------
B <- data.frame(m = factor(c("sep", "int"), levels = c("sep", "int")),
                rmse = c(S[["loo_rmse_sep"]], S[["loo_rmse_int"]]))
cap <- sprintf("paired ~Delta~'squared LOO error'~'%.4f [%.4f, %.4f], '*italic(P)*' = %.3f'",
               S[["loo_paired_mean"]], S[["loo_paired_lo"]], S[["loo_paired_hi"]], S[["loo_paired_p"]])
pC <- ggplot(B, aes(m, rmse, fill = m)) +
  geom_col(width = 0.5, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.4f", rmse)), vjust = -0.4, size = 2.2, colour = INK) +
  annotate("text", x = 1.5, y = S[["loo_rmse_sep"]] * 1.36, size = 2.1, colour = INK2, vjust = 1,
           label = sprintf("paired Δ squared LOO error\n%.4f [%.4f, %.4f], P = %.3f\nΔAIC = %.1f;  ti(T,c): P = %.1e",
                           S[["loo_paired_mean"]], S[["loo_paired_lo"]], S[["loo_paired_hi"]],
                           S[["loo_paired_p"]], S[["dAIC"]], S[["ti_p"]])) +
  scale_fill_manual(values = c(sep = "#C9C9C9", int = "#27519E"), guide = "none") +
  scale_x_discrete(labels = c(sep = "separable\nf(T) + h(c)", int = "interaction\n+ ti(T,c)")) +
  scale_y_continuous(limits = c(0, S[["loo_rmse_sep"]] * 1.38), expand = c(0, 0)) +
  labs(x = NULL, y = "Exact leave-one-out RMSE (log growth)") + theme_pub()

fig <- (pA | pB | pC) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"), plot.tag.position = c(0, 1.02))
dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(ROOT, "figures/FigureS11.png"), fig, width = 9.6, height = 3.3, dpi = 600, bg = "white")
ggsave(file.path(ROOT, "figures/FigureS11.pdf"), fig, width = 9.6, height = 3.3, device = cairo_pdf, bg = "white")
cat(sprintf("FigS11: LOO sep %.4f vs int %.4f | dAIC %.1f | ti P %.2e\n",
            S[["loo_rmse_sep"]], S[["loo_rmse_int"]], S[["dAIC"]], S[["ti_p"]]))
