## R_figureS4_magnitude.R — Supplementary Figure S4 (threshold-free magnitude,
## classification, shared programme), ggplot2 port
## RUN FROM the master "X0123 copy" folder:  Rscript revision/scripts_R/R_figureS4_magnitude.R
suppressMessages({library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(readr)})

TAB <- "tables/revision"; DE <- "tables/rnaseq"
OUT <- "figures"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source("scripts/fig_style.R")
TEMPCOL <- TEMP3
conds <- c("15_0","15_2","21_0","21_2","27_0","27_2")
th <- theme_pub(9) +
  theme(plot.title = element_text(size = 10, face = "plain", hjust = 0))

ld <- function(f) read_csv(file.path(DE, f), show_col_types = FALSE) |>
  select(gene_id, log2FoldChange, padj, significant)
d15 <- ld("DE_proth_at_15C_2vs0.csv"); d21 <- ld("DE_proth_at_21C_2vs0.csv")
d27 <- ld("DE_proth_at_27C_2vs0.csv"); t2715 <- ld("DE_temp_27vs15_proth0.csv")
L <- d15 |> rename(l15 = log2FoldChange, s15 = significant) |> select(gene_id, l15, s15) |>
  inner_join(d21 |> transmute(gene_id, l21 = log2FoldChange), by = "gene_id") |>
  inner_join(d27 |> transmute(gene_id, l27 = log2FoldChange), by = "gene_id")
conc <- read_csv(file.path(TAB, "magnitude/concordance_slopes.csv"), show_col_types = FALSE)
cls  <- read_csv(file.path(TAB, "magnitude/gene_classification_counts.csv"), show_col_types = FALSE)
Sm   <- read.csv(file.path(TAB, "magnitude/module_scores_by_condition.csv"), row.names = 1,
                 check.names = FALSE)

## ---- A: ECDF of |LFC| -------------------------------------------------------
ha <- L |> select(l15, l21, l27) |> pivot_longer(everything()) |>
  mutate(t = sub("l", "", name))
pA <- ggplot(ha, aes(abs(value), colour = t)) + stat_ecdf(linewidth = 0.6) +
  coord_cartesian(xlim = c(0, 2.5)) +
  scale_colour_manual(values = TEMPCOL, name = NULL, labels = paste0(c(15,21,27), " \u00b0C")) +
  labs(title = "A  Genome-wide drug effect size",
       x = expression("|log"[2]*"FC| (drug 2 vs 0, ashr)"), y = "ECDF") +
  th + theme(legend.position = c(0.8, 0.35))

## ---- B: retention scatter among 15C-significant -----------------------------
sig <- L |> filter(s15)
s21 <- conc |> filter(scope == "15C_significant", grepl("21C", comparison)) |> pull(slope_OLS_origin)
s27 <- conc |> filter(scope == "15C_significant", grepl("27C", comparison)) |> pull(slope_OLS_origin)
hb <- bind_rows(sig |> transmute(x = l15, y = l21, t = "21"),
                sig |> transmute(x = l15, y = l27, t = "27"))
pB <- ggplot(hb, aes(x, y, colour = t)) +
  geom_point(size = 0.4, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = s21, intercept = 0, colour = TEMPCOL["21"], linewidth = 0.6) +
  geom_abline(slope = s27, intercept = 0, colour = TEMPCOL["27"], linewidth = 0.6) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_colour_manual(values = TEMPCOL[c("21","27")], name = NULL,
                      labels = paste0(c(21, 27), " \u00b0C")) +
  labs(title = sprintf("B  Effect retention (slopes %.2f, %.2f)", s21, s27),
       x = "drug effect at 15 \u00b0C (log2FC)", y = "drug effect at 21 / 27 \u00b0C") +
  th + theme(legend.position = c(0.15, 0.9)) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

## ---- C: gene classes --------------------------------------------------------
lab_map <- c("fungicide-specific" = "fungicide-\nspecific",
             "shared, same direction" = "shared\n(same dir.)",
             "heat-specific" = "heat-\nspecific")
hc <- cls |> filter(class %in% names(lab_map)) |>
  mutate(lab = factor(lab_map[class], levels = unname(lab_map)))
pC <- ggplot(hc, aes(lab, n, fill = lab)) + geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("#6A3D9A", "#1B7837", "#D55E00"), guide = "none") +
  labs(title = "C  Gene classes (\u00a75.6)", x = NULL, y = "genes") + th

## ---- D: module scores heatmap ----------------------------------------------
hd <- Sm |> mutate(cond = rownames(Sm)) |>
  pivot_longer(-cond, names_to = "module", values_to = "z") |>
  mutate(cond = factor(cond, levels = conds),
         module = factor(module, levels = rev(colnames(Sm))))
pD <- ggplot(hd, aes(cond, module, fill = z)) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-0.8, 0.8), oob = scales::squish,
                       name = "module score\n(mean z)") +
  scale_x_discrete(labels = sub("_", "/", conds)) +
  labs(title = "D  Stress-programme module scores: 27 \u00b0C control already resembles drug-treated",
       x = NULL, y = NULL) + th + theme(axis.text.y = element_text(size = 7))

## ---- E: heat vs drug convergence -------------------------------------------
he <- t2715 |> transmute(gene_id, heat = log2FoldChange) |>
  inner_join(L |> select(gene_id, l15), by = "gene_id") |> filter(!is.na(heat), !is.na(l15))
r <- cor(he$heat, he$l15)
pE <- ggplot(he, aes(l15, heat)) + geom_point(size = 0.4, alpha = 0.25, colour = "#555555") +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  labs(title = sprintf("E  Heat vs drug convergence (r=%.2f)", r),
       x = "drug effect at 15 \u00b0C (log2FC)", y = "warming effect 27v15 \u00b0C, no drug") + th

layout <- (pA | pB | pC) / (pD | pE) + plot_layout(widths = c(1)) +
  plot_annotation(title = "Threshold-free magnitude, classification and shared-programme analyses",
                  theme = theme(plot.title = element_text(size = 11, hjust = 0.5)))
ggsave(file.path(OUT, "FigureS5.png"), layout, width = 13, height = 9, dpi = 200)
ggsave(file.path(OUT, "FigureS5.pdf"), layout, width = 13, height = 9, device = cairo_pdf)
cat("saved", file.path(OUT, "FigureS5.png"), "\n")
