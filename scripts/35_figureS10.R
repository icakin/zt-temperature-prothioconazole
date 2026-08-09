## R_figureS6_integration.R — Supplementary Figure S6 (transcriptome-physiology
## integration, membrane sweep, heat-control resemblance), ggplot2 port
## RUN FROM the master "X0123 copy" folder:  Rscript revision/scripts_R/R_figureS6_integration.R
suppressMessages({library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(readr)})

TAB <- "tables/revision"; OUT <- "figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source("scripts/fig_style.R")
TEMPCOL <- TEMP3
conds <- c("15_0","15_2","21_0","21_2","27_0","27_2")
th <- theme_pub(9) +
  theme(plot.title = element_text(size = 10, face = "plain", hjust = 0))

phys <- read_csv(file.path(TAB,"integration/physiology_condition_values.csv"), show_col_types=FALSE)
Sm   <- read.csv(file.path(TAB,"integration/module_scores_by_condition_extended.csv"),
                 row.names = 1, check.names = FALSE)
mem  <- read_csv(file.path(TAB,"integration/membrane_module_sweep.csv"), show_col_types=FALSE)
res  <- read_csv(file.path(TAB,"integration/heatctrl_resemblance_summary.csv"), show_col_types=FALSE)

D <- phys |> mutate(rib = Sm[cond, "Ribosome"], aox = Sm[cond, "AOX"],
                    unc = Sm[cond, "AOX"] - Sm[cond, "OXPHOS"],
                    t = as.character(T), dose = factor(dose))

mk <- function(df, xv, yv, xlab, ylab, letter, ttl) {
  r <- cor.test(df[[xv]], df[[yv]])
  ggplot(df, aes(.data[[xv]], .data[[yv]], colour = t, shape = dose)) +
    geom_point(size = 3, stroke = 0.7) +
    geom_text(aes(label = sub("_", "/", cond)), hjust = -0.3, vjust = -0.5, size = 2.6,
              colour = "black") +
    scale_colour_manual(values = TEMPCOL, guide = "none") +
    scale_shape_manual(values = c(`0` = 16, `2` = 15), guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    labs(title = sprintf("%s  %s (r=%.2f, p=%.3f)", letter, ttl, r$estimate, r$p.value),
         x = xlab, y = ylab) + th
}
pA <- mk(D, "rib", "growth", "Ribosome module score", expression("growth rate (h"^-1*")"),
         "A", "Growth vs ribosome")
pB <- mk(D, "aox", "respiration", "AOX score", expression("respiration (h"^-1*")"),
         "B", "Respiration vs AOX")
pC <- mk(D, "unc", "CUE", "AOX \u2212 OXPHOS score (uncoupling)", "model-derived CUE",
         "C", "CUE vs uncoupling")

## ---- D: membrane sweep heatmap ---------------------------------------------
keep <- c("Ergosterol","FA biosynthesis","FA elongation","Unsat. FA biosynthesis",
          "FA degradation","Glycerophospholipid","Glycerolipid","Ether lipid",
          "Sphingolipid","GPI anchor","Peroxisome","Starch/sucrose (storage)",
          "Cell wall (chitin/glucan)")
hd <- mem |> filter(module %in% keep) |>
  mutate(mlab = paste0(module, " (n=", n, ")")) |>
  select(mlab, drug_meanLFC_15C, drug_meanLFC_21C, drug_meanLFC_27C, heat_meanLFC_27v15) |>
  pivot_longer(-mlab, names_to = "contrast", values_to = "lfc") |>
  mutate(contrast = factor(contrast,
           levels = c("drug_meanLFC_15C","drug_meanLFC_21C","drug_meanLFC_27C","heat_meanLFC_27v15"),
           labels = c("drug 15\u00b0C","drug 21\u00b0C","drug 27\u00b0C","heat 27v15")))
hd$mlab <- factor(hd$mlab, levels = rev(unique(hd$mlab)))
pD <- ggplot(hd, aes(contrast, mlab, fill = lfc)) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-1, 1), oob = scales::squish, name = "mean log2FC") +
  labs(title = "D  Membrane / storage remodelling (mean shrunken log2FC)", x = NULL, y = NULL) +
  th + theme(axis.text.y = element_text(size = 7))

## ---- E: resemblance bars ----------------------------------------------------
res <- res |> mutate(lab = factor(c("vs 21\u00b0C untreated","vs 21\u00b0C drug",
                                    "vs 15\u00b0C untreated","vs 15\u00b0C drug"),
                     levels = rev(c("vs 21\u00b0C untreated","vs 21\u00b0C drug",
                                    "vs 15\u00b0C untreated","vs 15\u00b0C drug"))),
                     kind = c("untreated","drug","untreated","drug"))
pE <- ggplot(res, aes(n_DE_lfc1, lab, fill = kind)) + geom_col(width = 0.6) +
  scale_fill_manual(values = c(untreated = "grey60", drug = "#1B7837"), guide = "none") +
  labs(title = "E  What does 27 \u00b0C control resemble?",
       x = "genes different (padj<0.05, |LFC|\u22651)", y = NULL) + th

layout <- (pA | pB | pC) / (pD | pE) +
  plot_annotation(title = "Transcriptome-physiology integration (\u00a75.7) and membrane remodelling (\u00a75.4)",
                  theme = theme(plot.title = element_text(size = 11, hjust = 0.5)))
ggsave(file.path(OUT, "FigureS10.png"), layout, width = 13, height = 8.5, dpi = 200)
ggsave(file.path(OUT, "FigureS10.pdf"), layout, width = 13, height = 8.5, device = cairo_pdf)
cat("saved", file.path(OUT, "FigureS10.png"), "\n")
