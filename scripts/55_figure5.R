# =============================================================================
# 55_figure5.R -- Figure 5: prothioconazole and temperature converge on shared
# programmes that track whole-organism physiology.
# R port of 80_figure5_redesign.py (values unchanged; panel A reads the cached
# GSEA table, exactly as the Python does when that cache exists).
#   A  pooled drug NES vs warming NES, pathway level
#   B  gene-level alignment of the drug effect on the warming effect, per temperature
#   C  ribosome module score against growth rate, six condition means
# Also writes tables/revision/fig5/fig5_loo_sensitivity.csv (used by Fig. S13).
# Outputs: figures/Figure5.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork); library(ggrepel)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/fig_common.R")); anno_load(ROOT)

## ---- A: pathway-level convergence ------------------------------------------
CA <- read.csv(file.path(ROOT, "tables/revision/fig5/fig5A_gsea_cache.csv"), stringsAsFactors = FALSE)
CA$sig <- as.logical(CA$sig)
rA <- cor(CA$nes_drug, CA$nes_warm)
lim <- max(abs(c(CA$nes_drug, CA$nes_warm))) * 1.12
HI <- data.frame(
  kegg_id = c("map03010","map00190","map02010","map00982","map03050","map00020"),
  lab = c("Ribosome","OXPHOS","ABC transporters (efflux)","Drug metab. P450","Proteasome","TCA cycle"),
  col = c("#2F6BB3","#E08214","#7B3FA0","#D6604D","#11838E","#8C510A"), stringsAsFactors = FALSE)
HI <- merge(HI, CA[, c("kegg_id","nes_drug","nes_warm")], by = "kegg_id")

pA <- ggplot(CA, aes(nes_drug, nes_warm)) +
  geom_hline(yintercept = 0, colour = "#E6E6E6", linewidth = 0.28) +
  geom_vline(xintercept = 0, colour = "#E6E6E6", linewidth = 0.28) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#BBBBBB", linewidth = 0.35) +
  geom_point(data = filter(CA, !sig), colour = "#DADADA", size = 1.0, alpha = 0.6) +
  geom_point(data = filter(CA, sig), colour = "#4D4D4D", size = 1.5, alpha = 0.85) +
  geom_point(data = HI, aes(colour = lab), size = 2.1, show.legend = FALSE) +
  geom_text_repel(data = HI, aes(label = lab, colour = lab), size = 2.5, fontface = "bold",
                  segment.colour = "#999999", segment.size = 0.25, box.padding = 0.4,
                  max.overlaps = 20, show.legend = FALSE) +
  annotate("text", x = lim, y = -lim, hjust = 1, vjust = 0, size = 2.9, colour = "#333333",
           label = sprintf("r = %.2f (descriptive)", rA)) +
  scale_colour_manual(values = setNames(HI$col, HI$lab)) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "Common drug pathway enrichment\n(mean NES over 15 & 21 °C)",
       y = "Warming pathway enrichment (NES, 27 vs 15 °C)") +
  theme_classic(base_size = 9) + theme(axis.text = element_text(size = 7.5))

## ---- B: gene-level alignment ------------------------------------------------
warm <- load_de("temp_27vs15_proth0", ROOT)
drug <- lapply(c("15","21","27"), function(t) load_de(paste0("proth_at_", t, "C_2vs0"), ROOT))
names(drug) <- c("15","21","27")
common <- Reduce(intersect, c(list(warm$gene), lapply(drug, `[[`, "gene")))
xv <- warm$lfc[match(common, warm$gene)]
B <- bind_rows(lapply(names(drug), function(t)
  data.frame(x = xv, y = drug[[t]]$lfc[match(common, drug[[t]]$gene)], temp = t)))
FIT <- B %>% group_by(temp) %>% summarise(b = coef(lm(y ~ x))[2], a = coef(lm(y ~ x))[1],
                                          r = cor(x, y), .groups = "drop")
B$facet   <- factor(paste("Drug at", B$temp, "°C"), levels = paste("Drug at", c(15,21,27), "°C"))
FIT$facet <- factor(paste("Drug at", FIT$temp, "°C"), levels = levels(B$facet))
FITL <- do.call(rbind, lapply(seq_len(nrow(FIT)), function(i)
  data.frame(x = seq(-8, 8, length.out = 100),
             y = FIT$a[i] + FIT$b[i]*seq(-8, 8, length.out = 100),
             temp = FIT$temp[i], facet = FIT$facet[i])))
pB <- ggplot(B, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "#ECECEC", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "#ECECEC", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#E0E0E0", linewidth = 0.22) +
  geom_point(colour = "#9AA2AB", alpha = 0.12, size = 0.22) +
  geom_line(data = FITL, aes(x, y, colour = temp), linewidth = 0.95, lineend = "round") +
  geom_text(data = FIT, aes(-7.6, 6.6, label = sprintf("alignment coeff. = %+.2f\nr = %+.2f", b, r)),
            hjust = 0, vjust = 1, size = 2.5, colour = "#333333") +
  facet_wrap(~facet, nrow = 1) +
  scale_colour_manual(values = TEMP, guide = "none") +
  coord_cartesian(xlim = c(-8, 8), ylim = c(-7, 7)) +
  labs(x = expression("Warming effect, 27 vs 15 °C (log"[2]*"FC, no drug)"),
       y = expression("Prothioconazole effect (log"[2]*"FC)")) +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 7.5), panel.spacing = unit(2.2, "mm"))
# colour the facet titles by temperature is not possible in base ggplot; keep bold black

## ---- C: ribosome module vs growth ------------------------------------------
MS <- read.csv(file.path(ROOT, "tables/revision/integration/module_scores_by_condition_extended.csv"),
               check.names = FALSE)
PH <- read.csv(file.path(ROOT, "tables/revision/integration/physiology_condition_values.csv"))
D <- merge(PH, MS[, c("cond","Ribosome","AOX")], by = "cond")
D$temp <- factor(sub("_.*", "", D$cond), levels = c("15","21","27"))
D$dose <- sub(".*_", "", D$cond)
rC <- cor(D$Ribosome, D$growth)
pC <- ggplot(D, aes(Ribosome, growth)) +
  geom_line(aes(group = temp), colour = "#CFCFCF", linewidth = 0.4) +
  geom_point(aes(colour = temp, shape = dose), size = 2.6) +
  geom_text(aes(label = sprintf("%s °C, %s", temp, dose), colour = temp),
            size = 2.3, hjust = -0.15, vjust = -0.5, show.legend = FALSE) +
  annotate("text", x = min(D$Ribosome), y = max(D$growth), hjust = 0, vjust = 1, size = 2.7,
           colour = "#333333",
           label = sprintf("r = %+.2f\n(descriptive, n = 6\ncondition means)", rC)) +
  scale_colour_manual(values = TEMP, guide = "none") +
  scale_shape_manual(values = c(`0` = 16, `2` = 17), name = "Prothioconazole",
                     labels = c(expression(0~mg~L^-1), expression(2~mg~L^-1))) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.22))) +
  labs(x = "Ribosome module score (VST z, 117 genes)",
       y = expression("Growth rate (C C"^-1~h^-1*")")) +
  theme_classic(base_size = 9.5) +
  theme(axis.text = element_text(size = 8), legend.position = c(0.98, 0.04),
        legend.justification = c(1, 0), legend.background = element_blank(),
        legend.text = element_text(size = 7.4), legend.title = element_text(size = 7.8))

## ---- leave-one-condition-out sensitivity (feeds Fig. S13) ------------------
loo <- do.call(rbind, lapply(list(c("Ribosome","growth"), c("AOX","CUE")), function(mt) {
  x <- D[[mt[1]]]; y <- D[[mt[2]]]
  v <- sapply(seq_along(x), function(i) cor(x[-i], y[-i]))
  out <- data.frame(module = mt[1], trait = mt[2], r_full = round(cor(x, y), 3),
                    r_loo_min = round(min(v), 3), r_loo_max = round(max(v), 3))
  cbind(out, setNames(as.data.frame(t(round(v, 3))), paste0("r_wo_", D$cond))) }))
write.csv(loo, file.path(ROOT, "tables/revision/fig5/fig5_loo_sensitivity.csv"), row.names = FALSE)

fig <- (pA + pB + plot_layout(widths = c(1, 1.9))) / pC +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = list(c("A","B","C"))) &
  theme(plot.tag = element_text(size = 13, face = "bold"), plot.tag.position = c(0, 1.02))
dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(ROOT, "figures/Figure5.png"), fig, width = 10.3, height = 8.4, dpi = 400, bg = "white")
ggsave(file.path(ROOT, "figures/Figure5.pdf"), fig, width = 10.3, height = 8.4, device = cairo_pdf, bg = "white")
cat(sprintf("Figure5: pathway r = %.2f | alignment %s | ribosome-growth r = %+.2f\n",
    rA, paste(sprintf("%.2f", FIT$b), collapse = "/"), rC))
