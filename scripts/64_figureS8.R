# =============================================================================
# 64_figureS8.R -- Figure S8: gene-level convergence of the drug and warming
# programmes. R port of the supplementary panel in 28_figure5_S8.py.
# Each gene's PEAK prothioconazole effect across 15-27 C against its warming
# effect (27 vs 15 C, no drug). Genes significant in >=1 temperature are
# coloured by whether the two effects share a direction.
# Inputs : tables/rnaseq/DE_proth_at_{15,21,27}C_2vs0.csv, DE_temp_27vs15_proth0.csv
# Outputs: figures/FigureS8.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggrepel)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/fig_common.R")); anno_load(ROOT)

warm <- load_de("temp_27vs15_proth0", ROOT)
drug <- lapply(c("15","21","27"), function(t) load_de(paste0("proth_at_", t, "C_2vs0"), ROOT))
names(drug) <- c("15","21","27")
common <- Reduce(intersect, c(list(warm$gene), lapply(drug, `[[`, "gene")))

M <- sapply(drug, function(d) d$lfc[match(common, d$gene)])
Q <- sapply(drug, function(d) d$padj[match(common, d$gene)])
peak <- M[cbind(seq_len(nrow(M)), max.col(abs(M)))]           # largest |log2FC| across temperatures
sigany <- apply(Q < FDR_MAX & abs(M) >= LFC_MIN, 1, any)
D <- data.frame(gene = common, x = warm$lfc[match(common, warm$gene)], y = peak, sig = sigany)
D$grp <- ifelse(!D$sig, "ns", ifelse(sign(D$x) != sign(D$y), "divergent",
                                     ifelse(D$x > 0, "shared_up", "shared_dn")))
S <- filter(D, sig)
rS <- cor(S$x, S$y); conc <- 100 * mean(sign(S$x) == sign(S$y))

HL <- D %>% filter(gene %in% c(aoxg()[1], "Mycgr3G50464")) %>%
  mutate(lab = ifelse(gene == "Mycgr3G50464", "CDR1 (drug-specific)", "AOX (shared)"),
         col = ifelse(gene == "Mycgr3G50464", "#6A3D9A", "#B2182B"))

COLS <- c(ns = "#E8E8E8", shared_up = UP, shared_dn = DN, divergent = DIV)
fig <- ggplot(D, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "#E6E6E6", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#E6E6E6", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#BBBBBB", linewidth = 0.35) +
  geom_point(data = filter(D, !sig), colour = COLS[["ns"]], size = 0.3, alpha = 0.45) +
  geom_point(data = S, aes(colour = grp), size = 0.55, alpha = 0.55) +
  geom_point(data = HL, aes(fill = lab), shape = 21, size = 2.4, colour = "white", stroke = 0.35) +
  geom_text_repel(data = HL, aes(label = lab, colour = grp), size = 2.6, fontface = "bold",
                  segment.colour = "#999999", segment.size = 0.25, box.padding = 1.4, point.padding = 0.4, min.segment.length = 0,
                  show.legend = FALSE, colour = HL$col) +
  annotate("text", x = -7.7, y = 7.7, hjust = 0, vjust = 1, size = 3.1, colour = "#333333",
           label = sprintf("r = %.2f\n%.0f%% same direction", rS, conc)) +
  scale_colour_manual(values = COLS, breaks = c("shared_up","shared_dn","divergent"),
                      labels = c("shared ↑","shared ↓","divergent"), name = NULL) +
  scale_fill_manual(values = setNames(HL$col, HL$lab), guide = "none") +
  coord_fixed(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(x = expression("Warming effect, 27 vs 15 °C (log"[2]*"FC)"),
       y = expression("Peak prothioconazole effect across 15–27 °C (log"[2]*"FC)"),
       title = "Gene-level convergence of drug and warming") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_blank(), legend.text = element_text(size = 7.6),
        legend.key.height = unit(3.2, "mm")) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(ROOT, "figures/FigureS8.png"), fig, width = 5.4, height = 5.1, dpi = 600, bg = "white")
ggsave(file.path(ROOT, "figures/FigureS8.pdf"), fig, width = 5.4, height = 5.1, device = cairo_pdf, bg = "white")
cat(sprintf("FigS8: %d genes | r = %.2f | %.0f%% same direction\n", nrow(S), rS, conc))
