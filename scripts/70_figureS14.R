# =============================================================================
# 70_figureS14.R -- Figure S14: pathway-level temperature x drug interaction
# (GSEA on the DESeq2 interaction statistic, 27 vs 15 C contrast).
# Port of 82_figS_interaction.py.
#   positive NES = drug repression attenuated at 27 C (gene-expression machinery)
#   negative NES = drug induction attenuated at 27 C (detox / degradation)
# Inputs : tables/revision/inference/interaction_gsea_27v15.csv
#          tables/rnaseq/figure4C_pathways.csv
# Outputs: figures/FigureS14.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/fig_style.R"))

res <- read.csv(file.path(ROOT, "tables/revision/inference/interaction_gsea_27v15.csv"))
sig <- res %>% filter(FDR < 0.05) %>% arrange(NES)

DISP <- c("Metabolism of xenobiotics by cytochrome P450" = "Xenobiotic metabolism (P450)",
          "Ribosome biogenesis in eukaryotes"            = "Ribosome biogenesis",
          "Chloroalkane and chloroalkene degradation"    = "Chloroalkane degradation",
          "Pentose and glucuronate interconversions"     = "Pentose/glucuronate interconv.",
          "Valine, leucine and isoleucine degradation"   = "BCAA degradation")
nm  <- ifelse(is.na(sig$name) | sig$name == "", sig$kegg_id, sig$name)
sig$lab <- ifelse(nm %in% names(DISP), DISP[nm], nm)
sig$lab <- factor(sig$lab, levels = sig$lab)
sig$dir <- ifelse(sig$NES > 0, "up", "dn")

fig <- ggplot(sig, aes(NES, lab, fill = dir)) +
  geom_col(width = 0.66, alpha = 0.85) +
  geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.3) +
  annotate("text", x = max(sig$NES), y = 1.4, hjust = 1, vjust = 0, size = 2.5,
           colour = DIV_HIGH, label = "positive: drug repression\nattenuated at 27 °C\n(gene-expression machinery)") +
  annotate("text", x = min(sig$NES), y = nrow(sig) - 0.4, hjust = 0, vjust = 1, size = 2.5,
           colour = DIV_LOW, label = "negative: drug induction\nattenuated at 27 °C\n(detox / degradation)") +
  scale_x_continuous(expand = expansion(mult = 0.06)) +
  scale_fill_manual(values = c(up = DIV_HIGH, dn = DIV_LOW), guide = "none") +
  labs(x = "Interaction NES  (temperature × drug, 27 vs 15 °C contrast)", y = NULL,
       title = sprintf("Pathway-level temperature × drug interaction\n(GSEA, FDR < 0.05; n = %d of %d pathways)",
                       nrow(sig), nrow(res))) +
  theme_pub(base_size = 8) +
  theme(plot.title = element_text(size = 8.2, face = "bold", hjust = 0),
        axis.text.y = element_text(size = 7))

dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(ROOT, "figures/FigureS14.png"), fig, width = 6.6, height = 6.2, dpi = 600, bg = "white")
ggsave(file.path(ROOT, "figures/FigureS14.pdf"), fig, width = 6.6, height = 6.2, device = cairo_pdf, bg = "white")
cat(sprintf("FigS14: %d of %d pathways at FDR<0.05\n", nrow(sig), nrow(res)))
