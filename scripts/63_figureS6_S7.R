# =============================================================================
# 63_figureS6_S7.R -- Figures S6 and S7. R port of 33_figureS6_S7.py.
#   S6  full KEGG landscape of the prothioconazole response (GSEA FDR<0.05 in
#       >=1 temperature), cell = mean member-gene log2FC
#   S7  gene-level grid: genes with a temperature x drug interaction (LRT
#       FDR<0.05) and |log2FC| >= 4 in >=1 temperature
# Pathway selection and GSEA statistics are read from the cached tables written
# by the enrichment step, so no GSEA is re-run here.
# Inputs : tables/rnaseq/figureS1_pathways.csv, DE_proth_at_{15,21,27}C_2vs0.csv,
#          DE_interaction_LRT_any.csv, data/reference/gene_annotation.csv
# Outputs: figures/FigureS6.png/.pdf, figures/FigureS7.png/.pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/fig_common.R")); anno_load(ROOT)
TEMPS <- c("15","21","27")
stars <- function(q) ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "")))
dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)

## ============================ Figure S6 =====================================
PS <- read.csv(file.path(ROOT, "tables/rnaseq/figureS1_pathways.csv"), stringsAsFactors = FALSE)
PS <- PS[order(PS$meanLFC_15), ]
PS$lab <- factor(PS$kegg_name, levels = PS$kegg_name)
H6 <- PS %>% select(lab, meanLFC_15, meanLFC_21, meanLFC_27, FDR_15, FDR_21, FDR_27) %>%
  pivot_longer(starts_with("meanLFC_"), names_to = "k", values_to = "lfc") %>%
  mutate(temp = sub("meanLFC_", "", k),
         fdr = case_when(temp == "15" ~ FDR_15, temp == "21" ~ FDR_21, TRUE ~ FDR_27),
         txt = sprintf("%+.2f%s", lfc, stars(fdr)),
         temp = factor(paste0(temp, " °C"), levels = paste0(TEMPS, " °C")))
pS6 <- ggplot(H6, aes(temp, lab, fill = lfc)) +
  geom_tile() +
  geom_text(aes(label = txt, colour = abs(lfc) > 0.9), size = 1.9, show.legend = FALSE) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       limits = c(-1.5, 1.5), oob = scales::squish,
                       name = expression("mean member log"[2]*"FC")) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#222222")) +
  scale_y_discrete(limits = rev(levels(PS$lab))) +
  labs(x = expression(atop("Prothioconazole effect (2 vs 0 mg L"^-1*")",
                           "cell = mean member log"[2]*"FC; * FDR<0.05  ** <0.01  *** <0.001")), y = NULL,
       title = "KEGG pathways significant (GSEA FDR<0.05) in ≥1 temperature\nrestricted universe; rows ordered by effect direction") +
  theme_classic(base_size = 9) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7),
        plot.title = element_text(size = 8.5, face = "bold"),
        legend.key.width = unit(3, "mm"), legend.key.height = unit(9, "mm"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))
ggsave(file.path(ROOT, "figures/FigureS6.png"), pS6, width = 6.8, height = 12.5, dpi = 400, bg = "white", limitsize = FALSE)
ggsave(file.path(ROOT, "figures/FigureS6.pdf"), pS6, width = 6.8, height = 12.5, device = cairo_pdf, bg = "white", limitsize = FALSE)

## ============================ Figure S7 =====================================
DE <- lapply(TEMPS, function(t) load_de(paste0("proth_at_", t, "C_2vs0"), ROOT)); names(DE) <- TEMPS
LRT <- read.csv(file.path(ROOT, "tables/rnaseq/DE_interaction_LRT_any.csv"), stringsAsFactors = FALSE)
intp <- setNames(LRT$padj, LRT$gene_id)
L <- bind_rows(lapply(TEMPS, function(t) DE[[t]] %>% select(gene, lfc) %>% mutate(temp = t)))
W <- L %>% group_by(gene) %>% filter(n() == 3) %>%
  summarise(mx = max(abs(lfc)), dom = lfc[which.max(abs(lfc))], .groups = "drop") %>%
  filter(mx >= 4, !is.na(intp[gene]), intp[gene] < 0.05)
W <- W %>% mutate(up = dom > 0) %>% arrange(desc(up), desc(dom))

anno <- .anno_cache$anno; rownames(anno) <- anno$gene_id
UNIPROT <- c(Mycgr3G66678="luciferase domain (66678)", Mycgr3G18811="DUF6699 domain (18811)",
  Mycgr3G97077="secreted (97077)", Mycgr3G106329="secreted (106329)", Mycgr3G67799="secreted (67799)",
  Mycgr3G108482="secreted (108482)", Mycgr3G103091="secreted (103091)", Mycgr3G102617="secreted (102617)",
  Mycgr3G92048="secreted (92048)", Mycgr3G104082="membrane (104082)", Mycgr3G108147="membrane (108147)",
  Mycgr3G99182="membrane (99182)", Mycgr3G94648="coiled-coil (94648)")
shorten <- function(d) {
  for (x in c("Belongs to the ","Catalyzes ","Domain of unknown function","Partial ","Putative "))
    d <- gsub(x, "", d, fixed = TRUE)
  d <- gsub(" superfamily", "", d, fixed = TRUE); d <- gsub(" family", "", d, fixed = TRUE)
  substr(trimws(sub("\\.\\s.*$", "", d)), 1, 26)
}
wlabel <- function(g) {
  nm <- gsym(g)
  if (!is.na(nm)) return(nm)
  de <- anno[g, "description"]
  if (!is.na(de) && nzchar(de) && !grepl("uncharacter", tolower(de))) return(tolower(shorten(de)))
  if (g %in% names(UNIPROT)) return(unname(UNIPROT[g]))
  g
}
W$lab <- make.unique(vapply(W$gene, wlabel, ""))
W$lab <- factor(W$lab, levels = W$lab)
H7 <- L %>% filter(gene %in% W$gene) %>% left_join(W[, c("gene","lab")], by = "gene") %>%
  mutate(temp = factor(paste0(temp, " °C"), levels = paste0(TEMPS, " °C")))
nup <- sum(W$up)
pS7 <- ggplot(H7, aes(temp, lab, fill = lfc)) +
  geom_tile() +
  geom_hline(yintercept = nrow(W) + 0.5 - nup, colour = "black", linewidth = 0.35) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       limits = c(-8, 8), oob = scales::squish,
                       name = expression("Prothioconazole effect (log"[2]*"FC)")) +
  scale_y_discrete(limits = rev(levels(W$lab))) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 4.6), axis.text.x = element_text(size = 8),
        legend.key.width = unit(3, "mm"), legend.key.height = unit(9, "mm"),
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6))
hh <- max(9, nrow(W) * 0.21)
ggsave(file.path(ROOT, "figures/FigureS7.png"), pS7, width = 6.4, height = hh, dpi = 400, bg = "white", limitsize = FALSE)
ggsave(file.path(ROOT, "figures/FigureS7.pdf"), pS7, width = 6.4, height = hh, device = cairo_pdf, bg = "white", limitsize = FALSE)
cat(sprintf("FigS6: %d pathways | FigS7: %d genes (up=%d down=%d)\n", nrow(PS), nrow(W), nup, nrow(W) - nup))
