# =============================================================================
# 54_figure4.R -- Figure 4: the sole effect of temperature on the transcriptome.
# R port of 79_figure4_redesign.py (values unchanged; panel B reads the cached
# GSEA table, exactly as the Python does when that cache exists).
#   A  volcano panels for the three temperature steps, no prothioconazole
#   B  pathway heatmap, mean member-gene log2FC, * = cellwise GSEA FDR < 0.05
#   C  OXPHOS composite and the single AOX gene, control samples only
# Inputs : tables/rnaseq/DE_temp_{21vs15,27vs21,27vs15}_proth0.csv
#          tables/revision/fig4/fig4B_gsea_verify.csv
#          data/rnaseq/{sample_metadata,raw_count_matrix}.csv
# Outputs: figures/Figure4.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)
                                library(patchwork); library(ggrepel)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/fig_common.R")); anno_load(ROOT)
set.seed(1)

STEPS <- c(`15→21 °C` = "temp_21vs15_proth0", `21→27 °C` = "temp_27vs21_proth0",
           `15→27 °C` = "temp_27vs15_proth0")
FIXED <- c(Mycgr3G72918 = "AOX", Mycgr3G50464 = "CDR1", Mycgr3G108724 = "RPL7")

## ---- A: volcano panels ------------------------------------------------------
VB <- lapply(seq_along(STEPS), function(i)
  volcano_facet(load_de(STEPS[i], ROOT), extra = FIXED, xlim = c(-9, 9), label_per_cat = 0))
names(VB) <- names(STEPS)
DEV <- bind_rows(lapply(names(VB), function(t) VB[[t]]$de %>% mutate(step = t)))
LAB <- bind_rows(lapply(names(VB), function(t) VB[[t]]$lab %>% mutate(step = t)))
HDR <- bind_rows(lapply(names(VB), function(t) data.frame(step = t,
        hdr = sprintf("%s ↑   %s ↓", format(VB[[t]]$n_up, big.mark = ","), format(VB[[t]]$n_dn, big.mark = ",")),
        tot = sprintf("%s DEGs", format(VB[[t]]$n_tot, big.mark = ",")))))
for (d in c("DEV","LAB","HDR")) { x <- get(d); x$step <- factor(x$step, levels = names(STEPS)); assign(d, x) }
LAB_TXT <- LAB %>% filter(step != names(STEPS)[1])       # first step: rings only, no text
present <- unique(na.omit(DEV$cat[DEV$sig])); CC <- cat_cols()

pA <- ggplot() +
  geom_point(data = filter(DEV, !sig), aes(lfc, y), colour = "#ECECEC", size = 0.28, alpha = 0.55) +
  geom_point(data = filter(DEV, sig, is.na(cat)), aes(lfc, y), colour = "#BCC2C9", size = 0.32, alpha = 0.5) +
  geom_point(data = filter(DEV, sig, !is.na(cat)), aes(lfc, y, colour = cat), size = 0.75, alpha = 0.82) +
  geom_vline(xintercept = c(-1, 1), linetype = "33", colour = "#D4D4D4", linewidth = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = "33", colour = "#D4D4D4", linewidth = 0.25) +
  geom_point(data = LAB, aes(lfc, y, colour = cat), size = 1.5, shape = 21, fill = NA,
             stroke = 0.35, show.legend = FALSE) +
  geom_text_repel(data = LAB_TXT, aes(lfc, y, label = txt, colour = cat), size = 2.1,
                  fontface = "bold.italic", segment.colour = "#9A9A9A", segment.size = 0.2,
                  min.segment.length = 0, box.padding = 0.35, max.overlaps = 30, show.legend = FALSE) +
  geom_text(data = HDR, aes(0, 51, label = hdr), size = 2.3, colour = "#333333", vjust = 1) +
  geom_text(data = HDR, aes(0, 47, label = tot), size = 2.6, fontface = "bold", colour = "#111111", vjust = 1) +
  facet_wrap(~step, nrow = 1) +
  scale_colour_manual(values = CC, breaks = present, name = NULL, na.value = "#555555") +
  coord_cartesian(xlim = c(-9, 9), ylim = c(-2, 52)) +
  scale_x_continuous(breaks = c(-6, -3, 0, 3, 6)) +
  labs(x = expression("Temperature effect (log"[2]*"FC, no prothioconazole)"),
       y = expression(-log[10]~"adj "*italic(P))) +
  theme_classic(base_size = 8.5) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 7.5), legend.position = "bottom",
        legend.text = element_text(size = 6.2), legend.key.height = unit(2.6, "mm"),
        panel.spacing = unit(2.2, "mm")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6)))

## ---- B: pathway heatmap (cached GSEA selection, order preserved) ------------
VC <- read.csv(file.path(ROOT, "tables/revision/fig4/fig4B_gsea_verify.csv"), stringsAsFactors = FALSE)
VC$lab <- factor(VC$name, levels = VC$name)
H <- VC %>%
  select(lab, mean_21, mean_2721, mean_27, FDR_21, FDR_2721, FDR_27) %>%
  pivot_longer(starts_with("mean_"), names_to = "k", values_to = "lfc") %>%
  mutate(step = recode(sub("mean_", "", k), `21` = "15→21 °C", `2721` = "21→27 °C", `27` = "15→27 °C"),
         fdr = case_when(k == "mean_21" ~ FDR_21, k == "mean_2721" ~ FDR_2721, TRUE ~ FDR_27),
         txt = sprintf("%+.2f%s", lfc, ifelse(fdr < 0.05, "*", "")),
         step = factor(step, levels = names(STEPS)))
GIP <- c("map03010","map03020","map03008","map03040","map03013","map03030","map03018",
         "map03015","map03060","map03050","map03022")
bnd <- nrow(VC) + 1 - (which(diff(as.integer(!VC$kegg_id %in% GIP)) != 0) + 0.5)  # y scale is reversed
pB <- ggplot(H, aes(step, lab, fill = lfc)) +
  geom_tile() +
  geom_text(aes(label = txt, colour = abs(lfc) > 0.7), size = 1.9, show.legend = FALSE) +
  geom_hline(yintercept = bnd, colour = "white", linewidth = 0.85) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       limits = c(-1.3, 1.3), oob = scales::squish,
                       breaks = c(-1.3, 0, 1.3), labels = c("−1.3", "0", "+1.3"),
                       name = expression("mean member log"[2]*"FC")) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#222222")) +
  scale_y_discrete(limits = rev(levels(VC$lab))) +
  labs(x = "Sole effect of temperature; rows selected by GSEA FDR<0.001 at ≥1 step; * cellwise GSEA FDR<0.05",
       y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 6.6), axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 6.4),
        legend.key.width = unit(2.4, "mm"), legend.key.height = unit(7, "mm"),
        legend.title = element_text(size = 6.8), legend.text = element_text(size = 6.4))

## ---- C: control-only temperature series ------------------------------------
meta <- read.csv(file.path(ROOT, "data/rnaseq/sample_metadata.csv"), stringsAsFactors = FALSE)
meta <- meta[tolower(as.character(meta$is_control)) == "false", ]
meta$temperature <- as.character(as.integer(meta$temperature))
meta$prothioconazole <- as.character(as.integer(meta$prothioconazole))
cnt <- read.csv(file.path(ROOT, "data/rnaseq/raw_count_matrix.csv"), row.names = 1, check.names = FALSE)
colnames(cnt) <- sub("^.*\\((.*)\\)$", "\\1", colnames(cnt))
cnt <- as.matrix(cnt[, meta$sample]); norm <- norm_counts(cnt)
ctl <- meta$prothioconazole == "0"
tp <- function(genes, ttl, ylab) {
  g <- intersect(genes, rownames(norm))
  sc <- log2(colMeans(norm[g, , drop = FALSE]) + 1)
  D <- data.frame(v = as.numeric(sc)[ctl], temp = meta$temperature[ctl])
  S <- D %>% group_by(temp) %>% summarise(m = mean(v), e = ci95(v), .groups = "drop")
  ggplot(D, aes(temp, v, colour = temp)) +
    geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.8, alpha = 0.55) +
    geom_errorbar(data = S, aes(temp, ymin = m - e, ymax = m + e), width = 0.14,
                  linewidth = 0.42, inherit.aes = FALSE, colour = TEMP[S$temp]) +
    geom_point(data = S, aes(temp, m), shape = 95, size = 6, inherit.aes = FALSE, colour = TEMP[S$temp]) +
    scale_colour_manual(values = TEMP, guide = "none") +
    scale_x_discrete(labels = paste0(c(15, 21, 27), "°C")) +
    labs(x = NULL, y = ylab, title = ttl) +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 7.5))
}
pC1 <- tp(kg("map00190"), sprintf("OXPHOS composite\n(mean of %d genes)", length(kg("map00190"))),
          expression(atop("log"[2]*"(mean normalized", "expression + 1)")))
pC2 <- tp(aoxg(), "Alternative oxidase\n(AOX, single gene)",
          expression(atop("log"[2]*"(normalized", "expression + 1)")))

fig <- pA / (pB + (pC1 | pC2) + plot_layout(widths = c(1.15, 1))) +
  plot_layout(heights = c(1, 1.10)) +
  plot_annotation(tag_levels = list(c("A","B","C",""))) &
  theme(plot.tag = element_text(size = 13, face = "bold"), plot.tag.position = c(0, 1.02))

dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(ROOT, "figures/Figure4.png"), fig, width = 10.2, height = 7.35, dpi = 400, bg = "white")
ggsave(file.path(ROOT, "figures/Figure4.pdf"), fig, width = 10.2, height = 7.35, device = cairo_pdf, bg = "white")
cat(sprintf("Figure4: DEGs %s | %d pathways\n",
    paste(sapply(VB, function(v) v$n_tot), collapse = "/"), nrow(VC)))
