# =============================================================================
# 53_figure3.R -- Figure 3: the prothioconazole transcriptional response and its
# temperature dependence. R port of 78_figure3_redesign.py (values unchanged).
#   A  PCA of the 500 most variable genes
#   B  volcano panels at 15 / 21 / 27 C, KEGG mechanistic categories coloured
#   C  pathway heatmap, mean member-gene log2FC, * = cellwise GSEA FDR < 0.05
#   D  OXPHOS composite and the single AOX gene, per sample with 95% CI
# Inputs : data/rnaseq/{sample_metadata,raw_count_matrix}.csv
#          tables/rnaseq/DE_proth_at_{15,21,27}C_2vs0.csv, figure4C_pathways.csv
#          data/reference/gene_annotation.csv
# Outputs: figures/Figure3.png / .pdf
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)
                                library(patchwork); library(ggrepel)})
ARGS <- commandArgs(trailingOnly = TRUE); ROOT <- if (length(ARGS)) ARGS[1] else "."
source(file.path(ROOT, "scripts/fig_common.R")); anno_load(ROOT)
set.seed(1)

meta <- read.csv(file.path(ROOT, "data/rnaseq/sample_metadata.csv"), stringsAsFactors = FALSE)
meta <- meta[tolower(as.character(meta$is_control)) == "false", ]
meta$temperature <- as.character(as.integer(meta$temperature))
meta$prothioconazole <- as.character(as.integer(meta$prothioconazole))
cnt <- read.csv(file.path(ROOT, "data/rnaseq/raw_count_matrix.csv"), row.names = 1, check.names = FALSE)
colnames(cnt) <- sub("^.*\\((.*)\\)$", "\\1", colnames(cnt))
cnt <- as.matrix(cnt[, meta$sample])
norm <- norm_counts(cnt); logn <- log2(norm + 1)

## ---- A: PCA ----------------------------------------------------------------
X  <- t(logn)                                        # samples x genes
v  <- apply(X, 2, var); idx <- order(v, decreasing = TRUE)[1:500]
Xt <- scale(X[, idx], center = TRUE, scale = FALSE)
sv <- svd(Xt); pc <- sv$u[, 1:2] %*% diag(sv$d[1:2]); pv <- (sv$d^2 / sum(sv$d^2))[1:2] * 100
if (mean(pc[meta$temperature == "27", 1]) < mean(pc[meta$temperature == "15", 1])) pc[, 1] <- -pc[, 1]
if (mean(pc[meta$prothioconazole == "2", 2]) < mean(pc[meta$prothioconazole == "0", 2])) pc[, 2] <- -pc[, 2]
P <- data.frame(PC1 = pc[, 1], PC2 = pc[, 2], temp = meta$temperature, dose = meta$prothioconazole)

ellipse_df <- function(x, y, level = 2.30) {          # same chi-sq radius as the Python
  if (length(x) < 3) return(NULL)
  e <- eigen(cov(cbind(x, y))); r <- sqrt(e$values * level)
  th <- seq(0, 2*pi, length.out = 120)
  p <- e$vectors %*% rbind(r[1]*cos(th), r[2]*sin(th))
  data.frame(x = p[1, ] + mean(x), y = p[2, ] + mean(y))
}
EL <- do.call(rbind, lapply(split(P, list(P$temp, P$dose), drop = TRUE), function(g) {
  e <- ellipse_df(g$PC1, g$PC2); if (is.null(e)) return(NULL)
  cbind(e, temp = g$temp[1], dose = g$dose[1], grp = paste0(g$temp[1], g$dose[1])) }))

pA <- ggplot(P, aes(PC1, PC2)) +
  geom_path(data = EL, aes(x, y, colour = temp, group = grp, linetype = dose),
            linewidth = 0.32, alpha = 0.75, inherit.aes = FALSE) +
  geom_point(aes(colour = temp, shape = dose), size = 2.0, stroke = 0.25) +
  scale_colour_manual(values = TEMP, name = "Temperature", labels = paste0(c(15,21,27), " °C")) +
  scale_shape_manual(values = c(`0` = 16, `2` = 17), name = "Prothioconazole",
                     labels = c(expression(0~mg~L^-1), expression(2~mg~L^-1))) +
  scale_linetype_manual(values = c(`0` = "solid", `2` = "42"), guide = "none") +
  labs(x = sprintf("PC1 (%.0f%%) — temperature", pv[1]),
       y = sprintf("PC2 (%.0f%%) — prothioconazole", pv[2])) +
  theme_classic(base_size = 8.5) +
  theme(legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(0, "mm"),
        legend.margin = margin(0,0,0,0), legend.key.height = unit(3, "mm"),
        legend.text = element_text(size = 6.6), legend.title = element_text(size = 7.2),
        axis.text = element_text(size = 7.5))

## ---- B: volcano panels ------------------------------------------------------
FIXED <- c(Mycgr3G110231 = "CYP51/ERG11 (ns)", Mycgr3G72918 = "AOX",
           Mycgr3G50464 = "CDR1", Mycgr3G108724 = "RPL7")
VB <- lapply(c("15","21","27"), function(t) volcano_facet(load_de(paste0("proth_at_", t, "C_2vs0"), ROOT),
                                                          extra = FIXED, label_per_cat = 0))
names(VB) <- c("15","21","27")
DEV <- bind_rows(lapply(names(VB), function(t) VB[[t]]$de %>% mutate(temp = t)))
LAB <- bind_rows(lapply(names(VB), function(t) VB[[t]]$lab %>% mutate(temp = t)))
HDR <- bind_rows(lapply(names(VB), function(t) data.frame(temp = t,
          hdr = sprintf("%s ↑   %s ↓", format(VB[[t]]$n_up, big.mark = ","), format(VB[[t]]$n_dn, big.mark = ",")),
          tot = sprintf("%s DEGs", format(VB[[t]]$n_tot, big.mark = ",")))))
DEV$temp <- factor(DEV$temp, levels = c("15","21","27"), labels = paste0(c(15,21,27), " °C"))
LAB$temp <- factor(LAB$temp, levels = c("15","21","27"), labels = paste0(c(15,21,27), " °C"))
HDR$temp <- factor(HDR$temp, levels = c("15","21","27"), labels = paste0(c(15,21,27), " °C"))
present <- unique(na.omit(DEV$cat[DEV$sig]))
CC <- cat_cols()

pB <- ggplot() +
  geom_point(data = filter(DEV, !sig), aes(lfc, y), colour = "#ECECEC", size = 0.28, alpha = 0.55) +
  geom_point(data = filter(DEV, sig, is.na(cat)), aes(lfc, y), colour = "#BCC2C9", size = 0.32, alpha = 0.5) +
  geom_point(data = filter(DEV, sig, !is.na(cat)), aes(lfc, y, colour = cat), size = 0.75, alpha = 0.82) +
  geom_vline(xintercept = c(-1, 1), linetype = "33", colour = "#D4D4D4", linewidth = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = "33", colour = "#D4D4D4", linewidth = 0.25) +
  geom_point(data = LAB, aes(lfc, y, colour = cat), size = 1.5, shape = 21, fill = NA, stroke = 0.35, show.legend = FALSE) +
  geom_text_repel(data = LAB, aes(lfc, y, label = txt, colour = cat), size = 2.1, fontface = "bold.italic",
                  segment.colour = "#9A9A9A", segment.size = 0.2, min.segment.length = 0,
                  box.padding = 0.35, max.overlaps = 30, show.legend = FALSE) +
  geom_text(data = HDR, aes(0, 51, label = hdr), size = 2.3, colour = "#333333", vjust = 1) +
  geom_text(data = HDR, aes(0, 47, label = tot), size = 2.6, fontface = "bold", colour = "#111111", vjust = 1) +
  facet_wrap(~temp, nrow = 1) +
  scale_colour_manual(values = CC, breaks = present, name = NULL, na.value = "#555555") +
  coord_cartesian(xlim = c(-8, 8), ylim = c(-2, 52)) +
  scale_x_continuous(breaks = c(-6, -3, 0, 3, 6)) +
  labs(x = expression("Prothioconazole effect (log"[2]*"FC, 2 vs 0 mg L"^-1*")"),
       y = expression(-log[10]~"adj "*italic(P))) +
  theme_classic(base_size = 8.5) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 7.5), legend.position = "bottom",
        legend.text = element_text(size = 6.2), legend.key.height = unit(2.6, "mm"),
        panel.spacing = unit(2.2, "mm")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6)))

## ---- C: pathway heatmap -----------------------------------------------------
PC <- read.csv(file.path(ROOT, "tables/rnaseq/figure4C_pathways.csv"), stringsAsFactors = FALSE)
GRP <- c(map03010=0,map03020=0,map03008=0,map03040=0,map03013=0,map03030=0,
         map03018=0,map03015=0,map03060=0, map00020=1,map00190=1,map00240=1,map00230=1,
         map00980=2,map00071=2)
PC$grp <- GRP[PC$kegg_id]
PC <- PC %>% filter(!is.na(grp)) %>% arrange(grp, meanLFC_15)
DISP <- c("Metabolism of xenobiotics by cytochrome P450" = "Xenobiotic metabolism (P450)",
          "Ribosome biogenesis in eukaryotes" = "Ribosome biogenesis")
PC$lab <- ifelse(PC$kegg_name %in% names(DISP), DISP[PC$kegg_name], PC$kegg_name)
PC$lab <- factor(PC$lab, levels = PC$lab)
H <- PC %>% select(lab, grp, meanLFC_15, meanLFC_21, meanLFC_27, FDR_15, FDR_21, FDR_27) %>%
  pivot_longer(c(meanLFC_15, meanLFC_21, meanLFC_27), names_to = "k", values_to = "lfc") %>%
  mutate(temp = sub("meanLFC_", "", k),
         fdr = ifelse(temp == "15", FDR_15, ifelse(temp == "21", FDR_21, FDR_27)),
         txt = sprintf("%+.2f%s", lfc, ifelse(fdr < 0.05, "*", "")),
         temp = factor(paste0(temp, " °C"), levels = paste0(c(15,21,27), " °C")))
bnd <- nrow(PC) + 1 - (which(diff(PC$grp) != 0) + 0.5)   # y scale is reversed
pC <- ggplot(H, aes(temp, lab, fill = lfc)) +
  geom_tile() +
  geom_text(aes(label = txt, colour = abs(lfc) > 0.7), size = 1.95, show.legend = FALSE) +
  geom_hline(yintercept = bnd, colour = "white", linewidth = 0.85) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       limits = c(-1.3, 1.3), oob = scales::squish,
                       breaks = c(-1.3, 0, 1.3), labels = c("−1.3", "0", "+1.3"),
                       name = expression("mean member log"[2]*"FC")) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#222222")) +
  scale_y_discrete(limits = rev(levels(PC$lab))) +
  labs(x = expression(atop("Prothioconazole effect (2 vs 0 mg L"^-1*"), mean member-gene log"[2]*"FC",
                           "rows selected by GSEA FDR<0.001 at ">=1*" temperature; * cellwise GSEA FDR<0.05")),
       y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 6.6), axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 6.6),
        legend.key.width = unit(2.4, "mm"), legend.key.height = unit(7, "mm"),
        legend.title = element_text(size = 6.8), legend.text = element_text(size = 6.4))

## ---- D: OXPHOS composite and AOX -------------------------------------------
ox <- log2(colMeans(norm[intersect(kg("map00190"), rownames(norm)), , drop = FALSE]) + 1)
ax <- log2(norm[aoxg()[1], ] + 1)
mk <- function(vals, ttl, ylab) {
  D <- data.frame(v = as.numeric(vals), temp = meta$temperature, dose = meta$prothioconazole)
  S <- D %>% group_by(temp, dose) %>% summarise(m = mean(v), e = ci95(v), .groups = "drop")
  ggplot(D, aes(temp, v, colour = dose)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.22, dodge.width = 0.62),
               size = 0.75, alpha = 0.55) +
    geom_errorbar(data = S, aes(temp, ymin = m - e, ymax = m + e, colour = dose),
                  width = 0.16, linewidth = 0.42, position = position_dodge(0.62), inherit.aes = FALSE) +
    geom_point(data = S, aes(temp, m, colour = dose), shape = 95, size = 5,
               position = position_dodge(0.62), inherit.aes = FALSE) +
    scale_colour_manual(values = c(`0` = "#9A9A9A", `2` = "#433D84"), name = NULL,
                        labels = c("Control", "Prothioconazole")) +
    scale_x_discrete(labels = paste0(c(15,21,27), "°C")) +
    labs(x = NULL, y = ylab, title = ttl) +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 7.5), legend.position = "none")
}
pD1 <- mk(ox, sprintf("OXPHOS composite\n(mean of %d genes)", length(kg("map00190"))),
          expression(atop("log"[2]*"(mean normalized", "expression + 1)")))
pD2 <- mk(ax, "Alternative oxidase\n(AOX, single gene)",
          expression(atop("log"[2]*"(normalized", "expression + 1)"))) +
  theme(legend.position = c(0.03, 0.99), legend.justification = c(0, 1),
        legend.text = element_text(size = 6.0), legend.key.height = unit(2.6, "mm"),
        legend.key.width = unit(2.6, "mm"), legend.background = element_blank())

## ---- assemble ---------------------------------------------------------------
fig <- (pA + pB + plot_layout(widths = c(1, 1.95))) /
       (pC + (pD1 | pD2) + plot_layout(widths = c(1.15, 1))) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = list(c("A","B","C","D",""))) &
  theme(plot.tag = element_text(size = 13, face = "bold"), plot.tag.position = c(0, 1.02))

dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(ROOT, "figures/Figure3.png"), fig, width = 9.6, height = 7.35, dpi = 400, bg = "white")
ggsave(file.path(ROOT, "figures/Figure3.pdf"), fig, width = 9.6, height = 7.35, device = cairo_pdf, bg = "white")
cat(sprintf("Figure3: PC1 %.0f%% PC2 %.0f%% | DEGs %s\n", pv[1], pv[2],
    paste(sapply(VB, function(v) v$n_tot), collapse = "/")))
