## R_figureS3_ergosterol.R — Supplementary Figure S3 (ergosterol pathway), ggplot2 port
## RUN FROM the master "X0123 copy" folder:  Rscript revision/scripts_R/R_figureS3_ergosterol.R
suppressMessages({library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(readr)})

TAB <- "tables/revision"; OUT <- "figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source("scripts/fig_style.R")
TEMPCOL <- TEMP3
conds <- c("15_0","15_2","21_0","21_2","27_0","27_2")
th <- theme_pub(9) +
  theme(plot.title = element_text(size = 10, face = "plain", hjust = 0))

erg <- read_csv(file.path(TAB,"ergosterol/erg_pathway_gene_table.csv"), show_col_types = FALSE)
sc  <- read_csv(file.path(TAB,"ergosterol/erg_pathway_scores_per_sample.csv"), show_col_types = FALSE)
vsd <- read.csv(file.path(TAB,"expression/vsd_matrix.csv"), row.names = 1, check.names = FALSE)
cd  <- read_csv(file.path(TAB,"expression/coldata.csv"), show_col_types = FALSE) |>
  mutate(cond = paste(as.integer(temperature), as.integer(prothioconazole), sep = "_"))

core <- erg |> filter(tier == "core enzyme", gene_id %in% rownames(vsd))
core$label <- ifelse(is.na(core$gene_name) | core$gene_name == "", core$gene_id, core$gene_name)
core$label <- paste0(core$label, " (", sub("Mycgr3G", "", core$gene_id), ")")
core$label <- factor(core$label, levels = rev(core$label))   # keep pathway order

## ---- A: group-mean VST z heatmap -------------------------------------------
m <- as.matrix(vsd[core$gene_id, ])
z <- t(scale(t(m)))
gm <- sapply(conds, function(cc) rowMeans(z[, cd$sample[cd$cond == cc], drop = FALSE]))
ha <- as.data.frame(gm) |> mutate(gene = core$label) |>
  pivot_longer(-gene, names_to = "cond", values_to = "z") |>
  mutate(cond = factor(cond, levels = conds))
pA <- ggplot(ha, aes(cond, gene, fill = z)) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1.6, 1.6),
                       oob = scales::squish) +
  scale_x_discrete(labels = c("15\u00b0C\n0","15\u00b0C\n2","21\u00b0C\n0","21\u00b0C\n2","27\u00b0C\n0","27\u00b0C\n2")) +
  labs(title = "A  Ergosterol-pathway expression (VST z, group mean)", x = NULL, y = NULL) +
  th + theme(axis.text.y = element_text(size = 5.5))

## ---- B: drug LFC per gene ---------------------------------------------------
hb <- core |> select(label, lfc_drug15, lfc_drug21, lfc_drug27) |>
  pivot_longer(-label, names_to = "t", values_to = "lfc") |>
  mutate(t = sub("lfc_drug", "", t))
pB <- ggplot(hb, aes(lfc, label, colour = t)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 1.1, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = TEMPCOL, name = NULL, labels = paste0(c(15,21,27), " \u00b0C")) +
  labs(title = "B  Drug effect per gene",
       x = expression("log"[2]*"FC (2 vs 0 mg L"^-1*", ashr)"), y = NULL) +
  th + theme(axis.text.y = element_blank(), legend.position = c(0.85, 0.9))

## ---- C: pathway score per sample -------------------------------------------
sc2 <- sc |> mutate(t = sub("_.*", "", cond), cond = factor(cond, levels = conds))
pC <- ggplot(sc2, aes(cond, erg_core, colour = t)) +
  geom_jitter(width = 0.08, size = 1.6) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.35, colour = "black") +
  scale_colour_manual(values = TEMPCOL, guide = "none") +
  scale_x_discrete(labels = sub("_", "/", conds)) +
  labs(title = "C  Pathway activity score per sample", x = NULL,
       y = "Ergosterol-pathway score (mean z)") + th

## ---- D: CYP51 vs AOX --------------------------------------------------------
picks <- c(`CYP51/ERG11` = "Mycgr3G110231", AOX = "Mycgr3G72918")
hd <- lapply(names(picks), function(nm) {
  g <- picks[[nm]]
  sapply(conds, function(cc) {
    v <- as.numeric(vsd[g, cd$sample[cd$cond == cc]])
    c(mean = mean(v), sem = sd(v)/sqrt(length(v)))
  }) |> t() |> as.data.frame() |> mutate(cond = conds, gene = nm)
}) |> bind_rows() |> mutate(cond = factor(cond, levels = conds))
pD <- ggplot(hd, aes(cond, mean, fill = gene)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = c(`CYP51/ERG11` = "#1B7837", AOX = "#E6820D"), name = NULL) +
  scale_x_discrete(labels = sub("_", "/", conds)) +
  labs(title = "D  CYP51 vs AOX", x = NULL, y = "VST expression") +
  th + theme(legend.position = c(0.25, 0.95))

## ---- E: interaction coefficient per gene ------------------------------------
he <- core |> mutate(sig = !is.na(padj_int27v15) & padj_int27v15 < 0.05)
pE <- ggplot(he, aes(lfc_int27v15_MLE, label, fill = sig)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"), guide = "none") +
  labs(title = "E  Temperature \u00d7 drug interaction",
       x = "Interaction LFC (27\u221215), red = padj<0.05", y = NULL) +
  th + theme(axis.text.y = element_text(size = 5))

layout <- (pA | pB) / (pC | pD | pE) + plot_layout(heights = c(1.35, 1)) +
  plot_annotation(title = "Ergosterol-pathway engagement across temperature \u00d7 prothioconazole",
                  theme = theme(plot.title = element_text(size = 11, hjust = 0.5)))
ggsave(file.path(OUT, "FigureS4.png"), layout, width = 13, height = 9, dpi = 200)
ggsave(file.path(OUT, "FigureS4.pdf"), layout, width = 13, height = 9, device = cairo_pdf)
cat("saved", file.path(OUT, "FigureS4.png"), "\n")
