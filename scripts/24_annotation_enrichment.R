#!/usr/bin/env Rscript
###############################################################################
## annotate_and_enrich.R
##   (1) merge eggNOG functional annotation (name, description, KEGG, Pfam, EC)
##       into every DE_*.csv table
##   (2) hypergeometric enrichment per comparison for BOTH:
##         - GO terms
##         - KEGG pathways   (highlights ergosterol/steroid biosynthesis,
##                            the prothioconazole target pathway)
##   Self-contained: base R + ggplot2 only (no clusterProfiler).
##
##   Prereqs: run rnaseq_analysis.R (DE tables) and run_eggnog.sh (annotation).
##   Run:  conda activate rnaseq-r ; Rscript annotate_and_enrich.R
###############################################################################
suppressPackageStartupMessages(library(ggplot2))

cfg <- list(
  dir_de   = "tables/rnaseq",
  dir_fig  = "figures/diagnostics",
  anno_csv = "data/reference/gene_annotation.csv",
  go_csv   = "data/reference/go_terms.csv",
  kegg_csv = "data/reference/kegg_pathways.csv",
  alpha = 0.05, lfc = 1,          # DE-calling thresholds
  go_alpha = 0.05,                # enrichment significance
  min_set = 5, min_sig = 10, top_n = 20,
  ergosterol = c("map00100")      # steroid/ergosterol biosynthesis (highlight)
)
enr_csv_dir <- file.path(cfg$dir_de, "enrichment")
enr_fig_dir <- file.path(cfg$dir_fig, "enrichment")
dir.create(enr_csv_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(enr_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- load annotation ------------------------------------------------------
stopifnot(file.exists(cfg$anno_csv))
anno <- read.csv(cfg$anno_csv, stringsAsFactors = FALSE, colClasses = "character")
for (cc in c("gene_name","description","GO","KEGG","EC","PFAM"))
  if (!cc %in% names(anno)) anno[[cc]] <- ""
anno[is.na(anno)] <- ""

mapcol <- function(col) setNames(anno[[col]], anno$gene_id)
name_map <- mapcol("gene_name"); desc_map <- mapcol("description")
kegg_col <- mapcol("KEGG"); pfam_col <- mapcol("PFAM"); ec_col <- mapcol("EC")

split_set <- function(x) { x <- strsplit(ifelse(is.na(x), "", x), ";", fixed = TRUE)
                           lapply(x, function(v) v[nzchar(v)]) }
gene2go   <- setNames(split_set(anno$GO),   anno$gene_id)
gene2kegg <- setNames(split_set(anno$KEGG), anno$gene_id)

read_names <- function(path, idc, nmc) {
  if (!file.exists(path)) return(setNames(character(), character()))
  d <- read.csv(path, stringsAsFactors = FALSE, colClasses = "character")
  setNames(d[[nmc]], d[[idc]])
}
go_name   <- read_names(cfg$go_csv,   "go_id",   "go_name")
go_ns     <- read_names(cfg$go_csv,   "go_id",   "go_namespace")
kegg_name <- read_names(cfg$kegg_csv, "kegg_id", "kegg_name")

message("Annotation: ", nrow(anno), " genes; ",
        sum(lengths(gene2go) > 0), " with GO; ",
        sum(lengths(gene2kegg) > 0), " with KEGG.")

## ---- generic hypergeometric enrichment ------------------------------------
enrich <- function(sig, bg, gene2term, term_name, term_ns = NULL) {
  bg  <- bg[lengths(gene2term[bg]) > 0]
  sig <- sig[lengths(gene2term[sig]) > 0]
  if (length(sig) < cfg$min_sig) return(NULL)
  N <- length(bg); draw <- length(sig)
  bgT  <- table(unlist(gene2term[bg]))
  sigT <- table(unlist(gene2term[sig]))
  test <- names(bgT)[bgT >= cfg$min_set]
  res <- lapply(test, function(t) {
    K <- as.integer(bgT[[t]]); x <- if (t %in% names(sigT)) as.integer(sigT[[t]]) else 0L
    if (x == 0L) return(NULL)
    hits <- sig[vapply(gene2term[sig], function(g) t %in% g, logical(1))]
    data.frame(term = t,
               name = term_name[t] %||% "",
               namespace = if (is.null(term_ns)) "" else (term_ns[t] %||% ""),
               n_sig = x, n_bg = K, gene_ratio = x/draw,
               enrichment = (x/draw)/(K/N),
               pvalue = phyper(x-1L, K, N-K, draw, lower.tail = FALSE),
               genes = paste(hits, collapse = ";"), stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  if (is.null(res) || !nrow(res)) return(NULL)
  res$padj <- p.adjust(res$pvalue, "BH")
  res[order(res$padj, res$pvalue), ]
}
`%||%` <- function(a, b) if (is.null(a) || length(a)==0 || is.na(a) || a=="") b else a

dotplot <- function(sig_res, title, file) {
  pl <- head(sig_res[sig_res$name != "", ], cfg$top_n)
  if (!nrow(pl)) return(invisible())
  pl$label <- factor(pl$name, levels = rev(unique(pl$name)))
  p <- ggplot(pl, aes(gene_ratio, label, size = n_sig, color = padj)) +
    geom_point() +
    scale_color_gradient(low = "#D55E00", high = "#56B4E9", name = "padj") +
    scale_size_continuous(name = "DE genes") +
    labs(x = "Gene ratio (DE in set / DE total)", y = NULL, title = title) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  h <- min(9, 1.5 + 0.32 * nrow(pl))
  ggsave(paste0(file, ".pdf"), p, width = 8.5, height = h)
  ggsave(paste0(file, ".png"), p, width = 8.5, height = h, dpi = 300)
}

## ---- iterate over DE tables -----------------------------------------------
de_files <- list.files(cfg$dir_de, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("DE_summary", de_files)]
message("Found ", length(de_files), " DE tables.")

ergo_track <- list()   # collect ergosterol-pathway result across comparisons

for (f in de_files) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  comp <- if ("comparison" %in% names(df)) df$comparison[1] else sub("^DE_|\\.csv$","",basename(f))

  ## (1) annotate DE table in place
  df$gene_name   <- name_map[df$gene_id]
  df$description <- desc_map[df$gene_id]
  df$KEGG <- kegg_col[df$gene_id]; df$PFAM <- pfam_col[df$gene_id]; df$EC <- ec_col[df$gene_id]
  front <- c("gene_id","gene_name","description")
  df <- df[, c(front, setdiff(names(df), front))]
  write.csv(df, f, row.names = FALSE)

  bg  <- df$gene_id[!is.na(df$padj)]
  sig <- df$gene_id[!is.na(df$padj) & df$padj < cfg$alpha & abs(df$log2FoldChange) >= cfg$lfc]

  ## (2a) GO enrichment
  go_res <- enrich(sig, bg, gene2go, go_name, go_ns)
  if (!is.null(go_res)) {
    write.csv(go_res, file.path(enr_csv_dir, paste0("GO_", comp, ".csv")), row.names = FALSE)
    dotplot(go_res[go_res$padj < cfg$go_alpha, ], paste0("GO enrichment: ", comp),
            file.path(enr_fig_dir, paste0("GO_", comp)))
  }
  ## (2b) KEGG enrichment
  kg_res <- enrich(sig, bg, gene2kegg, kegg_name)
  if (!is.null(kg_res)) {
    write.csv(kg_res, file.path(enr_csv_dir, paste0("KEGG_", comp, ".csv")), row.names = FALSE)
    dotplot(kg_res[kg_res$padj < cfg$go_alpha, ], paste0("KEGG enrichment: ", comp),
            file.path(enr_fig_dir, paste0("KEGG_", comp)))
    er <- kg_res[kg_res$term %in% cfg$ergosterol, ]
    if (nrow(er)) ergo_track[[comp]] <- data.frame(comparison = comp, er[, c("name","n_sig","n_bg","enrichment","padj")])
  }
  ng <- sum(!is.na(df$padj) & df$padj < cfg$alpha & abs(df$log2FoldChange) >= cfg$lfc)
  message(sprintf("  %-34s  %d DE | GO sig: %s | KEGG sig: %s", comp, ng,
                  if (is.null(go_res)) 0 else sum(go_res$padj < cfg$go_alpha),
                  if (is.null(kg_res)) 0 else sum(kg_res$padj < cfg$go_alpha)))
}

## ---- ergosterol pathway across comparisons (the headline) -----------------
if (length(ergo_track)) {
  et <- do.call(rbind, ergo_track)
  write.csv(et, file.path(enr_csv_dir, "ergosterol_pathway_across_comparisons.csv"), row.names = FALSE)
  message("\nErgosterol/steroid-biosynthesis (KEGG map00100) across comparisons:")
  print(et, row.names = FALSE)
}
message("\nDONE.")
message("  Annotated DE tables : ", cfg$dir_de, "/DE_*.csv")
message("  Enrichment tables   : ", enr_csv_dir, " (GO_*.csv, KEGG_*.csv)")
message("  Enrichment figures  : ", enr_fig_dir)
