# =============================================================================
# fig_common.R -- shared annotation, KEGG categories, DE loaders and volcano
# helper for the transcriptomics figures (R port of fig_common.py).
# Values (categories, colours, thresholds) are identical to the Python version.
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggrepel)})

DE_DIR <- "tables/rnaseq"
TEMP   <- c(`15` = "#27519E", `21` = "#64A5DE", `27` = "#D64B21")
UP <- "#B2182B"; DN <- "#2166AC"; DIV <- "#E08214"
LFC_MIN <- 1; FDR_MAX <- 0.05

.anno_cache <- new.env(parent = emptyenv())

anno_load <- function(root = ".") {
  if (!is.null(.anno_cache$anno)) return(invisible(NULL))
  a <- read.csv(file.path(root, "data/reference/gene_annotation.csv"),
                stringsAsFactors = FALSE, colClasses = "character")
  .anno_cache$anno <- a
  kk <- strsplit(ifelse(is.na(a$KEGG), "", a$KEGG), ";", fixed = TRUE)
  names(kk) <- a$gene_id
  .anno_cache$gk <- lapply(kk, function(v) v[nzchar(v)])
  # kegg id -> member genes
  flat <- data.frame(gene = rep(a$gene_id, lengths(.anno_cache$gk)),
                     kid  = unlist(.anno_cache$gk, use.names = FALSE),
                     stringsAsFactors = FALSE)
  .anno_cache$kegg_sets <- split(flat$gene, flat$kid)
  .anno_cache$aoxg <- a$gene_id[grepl("alternative oxidase", tolower(a$description))]
  .anno_cache$sym  <- setNames(trimws(ifelse(is.na(a$gene_name), "", a$gene_name)), a$gene_id)
  invisible(NULL)
}
kg     <- function(mid) { anno_load(); s <- .anno_cache$kegg_sets[[mid]]; if (is.null(s)) character(0) else s }
aoxg   <- function()    { anno_load(); .anno_cache$aoxg }
gsym   <- function(g)   { anno_load(); s <- .anno_cache$sym[g]
                          ifelse(!is.na(s) & nzchar(s) & !toupper(s) %in% c("NA","NAN"), s, NA_character_) }

cats <- function() {
  anno_load()
  list(
    list(lab = "Efflux (ABC transporters)",   col = "#7B3FA0", set = kg("map02010")),
    list(lab = "Ergosterol / sterol",         col = "#1B7837", set = kg("map00100")),
    list(lab = "Cytochrome-P450 detox",       col = "#D6604D", set = union(kg("map00982"), kg("map00980"))),
    list(lab = "Respiration (OXPHOS / AOX)",  col = "#E08214", set = union(kg("map00190"), aoxg())),
    list(lab = "Translation (ribosome)",      col = "#2F6BB3", set = kg("map03010")),
    list(lab = "Proteasome",                  col = "#11838E", set = kg("map03050")),
    list(lab = "TCA cycle",                   col = "#8C510A", set = kg("map00020")))
}
cat_of <- function(genes) {                 # first matching category label, else NA
  CA <- cats(); out <- rep(NA_character_, length(genes))
  for (k in CA) { hit <- is.na(out) & genes %in% k$set; out[hit] <- k$lab }
  out
}
cat_cols <- function() { CA <- cats(); setNames(vapply(CA, `[[`, "", "col"), vapply(CA, `[[`, "", "lab")) }

ci95 <- function(v) { v <- as.numeric(v); n <- length(v)
  if (n > 1) sd(v)/sqrt(n) * qt(0.975, n - 1) else 0 }

load_de <- function(name, root = ".") {
  d <- read.csv(file.path(root, DE_DIR, paste0("DE_", name, ".csv")), stringsAsFactors = FALSE)
  d <- d[!is.na(d$log2FoldChange) & !is.na(d$padj), c("gene_id","log2FoldChange","lfcSE","padj")]
  names(d) <- c("gene","lfc","se","padj"); d
}

## median-of-ratios normalisation, matching the Python implementation ---------
norm_counts <- function(cnt) {
  logc <- log(replace(cnt, cnt == 0, NA))
  gm <- rowMeans(logc, na.rm = TRUE)   # pandas .mean() skips NaN
  ok <- is.finite(gm)
  sf <- exp(apply(sweep(logc[ok, , drop = FALSE], 1, gm[ok], "-"), 2, median, na.rm = TRUE))
  sweep(cnt, 2, sf, "/")
}

## one volcano facet ----------------------------------------------------------
volcano_facet <- function(de, extra = NULL, xlim = c(-8, 8), ylim = c(-2, 52),
                          label_per_cat = 1) {
  de <- de %>% mutate(y = -log10(pmax(padj, 1e-50)),
                      sig = padj < FDR_MAX & abs(lfc) >= LFC_MIN,
                      cat = cat_of(gene))
  lab <- NULL
  if (label_per_cat > 0)
    lab <- de %>% filter(sig, !is.na(cat), !is.na(gsym(gene))) %>%
      group_by(cat) %>% slice_max(abs(lfc), n = label_per_cat) %>% ungroup() %>%
      mutate(txt = gsym(gene))
  if (!is.null(extra)) {
    ex <- de %>% filter(gene %in% names(extra)) %>% mutate(txt = unname(extra[gene]))
    lab <- bind_rows(lab, ex)
  }
  list(de = de, lab = lab,
       n_up = sum(de$sig & de$lfc > 0), n_dn = sum(de$sig & de$lfc < 0),
       n_tot = sum(de$sig))
}
