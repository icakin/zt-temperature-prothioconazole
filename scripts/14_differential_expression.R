#!/usr/bin/env Rscript
###############################################################################
## RNA-seq analysis workflow  --  Project 12129 (run X0123)
##
## Experiment: Zymoseptoria tritici, 3 temperatures x 2 prothioconazole doses
##   Temperature        : 15, 21, 27  (degrees C)
##   Prothioconazole     : 0, 2        (micrograms / mL)
##   Replicates          : up to 8 per group (some groups have 7)
##   Library             : paired-end (_R1_001 / _R2_001)
##
## This single script covers:
##   (1) project setup & directory creation
##   (2) sample-metadata construction from FASTQ filenames
##   (3) QC wrappers (FastQC / FastQ Screen / MultiQC) + R QC plots
##   (4) read quantification  -- Salmon (recommended) OR STAR+featureCounts
##                              OR start directly from an existing count matrix
##   (5) DESeq2 differential expression (interaction design + group factor)
##   (6) low-count filtering, thresholds, annotation, CSV export
##   (7) journal-quality figures (PCA, distance heatmap, volcano, MA,
##       DEG heatmaps, interaction heatmap, gene expression panels,
##       UpSet/Venn, optional GO/KEGG)
##   (8) reproducibility (package checks, sessionInfo)
##
## HOW TO RUN
##   Run from the PROJECT ROOT (the folder that contains 01_reads/):
##       Rscript rnaseq_analysis.R
##   or open in RStudio with the project root as the working directory and
##   source it section by section. Edit ONLY the CONFIG block below.
##
## WHAT YOU MUST PROVIDE MANUALLY  (see README_workflow.md for details)
##   - A reference transcriptome FASTA (for Salmon)            [REF_TRANSCRIPTOME]
##       e.g. Zymoseptoria tritici IPO323 cDNA from Ensembl Fungi
##   - OR a genome FASTA + GTF/GFF annotation (for STAR+featureCounts)
##       [REF_GENOME_FASTA], [REF_GTF]
##   - Optional: a gene annotation table (gene_id -> name/description/GO)
##       [GENE_ANNOTATION_CSV]
##   The heavy external steps (Salmon/STAR) normally run on a compute cluster;
##   this script can either launch them (RUN_QUANT=TRUE) or just import results.
###############################################################################


## ===========================================================================
## 0. CONFIG  --  EDIT THIS BLOCK, then leave the rest of the script untouched
## ===========================================================================

config <- list(

  ## --- Project root & inputs -------------------------------------------------
  ## Default: the current working directory. Override with an absolute path.
  project_root = normalizePath(getwd(), mustWork = FALSE),
  reads_dir    = "01_reads",          # relative to project_root
  fastq_glob   = "\\.fastq\\.gz$",    # regex used to detect read files

  ## Filenames look like:  12129_15_0_R1_S9_R1_001.fastq.gz
  ##   <project>_<temp>_<proth>_R<replicate>_S<sampleID>_<readpair>_001.fastq.gz
  ## NOTE the deliberate ambiguity: R<replicate> (R1..R8) is the biological
  ## replicate, while the trailing _R1_/_R2_ is the paired-end read tag.
  read1_tag = "_R1_001",
  read2_tag = "_R2_001",

  ## --- Output layout (created automatically) ---------------------------------
  out_root   = "results",
  dir_qc     = "data/qc/salmon_qc",
  dir_quant  = "data/rnaseq/salmon_quant",
  dir_counts = "data/rnaseq",
  dir_de     = "tables/rnaseq",
  dir_fig    = "figures/diagnostics",
  dir_meta   = "data/rnaseq",

  ## --- Quantification --------------------------------------------------------
  ## One of: "salmon"  (recommended; transcriptome, alignment-free)
  ##         "featurecounts" (genome alignment + GTF gene counting)
  ##         "count_matrix"  (skip quantification; load an existing matrix)
  quant_method = "salmon",

  ## Set RUN_QUANT = TRUE to actually launch Salmon/STAR from R (needs the
  ## tools on PATH and a lot of CPU/RAM). FALSE = assume quant already done and
  ## only IMPORT the results. Default FALSE: safest on a laptop.
  run_quant = FALSE,
  threads   = 8,

  ## ---- USER-PROVIDED REFERENCE FILES  (fill these in) -----------------------
  ##   Salmon route:  (paths set by reference_setup.sh)
  ref_transcriptome = "data/reference/Zymoseptoria_tritici.MG2.cdna.all.fa.gz",
  salmon_index      = "data/reference/salmon_index",   # built if missing
  tx2gene_csv       = "data/reference/tx2gene.csv",     # 2-col (transcript_id, gene_id)
  ##   STAR + featureCounts route:
  ref_genome_fasta  = "data/reference/Zymoseptoria_tritici.MG2.dna.toplevel.fa.gz",
  ref_gtf           = "data/reference/Zymoseptoria_tritici.MG2.63.gtf.gz",
  star_index        = "data/reference/star_index",
  ##   Existing count-matrix route:
  count_matrix_csv  = "",   # genes (rows) x samples (cols), first col = gene_id

  ## --- Optional gene annotation ---------------------------------------------
  ## 2+ column CSV keyed on gene_id, e.g. columns: gene_id,gene_name,description
  gene_annotation_csv = "",

  ## --- DESeq2 / thresholds ---------------------------------------------------
  min_count        = 10,    # gene kept if it has >= min_count in ...
  min_samples      = 7,     # ... at least this many samples (smallest group n)
  alpha            = 0.05,  # adjusted-p significance cutoff
  lfc_threshold    = 1,     # |log2FC| cutoff for "DE" calls and figures
  reference_temp   = "15",  # baseline temperature level
  reference_proth  = "0",   # baseline prothioconazole level

  ## --- Figures ---------------------------------------------------------------
  fig_dpi   = 300,
  fig_width = 7,            # inches (default; individual plots may override)
  fig_height= 6,
  n_top_genes = 40,         # genes shown in the "top DEG" heatmap

  ## --- Control samples to EXCLUDE from modelling -----------------------------
  ## Anything whose parsed group is one of these is treated as a QC control.
  control_labels = c("BLANK", "Nc")
)


## ===========================================================================
## 1. PACKAGES  --  install checks (CRAN + Bioconductor), then load
## ===========================================================================

## Comment out install_missing() if your environment is already provisioned
## (e.g. a conda/renv lockfile). Loading is always attempted.

cran_pkgs <- c("ggplot2", "data.table", "RColorBrewer", "pheatmap",
               "ggrepel", "scales", "UpSetR", "matrixStats", "cowplot")
bioc_pkgs <- c("DESeq2", "tximport", "tximeta", "apeglm", "ashr", "vsn",
               "ComplexHeatmap", "EnhancedVolcano", "Rsubread",
               "AnnotationDbi", "clusterProfiler")

install_missing <- function(cran = cran_pkgs, bioc = bioc_pkgs) {
  miss_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_cran)) install.packages(miss_cran, repos = "https://cloud.r-project.org")
  if (length(bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    miss_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss_bioc)) BiocManager::install(miss_bioc, update = FALSE, ask = FALSE)
  }
}
## --> Uncomment to auto-install on first run:
# install_missing()

## Load what we can; warn (do not stop) for optional packages so the metadata
## and QC sections still run on a minimal setup.
suppressPackageStartupMessages({
  .need <- function(p, optional = FALSE) {
    ok <- requireNamespace(p, quietly = TRUE)
    if (ok) library(p, character.only = TRUE)
    else if (optional) message("[optional] package not installed: ", p)
    else warning("Required package not installed: ", p,
                 " -- some sections will be skipped.", call. = FALSE)
    invisible(ok)
  }
  for (p in c("ggplot2", "data.table", "RColorBrewer", "scales")) .need(p)
  for (p in c("DESeq2", "tximport", "pheatmap", "ggrepel", "UpSetR",
              "EnhancedVolcano", "ComplexHeatmap", "matrixStats", "cowplot",
              "apeglm", "Rsubread", "clusterProfiler")) .need(p, optional = TRUE)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a


## ===========================================================================
## 2. PROJECT SETUP  --  resolve paths, create output folders, detect reads
## ===========================================================================

setwd(config$project_root)
message("Project root: ", config$project_root)

out_dirs <- c(config$out_root, config$dir_meta, config$dir_qc, config$dir_quant,
              config$dir_counts, config$dir_de, config$dir_fig)
for (d in out_dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

reads_path <- file.path(config$project_root, config$reads_dir)
if (!dir.exists(reads_path))
  warning("Reads directory not found: ", reads_path,
          " -- detection will return nothing.", call. = FALSE)

## Detect ALL fastq.gz files, then keep R1 files as the per-sample anchor.
all_fastq <- list.files(reads_path, pattern = config$fastq_glob, full.names = TRUE)
message("Detected ", length(all_fastq), " FASTQ files in ", config$reads_dir)

r1_files <- all_fastq[grepl(config$read1_tag, basename(all_fastq), fixed = TRUE)]
r2_files <- all_fastq[grepl(config$read2_tag, basename(all_fastq), fixed = TRUE)]

## ---- Detect paired vs single-end ------------------------------------------
is_paired <- length(r2_files) > 0 &&
             length(r2_files) >= 0.5 * length(r1_files)
message("Library layout detected as: ",
        if (is_paired) "PAIRED-END" else "SINGLE-END",
        "  (R1=", length(r1_files), ", R2=", length(r2_files), ")")
if (length(r1_files) == 0 && length(all_fastq) > 0) {
  warning("No files matched read1_tag='", config$read1_tag,
          "'. Treating every fastq as a single-end sample.", call. = FALSE)
  r1_files <- all_fastq
  is_paired <- FALSE
}


## ===========================================================================
## 3. SAMPLE METADATA  --  parse filenames, build & save the table
## ===========================================================================
## Strategy: strip the read-pair tag + extension to get a clean sample id, then
## match the treatment pattern. Anything that does NOT match (BLANK, Nc, ...)
## is recorded as a control and excluded from the model.

strip_to_sample <- function(x) {
  b <- basename(x)
  b <- sub(paste0(config$read1_tag, "\\.fastq\\.gz$"), "", b)
  b <- sub(paste0(config$read2_tag, "\\.fastq\\.gz$"), "", b)
  b <- sub("\\.fastq\\.gz$", "", b)
  b
}

sample_ids <- strip_to_sample(r1_files)

## Treatment pattern:  <project>_<temp>_<proth>_R<rep>_S<id>
trt_pat <- "^(\\d+)_([0-9]+)_([0-9]+)_R([0-9]+)_S([0-9]+)$"

parse_one <- function(s) {
  m <- regmatches(s, regexec(trt_pat, s))[[1]]
  if (length(m) == 6) {
    data.frame(sample        = s,
               project       = m[2],
               temperature   = m[3],
               prothioconazole = m[4],
               replicate     = as.integer(m[5]),
               sample_id     = paste0("S", m[6]),
               is_control    = FALSE,
               stringsAsFactors = FALSE)
  } else {
    ## Control / non-conforming name: try to grab a trailing S<id>
    sid <- sub(".*_(S[0-9]+)$", "\\1", s)
    lab <- sub("^\\d+_", "", s)                 # drop project prefix
    lab <- sub("_S[0-9]+$", "", lab)            # -> "BLANK" / "Nc"
    data.frame(sample = s, project = sub("_.*", "", s),
               temperature = NA_character_, prothioconazole = NA_character_,
               replicate = NA_integer_,
               sample_id = ifelse(grepl("^S[0-9]+$", sid), sid, NA_character_),
               is_control = TRUE, control_label = lab,
               stringsAsFactors = FALSE)
  }
}

meta_list <- lapply(sample_ids, parse_one)
## rbind with fill (control rows have an extra column)
all_cols <- unique(unlist(lapply(meta_list, names)))
meta_list <- lapply(meta_list, function(d) {
  for (cc in setdiff(all_cols, names(d))) d[[cc]] <- NA
  d[all_cols]
})
metadata <- do.call(rbind, meta_list)

## Force-flag known control labels even if they happened to parse oddly
if ("control_label" %in% names(metadata)) {
  metadata$is_control <- metadata$is_control |
    metadata$control_label %in% config$control_labels
}

## Attach file paths
metadata$fastq_1 <- r1_files[match(metadata$sample, sample_ids)]
if (is_paired) {
  r2_ids <- strip_to_sample(r2_files)
  metadata$fastq_2 <- r2_files[match(metadata$sample, r2_ids)]
} else {
  metadata$fastq_2 <- NA_character_
}

## Human-readable, ordered factors + a combined treatment group
temp_levels  <- c("15", "21", "27")
proth_levels <- c("0", "2")
metadata$temperature     <- factor(metadata$temperature,     levels = temp_levels)
metadata$prothioconazole <- factor(metadata$prothioconazole, levels = proth_levels)
metadata$temp_label  <- ifelse(is.na(metadata$temperature), NA,
                               paste0(as.character(metadata$temperature), "C"))
metadata$proth_label <- ifelse(is.na(metadata$prothioconazole), NA,
                               paste0(as.character(metadata$prothioconazole), "ug/mL"))
## Combined group, e.g. "15_0", used for simple pairwise contrasts
metadata$group <- ifelse(metadata$is_control, "control",
                         paste(metadata$temperature, metadata$prothioconazole, sep = "_"))
group_levels <- as.vector(t(outer(temp_levels, proth_levels, paste, sep = "_")))
metadata$group <- factor(metadata$group,
                         levels = c(group_levels, "control"))

## Order rows for readability
metadata <- metadata[order(metadata$is_control, metadata$temperature,
                           metadata$prothioconazole, metadata$replicate), ]
rownames(metadata) <- metadata$sample

## ---- Report & save --------------------------------------------------------
n_trt  <- sum(!metadata$is_control)
n_ctrl <- sum(metadata$is_control)
message("Parsed ", nrow(metadata), " samples: ", n_trt, " treatment + ",
        n_ctrl, " control.")
message("Replicates per treatment group:")
print(table(droplevels(metadata$group[!metadata$is_control])))

meta_out <- file.path(config$dir_meta, "sample_metadata.csv")
write.csv(metadata, meta_out, row.names = FALSE)
message("Saved metadata -> ", meta_out)

## Modelling subset (no controls, complete design)
coldata <- droplevels(metadata[!metadata$is_control, ])
coldata$temperature     <- relevel(coldata$temperature,     ref = config$reference_temp)
coldata$prothioconazole <- relevel(coldata$prothioconazole, ref = config$reference_proth)


## ===========================================================================
## 4. QUALITY CONTROL
## ===========================================================================
## (a) Shell wrappers for FastQC / FastQ Screen / MultiQC. This project ALREADY
##     has 09_QC_reports/ (FastQC v0.12.1, FastQ Screen v0.16.0, MultiQC v1.35),
##     so by default we REUSE those and only (re)build the MultiQC summary.
##     Set rerun_fastqc = TRUE to regenerate from scratch.
## (b) R-side QC plots come later from the DESeq2 object (library sizes, PCA,
##     sample-distance heatmap).

run_fastqc_multiqc <- function(fastq = metadata$fastq_1[!is.na(metadata$fastq_1)],
                               out = config$dir_qc,
                               existing_qc = "09_QC_reports",
                               rerun_fastqc = FALSE,
                               threads = config$threads) {
  has <- function(tool) nzchar(Sys.which(tool))
  if (rerun_fastqc) {
    if (!has("fastqc")) { message("fastqc not on PATH; skipping."); }
    else {
      message("Running FastQC on ", length(fastq), " files ...")
      system2("fastqc", c("-t", threads, "-o", shQuote(out), shQuote(fastq)))
    }
  } else {
    message("Reusing existing FastQC/FastQ Screen output in ", existing_qc, "/")
  }
  ## MultiQC aggregates everything it finds under the given search paths.
  if (has("multiqc")) {
    message("Aggregating QC with MultiQC ...")
    system2("multiqc", c(shQuote(existing_qc), shQuote(out),
                         "-o", shQuote(out), "-n", "multiqc_report", "-f"))
  } else {
    message("multiqc not on PATH; open ", existing_qc,
            "/ reports manually or `pip install multiqc`.")
  }
}
## --> Uncomment to (re)build the MultiQC summary:
# run_fastqc_multiqc()

## R-based QC plot: parse the existing FastQ Screen *_screen.txt files (if any)
## into a single contamination summary figure. Safe no-op if absent.
qc_fastq_screen_plot <- function(screen_dir = "09_QC_reports/fastq_screen",
                                 out = config$dir_qc) {
  txts <- list.files(screen_dir, pattern = "_screen\\.txt$", full.names = TRUE)
  if (!length(txts) || !requireNamespace("ggplot2", quietly = TRUE)) {
    message("No FastQ Screen txt files (or ggplot2 missing); skipping QC plot.")
    return(invisible(NULL))
  }
  rd <- function(f) {
    x <- tryCatch(read.delim(f, skip = 1, header = TRUE, comment.char = "%",
                             stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(x) || !"Genome" %in% names(x)) return(NULL)
    x$sample <- sub("_R1_001_screen\\.txt$", "", basename(f))
    x
  }
  df <- do.call(rbind, Filter(Negate(is.null), lapply(txts, rd)))
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  pcol <- grep("X.One_hit_one_genome|.*one.*genome", names(df),
               ignore.case = TRUE, value = TRUE)[1]
  if (is.na(pcol)) return(invisible(NULL))
  df$pct <- suppressWarnings(as.numeric(df[[pcol]]))
  p <- ggplot2::ggplot(df, ggplot2::aes(Genome, sample, fill = pct)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "% mapped\n(1 genome)", option = "C") +
    ggplot2::labs(title = "FastQ Screen: contamination overview",
                  x = "Reference genome", y = "Sample") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_fig(p, "qc_fastq_screen_overview",
           width = 8, height = max(4, 0.18 * length(unique(df$sample))))
}


## ===========================================================================
## 5. READ QUANTIFICATION
## ===========================================================================
## Three mutually exclusive routes selected by config$quant_method.
## All routes converge on a DESeqDataSet `dds`.
##
##  REQUIRED USER INPUTS (left blank in CONFIG):
##    salmon         -> config$ref_transcriptome (+ tx2gene_csv)
##    featurecounts  -> config$ref_genome_fasta + config$ref_gtf
##    count_matrix   -> config$count_matrix_csv
## ---------------------------------------------------------------------------

stopifnot(requireNamespace("DESeq2", quietly = TRUE))
design_formula <- ~ temperature + prothioconazole + temperature:prothioconazole

## ---- 5A. SALMON (recommended) ---------------------------------------------
build_salmon <- function() {
  if (!nzchar(config$ref_transcriptome) && config$run_quant)
    stop("config$ref_transcriptome is empty -- provide a transcriptome FASTA.")
  idx <- config$salmon_index
  ## (i) Build index + quantify, only if RUN_QUANT and salmon is available.
  if (config$run_quant && nzchar(Sys.which("salmon"))) {
    if (!dir.exists(idx)) {
      message("Building Salmon index ...")
      system2("salmon", c("index", "-t", shQuote(config$ref_transcriptome),
                          "-i", shQuote(idx), "-k", "31", "-p", config$threads))
    }
    for (i in which(!coldata$is_control)) {
      s  <- coldata$sample[i]
      od <- file.path(config$dir_quant, s)
      if (file.exists(file.path(od, "quant.sf"))) next
      args <- c("quant", "-i", shQuote(idx), "-l", "A",
                "-p", config$threads,
                "--gcBias", "--seqBias", "-o", shQuote(od))
      if (is_paired) args <- c(args, "-1", shQuote(coldata$fastq_1[i]),
                                      "-2", shQuote(coldata$fastq_2[i]))
      else           args <- c(args, "-r", shQuote(coldata$fastq_1[i]))
      message("Salmon quant: ", s)
      system2("salmon", args)
    }
  } else if (config$run_quant) {
    stop("RUN_QUANT=TRUE but salmon is not on PATH.")
  } else {
    message("RUN_QUANT=FALSE: importing pre-computed Salmon quant.sf files.")
  }
  ## (ii) Import with tximport.
  files <- file.path(config$dir_quant, coldata$sample, "quant.sf")
  names(files) <- coldata$sample
  ok <- file.exists(files)
  if (!all(ok)) {
    warning("Missing quant.sf for: ", paste(coldata$sample[!ok], collapse = ", "),
            "\n  -> run Salmon first (see README) or set quant_method='count_matrix'.",
            call. = FALSE)
    return(NULL)
  }
  if (!nzchar(config$tx2gene_csv))
    stop("config$tx2gene_csv is empty -- a transcript->gene map is required for gene-level DE.")
  tx2gene <- read.csv(config$tx2gene_csv, header = TRUE, stringsAsFactors = FALSE)
  txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2],
                            dropInfReps = TRUE)
  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = coldata,
                                          design = design_formula)
  list(dds = dds, txi = txi)
}

## ---- 5B. STAR + featureCounts (genome route) ------------------------------
build_featurecounts <- function() {
  if (config$run_quant) {
    if (!nzchar(config$ref_genome_fasta) || !nzchar(config$ref_gtf))
      stop("featureCounts route needs ref_genome_fasta + ref_gtf.")
    ## Align with STAR (external). Build index once, then map each sample.
    if (nzchar(Sys.which("STAR"))) {
      if (!dir.exists(config$star_index)) {
        dir.create(config$star_index, recursive = TRUE, showWarnings = FALSE)
        system2("STAR", c("--runMode", "genomeGenerate",
                          "--genomeDir", shQuote(config$star_index),
                          "--genomeFastaFiles", shQuote(config$ref_genome_fasta),
                          "--sjdbGTFfile", shQuote(config$ref_gtf),
                          "--runThreadN", config$threads))
      }
      for (i in seq_len(nrow(coldata))) {
        s <- coldata$sample[i]
        pref <- file.path(config$dir_quant, paste0(s, "."))
        if (file.exists(paste0(pref, "Aligned.sortedByCoord.out.bam"))) next
        rf <- c(coldata$fastq_1[i], if (is_paired) coldata$fastq_2[i])
        system2("STAR", c("--genomeDir", shQuote(config$star_index),
                          "--readFilesIn", paste(shQuote(rf), collapse = " "),
                          "--readFilesCommand", "zcat",
                          "--runThreadN", config$threads,
                          "--outSAMtype", "BAM", "SortedByCoordinate",
                          "--outFileNamePrefix", shQuote(pref)))
      }
    } else stop("RUN_QUANT=TRUE but STAR not on PATH.")
  }
  ## Count with Rsubread::featureCounts (R-native; no external featureCounts).
  bams <- file.path(config$dir_quant,
                    paste0(coldata$sample, ".Aligned.sortedByCoord.out.bam"))
  names(bams) <- coldata$sample
  if (!all(file.exists(bams))) {
    warning("Missing BAMs; run STAR or use quant_method='count_matrix'.", call. = FALSE)
    return(NULL)
  }
  fc <- Rsubread::featureCounts(files = bams, annot.ext = config$ref_gtf,
                                isGTFAnnotation = TRUE, GTF.featureType = "exon",
                                GTF.attrType = "gene_id",
                                isPairedEnd = is_paired, nthreads = config$threads)
  mat <- fc$counts; colnames(mat) <- coldata$sample
  dds <- DESeq2::DESeqDataSetFromMatrix(mat, colData = coldata,
                                        design = design_formula)
  list(dds = dds, fc = fc)
}

## ---- 5C. START FROM AN EXISTING COUNT MATRIX ------------------------------
build_from_matrix <- function() {
  if (!nzchar(config$count_matrix_csv) || !file.exists(config$count_matrix_csv))
    stop("quant_method='count_matrix' but config$count_matrix_csv is missing.")
  cm <- read.csv(config$count_matrix_csv, row.names = 1, check.names = FALSE)
  cm <- cm[, intersect(coldata$sample, colnames(cm)), drop = FALSE]
  cd <- coldata[match(colnames(cm), coldata$sample), ]
  mat <- as.matrix(round(cm))
  dds <- DESeq2::DESeqDataSetFromMatrix(mat, colData = cd, design = design_formula)
  list(dds = dds)
}

## ---- Dispatch -------------------------------------------------------------
quant <- switch(config$quant_method,
                salmon        = build_salmon(),
                featurecounts = build_featurecounts(),
                count_matrix  = build_from_matrix(),
                stop("Unknown quant_method: ", config$quant_method))

if (is.null(quant)) {
  message("\n*** Quantification not available yet. ***\n",
          "Provide the reference inputs in CONFIG and either set run_quant=TRUE\n",
          "or supply pre-computed results, then re-run. Metadata & QC steps above\n",
          "have completed and been saved. Stopping cleanly here.")
  ## Save a session log even on early exit, then stop without error.
  writeLines(capture.output(sessionInfo()),
             file.path(config$out_root, "sessionInfo.txt"))
  quit(save = "no", status = 0)
}
dds <- quant$dds

## Persist the raw count matrix
raw_counts <- DESeq2::counts(dds)
write.csv(data.frame(gene_id = rownames(raw_counts), raw_counts,
                     check.names = FALSE),
          file.path(config$dir_counts, "raw_count_matrix.csv"), row.names = FALSE)


## ===========================================================================
## 6. DIFFERENTIAL EXPRESSION  (DESeq2)
## ===========================================================================

## ---- 6.1 Low-count filtering ----------------------------------------------
keep <- rowSums(DESeq2::counts(dds) >= config$min_count) >= config$min_samples
message("Filtering: kept ", sum(keep), " / ", length(keep),
        " genes (>= ", config$min_count, " counts in >= ",
        config$min_samples, " samples).")
dds <- dds[keep, ]

## ---- 6.2 Fit the interaction model ----------------------------------------
## design = ~ temperature + prothioconazole + temperature:prothioconazole
DESeq2::design(dds) <- design_formula
dds <- DESeq2::DESeq(dds)
message("Model coefficients available: ",
        paste(DESeq2::resultsNames(dds), collapse = ", "))

## ---- 6.3 A parallel GROUP model for clean pairwise contrasts --------------
## Using a single combined factor makes the requested pairwise comparisons
## (e.g. 15_2 vs 15_0) trivial and unambiguous.
dds_grp <- dds
dds_grp$group <- droplevels(dds_grp$group)
DESeq2::design(dds_grp) <- ~ group
dds_grp <- DESeq2::DESeq(dds_grp)

## Variance-stabilised data for all downstream visualisation
vsd <- DESeq2::vst(dds, blind = FALSE)
saveRDS(dds,  file.path(config$dir_de, "dds.rds"))
saveRDS(vsd,  file.path(config$dir_de, "vsd.rds"))

## ---- 6.4 Helper: annotate + threshold + save a results table --------------
gene_anno <- NULL
if (nzchar(config$gene_annotation_csv) && file.exists(config$gene_annotation_csv)) {
  gene_anno <- read.csv(config$gene_annotation_csv, stringsAsFactors = FALSE)
  message("Loaded gene annotation: ", nrow(gene_anno), " rows.")
}

tidy_results <- function(res, name, lfc_shrink = NULL) {
  r <- as.data.frame(if (is.null(lfc_shrink)) res else lfc_shrink)
  r$gene_id <- rownames(r)
  ## carry the unshrunken p-values/padj if shrinkage dropped them
  if (!is.null(lfc_shrink)) {
    base <- as.data.frame(res)
    r$pvalue <- base$pvalue; r$padj <- base$padj
  }
  r$comparison <- name
  r$significant <- !is.na(r$padj) & r$padj < config$alpha &
                   abs(r$log2FoldChange) >= config$lfc_threshold
  r$direction <- ifelse(!r$significant, "ns",
                        ifelse(r$log2FoldChange > 0, "up", "down"))
  if (!is.null(gene_anno) && "gene_id" %in% names(gene_anno))
    r <- merge(r, gene_anno, by = "gene_id", all.x = TRUE, sort = FALSE)
  r <- r[order(r$padj), ]
  front <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE",
             "pvalue", "padj", "significant", "direction", "comparison")
  r <- r[, c(intersect(front, names(r)), setdiff(names(r), front))]
  write.csv(r, file.path(config$dir_de, paste0("DE_", name, ".csv")),
            row.names = FALSE)
  message(sprintf("  %-22s  %d significant (padj<%.2g, |lfc|>=%g)",
                  name, sum(r$significant, na.rm = TRUE),
                  config$alpha, config$lfc_threshold))
  r
}

## Convenience: pairwise contrast on the GROUP model
grp_contrast <- function(a, b, name) {
  res <- DESeq2::results(dds_grp, contrast = c("group", a, b),
                         alpha = config$alpha)
  shr <- tryCatch(
    DESeq2::lfcShrink(dds_grp, contrast = c("group", a, b), res = res,
                      type = "ashr"),
    error = function(e) NULL)
  tidy_results(res, name, lfc_shrink = shr)
}

de <- list()

## ---- 6.5 Prothioconazole effect AT EACH temperature -----------------------
de[["proth_at_15"]] <- grp_contrast("15_2", "15_0", "proth_at_15C_2vs0")
de[["proth_at_21"]] <- grp_contrast("21_2", "21_0", "proth_at_21C_2vs0")
de[["proth_at_27"]] <- grp_contrast("27_2", "27_0", "proth_at_27C_2vs0")

## ---- 6.6 Temperature effects WITHOUT prothioconazole ----------------------
de[["temp_21v15_0"]] <- grp_contrast("21_0", "15_0", "temp_21vs15_proth0")
de[["temp_27v15_0"]] <- grp_contrast("27_0", "15_0", "temp_27vs15_proth0")
de[["temp_27v21_0"]] <- grp_contrast("27_0", "21_0", "temp_27vs21_proth0")

## ---- 6.7 Temperature effects WITH prothioconazole -------------------------
de[["temp_21v15_2"]] <- grp_contrast("21_2", "15_2", "temp_21vs15_proth2")
de[["temp_27v15_2"]] <- grp_contrast("27_2", "15_2", "temp_27vs15_proth2")
de[["temp_27v21_2"]] <- grp_contrast("27_2", "21_2", "temp_27vs21_proth2")

## ---- 6.8 INTERACTION effects (temperature x prothioconazole) --------------
## Interaction coefficients answer: "does the prothioconazole response differ
## between temperatures?" Pull every interaction term from the full model.
inter_terms <- grep("^temperature.*\\.prothioconazole", DESeq2::resultsNames(dds),
                    value = TRUE)
for (term in inter_terms) {
  res <- DESeq2::results(dds, name = term, alpha = config$alpha)
  nm  <- paste0("interaction_", gsub("[^A-Za-z0-9]+", "_", term))
  de[[nm]] <- tidy_results(res, nm)
}
## Combined LRT: ANY interaction (drops the whole interaction block)
dds_lrt <- DESeq2::DESeq(dds, test = "LRT",
                         reduced = ~ temperature + prothioconazole)
res_lrt <- DESeq2::results(dds_lrt, alpha = config$alpha)
de[["interaction_LRT_any"]] <- tidy_results(res_lrt, "interaction_LRT_any")

## ---- 6.9 Combined summary of all comparisons ------------------------------
de_summary <- do.call(rbind, lapply(names(de), function(n) {
  r <- de[[n]]
  data.frame(comparison = n,
             n_tested = sum(!is.na(r$padj)),
             n_sig    = sum(r$significant, na.rm = TRUE),
             n_up     = sum(r$direction == "up",   na.rm = TRUE),
             n_down   = sum(r$direction == "down", na.rm = TRUE))
}))
write.csv(de_summary, file.path(config$dir_de, "DE_summary.csv"), row.names = FALSE)
message("\nDE summary:"); print(de_summary)


## ===========================================================================
## 7. JOURNAL-QUALITY FIGURES
## ===========================================================================
## House style: colorblind-friendly palettes, clean theme, PDF + PNG export.

## Colorblind-safe palettes (Okabe-Ito)
ok_ito  <- c("#000000","#E69F00","#56B4E9","#009E73",
             "#F0E442","#0072B2","#D55E00","#CC79A7")
temp_pal  <- c("15" = "#0072B2", "21" = "#E69F00", "27" = "#D55E00")
proth_pal <- c("0"  = "#999999", "2"  = "#009E73")
proth_shape <- c("0" = 16, "2" = 17)

theme_pub <- function(base = 11) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  ggplot2::theme_bw(base_size = base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(linewidth = 0.25),
                   plot.title  = ggplot2::element_text(face = "bold", size = base + 1),
                   legend.key  = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey92", colour = NA))
}

## Save any ggplot / grid object as BOTH pdf + png
save_fig <- function(plot, name, width = config$fig_width,
                     height = config$fig_height, dir = config$dir_fig) {
  pdf_f <- file.path(dir, paste0(name, ".pdf"))
  png_f <- file.path(dir, paste0(name, ".png"))
  is_gg <- inherits(plot, "ggplot")
  if (is_gg) {
    ggplot2::ggsave(pdf_f, plot, width = width, height = height, device = "pdf")
    ggplot2::ggsave(png_f, plot, width = width, height = height,
                    dpi = config$fig_dpi, device = "png")
  } else {
    ## grid/pheatmap/ComplexHeatmap objects
    grDevices::pdf(pdf_f, width = width, height = height)
      draw_obj(plot); grDevices::dev.off()
    grDevices::png(png_f, width = width, height = height, units = "in",
                   res = config$fig_dpi)
      draw_obj(plot); grDevices::dev.off()
  }
  message("  figure -> ", name, " (.pdf/.png)")
}
draw_obj <- function(x) {
  if (inherits(x, "Heatmap") || inherits(x, "HeatmapList"))
    ComplexHeatmap::draw(x)
  else if (inherits(x, "pheatmap")) grid::grid.draw(x$gtable)
  else if (inherits(x, "gtable"))   grid::grid.draw(x)
  else print(x)
}

## ---- 7.1 PCA (color = temperature, shape = prothioconazole) ---------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  pca <- DESeq2::plotPCA(vsd, intgroup = c("temperature", "prothioconazole"),
                         returnData = TRUE)
  pv  <- round(100 * attr(pca, "percentVar"))
  p_pca <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2,
              color = temperature, shape = prothioconazole)) +
    ggplot2::geom_point(size = 3.2, alpha = 0.9) +
    ggplot2::scale_color_manual(values = temp_pal, name = "Temperature (°C)") +
    ggplot2::scale_shape_manual(values = proth_shape,
              name = "Prothioconazole (µg/mL)") +
    ggplot2::labs(x = paste0("PC1 (", pv[1], "%)"),
                  y = paste0("PC2 (", pv[2], "%)"),
                  title = "Sample PCA (VST)") +
    theme_pub()
  if (requireNamespace("ggrepel", quietly = TRUE))
    p_pca <- p_pca + ggrepel::geom_text_repel(
      ggplot2::aes(label = name), size = 2.4, max.overlaps = 20,
      show.legend = FALSE)
  save_fig(p_pca, "PCA_temperature_prothioconazole", width = 7, height = 5.5)
}

## ---- 7.2 Sample-to-sample distance heatmap --------------------------------
if (requireNamespace("pheatmap", quietly = TRUE)) {
  sampleDist <- dist(t(SummarizedExperiment::assay(vsd)))
  dm <- as.matrix(sampleDist)
  ann <- data.frame(Temperature = vsd$temperature,
                    Prothioconazole = vsd$prothioconazole,
                    row.names = colnames(vsd))
  ann_cols <- list(Temperature = temp_pal, Prothioconazole = proth_pal)
  ph <- pheatmap::pheatmap(dm, clustering_distance_rows = sampleDist,
            clustering_distance_cols = sampleDist,
            annotation_col = ann, annotation_colors = ann_cols,
            col = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255),
            main = "Sample-to-sample distances (VST)", silent = TRUE)
  save_fig(ph, "sample_distance_heatmap", width = 8, height = 7)
}

## ---- 7.3 Volcano plots (EnhancedVolcano) for key comparisons --------------
volcano_set <- c("proth_at_15C_2vs0", "proth_at_21C_2vs0", "proth_at_27C_2vs0",
                 "temp_27vs15_proth0", "temp_27vs15_proth2", "interaction_LRT_any")
if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  for (nm in names(de)) {
    r <- de[[nm]]
    cmp <- unique(r$comparison)
    if (!(cmp %in% volcano_set)) next
    p <- EnhancedVolcano::EnhancedVolcano(r, lab = r$gene_id,
            x = "log2FoldChange", y = "padj",
            pCutoff = config$alpha, FCcutoff = config$lfc_threshold,
            title = cmp, subtitle = NULL, caption = NULL,
            pointSize = 1.6, labSize = 3.0,
            col = c("grey70", "#56B4E9", "#E69F00", "#D55E00"),
            colAlpha = 0.8, legendPosition = "right")
    save_fig(p, paste0("volcano_", cmp), width = 7.5, height = 6.5)
  }
}

## ---- 7.4 MA plots ---------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  for (nm in names(de)) {
    r <- de[[nm]]; cmp <- unique(r$comparison)
    if (!(cmp %in% volcano_set)) next
    r$col <- ifelse(r$significant, r$direction, "ns")
    p <- ggplot2::ggplot(r, ggplot2::aes(baseMean, log2FoldChange, color = col)) +
      ggplot2::geom_point(size = 0.7, alpha = 0.6) +
      ggplot2::scale_x_log10(labels = scales::comma) +
      ggplot2::scale_color_manual(values = c(ns = "grey75",
                up = "#D55E00", down = "#0072B2"), name = NULL) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = c(-1, 1) * config$lfc_threshold,
                lty = 2, linewidth = 0.3, color = "grey40") +
      ggplot2::labs(title = paste0("MA plot: ", cmp),
                x = "Mean of normalized counts", y = "log2 fold change") +
      theme_pub()
    save_fig(p, paste0("MA_", cmp), width = 6.5, height = 5)
  }
}

## ---- 7.5 Heatmap of top DE genes ------------------------------------------
## Union of the most significant genes across the pairwise comparisons.
top_genes <- unique(unlist(lapply(de, function(r) {
  rr <- r[r$significant %in% TRUE, ]
  rr <- rr[order(rr$padj), ]
  head(rr$gene_id, 15)
})))
top_genes <- head(top_genes, config$n_top_genes)
if (length(top_genes) >= 2 && requireNamespace("pheatmap", quietly = TRUE)) {
  mat <- SummarizedExperiment::assay(vsd)[top_genes, , drop = FALSE]
  mat <- mat - rowMeans(mat)                      # center per gene (z-ish)
  ann <- data.frame(Temperature = vsd$temperature,
                    Prothioconazole = vsd$prothioconazole,
                    row.names = colnames(vsd))
  col_ord <- order(vsd$temperature, vsd$prothioconazole)
  ph <- pheatmap::pheatmap(mat[, col_ord], cluster_cols = FALSE,
            annotation_col = ann,
            annotation_colors = list(Temperature = temp_pal,
                                     Prothioconazole = proth_pal),
            color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(255),
            show_colnames = FALSE, fontsize_row = 7,
            main = "Top differentially expressed genes (centered VST)",
            silent = TRUE)
  save_fig(ph, "heatmap_top_DE_genes", width = 8, height = 9)
}

## ---- 7.6 Heatmap of significant INTERACTION genes -------------------------
int_tab <- de[["interaction_LRT_any"]]
if (!is.null(int_tab)) {
  int_genes <- head(int_tab$gene_id[int_tab$significant %in% TRUE], config$n_top_genes)
  if (length(int_genes) >= 2 && requireNamespace("pheatmap", quietly = TRUE)) {
    mat <- SummarizedExperiment::assay(vsd)[int_genes, , drop = FALSE]
    mat <- mat - rowMeans(mat)
    ann <- data.frame(Temperature = vsd$temperature,
                      Prothioconazole = vsd$prothioconazole,
                      row.names = colnames(vsd))
    col_ord <- order(vsd$temperature, vsd$prothioconazole)
    ph <- pheatmap::pheatmap(mat[, col_ord], cluster_cols = FALSE,
              annotation_col = ann,
              annotation_colors = list(Temperature = temp_pal,
                                       Prothioconazole = proth_pal),
              color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "PuOr")))(255),
              show_colnames = FALSE, fontsize_row = 7,
              main = "Genes with significant temperature x prothioconazole interaction",
              silent = TRUE)
    save_fig(ph, "heatmap_interaction_genes", width = 8, height = 9)
  }
}

## ---- 7.7 Expression panels for selected top genes (all 6 groups) ----------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  sel <- head(top_genes, 6)
  if (length(sel) >= 1) {
    norm <- DESeq2::counts(dds, normalized = TRUE)[sel, , drop = FALSE]
    long <- do.call(rbind, lapply(sel, function(g) data.frame(
      gene = g, count = norm[g, ],
      temperature = dds$temperature, prothioconazole = dds$prothioconazole)))
    p <- ggplot2::ggplot(long, ggplot2::aes(temperature, count,
              color = prothioconazole, group = prothioconazole)) +
      ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 0.7) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(
              jitter.width = 0.12, dodge.width = 0.4), size = 1.4, alpha = 0.8) +
      ggplot2::scale_color_manual(values = proth_pal,
              name = "Prothioconazole (µg/mL)") +
      ggplot2::facet_wrap(~ gene, scales = "free_y") +
      ggplot2::labs(x = "Temperature (°C)",
                y = "Normalized count", title = "Top gene expression across treatments") +
      theme_pub()
    save_fig(p, "expression_top_genes_by_treatment", width = 9, height = 6)
  }
}

## ---- 7.8 UpSet / Venn of DE-gene overlap ----------------------------------
sig_sets <- lapply(de, function(r) r$gene_id[r$significant %in% TRUE])
names(sig_sets) <- vapply(de, function(r) unique(r$comparison), character(1))
sig_sets <- sig_sets[vapply(sig_sets, length, integer(1)) > 0]
if (length(sig_sets) >= 2 && requireNamespace("UpSetR", quietly = TRUE)) {
  grDevices::pdf(file.path(config$dir_fig, "upset_DE_overlap.pdf"),
                 width = 10, height = 6)
    print(UpSetR::upset(UpSetR::fromList(sig_sets), nsets = length(sig_sets),
          order.by = "freq", mb.ratio = c(0.6, 0.4),
          main.bar.color = "#0072B2", sets.bar.color = "#D55E00"))
  grDevices::dev.off()
  grDevices::png(file.path(config$dir_fig, "upset_DE_overlap.png"),
                 width = 10, height = 6, units = "in", res = config$fig_dpi)
    print(UpSetR::upset(UpSetR::fromList(sig_sets), nsets = length(sig_sets),
          order.by = "freq", mb.ratio = c(0.6, 0.4),
          main.bar.color = "#0072B2", sets.bar.color = "#D55E00"))
  grDevices::dev.off()
  message("  figure -> upset_DE_overlap (.pdf/.png)")
}

## ---- 7.9 OPTIONAL: GO / KEGG enrichment -----------------------------------
## Requires a gene->GO/KEGG mapping (config$gene_annotation_csv with a 'GO'
## column, OR an OrgDb). Zymoseptoria tritici has no standard Bioconductor
## OrgDb, so the most reproducible route is clusterProfiler::enricher() with a
## custom TERM2GENE built from the annotation (see README). Runs only if the
## annotation provides GO terms.
run_enrichment <- function() {
  if (is.null(gene_anno) || !"GO" %in% names(gene_anno) ||
      !requireNamespace("clusterProfiler", quietly = TRUE)) {
    message("Skipping enrichment (no GO annotation / clusterProfiler).")
    return(invisible(NULL))
  }
  ## Build TERM2GENE: one row per (GO term, gene). 'GO' may be ';'-separated.
  t2g <- do.call(rbind, lapply(seq_len(nrow(gene_anno)), function(i) {
    gos <- strsplit(as.character(gene_anno$GO[i]), "[;,| ]+")[[1]]
    gos <- gos[nzchar(gos)]
    if (!length(gos)) return(NULL)
    data.frame(term = gos, gene = gene_anno$gene_id[i])
  }))
  for (nm in names(de)) {
    r <- de[[nm]]; cmp <- unique(r$comparison)
    sig <- r$gene_id[r$significant %in% TRUE]
    if (length(sig) < 10) next
    eg <- tryCatch(clusterProfiler::enricher(sig, TERM2GENE = t2g,
                   pvalueCutoff = 0.05), error = function(e) NULL)
    if (is.null(eg) || nrow(as.data.frame(eg)) == 0) next
    write.csv(as.data.frame(eg),
              file.path(config$dir_de, paste0("enrichment_", cmp, ".csv")),
              row.names = FALSE)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- clusterProfiler::dotplot(eg, showCategory = 20) +
           ggplot2::ggtitle(paste0("Enrichment: ", cmp))
      save_fig(p, paste0("enrichment_", cmp), width = 8, height = 7)
    }
  }
}
run_enrichment()


## ===========================================================================
## 8. REPRODUCIBILITY
## ===========================================================================
writeLines(capture.output(sessionInfo()),
           file.path(config$out_root, "sessionInfo.txt"))
saveRDS(list(de = de, summary = de_summary),
        file.path(config$dir_de, "all_DE_results.rds"))

message("\n==============================================================")
message("DONE. Outputs written under: ", config$out_root, "/")
message("  metadata : ", config$dir_meta)
message("  QC       : ", config$dir_qc)
message("  counts   : ", config$dir_counts)
message("  DE tables: ", config$dir_de)
message("  figures  : ", config$dir_fig)
message("==============================================================")
###############################################################################
## END
###############################################################################
