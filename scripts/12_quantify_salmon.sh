#!/usr/bin/env bash
###############################################################################
## run_salmon.sh  --  build a Salmon index + quantify all samples (Project 12129)
##
## Heavy step: run this on a cluster / workstation with enough CPU & RAM, NOT in
## the R session. The R script (scripts/03_differential_expression.R) will then IMPORT the quant.sf
## files with run_quant=FALSE.
##
## EDIT the three variables below, then:   bash run_salmon.sh
###############################################################################
set -euo pipefail

# ---- USER INPUTS -----------------------------------------------------------
TRANSCRIPTOME="data/reference/Zymoseptoria_tritici.MG2.cdna.all.fa.gz"
READS_DIR="01_reads"
OUT_DIR="data/rnaseq/salmon_quant"
THREADS=4          # 4 is comfortable on a MacBook Air; raise if on a server
# ----------------------------------------------------------------------------

INDEX="data/reference/salmon_index"
mkdir -p "$OUT_DIR"

# Build index once (decoy-aware indexing is even better if you add a genome decoy)
if [ ! -d "$INDEX" ]; then
  echo "[salmon] building index from $TRANSCRIPTOME"
  salmon index -t "$TRANSCRIPTOME" -i "$INDEX" -k 31 -p "$THREADS"
fi

# Quantify every paired-end sample. Skips control/blank by requiring the
# <temp>_<proth>_R<rep> pattern in the name.
shopt -s nullglob
for R1 in "$READS_DIR"/*_R1_001.fastq.gz; do
  base=$(basename "$R1" _R1_001.fastq.gz)
  # only treatment samples: <project>_<temp>_<proth>_R<rep>_S<id>
  if [[ ! "$base" =~ ^[0-9]+_[0-9]+_[0-9]+_R[0-9]+_S[0-9]+$ ]]; then
    echo "[salmon] skipping non-treatment sample: $base"
    continue
  fi
  R2="${R1/_R1_001/_R2_001}"
  out="$OUT_DIR/$base"
  if [ -f "$out/quant.sf" ]; then echo "[salmon] done: $base"; continue; fi
  echo "[salmon] quant: $base"
  salmon quant -i "$INDEX" -l A \
      -1 "$R1" -2 "$R2" \
      -p "$THREADS" --gcBias --seqBias \
      -o "$out"
done

echo "[salmon] all done -> $OUT_DIR"
echo "Next: in R set config\$quant_method='salmon', run_quant=FALSE, then run scripts/03_differential_expression.R"
