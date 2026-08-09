#!/usr/bin/env bash
###############################################################################
## reference_setup.sh
##   Download the Zymoseptoria tritici IPO323 reference (assembly MG2) from
##   Ensembl Fungi and build the tx2gene map needed by the Salmon route.
##
##   IPO323 IS the reference strain for this species, so the "reference genome"
##   is just the canonical Ensembl assembly:
##       strain   : IPO323
##       assembly : MG2   (= JGI MYCGR v2.0, NCBI GCA_000219625.1)
##
##   This script:
##     1. finds the latest Ensembl Fungi release automatically
##     2. downloads -> data/reference/ :
##          - transcriptome  (cDNA FASTA)   for Salmon
##          - genome FASTA    (optional, for STAR/featureCounts or decoy-aware index)
##          - GTF annotation  (to build tx2gene + for the genome route)
##     3. builds data/reference/tx2gene.csv  (transcript_id,gene_id) with awk
##
##   Run from the project root:   bash reference_setup.sh
##   Needs internet + wget (or curl). ~tens of MB, takes under a minute.
###############################################################################
set -euo pipefail

SPECIES="zymoseptoria_tritici"      # FTP path is lowercase
SPECIES_CAP="Zymoseptoria_tritici"  # file names are Capitalised
OUTDIR="reference"
mkdir -p "$OUTDIR"

# --- pick a downloader -------------------------------------------------------
if command -v wget >/dev/null 2>&1; then
  DL() { wget -q -O "$2" "$1"; }
  LIST() { wget -q -O - "$1"; }
elif command -v curl >/dev/null 2>&1; then
  DL() { curl -fsSL -o "$2" "$1"; }
  LIST() { curl -fsSL "$1"; }
else
  echo "ERROR: need wget or curl installed." >&2; exit 1
fi

# --- mirrors (tried in order) ------------------------------------------------
MIRRORS=(
  "https://ftp.ensemblgenomes.org/pub/fungi"
  "https://ftp.ebi.ac.uk/ensemblgenomes/pub/fungi"
  "http://ftp.ensemblgenomes.org/pub/fungi"
)

# --- find the latest release-NN directory on the first reachable mirror ------
BASE=""; RELEASE=""
for m in "${MIRRORS[@]}"; do
  echo "[ref] probing mirror: $m"
  listing="$(LIST "$m/" 2>/dev/null || true)"
  rel="$(printf '%s\n' "$listing" | grep -oE 'release-[0-9]+' | sort -t- -k2 -n | tail -1 || true)"
  if [ -n "$rel" ]; then BASE="$m"; RELEASE="$rel"; break; fi
done
[ -n "$BASE" ] || { echo "ERROR: could not reach any Ensembl Fungi mirror." >&2; exit 1; }
echo "[ref] using $BASE / $RELEASE"

FASTA_CDNA_DIR="$BASE/$RELEASE/fasta/$SPECIES/cdna"
FASTA_DNA_DIR="$BASE/$RELEASE/fasta/$SPECIES/dna"
GTF_DIR="$BASE/$RELEASE/gtf/$SPECIES"

# --- helper: resolve an exact filename from a directory listing --------------
resolve() {  # $1 = dir url, $2 = grep pattern
  LIST "$1/" 2>/dev/null | grep -oE "$2" | sort -u | tail -1
}

CDNA_FILE="$(resolve "$FASTA_CDNA_DIR" "${SPECIES_CAP}\.[^\"]*cdna\.all\.fa\.gz")"
DNA_FILE="$(resolve  "$FASTA_DNA_DIR"  "${SPECIES_CAP}\.[^\"]*dna\.toplevel\.fa\.gz")"
GTF_FILE="$(resolve  "$GTF_DIR"        "${SPECIES_CAP}\.[^\"]*\.[0-9]+\.gtf\.gz")"

echo "[ref] transcriptome : ${CDNA_FILE:-NOT FOUND}"
echo "[ref] genome        : ${DNA_FILE:-NOT FOUND}"
echo "[ref] annotation    : ${GTF_FILE:-NOT FOUND}"

# --- download ----------------------------------------------------------------
[ -n "$CDNA_FILE" ] && { echo "[ref] downloading cDNA";   DL "$FASTA_CDNA_DIR/$CDNA_FILE" "$OUTDIR/$CDNA_FILE"; }
[ -n "$DNA_FILE"  ] && { echo "[ref] downloading genome"; DL "$FASTA_DNA_DIR/$DNA_FILE"  "$OUTDIR/$DNA_FILE";  }
[ -n "$GTF_FILE"  ] && { echo "[ref] downloading GTF";    DL "$GTF_DIR/$GTF_FILE"        "$OUTDIR/$GTF_FILE";  }

# --- build tx2gene.csv from the GTF (no R needed) ----------------------------
if [ -n "${GTF_FILE:-}" ] && [ -f "$OUTDIR/$GTF_FILE" ]; then
  echo "[ref] building tx2gene.csv"
  # Portable: gzip -dc works on macOS + Linux (zcat on macOS wants .Z).
  # Pure-awk extraction (no GNU sed) handles either attribute order. Every
  # transcript/exon/CDS line in an Ensembl GTF carries both ids, so we read all
  # lines that have a transcript_id and de-duplicate.
  {
    echo "transcript_id,gene_id"
    gzip -dc "$OUTDIR/$GTF_FILE" | awk -F'\t' '
      $9 ~ /transcript_id/ {
        g=""; t="";
        if (match($9, /gene_id "[^"]+"/))      { s=substr($9,RSTART,RLENGTH); gsub(/gene_id "|"/,"",s);      g=s }
        if (match($9, /transcript_id "[^"]+"/)) { s=substr($9,RSTART,RLENGTH); gsub(/transcript_id "|"/,"",s); t=s }
        if (g!="" && t!="") print t","g
      }' | sort -u
  } > "$OUTDIR/tx2gene.csv"
  echo "[ref] tx2gene rows: $(( $(wc -l < "$OUTDIR/tx2gene.csv") - 1 ))"
fi

echo ""
echo "[ref] DONE. Files in $OUTDIR/:"
ls -lh "$OUTDIR"
echo ""
echo "Now set in scripts/03_differential_expression.R CONFIG:"
echo "  ref_transcriptome = \"$OUTDIR/${CDNA_FILE:-<cdna.fa.gz>}\""
echo "  tx2gene_csv       = \"$OUTDIR/tx2gene.csv\""
echo "  ref_genome_fasta  = \"$OUTDIR/${DNA_FILE:-<genome.fa.gz>}\"   # only for STAR route"
echo "  ref_gtf           = \"$OUTDIR/${GTF_FILE:-<annotation.gtf.gz>}\""
