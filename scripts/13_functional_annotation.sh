#!/usr/bin/env bash
###############################################################################
## run_eggnog.sh  --  state-of-the-art functional annotation for Z. tritici
##
##   Runs eggNOG-mapper on the Zymoseptoria tritici (MG2) proteome and builds a
##   rich, gene-level annotation table with GO + KEGG pathways + EC + Pfam +
##   COG + functional descriptions, keyed on Mycgr3G gene IDs so it joins to
##   your DESeq2 results.  Also fetches GO and KEGG term-name tables for plots.
##
##   Outputs:
##     data/reference/gene_annotation.csv   gene_id,gene_name,description,COG,GO,KEGG,EC,PFAM
##     data/reference/go_terms.csv          go_id,go_name,go_namespace      (from go-basic.obo)
##     data/reference/kegg_pathways.csv     kegg_id,kegg_name               (from KEGG REST)
##     tables/rnaseq/annotation/ztritici.emapper.annotations   (raw eggNOG output)
##
##   Run from project root:   bash run_eggnog.sh
##   Needs: conda, internet, curl, python3, ~15 GB disk for the eggNOG DB.
###############################################################################
set -euo pipefail

OUTDIR="reference"
ANNO_OUT="tables/rnaseq/annotation"
EGGDB="$OUTDIR/eggnog_db"
THREADS=4
mkdir -p "$OUTDIR" "$ANNO_OUT" "$EGGDB"

# --- 0. disk-space sanity check ---------------------------------------------
avail_kb=$(df -k . | awk 'NR==2{print $4}')
if [ "${avail_kb:-0}" -lt 60000000 ]; then
  echo "WARNING: less than ~60 GB free here. The eggNOG v5 DB needs ~50 GB"
  echo "         (eggnog.db ~48 GB + diamond db ~4.5 GB). If space is tight,"
  echo "         tell me and we'll use the lighter Ensembl-GO + KEGG-REST route."
fi

# --- 1. download the proteome (pep) from Ensembl Fungi ----------------------
# auto-detect mirror + latest release (same logic as reference_setup.sh)
MIRRORS=( "https://ftp.ensemblgenomes.org/pub/fungi" "https://ftp.ebi.ac.uk/ensemblgenomes/pub/fungi" )
BASE=""; REL=""
for m in "${MIRRORS[@]}"; do
  rel=$(curl -s "$m/" 2>/dev/null | grep -oE 'release-[0-9]+' | sort -t- -k2 -n | tail -1 || true)
  [ -n "$rel" ] && { BASE="$m"; REL="$rel"; break; }
done
[ -n "$BASE" ] || { echo "ERROR: cannot reach Ensembl Fungi FTP." >&2; exit 1; }
PEP_DIR="$BASE/$REL/fasta/zymoseptoria_tritici/pep"
PEP_FILE=$(curl -s "$PEP_DIR/" | grep -oE 'Zymoseptoria_tritici\.[^"]*pep\.all\.fa\.gz' | sort -u | tail -1)
[ -n "$PEP_FILE" ] || { echo "ERROR: could not find proteome at $PEP_DIR" >&2; exit 1; }
if [ ! -f "$OUTDIR/$PEP_FILE" ]; then
  echo "[egg] downloading proteome $PEP_FILE"
  curl -s "$PEP_DIR/$PEP_FILE" -o "$OUTDIR/$PEP_FILE"
fi
PEP_FA="$OUTDIR/proteome.fa"
gzip -dc "$OUTDIR/$PEP_FILE" > "$PEP_FA"
echo "[egg] proteins: $(grep -c '^>' "$PEP_FA")"

# --- 2. eggNOG-mapper conda env ---------------------------------------------
if ! conda env list | grep -q '^eggnog '; then
  echo "[egg] creating conda env 'eggnog' (eggnog-mapper + diamond) ..."
  if ! conda create -n eggnog -y -c conda-forge -c bioconda eggnog-mapper; then
    echo "[egg] arm64 solve failed; retrying with osx-64 (Rosetta) ..."
    CONDA_SUBDIR=osx-64 conda create -n eggnog -y -c conda-forge -c bioconda eggnog-mapper
  fi
fi
run() { conda run -n eggnog "$@"; }

# --- 3. download eggNOG database (one-time) ---------------------------------
# NOTE: the tool's built-in URL (eggnogdb.embl.de) is dead; the working mirror
# is eggnog5.embl.de. We download the 3 files needed for diamond mode directly.
# Uncompressed size is ~50 GB (eggnog.db is the bulk). curl -C - resumes.
DBHOST="http://eggnog5.embl.de/download/emapperdb-5.0.2"
# resumable download that auto-retries through connection drops (curl error 18)
fetch() {  # url, dest
  local url="$1" dest="$2" n=0
  until curl -fL -C - --retry 5 --retry-delay 5 --connect-timeout 30 -o "$dest" "$url"; do
    n=$((n+1)); [ "$n" -ge 30 ] && { echo "[egg] giving up after 30 tries: $url" >&2; return 1; }
    echo "[egg] connection dropped; resuming (attempt $n) ..."; sleep 5
  done
}
dl_gz() {  # url, final_uncompressed_path
  local url="$1" final="$2" gz="$2.gz"
  [ -f "$final" ] && { echo "[egg] present: $(basename "$final")"; return; }
  echo "[egg] downloading $(basename "$url") (resumable) ..."
  fetch "$url" "$gz" && echo "[egg] decompressing $(basename "$gz") ..." && gunzip "$gz"
}
if [ ! -f "$EGGDB/eggnog.db" ]; then dl_gz "$DBHOST/eggnog.db.gz" "$EGGDB/eggnog.db"; fi
if [ ! -f "$EGGDB/eggnog_proteins.dmnd" ]; then dl_gz "$DBHOST/eggnog_proteins.dmnd.gz" "$EGGDB/eggnog_proteins.dmnd"; fi
if [ ! -f "$EGGDB/eggnog.taxa.db" ]; then
  echo "[egg] downloading eggnog.taxa.tar.gz (resumable) ..."
  fetch "$DBHOST/eggnog.taxa.tar.gz" "$EGGDB/eggnog.taxa.tar.gz" \
    && tar -zxf "$EGGDB/eggnog.taxa.tar.gz" -C "$EGGDB" && rm -f "$EGGDB/eggnog.taxa.tar.gz"
fi

# --- 4. run emapper (diamond mode) ------------------------------------------
if [ ! -f "$ANNO_OUT/ztritici.emapper.annotations" ]; then
  echo "[egg] running eggNOG-mapper (diamond) ..."
  run emapper.py -i "$PEP_FA" --itype proteins -m diamond \
      --data_dir "$EGGDB" --output_dir "$ANNO_OUT" -o ztritici \
      --cpu "$THREADS" --override
fi
echo "[egg] annotations: $ANNO_OUT/ztritici.emapper.annotations"

# --- 5. GO + KEGG term-name tables ------------------------------------------
if [ ! -f "$OUTDIR/go_terms.csv" ]; then
  echo "[egg] fetching GO term names (go-basic.obo) ..."
  curl -s "https://purl.obolibrary.org/obo/go/go-basic.obo" -o "$OUTDIR/_go.obo"
fi
if [ ! -f "$OUTDIR/kegg_pathways.csv" ]; then
  echo "[egg] fetching KEGG pathway names ..."
  curl -s "https://rest.kegg.jp/list/pathway" -o "$OUTDIR/_kegg.tsv"
fi

# --- 6. assemble gene-level annotation (python) -----------------------------
echo "[egg] building data/reference/gene_annotation.csv ..."
python3 - "$OUTDIR" "$ANNO_OUT" "$PEP_FA" <<'PY'
import sys, os, csv, collections, re
outdir, anno_out, pep_fa = sys.argv[1], sys.argv[2], sys.argv[3]

# protein -> gene from pep FASTA headers ( ... gene:Mycgr3G..... )
prot2gene = {}
with open(pep_fa) as f:
    for line in f:
        if not line.startswith(">"): continue
        pid = line[1:].split()[0]
        m = re.search(r"gene:(\S+)", line)
        prot2gene[pid] = m.group(1) if m else pid

# parse emapper annotations (header line starts with '#query')
ann_file = os.path.join(anno_out, "ztritici.emapper.annotations")
rows = []; header = None
with open(ann_file) as f:
    for line in f:
        if line.startswith("#query") or line.startswith("query"):
            header = line.lstrip("#").rstrip("\n").split("\t"); continue
        if line.startswith("#") or not line.strip(): continue
        if header is None: continue
        rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
print(f"[egg] annotated proteins: {len(rows)}")

def clean(v): return "" if v in (None, "-", "") else v
def split_terms(v):
    v = clean(v)
    return [x for x in re.split(r"[,;]", v) if x and x != "-"]

# aggregate per gene
g = collections.defaultdict(lambda: dict(name="", desc="", cog=set(), go=set(),
                                         kegg=set(), ec=set(), pfam=set()))
for r in rows:
    gene = prot2gene.get(r.get("#query", r.get("query","")), r.get("#query",""))
    rec = g[gene]
    if not rec["name"]: rec["name"] = clean(r.get("Preferred_name",""))
    if not rec["desc"]: rec["desc"] = clean(r.get("Description",""))
    rec["cog"] |= set(clean(r.get("COG_category","")))
    rec["go"]  |= set(split_terms(r.get("GOs","")))
    # KEGG_Pathway holds ko##### and map#####; keep map##### (reference pathways)
    rec["kegg"]|= set(x for x in split_terms(r.get("KEGG_Pathway","")) if x.startswith("map"))
    rec["ec"]  |= set(split_terms(r.get("EC","")))
    rec["pfam"]|= set(split_terms(r.get("PFAMs","")))

with open(os.path.join(outdir, "gene_annotation.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene_id","gene_name","description","COG","GO","KEGG","EC","PFAM"])
    for gene in sorted(g):
        r = g[gene]
        w.writerow([gene, r["name"], r["desc"], "".join(sorted(r["cog"])),
                    ";".join(sorted(r["go"])), ";".join(sorted(r["kegg"])),
                    ";".join(sorted(r["ec"])), ";".join(sorted(r["pfam"]))])
ng = len(g)
print(f"[egg] genes annotated: {ng}  with GO: {sum(1 for x in g.values() if x['go'])}"
      f"  with KEGG: {sum(1 for x in g.values() if x['kegg'])}")

# GO term names from obo
obo = os.path.join(outdir, "_go.obo")
if os.path.exists(obo):
    with open(os.path.join(outdir,"go_terms.csv"),"w",newline="") as f:
        w=csv.writer(f); w.writerow(["go_id","go_name","go_namespace"])
        cur={}
        for line in open(obo):
            line=line.rstrip("\n")
            if line=="[Term]": cur={}
            elif line.startswith("id: GO:"): cur["id"]=line[4:]
            elif line.startswith("name: "): cur["name"]=line[6:]
            elif line.startswith("namespace: "):
                cur["ns"]=line[11:]
                if cur.get("id"): w.writerow([cur["id"],cur.get("name",""),cur["ns"]])
    os.remove(obo)

# KEGG pathway names from REST list
kegg=os.path.join(outdir,"_kegg.tsv")
if os.path.exists(kegg):
    with open(os.path.join(outdir,"kegg_pathways.csv"),"w",newline="") as f:
        w=csv.writer(f); w.writerow(["kegg_id","kegg_name"])
        for line in open(kegg):
            p=line.rstrip("\n").split("\t")
            if len(p)<2: continue
            kid=p[0].replace("path:","")          # e.g. map00100
            w.writerow([kid, p[1]])
    os.remove(kegg)
print("[egg] wrote go_terms.csv + kegg_pathways.csv")
PY

rm -f "$PEP_FA"
echo "[egg] DONE."
echo "  -> data/reference/gene_annotation.csv  (GO + KEGG + EC + Pfam + COG + description)"
echo "  -> data/reference/go_terms.csv , data/reference/kegg_pathways.csv"
echo "Next: conda activate rnaseq-r ; Rscript scripts/04_annotation_enrichment.R"
