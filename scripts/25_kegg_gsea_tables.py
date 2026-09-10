#!/usr/bin/env python3
# =============================================================================
# 06_module_enrichment.py   —   DATA-DRIVEN pathway & gene selection for Fig 4
#
# No hand-picked pathways. Selection is entirely quantitative:
#
#   UNIVERSE (declared a priori, category-level):
#     KEGG Metabolism (map00xxx) + Genetic Information Processing (map03xxx)
#     + a short whitelist of degradation/transport cellular pathways
#     (peroxisome, autophagy, mitophagy, protein processing in ER,
#      ubiquitin-mediated proteolysis, ABC transporters).
#     EXCLUDED: Human Diseases (map05xxx), Organismal/Environmental systems
#     (other map04xxx, map02xxx), and KEGG global/overview maps (map01xxx).
#
#   PATHWAYS: GSEA (gseapy prerank, genes ranked by DESeq2 Wald = lfc/lfcSE)
#     over the restricted universe, per temperature. Selection rule:
#     pathways with BH-FDR < 0.05 in >=1 temperature, top 15 by max |NES|.
#
#   GENES: significant at 15 C (padj<0.05), top 12 up + 12 down by log2FC.
#
# Outputs (results/...):
#   04_differential_expression/kegg_gsea_all_<T>C.csv     full restricted GSEA
#   04_differential_expression/figure4C_pathways.csv      selected pathways
#   04_differential_expression/figure4_top_genes.csv      selected genes
#   05_figures/Figure4_datadriven.png / .pdf              2-panel figure
# =============================================================================
import csv, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import TwoSlopeNorm
import gseapy as gp

DE_DIR="tables/rnaseq"; FIG_DIR="figures/diagnostics"
TEMPS=["15","21","27"]; SEED,PERM,MIN_SIZE,MAX_SIZE=42,1000,10,500

# ---- restricted-universe rule ------------------------------------------------
GIP={"map03020","map03022","map03040","map03010","map03013","map03015","map03008",
     "map03060","map03050","map03018","map03030","map03410","map03420","map03430",
     "map03440","map03450","map03460"}            # Genetic Information Processing
CELL={"map04146","map04138","map04139","map04141","map04120","map04122","map04130","map02010"}
def in_universe(mid):
    if mid.startswith("map00"): return True         # Metabolism (excl. global map01 overviews)
    if mid in GIP:               return True         # Genetic Information Processing (explicit)
    if mid in CELL:              return True         # selected cellular / transport
    return False                                    # excl. map01 overview, organismal, diseases

# ---- annotation --------------------------------------------------------------
anno={r["gene_id"]:r for r in csv.DictReader(open("data/reference/gene_annotation.csv"))}
gene_kegg={g:set((anno[g].get("KEGG") or "").split(";"))-{""} for g in anno}
kegg_name={}
try:
    for r in csv.DictReader(open("data/reference/kegg_pathways.csv")): kegg_name[r["kegg_id"]]=r["kegg_name"]
except FileNotFoundError: pass
kegg_sets={}
for g,ks in gene_kegg.items():
    for k in ks:
        if in_universe(k): kegg_sets.setdefault(k,[]).append(g)
def kg(mid): return kegg_sets.get(mid,[])

def load_de(name):
    d={}
    for r in csv.DictReader(open(f"{DE_DIR}/DE_{name}.csv")):
        try:
            _nm=(r.get("gene_name") or "").strip(); _ds=(r.get("description") or "").strip()
            if _nm.upper() in ("NA","NAN","."): _nm=""
            if _ds.upper() in ("NA","NAN","."): _ds=""
            d[r["gene_id"]]=dict(lfc=float(r["log2FoldChange"]),se=float(r["lfcSE"]),
                                 padj=float(r["padj"]),name=_nm,desc=_ds)
        except (ValueError,KeyError): pass
    return d
DE={t:load_de(f"proth_at_{t}C_2vs0") for t in TEMPS}

def ranking(de):
    rows=[(g,v["lfc"]/v["se"]) for g,v in de.items() if v["se"] and np.isfinite(v["lfc"]/v["se"])]
    return pd.DataFrame(rows,columns=["gene","stat"]).sort_values("stat",ascending=False).reset_index(drop=True)

def gsea_cols(res):
    c=list(res.columns)
    return ("Term" if "Term" in c else c[0],
            next(x for x in c if x.upper()=="NES"),
            next(x for x in c if "FDR" in x.upper()))

# ============================== GSEA (restricted) =============================
NES={}; FDR={}
for t in TEMPS:
    pre=gp.prerank(rnk=ranking(DE[t]),gene_sets=kegg_sets,min_size=MIN_SIZE,max_size=MAX_SIZE,
                   permutation_num=PERM,seed=SEED,outdir=None,no_plot=True,verbose=False)
    res=pre.res2d.copy(); term,nes,q=gsea_cols(res)
    res=res.rename(columns={term:"kegg_id",nes:"NES",q:"FDR_q"})
    res["NES"]=pd.to_numeric(res["NES"],errors="coerce"); res["FDR_q"]=pd.to_numeric(res["FDR_q"],errors="coerce")
    res["kegg_name"]=res["kegg_id"].map(lambda k:kegg_name.get(k,"")); res["n_genes"]=res["kegg_id"].map(lambda k:len(kg(k)))
    res.sort_values("FDR_q")[["kegg_id","kegg_name","n_genes","NES","FDR_q"]].to_csv(f"{DE_DIR}/kegg_gsea_all_{t}C.csv",index=False)
    for r in res.itertuples(): NES[(r.kegg_id,t)]=r.NES; FDR[(r.kegg_id,t)]=r.FDR_q
    print(f"[GSEA restricted] {t} C: {len(res)} pathways, {int((res.FDR_q<0.05).sum())} sig (FDR<0.05)")

# --------- pathway selection at two FDR cutoffs ------------------------------
def mean_lfc(k,t):
    v=[DE[t][g]["lfc"] for g in kg(k) if g in DE[t]]; return float(np.mean(v)) if v else np.nan
allp={k for (k,t) in NES}
def build(cut):
    sel=sorted([k for k in allp if any(FDR.get((k,t),1)<cut for t in TEMPS)], key=lambda k:mean_lfc(k,"15"))
    rows=[]
    for k in sel:
        row=dict(kegg_id=k,kegg_name=kegg_name.get(k,k),n_genes=len(kg(k)))
        for t in TEMPS: row[f"meanLFC_{t}"]=mean_lfc(k,t); row[f"NES_{t}"]=NES.get((k,t),np.nan); row[f"FDR_{t}"]=FDR.get((k,t),np.nan)
        rows.append(row)
    return pd.DataFrame(rows)
pathdf    = build(0.001)   # main Figure 4C (strict)
pathdf_s1 = build(0.05)    # Supplementary Figure S1 (fuller)
pathdf.to_csv(f"{DE_DIR}/figure4C_pathways.csv",index=False)
pathdf_s1.to_csv(f"{DE_DIR}/figureS1_pathways.csv",index=False)

# --------- select top genes by |log2FC| (sig at 15 C) -------------------------
LFC_CUT=4.0   # threshold-driven: padj<0.05 AND |log2FC|>=4 at 15 C (no top-N)
sig15=[g for g,v in DE["15"].items() if v["padj"]<0.05 and abs(v["lfc"])>=LFC_CUT]
up=sorted([g for g in sig15 if DE["15"][g]["lfc"]>0],key=lambda g:-DE["15"][g]["lfc"])
dn=sorted([g for g in sig15 if DE["15"][g]["lfc"]<0],key=lambda g:-DE["15"][g]["lfc"])
genes=up+dn
def glabel(g):
    v=DE["15"][g]; lab=v["name"] or (v["desc"][:28] if v["desc"] else g)
    return lab if len(lab)<=30 else lab[:27]+"..."
grows=[]
for g in genes:
    r=dict(gene_id=g,label=glabel(g))
    for t in TEMPS: r[f"lfc_{t}"]=DE[t][g]["lfc"] if g in DE[t] else np.nan; r[f"padj_{t}"]=DE[t][g]["padj"] if g in DE[t] else np.nan
    grows.append(r)
genedf=pd.DataFrame(grows); genedf.to_csv(f"{DE_DIR}/figure4_top_genes.csv",index=False)

# ---------------------------------------------------------------------------
# NOTE: this script now stops at the tables. The figures it used to draw
# (Figure4_datadriven, FigureS2_genes_wide) are superseded by the R figure
# scripts, which read the tables written above. Every manuscript figure is
# generated by R; see scripts/make_main_figures.R and make_supp_figures.R.
# ---------------------------------------------------------------------------
print("KEGG GSEA tables written:",
      "kegg_gsea_all_{15,21,27}C.csv, figure4C_pathways.csv,",
      "figureS1_pathways.csv, figure4_top_genes.csv")
