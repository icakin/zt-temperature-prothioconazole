#!/usr/bin/env python3
# 30_verify_ribosome.py — hardening tests for the ribosome-sector mechanism.
#   V1 per-gene temperature modulation: is ribosomal repression significantly
#      stronger at 15 than 27 across the 153 genes? (paired Wilcoxon)
#   V2 specificity: are ribosomal genes more repressed than the genome-wide
#      background at each temperature? (Mann-Whitney vs all genes)
#   V3 growth-law: does observed log(growth drug/ctrl) track log(ribosome
#      drug/ctrl) with slope ~1? bootstrap CI over genes & vials.
#   V4 robustness: does the temperature pattern hold across 3 gene-set definitions?
#   INT integrated per-condition prediction (etcGEM temperature x ribosome sector).
import os, sys, json, logging, re
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__)); logging.getLogger("cobra").setLevel(logging.ERROR)
OUT="tables/revision/v3build"; LOG=open(os.path.join(OUT,"72_log.txt"),"w")
def log(*a): s=" ".join(str(x) for x in a); print(s); LOG.write(s+"\n"); LOG.flush()

# per-sample counts -> per-condition replicate sums for ribosomal genes
counts=pd.read_csv("tables/rnaseq_stage/normalized_counts.csv",index_col=0)
de=pd.read_csv("tables/rnaseq_stage/DE_proth_at_21C_2vs0.csv")
desc=de["description"].fillna("").str.lower(); nm=de["gene_name"].fillna("").str.lower()
kegg=de["KEGG"].fillna("").astype(str); pfam=de["PFAM"].fillna("").str.lower()
setmap={
 "KEGG_ko03010": de[kegg.str.contains("map03010")]["gene_id"].tolist(),
 "PFAM_ribosomal": de[pfam.str.contains("ribosomal_")]["gene_id"].tolist(),
 "name_rpl/rps": de[nm.str.match(r"rp[ls]")]["gene_id"].tolist(),
 "union": de[(desc.str.contains("ribosomal protein"))|(kegg.str.contains("map03010"))|
             (pfam.str.contains("ribosomal_"))|(nm.str.match(r"rp[ls]"))]["gene_id"].tolist(),
}
ribo=[gname for gname in setmap["union"] if gname in counts.index]
log(f"ribosomal genes (union, in counts): {len(ribo)}")
def cols(T,D): return [c for c in counts.columns if re.match(fr"12129_{T}_{D}_R",c)]

# ---- V1: per-gene repression by temperature + paired test 15 vs 27 ----
def gene_lfc(T):  # log2 drug/ctrl per ribosomal gene (condition means)
    c0=counts.loc[ribo,cols(T,0)].mean(axis=1); c2=counts.loc[ribo,cols(T,2)].mean(axis=1)
    return np.log2((c2+1)/(c0+1))
l15,l21,l27=gene_lfc(15),gene_lfc(21),gene_lfc(27)
log("\n[V1] median ribosomal log2FC (drug/ctrl): 15C=%.2f 21C=%.2f 27C=%.2f"%(l15.median(),l21.median(),l27.median()))
w=stats.wilcoxon(l15,l27); log(f"[V1] paired Wilcoxon 15 vs 27: stat={w.statistic:.0f} p={w.pvalue:.2e} "
    f"(more repressed at 15 in {int((l15<l27).sum())}/{len(ribo)} genes)")

# ---- V2: specificity vs genome-wide background ----
allg=[gname for gname in counts.index]
def all_lfc(T):
    c0=counts.loc[allg,cols(T,0)].mean(axis=1); c2=counts.loc[allg,cols(T,2)].mean(axis=1)
    return np.log2((c2+1)/(c0+1))
for T,lr in [(15,l15),(21,l21),(27,l27)]:
    bg=all_lfc(T); mw=stats.mannwhitneyu(lr,bg,alternative="less")
    log(f"[V2] {T}C ribosomal median {lr.median():.2f} vs background {bg.median():.2f}  MWU p={mw.pvalue:.2e}")

# ---- V3: growth-law slope (bootstrap) ----
d=pd.read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv")
d["conc"]=d["Dose"].replace("Control","0").astype(float); d["Tc"]=d["T"].astype(int)
d=d[d.growth_C_per_C_h>0]
def ribo_ratio(T): return counts.loc[ribo,cols(T,2)].sum().mean()/counts.loc[ribo,cols(T,0)].sum().mean()
def growth_ratio(T):
    g0=d[(d["Tc"]==T)&(d["conc"]==0)]["growth_C_per_C_h"].mean()
    g2=d[(d["Tc"]==T)&(d["conc"]==2)]["growth_C_per_C_h"].mean(); return g2/g0
xs=np.array([np.log(ribo_ratio(T)) for T in [15,21,27]]); ys=np.array([np.log(growth_ratio(T)) for T in [15,21,27]])
slope=np.sum(xs*ys)/np.sum(xs*xs)  # through-origin (growth law)
log(f"\n[V3] growth-law through-origin slope beta={slope:.2f} (1.0 = growth proportional to ribosomes)")
log(f"[V3] ribosome ratio vs growth ratio by T: "+", ".join(f"{T}: R={np.exp(x):.2f} g={np.exp(y):.2f}" for T,x,y in zip([15,21,27],xs,ys)))
# bootstrap beta over ribosomal genes
rng=np.random.default_rng(0); bs=[]
for _ in range(1000):
    rb=list(rng.choice(ribo,len(ribo),replace=True))
    xb=np.array([np.log(counts.loc[rb,cols(T,2)].sum().mean()/counts.loc[rb,cols(T,0)].sum().mean()) for T in [15,21,27]])
    bs.append(np.sum(xb*ys)/np.sum(xb*xb))
log(f"[V3] beta 95%% CI over gene bootstrap: [{np.percentile(bs,2.5):.2f}, {np.percentile(bs,97.5):.2f}]")

# ---- V4: robustness across gene-set definitions ----
log("\n[V4] temperature pattern across gene-set definitions (repression %):")
for name,gs in setmap.items():
    gs=[gname for gname in gs if gname in counts.index]
    if len(gs)<10: continue
    rep={T:100*(1-counts.loc[gs,cols(T,2)].sum().mean()/counts.loc[gs,cols(T,0)].sum().mean()) for T in [15,21,27]}
    log(f"   {name:16s} (n={len(gs):3d}): 15C={rep[15]:.0f}%  21C={rep[21]:.0f}%  27C={rep[27]:.0f}%  hold(15>27):{rep[15]>rep[27]}")

json.dump(dict(n_ribo=len(ribo),wilcox_p=float(w.pvalue),beta=float(slope),
               beta_ci=[float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))]),
          open(os.path.join(OUT,"72_verify.json"),"w"),indent=2)
log("\nwrote 72_verify.json"); LOG.close()
