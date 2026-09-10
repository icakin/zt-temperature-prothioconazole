#!/usr/bin/env python3
# =============================================================================
# 32_inference_upgrades.py — inferential upgrades (existing data only):
#  U1 sample-level temperature x drug interaction on the RIBOSOME module score
#     (n=45 biological samples; replaces gene-level pseudoreplicated tests) + AOX.
#  U2 pathway-level enrichment on the EXISTING factorial interaction statistic
#     (DESeq2 interaction contrast 27 vs 15; GSEA prerank; correlation caveat noted).
#  U3 state-convergence S(27): sample-resampling bootstrap CI + sample-label
#     permutation null (preserves gene-gene correlation).
#  U4 posterior dose-vs-control contrasts for E_a and for the trait gap.
# Outputs -> tables/revision/inference/.
# =============================================================================
import os, sys, types, re, json
import numpy as np, pandas as pd
sys.modules.setdefault("gseapy_placeholder", None)
import statsmodels.formula.api as smf
OUT="tables/revision/inference"; os.makedirs(OUT,exist_ok=True)
LOG=open(f"{OUT}/81_log.txt","w")
def log(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.write(s+"\n"); LOG.flush()
rng=np.random.default_rng(7)

# ---- shared data ----
cnt=pd.read_csv("tables/rnaseq_stage/normalized_counts.csv",index_col=0)
de21=pd.read_csv("tables/rnaseq_stage/DE_proth_at_21C_2vs0.csv")
desc=de21["description"].fillna("").str.lower(); nm=de21["gene_name"].fillna("").str.lower()
kegg=de21["KEGG"].fillna("").astype(str); pfam=de21["PFAM"].fillna("").str.lower()
# ribosome module = KEGG ribosome pathway (map03010), identical to Fig 5C / module-trait analysis
ribo=[g for g in de21[kegg.str.contains("map03010")]["gene_id"] if g in cnt.index]
meta=pd.DataFrame([re.match(r"12129_(\d+)_(\d+)_R(\d+)",c).groups() for c in cnt.columns],
                  columns=["T","D","rep"],index=cnt.columns)
logn=np.log2(cnt+1)

# ===== U1: sample-level interaction on module scores =====
log("=== U1: sample-level temperature x drug interaction (module scores) ===")
def interaction_test(series,label):
    df=meta.copy(); df["y"]=series
    m=smf.ols("y ~ C(T)*C(D)",data=df).fit()
    m0=smf.ols("y ~ C(T)+C(D)",data=df).fit()
    from statsmodels.stats.anova import anova_lm
    a=anova_lm(m0,m)
    F=a["F"].iloc[1]; p=a["Pr(>F)"].iloc[1]
    log(f"[{label}] interaction F({int(a['df_resid'].iloc[0]-a['df_resid'].iloc[1])},{int(a['df_resid'].iloc[1])}) = {F:.2f}, p = {p:.2e}")
    # planned contrast: drug effect at 27 minus drug effect at 15 (and vs 21)
    for tt in ["21","27"]:
        c=m.params.get(f"C(T)[T.{tt}]:C(D)[T.2]",np.nan)
        se=m.bse.get(f"C(T)[T.{tt}]:C(D)[T.2]",np.nan)
        log(f"[{label}] (drug effect @{tt}) - (drug effect @15) = {c:+.3f} ± {1.96*se:.3f} (95% CI)")
    return m
ribo_score=logn.loc[ribo].mean(axis=0)
interaction_test(ribo_score,f"Ribosome module ({len(ribo)} genes, per-sample mean log2)")
aox=logn.loc["Mycgr3G72918"]
interaction_test(aox,"AOX (single gene, log2)")

# ===== U2: pathway enrichment on the existing interaction statistic =====
log("\n=== U2: GSEA on DESeq2 interaction statistic (27 vs 15 contrast) ===")
import gseapy as gp
inter=pd.read_csv("tables/rnaseq/DE_interaction_temperature27_prothioconazole2.csv")
inter=inter[inter["stat"].notna()]
rnk=inter[["gene_id","stat"]].rename(columns={"gene_id":"gene"}).sort_values("stat",ascending=False)
# same restricted universe as fig_common
import csv
anno={r["gene_id"]:r for r in csv.DictReader(open("data/reference/gene_annotation.csv"))}
gk={g:set((anno[g].get("KEGG") or "").split(";"))-{""} for g in anno}
GIP={"map03020","map03022","map03040","map03010","map03013","map03015","map03008","map03060",
     "map03050","map03018","map03030","map03410","map03420","map03430","map03440","map03450","map03460"}
CELL={"map04146","map04138","map04139","map04141","map04120","map04122","map04130","map02010"}
inuni=lambda m: m.startswith("map00") or m in GIP or m in CELL
ksets={}
for g,ks in gk.items():
    for k in ks:
        if inuni(k): ksets.setdefault(k,[]).append(g)
pre=gp.prerank(rnk=rnk,gene_sets=ksets,min_size=10,max_size=500,permutation_num=1000,
               seed=42,outdir=None,no_plot=True,verbose=False)
res=pre.res2d.rename(columns=lambda c:c.strip())
term="Term" if "Term" in res.columns else [c for c in res.columns if c.lower()=="term"][0]
nes=[c for c in res.columns if c.upper()=="NES"][0]; fdr=[c for c in res.columns if "FDR" in c.upper()][0]
res=res[[term,nes,fdr]].rename(columns={term:"kegg_id",nes:"NES",fdr:"FDR"})
res["NES"]=pd.to_numeric(res.NES,errors="coerce"); res["FDR"]=pd.to_numeric(res.FDR,errors="coerce")
res=res.sort_values("NES")
res.to_csv(f"{OUT}/interaction_gsea_27v15.csv",index=False)
top=res[res.FDR<0.05]
log(f"pathways with interaction FDR<0.05: {len(top)}")
for r_ in top.itertuples():
    log(f"  {r_.kegg_id}  NES={r_.NES:+.2f}  FDR={r_.FDR:.4f}")
rb=res[res.kegg_id=="map03010"]
if len(rb): log(f"RIBOSOME interaction: NES={rb.NES.iloc[0]:+.2f}, FDR={rb.FDR.iloc[0]:.4f} "
                "(positive NES = drug repression attenuated at 27 C)")

# ===== U3: state projection S(27) - sample bootstrap CI + label permutation =====
log("\n=== U3: state-convergence S(27) with sample-level statistics ===")
def cols(T,D): return [c for c in cnt.columns if re.match(fr"12129_{T}_{D}_R",c)]
def unit(v): n=np.linalg.norm(v); return v/n if n>0 else v
genes=[g for g in de21.gene_id if g in cnt.index]
bm=cnt.loc[genes].min(axis=1)
genes=[g for g in genes if True]
L=logn.loc[genes]
def signature(c15_0,c15_2,c21_0,c21_2):
    d15=L[c15_2].mean(1)-L[c15_0].mean(1); d21=L[c21_2].mean(1)-L[c21_0].mean(1)
    return unit(unit(d15.values)+unit(d21.values))
def S27(d, c15,c27):
    h=L[c27].mean(1).values-L[c15].mean(1).values
    return float(np.dot(d,h)/np.dot(d,d))
c15_0,c15_2,c21_0,c21_2,c27_0=cols("15","0"),cols("15","2"),cols("21","0"),cols("21","2"),cols("27","0")
d=signature(c15_0,c15_2,c21_0,c21_2)
S=S27(d,c15_0,c27_0)
# scale to % of heat displacement (cos-based variance explained as before)
h=L[c27_0].mean(1).values-L[c15_0].mean(1).values
cos=np.dot(d,h)/(np.linalg.norm(d)*np.linalg.norm(h))
log(f"S(27) = {S:.2f}; cos(drug,heat) = {cos:.3f}; variance explained = {cos**2:.3f}")
# bootstrap over samples (resample within each condition; recompute signature+projection)
B=500; bs=[]
for b in range(B):
    r15_0=list(rng.choice(c15_0,len(c15_0),replace=True)); r15_2=list(rng.choice(c15_2,len(c15_2),replace=True))
    r21_0=list(rng.choice(c21_0,len(c21_0),replace=True)); r21_2=list(rng.choice(c21_2,len(c21_2),replace=True))
    r27_0=list(rng.choice(c27_0,len(c27_0),replace=True))
    db=signature(r15_0,r15_2,r21_0,r21_2)
    hb=L[r27_0].mean(1).values-L[r15_0].mean(1).values
    bs.append(np.dot(db,hb)/(np.linalg.norm(db)*np.linalg.norm(hb)))
log(f"cos bootstrap 95% CI (samples resampled within condition): [{np.percentile(bs,2.5):.3f}, {np.percentile(bs,97.5):.3f}]")
# sample-label permutation: shuffle which control samples are '15' vs '27' (preserves gene-gene corr)
ctrl=c15_0+c27_0; n15=len(c15_0); perm=[]
for b in range(1000):
    idx=rng.permutation(len(ctrl))
    p15=[ctrl[i] for i in idx[:n15]]; p27=[ctrl[i] for i in idx[n15:]]
    hp=L[p27].mean(1).values-L[p15].mean(1).values
    perm.append(np.dot(d,hp)/(np.linalg.norm(d)*np.linalg.norm(hp)))
pval=float(np.mean(np.array(perm)>=cos))
log(f"sample-label permutation p (cos >= observed): {pval:.4f} (perm 95th pct {np.percentile(perm,95):.3f})")

# ===== U4: posterior E_a contrasts =====
log("\n=== U4: E_a dose-vs-control and trait-gap contrasts (posterior draws) ===")
dg=pd.read_csv("tables/physiology/posterior_sharpe_schoolfield_tpc_growth_C_per_C_h_by_dose.csv")
dr=pd.read_csv("tables/physiology/posterior_sharpe_schoolfield_tpc_respiration_C_per_C_h_by_dose.csv")
def norm_dose(x):
    x=str(x); return "Control" if x in ("0","Control") else x.rstrip("0").rstrip(".") if "." in x else x
for df in (dg,dr): df["Dose"]=df["Dose"].map(norm_dose)
piv_g=dg.pivot_table(index=".draw",columns="Dose",values="E")
piv_r=dr.pivot_table(index=".draw",columns="Dose",values="E")
rows=[]
for dcol in [c for c in piv_g.columns if c!="Control"]:
    dgc=piv_g[dcol]-piv_g["Control"]; drc=piv_r[dcol]-piv_r["Control"]
    gap=(piv_r[dcol]-piv_g[dcol])-(piv_r["Control"]-piv_g["Control"])
    rows.append(dict(dose=dcol,
        dEa_growth=f"{dgc.median():+.2f} [{dgc.quantile(.025):+.2f},{dgc.quantile(.975):+.2f}]",
        dEa_resp=f"{drc.median():+.2f} [{drc.quantile(.025):+.2f},{drc.quantile(.975):+.2f}]",
        dGap=f"{gap.median():+.2f} [{gap.quantile(.025):+.2f},{gap.quantile(.975):+.2f}]"))
T=pd.DataFrame(rows); T.to_csv(f"{OUT}/Ea_dose_contrasts.csv",index=False)
log(T.to_string(index=False))
log("\nwrote -> tables/revision/inference/")
LOG.close()
