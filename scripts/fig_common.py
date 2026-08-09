#!/usr/bin/env python3
# =============================================================================
# fig_common.py — shared style, data loaders, GSEA and plotting helpers used by
# Figures 4, 5 and 6. Keeping these in one place guarantees the figures stay
# consistent (identical KEGG categories, identical GSEA universe, identical
# volcano styling) and removes the duplication that had accumulated.
# =============================================================================
import csv, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as _pe
try:
    from adjustText import adjust_text; HAVE_AT=True
except Exception:
    HAVE_AT=False
import gseapy as gp

DE_DIR="tables/rnaseq"; FIG_DIR="figures/diagnostics"
SEED,PERM,MIN_SIZE,MAX_SIZE=42,1000,10,500
TEMP={"15":"#27519E","21":"#64A5DE","27":"#D64B21"}      # temperature palette (= fig_style ramp members)
UP,DN,DIV="#B2182B","#2166AC","#E08214"                  # shared-up / shared-down / divergent
STROKE=[_pe.withStroke(linewidth=1.7,foreground="white")]

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["DejaVu Sans"],
                     "font.size":9,"axes.linewidth":0.6,"pdf.fonttype":42,"ps.fonttype":42})

# ---- annotation -------------------------------------------------------------
anno={r["gene_id"]:r for r in csv.DictReader(open("data/reference/gene_annotation.csv"))}
gk={g:set((anno[g].get("KEGG") or "").split(";"))-{""} for g in anno}
kg=lambda mid:[g for g,k in gk.items() if mid in k]
aoxg=[g for g in anno if "alternative oxidase" in (anno[g].get("description","") or "").lower()]
kegg_name={}
if os.path.exists("data/reference/kegg_pathways.csv"):
    for r in csv.DictReader(open("data/reference/kegg_pathways.csv")): kegg_name[r["kegg_id"]]=r["kegg_name"]

# ---- shared KEGG mechanistic categories (identical in Fig 4 & Fig 5) --------
CATS=[("Efflux (ABC transporters)","#7B3FA0",set(kg("map02010"))),
      ("Ergosterol / sterol","#1B7837",set(kg("map00100"))),
      ("Cytochrome-P450 detox","#D6604D",set(kg("map00982"))|set(kg("map00980"))),
      ("Respiration (OXPHOS / AOX)","#E08214",set(kg("map00190"))|set(aoxg)),
      ("Translation (ribosome)","#2F6BB3",set(kg("map03010"))),
      ("Proteasome","#11838E",set(kg("map03050"))),
      ("TCA cycle","#8C510A",set(kg("map00020")))]
CATSET={g for _,_,s in CATS for g in s}

def sym(g):
    """Real gene symbol only (no orphan IDs, no descriptions); None if absent."""
    n=(anno.get(g,{}).get("gene_name") or "").strip()
    return n if n and n.upper() not in ("NA","NAN") else None

def stars(q): return "***" if q<0.001 else "**" if q<0.01 else "*" if q<0.05 else ""

from scipy import stats as _stats
def ci95(v):
    """Half-width of the 95% confidence interval of the mean (t-based)."""
    v=np.asarray(v,dtype=float); n=len(v)
    return float(v.std(ddof=1)/np.sqrt(n)*_stats.t.ppf(0.975,n-1)) if n>1 else 0.0

# ---- DE loaders -------------------------------------------------------------
def load_de(name):    # {gene:(lfc,padj)}
    d={}
    for r in csv.DictReader(open(f"{DE_DIR}/DE_{name}.csv")):
        try: d[r["gene_id"]]=(float(r["log2FoldChange"]),float(r["padj"]))
        except (ValueError,KeyError): pass
    return d
def load_full(name):  # {gene:{lfc,padj,se}}
    d={}
    for r in csv.DictReader(open(f"{DE_DIR}/DE_{name}.csv")):
        try: d[r["gene_id"]]=dict(lfc=float(r["log2FoldChange"]),padj=float(r["padj"]),se=float(r["lfcSE"]))
        except (ValueError,KeyError): pass
    return d

# ---- restricted KEGG universe + GSEA (identical across figures) -------------
_GIP={"map03020","map03022","map03040","map03010","map03013","map03015","map03008",
      "map03060","map03050","map03018","map03030","map03410","map03420","map03430",
      "map03440","map03450","map03460"}                         # Genetic Information Processing
_CELL={"map04146","map04138","map04139","map04141","map04120","map04122","map04130","map02010"}
in_universe=lambda m: m.startswith("map00") or m in _GIP or m in _CELL
kegg_sets={}
for _g,_ks in gk.items():
    for _k in _ks:
        if in_universe(_k): kegg_sets.setdefault(_k,[]).append(_g)

def _ranking(de_full):
    rows=[(g,v["lfc"]/v["se"]) for g,v in de_full.items() if v["se"] and np.isfinite(v["lfc"]/v["se"])]
    return pd.DataFrame(rows,columns=["gene","stat"]).sort_values("stat",ascending=False).reset_index(drop=True)
def run_gsea(de_full):
    """GSEA prerank over the restricted universe -> {kegg_id:(NES,FDR)}."""
    pre=gp.prerank(rnk=_ranking(de_full),gene_sets=kegg_sets,min_size=MIN_SIZE,max_size=MAX_SIZE,
                   permutation_num=PERM,seed=SEED,outdir=None,no_plot=True,verbose=False)
    res=pre.res2d.copy(); c=list(res.columns)
    term="Term" if "Term" in c else c[0]
    nes=next(x for x in c if x.upper()=="NES"); q=next(x for x in c if "FDR" in x.upper())
    res=res.rename(columns={term:"kegg_id",nes:"NES",q:"FDR_q"})
    res["NES"]=pd.to_numeric(res["NES"],errors="coerce"); res["FDR_q"]=pd.to_numeric(res["FDR_q"],errors="coerce")
    return {r.kegg_id:(r.NES,r.FDR_q) for r in res.itertuples()}

def mean_member_lfc(kegg_id,de):
    """Mean log2FC over a pathway's members present in a DE table {g:{lfc,..}}."""
    v=[de[g]["lfc"] for g in kegg_sets.get(kegg_id,[]) if g in de]
    return float(np.mean(v)) if v else np.nan

# ---- shared volcano facet ---------------------------------------------------
def volcano_facet(ax,lp,title,xlim=(-8,8),ylim=(-2,52),show_ylabel=False,label_per_cat=1,
                  fs_title=10.5,fs_count=7.2,fs_tick=7,fs_ylabel=7.8,fs_label=6.6,
                  cat_s=12,bg_s=5,sig_s=6,extra=None):
    """Draw one volcano panel. lp = {gene:(lfc,padj)}. KEGG categories coloured;
    significant (padj<0.05 & |log2FC|>=1) genes highlighted; top symbol-named gene
    per category labelled in italic, white-outlined, repelled text. Font sizes and
    marker sizes are parameters so the panel works at any physical figure size."""
    gl=list(lp)
    xs=np.array([lp[g][0] for g in gl]); ys=np.array([-np.log10(max(lp[g][1],1e-50)) for g in gl])
    sig=np.array([(lp[g][1]<0.05 and abs(lp[g][0])>=1) for g in gl])
    incat=np.array([g in CATSET for g in gl])
    ax.scatter(xs[~sig],ys[~sig],s=bg_s,color="#ececec",alpha=0.55,ec="none",zorder=1,rasterized=True)
    ax.scatter(xs[sig&~incat],ys[sig&~incat],s=sig_s,color="#bcc2c9",alpha=0.5,ec="none",zorder=2,rasterized=True)
    for _,c,s in CATS:
        m=np.array([g in s for g in gl])&sig
        if m.any(): ax.scatter(xs[m],ys[m],s=cat_s,color=c,alpha=0.82,ec="white",lw=0.2,zorder=4)
    for xv in (1,-1): ax.axvline(xv,ls=(0,(3,3)),color="#d4d4d4",lw=0.7)
    ax.axhline(-np.log10(0.05),ls=(0,(3,3)),color="#d4d4d4",lw=0.7)
    nu=int((sig&(xs>0)).sum()); nd=int((sig&(xs<0)).sum())
    ax.set_title(title,fontsize=fs_title,fontweight="bold",pad=6)
    ax.text(0.5,0.995,f"{nu:,} ↑   {nd:,} ↓",transform=ax.transAxes,ha="center",va="top",fontsize=fs_count,color="#333")
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_xticks([-6,-3,0,3,6])
    ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=fs_tick)
    if show_ylabel: ax.set_ylabel("$-$log$_{10}$ adj P",fontsize=fs_ylabel)
    else: ax.set_yticklabels([])
    sigset={g for g,sg in zip(gl,sig) if sg}
    cand=[g for g in gl if g in sigset and sym(g)]
    labels={}                                   # gene -> text
    for _,c,s in CATS:
        for g in sorted([g for g in cand if g in s],key=lambda g:-abs(lp[g][0]))[:label_per_cat]:
            labels.setdefault(g,sym(g))
    if extra:                                   # force-label specific genes (e.g. the drug target)
        for g,txt in extra.items():
            if g in lp:
                labels[g]=txt
                c=next((c for _,c,s in CATS if g in s),"#555")
                ax.scatter([lp[g][0]],[-np.log10(max(lp[g][1],1e-50))],s=cat_s*2.4,color=c,
                           ec="black",lw=0.8,zorder=5)   # mark it so the point is visible
    tx=[]
    for g,txt in labels.items():
        c=next((c for _,c,s in CATS if g in s),"#555")
        tt=ax.text(lp[g][0],-np.log10(max(lp[g][1],1e-50)),txt,fontsize=fs_label,color=c,
                   fontweight="bold",fontstyle="italic",zorder=6); tt.set_path_effects(STROKE); tx.append(tt)
    if HAVE_AT and tx:
        adjust_text(tx,ax=ax,arrowprops=dict(arrowstyle="-",color="#9a9a9a",lw=0.45),
                    expand=(2.2,2.8),force_text=(1.3,2.0))

def cat_handles(only=None):
    """Legend handles for the categories; pass `only` (set of labels) to show
    just the categories that actually appear in a given figure."""
    return [Line2D([0],[0],marker="o",ls="",mfc=c,mec="white",ms=6,label=l)
            for l,c,_ in CATS if (only is None or l in only)]

def panel_letters(fig,items,fs=16,dy=0.014):
    """items: list of (ax, 'A', dx). Places bold letters at each axis' top-left."""
    fig.canvas.draw()
    for ax,L,dx in items:
        p=ax.get_position(); fig.text(p.x0-dx,p.y1+dy,L,fontsize=fs,fontweight="bold",va="bottom",ha="left")
