#!/usr/bin/env python3
# =============================================================================
# 05_figure4_drug_mechanism.py — Figure 4: prothioconazole mechanism.
#   A  PCA of all RNA-seq samples (6 groups: temperature × dose)
#   B  Inclusive volcano of the drug response at 15/21/27 C (shared categories)
#   C  Data-driven KEGG pathway heatmap (GSEA, FDR<0.001; figure4C_pathways.csv)
#   D  Classical respiration (OXPHOS) vs alternative oxidase (AOX), control vs drug
# =============================================================================
from fig_common import *
from matplotlib import gridspec
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Ellipse

# ---- normalized counts + metadata (PCA and panel D) -------------------------
meta=pd.read_csv("data/rnaseq/sample_metadata.csv"); meta=meta[meta.is_control==False].copy()
meta["temperature"]=meta.temperature.astype(int).astype(str)
meta["prothioconazole"]=meta.prothioconazole.astype(int).astype(str); meta=meta.set_index("sample")
cnt=pd.read_csv("data/rnaseq/raw_count_matrix.csv",index_col=0); cnt.columns=[c.split("(")[1].rstrip(")") if "(" in c else c for c in cnt.columns]; cnt=cnt[meta.index]
logc=np.log(cnt.replace(0,np.nan)); gm=logc.mean(1); ok=gm.notna()
sf=np.exp((logc[ok.values].sub(gm[ok],axis=0)).median(0)); norm=cnt.div(sf,1); logn=np.log2(norm+1)
D={k:load_de(f"proth_at_{k}C_2vs0") for k in ["15","21","27"]}

fig=plt.figure(figsize=(9.6,7.35))
outer=gridspec.GridSpec(2,1,figure=fig,height_ratios=[1,1.05],hspace=0.42,left=0.07,right=0.95,top=0.93,bottom=0.07)
row0=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=outer[0],width_ratios=[1,1.95],wspace=0.26)
row1=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=outer[1],width_ratios=[1.15,1],wspace=0.34)

# ===== A PCA (6 groups: temperature × dose) ==================================
axA=fig.add_subplot(row0[0])
X=logn.T.values; idx=np.argsort(X.var(0))[::-1][:500]; Xt=X[:,idx]-X[:,idx].mean(0)
U,S,_=np.linalg.svd(Xt,full_matrices=False); pc=U[:,:2]*S[:2]; pv=(S**2/(S**2).sum())[:2]*100
if pc[meta.temperature.eq('27').values,0].mean()<pc[meta.temperature.eq('15').values,0].mean(): pc[:,0]*=-1
if pc[meta.prothioconazole.eq('2').values,1].mean()<pc[meta.prothioconazole.eq('0').values,1].mean(): pc[:,1]*=-1
def ell(ax,x,y,c,ls):                                      # outline-only (no fill)
    if len(x)<3: return
    cov=np.cov(x,y); val,vec=np.linalg.eigh(cov); o=val.argsort()[::-1]; val,vec=val[o],vec[:,o]
    ang=np.degrees(np.arctan2(vec[1,0],vec[0,0])); w,h=2*np.sqrt(val*2.30)
    ax.add_patch(Ellipse((x.mean(),y.mean()),w,h,angle=ang,facecolor="none",edgecolor=c,lw=1.2,ls=ls,zorder=1,alpha=0.9))
shp={"0":"o","2":"^"}; lsd={"0":"-","2":(0,(4,2))}
for t in ["15","21","27"]:                                 # 6 groups
    for p in ["0","2"]:
        m=((meta.temperature==t)&(meta.prothioconazole==p)).values
        ell(axA,pc[m,0],pc[m,1],TEMP[t],lsd[p])
        axA.scatter(pc[m,0],pc[m,1],c=TEMP[t],marker=shp[p],s=40,edgecolor="white",lw=0.6,zorder=3)
axA.set_xlabel(f"PC1 ({pv[0]:.0f}%) — temperature",fontsize=8.5); axA.set_ylabel(f"PC2 ({pv[1]:.0f}%) — Prothioconazole",fontsize=8.5)
axA.tick_params(labelsize=7.5); axA.spines[["top","right"]].set_visible(False)
_y0,_y1=axA.get_ylim(); axA.set_ylim(_y0,_y1+0.22*(_y1-_y0))   # headroom for legend
l1=axA.legend(handles=[Line2D([0],[0],marker="s",color="w",markerfacecolor=TEMP[t],markersize=8,label=f"{t} °C") for t in ["15","21","27"]],
              title="Temperature",loc="upper left",frameon=False,fontsize=7.5,title_fontsize=8)
axA.add_artist(l1)
axA.legend(handles=[Line2D([0],[0],marker="o",color="#555",ls="-",label="0 mg L$^{-1}$"),
                    Line2D([0],[0],marker="^",color="#555",ls=(0,(4,2)),label="2 mg L$^{-1}$")],
           title="Prothioconazole",loc="upper right",frameon=False,fontsize=7.5,title_fontsize=8)

# ===== B inclusive volcano (shared helper) ===================================
gsB=gridspec.GridSpecFromSubplotSpec(1,3,subplot_spec=row0[1],wspace=0.13)
axBs=[fig.add_subplot(gsB[i]) for i in range(3)]; axB0=axBs[0]
for i,t in enumerate(["15","21","27"]):
    volcano_facet(axBs[i],D[t],f"{t} °C",xlim=(-8,8),show_ylabel=(i==0),
                  fs_title=11,fs_count=7.5,fs_tick=7.5,fs_ylabel=8.5,fs_label=7,
                  extra=({"Mycgr3G110231":"CYP51/ERG11"} if i==0 else None))   # azole target, 15 C panel
fig.text((axB0.get_position().x0+0.95)/2,axB0.get_position().y0-0.05,
         "Prothioconazole effect (log$_2$FC, 2 vs 0 mg L$^{-1}$)",ha="center",fontsize=9)
present={lab for t in ["15","21","27"] for lab,c,s in CATS
         if any((g in D[t] and D[t][g][1]<0.05 and abs(D[t][g][0])>=1) for g in s)}
axB0.legend(handles=cat_handles(only=present),fontsize=6.5,loc="upper left",bbox_to_anchor=(0.0,-0.155),
            ncol=4,frameon=False,handletextpad=0.2,columnspacing=0.8)

# ===== C data-driven pathway heatmap (FDR<0.001) =============================
axC=fig.add_subplot(row1[0])
_pc=pd.read_csv(f"{DE_DIR}/figure4C_pathways.csv")
M=_pc[["meanLFC_15","meanLFC_21","meanLFC_27"]].values
im=axC.imshow(M,cmap="RdBu_r",norm=TwoSlopeNorm(vcenter=0,vmin=-1.3,vmax=1.3),aspect="auto")
axC.set_xticks(range(3)); axC.set_xticklabels(["15 °C","21 °C","27 °C"],fontsize=8)
_disp={"Metabolism of xenobiotics by cytochrome P450":"Xenobiotic metabolism (P450)",
       "Ribosome biogenesis in eukaryotes":"Ribosome biogenesis"}
axC.set_yticks(range(len(_pc))); axC.set_yticklabels([_disp.get(n,n) for n in _pc["kegg_name"]],fontsize=7)
for i in range(len(_pc)):
    for j,k in enumerate(["15","21","27"]):
        if not np.isnan(M[i,j]):
            axC.text(j,i,f"{M[i,j]:+.2f}{stars(_pc.iloc[i][f'FDR_{k}'])}",ha="center",va="center",
                     fontsize=6.3,color="white" if abs(M[i,j])>0.7 else "#222")
axC.set_xlabel("Prothioconazole effect (2 vs 0 mg L$^{-1}$)\nGSEA pathways (FDR<0.001);  * FDR<0.05  ** <0.01  *** <0.001",fontsize=7)
for s in axC.spines.values(): s.set_visible(False)
axC.tick_params(length=0)
cb=fig.colorbar(im,ax=axC,fraction=0.045,pad=0.06); cb.set_ticks([-1.3,0,1.3]); cb.set_ticklabels(["−1.3","0","+1.3"]); cb.ax.tick_params(labelsize=7)
cb.set_label("mean member log$_2$FC",fontsize=7.5,labelpad=2)
cb.ax.text(0.5,1.04,"↑ induced",transform=cb.ax.transAxes,ha="center",va="bottom",fontsize=6.5,color="#B2182B",fontweight="bold")
cb.ax.text(0.5,-0.04,"↓ repressed",transform=cb.ax.transAxes,ha="center",va="top",fontsize=6.5,color="#2166AC",fontweight="bold")

# ===== D classical (OXPHOS) vs alternative (AOX) respiration =================
gsD=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=row1[1],wspace=0.5)
ox_score=norm.loc[kg("map00190")].mean(axis=0); aox_expr=norm.loc[aoxg[0]]; temps=["15","21","27"]
def grp_bars(ax,series,title,ylab,legend=False):
    x=np.arange(3); w=0.36
    for off,p,col,lab in [(-w/2,"0","#9A9A9A","Control (0 mg L$^{-1}$)"),(w/2,"2","#433D84","Prothioconazole (2 mg L$^{-1}$)")]:
        vals=[series[((meta.temperature==t)&(meta.prothioconazole==p)).values] for t in temps]
        means=[float(v.mean()) for v in vals]; errs=[ci95(v) for v in vals]
        ax.bar(x+off,means,w,yerr=errs,color=col,ec="white",lw=0.5,label=lab,zorder=2,
               error_kw=dict(ecolor="#333",elinewidth=0.8,capsize=2.5))
    ax.set_xticks(x); ax.set_xticklabels([f"{t}°C" for t in temps],fontsize=8); ax.tick_params(labelsize=7.5)
    ax.set_title(title,fontsize=9,fontweight="bold"); ax.set_ylabel(ylab,fontsize=8)
    ax.set_ylim(top=ax.get_ylim()[1]*1.20)                 # headroom so legend clears the bars
    ax.spines[["top","right"]].set_visible(False)
    if legend: ax.legend(frameon=False,fontsize=7,loc="upper left")
axd1=fig.add_subplot(gsD[0]); grp_bars(axd1,ox_score,"Classical respiration\n(OXPHOS, %d genes)"%len(kg("map00190")),"mean normalized\nexpression")
axd2=fig.add_subplot(gsD[1]); grp_bars(axd2,aox_expr,"Alternative oxidase\n(AOX)","normalized\nexpression",legend=True)

panel_letters(fig,[(axA,"A",0.060),(axd1,"D",0.040),(axB0,"B",0.055),(axC,"C",0.065)],fs=15)
fig.savefig(f"{FIG_DIR}/Figure4_drug_mechanism.png",dpi=400,bbox_inches="tight")
fig.savefig(f"{FIG_DIR}/Figure4_drug_mechanism.pdf",bbox_inches="tight")
print("Figure 4 done")
