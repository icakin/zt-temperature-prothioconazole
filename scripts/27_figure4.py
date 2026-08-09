#!/usr/bin/env python3
# =============================================================================
# 10_figure5_temperature.py — Figure 5: SOLE effect of temperature (no drug).
#   A  Inclusive volcano of the pure thermal response (15→21, 21→27, 15→27 °C)
#   B  Pathway GSEA heatmap (restricted universe; mean member log2FC; FDR<0.001)
#   C  Classical respiration (OXPHOS) vs alternative oxidase (AOX), control only
# Shared style/loaders/GSEA/volcano come from fig_common.
# =============================================================================
from fig_common import *
from matplotlib import gridspec
from matplotlib.colors import TwoSlopeNorm

# warming contrasts (prothioconazole = 0)
C21=load_full("temp_21vs15_proth0"); C2721=load_full("temp_27vs21_proth0"); C27=load_full("temp_27vs15_proth0")
de_by={"21":C21,"2721":C2721,"27":C27}
STEPS=[("15→21 °C",C21),("21→27 °C",C2721),("15→27 °C",C27)]

# GSEA per warming step, then select pathways at FDR<0.001 in >=1 step
G={key:run_gsea(de) for key,de in [("21",C21),("2721",C2721),("27",C27)]}
CUT=0.001
allp={k for key in G for k in G[key]}
selp=sorted([k for k in allp if any(G[key].get(k,(0,1))[1]<CUT for key in G)],
            key=lambda k:np.nan_to_num(mean_member_lfc(k,C27)))

# normalized counts for panel C (control / no-drug samples only)
meta=pd.read_csv("data/rnaseq/sample_metadata.csv"); meta=meta[meta.is_control==False].copy()
meta["temperature"]=meta.temperature.astype(int).astype(str)
meta["prothioconazole"]=meta.prothioconazole.astype(int).astype(str); meta=meta.set_index("sample")
cnt=pd.read_csv("data/rnaseq/raw_count_matrix.csv",index_col=0); cnt.columns=[c.split("(")[1].rstrip(")") if "(" in c else c for c in cnt.columns]; cnt=cnt[meta.index]
logc=np.log(cnt.replace(0,np.nan)); gm=logc.mean(1); ok=gm.notna()
sf=np.exp((logc[ok.values].sub(gm[ok],axis=0)).median(0)); norm=cnt.div(sf,1)
ctrl=meta[meta.prothioconazole=="0"]

fig=plt.figure(figsize=(10.2,7.35))
outer=gridspec.GridSpec(2,1,figure=fig,height_ratios=[1,1.10],hspace=0.40,left=0.07,right=0.95,top=0.92,bottom=0.07)
row1=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=outer[1],width_ratios=[1.15,1],wspace=0.30)

# ===== A volcano of the pure thermal response (full-width, shared helper) =====
gsB=gridspec.GridSpecFromSubplotSpec(1,3,subplot_spec=outer[0],wspace=0.12)
axBs=[fig.add_subplot(gsB[i]) for i in range(3)]; axB0=axBs[0]
for i,(lab,de) in enumerate(STEPS):
    lp={g:(v["lfc"],v["padj"]) for g,v in de.items()}
    ex={"Mycgr3G67342":"ADH"} if i==0 else None   # label the lone P450 point in 15->21
    volcano_facet(axBs[i],lp,lab,xlim=(-9,9),show_ylabel=(i==0),extra=ex)
fig.text((axB0.get_position().x0+0.95)/2,axB0.get_position().y0-0.052,
         "Temperature effect (log$_2$FC, no prothioconazole)",ha="center",fontsize=8.5)
present={lab for _,de in STEPS for lab,c,s in CATS
         if any((g in de and de[g]["padj"]<0.05 and abs(de[g]["lfc"])>=1) for g in s)}
axB0.legend(handles=cat_handles(only=present),fontsize=6.5,loc="upper left",bbox_to_anchor=(0.0,-0.16),
            ncol=4,frameon=False,handletextpad=0.2,columnspacing=0.8)

# ===== B pathway GSEA heatmap (mean member log2FC; same format as Fig 4C) =====
axC=fig.add_subplot(row1[0])
COLS=[("21","15→21 °C"),("2721","21→27 °C"),("27","15→27 °C")]
MC=np.array([[mean_member_lfc(k,de_by[t]) for t,_ in COLS] for k in selp])
vmaxC=float(np.nanmax(np.abs(MC))); normC=TwoSlopeNorm(vcenter=0,vmin=-vmaxC,vmax=vmaxC)
imC=axC.imshow(MC,cmap="RdBu_r",norm=normC,aspect="auto")
axC.set_xticks(range(len(COLS))); axC.set_xticklabels([l for _,l in COLS],fontsize=8)
axC.set_yticks(range(len(selp))); axC.set_yticklabels([kegg_name.get(k,k) for k in selp],fontsize=7)
for i,k in enumerate(selp):
    for j,(t,_) in enumerate(COLS):
        v=MC[i,j]
        if np.isfinite(v):
            axC.text(j,i,f"{v:+.2f}{stars(G[t].get(k,(0,1))[1])}",ha="center",va="center",
                     fontsize=6.2,color="white" if abs(v)>vmaxC*0.62 else "#222")
for s in axC.spines.values(): s.set_visible(False)
axC.tick_params(length=0)
cbC=fig.colorbar(imC,ax=axC,fraction=0.05,pad=0.03); cbC.set_label("cell = mean member log$_2$FC",fontsize=7.5); cbC.ax.tick_params(labelsize=7)
cbC.ax.text(0.5,1.03,"↑ up with warming",transform=cbC.ax.transAxes,ha="center",va="bottom",fontsize=6,color="#B2182B",fontweight="bold")
cbC.ax.text(0.5,-0.03,"↓ down with warming",transform=cbC.ax.transAxes,ha="center",va="top",fontsize=6,color="#2166AC",fontweight="bold")
axC.set_xlabel(f"Sole effect of temperature; GSEA FDR<{CUT}   * <0.05  ** <0.01  *** <0.001",fontsize=6.8)

# ===== C OXPHOS vs AOX mean expression across temperature (control only) =====
gsD=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=row1[1],wspace=0.5)
def tbars(ax,genes,title):
    sc=norm.loc[[g for g in genes if g in norm.index]].mean(0)
    for i,t in enumerate(["15","21","27"]):
        v=sc[ctrl[ctrl.temperature==t].index]
        ax.bar(i,v.mean(),0.64,yerr=ci95(v),color=TEMP[t],ec="white",lw=0.5,
               error_kw=dict(ecolor="#333",elinewidth=0.8,capsize=2.5))
    ax.set_xticks(range(3)); ax.set_xticklabels(["15°C","21°C","27°C"],fontsize=8.5)
    ax.set_title(title,fontsize=8.7,fontweight="bold"); ax.set_ylabel("mean normalized\nexpression",fontsize=8)
    ax.spines[["top","right"]].set_visible(False)
axd1=fig.add_subplot(gsD[0]); tbars(axd1,kg("map00190"),"Classical respiration\n(OXPHOS, %d genes)"%len(kg("map00190")))
axd2=fig.add_subplot(gsD[1]); tbars(axd2,aoxg,"Alternative oxidase\n(AOX)")

panel_letters(fig,[(axB0,"A",0.045),(axC,"B",0.058),(axd1,"C",0.030)])
fig.savefig(f"{FIG_DIR}/Figure5_temperature.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{FIG_DIR}/Figure5_temperature.pdf",bbox_inches="tight")
print(f"Figure 5 done — {len(selp)} pathways; control n/temp:",{t:int((ctrl.temperature==t).sum()) for t in ['15','21','27']})
