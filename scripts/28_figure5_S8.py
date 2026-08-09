#!/usr/bin/env python3
# =============================================================================
# 11_figure6_interaction.py — Figure 6: temperature x drug COLLISION (explains
# the Bliss interaction in Fig 3). Both panels span all three temperatures.
#   A  Pathway convergence: drug NES (pooled over engaged 15+21 C) vs warming NES
#   B  Drug effect at 15/21/27 C vs warming -> alignment collapses with temperature
# Supplementary FigureS3: gene-level convergence (peak drug effect vs warming).
# Shared style/loaders/GSEA come from fig_common.
# =============================================================================
from fig_common import *
from matplotlib import gridspec

drug={t:load_full(f"proth_at_{t}C_2vs0") for t in ["15","21","27"]}
warm=load_full("temp_27vs15_proth0")
G15=run_gsea(drug["15"]); G21=run_gsea(drug["21"]); Gtemp=run_gsea(warm)
common=[g for g in warm if all(g in drug[t] for t in ["15","21","27"])]
xv=np.array([warm[g]["lfc"] for g in common])

# ===================== MAIN FIGURE 6 (two panels) ============================
fig=plt.figure(figsize=(10.3,8.25))
gs=gridspec.GridSpec(2,2,figure=fig,width_ratios=[1,1.32],height_ratios=[1.35,1],
                     wspace=0.24,hspace=0.34,left=0.07,right=0.96,top=0.94,bottom=0.07)

# ---- A pathway convergence (drug pooled 15+21 vs warming) ----
axA=fig.add_subplot(gs[0])
paths=[k for k in kegg_sets if k in G15 and k in G21 and k in Gtemp]
ndrug={k:0.5*(G15[k][0]+G21[k][0]) for k in paths}; fdrug={k:min(G15[k][1],G21[k][1]) for k in paths}
nd=np.array([ndrug[k] for k in paths]); nt=np.array([Gtemp[k][0] for k in paths])
sigB=np.array([(fdrug[k]<0.05 or Gtemp[k][1]<0.05) for k in paths])
lim=np.nanmax(np.abs(np.concatenate([nd,nt])))*1.12
axA.axhline(0,color="#e6e6e6",lw=0.8); axA.axvline(0,color="#e6e6e6",lw=0.8)
axA.plot([-lim,lim],[-lim,lim],ls="--",color="#bbb",lw=1)
axA.scatter(nd[~sigB],nt[~sigB],s=18,color="#dadada",alpha=0.6,ec="none",zorder=2)
axA.scatter(nd[sigB],nt[sigB],s=38,color="#4D4D4D",alpha=0.85,ec="white",lw=0.4,zorder=3)
rB=np.corrcoef(nd,nt)[0,1]; axA.text(0.97,0.04,f"r = {rB:.2f}",transform=axA.transAxes,va="bottom",ha="right",fontsize=10,color="#333")
# label key pathways, coloured by their KEGG category (consistent with Figs 4/5)
HILITE={"map03010":("Ribosome","#2F6BB3"),"map00190":("OXPHOS","#E08214"),
        "map02010":("ABC transporters (efflux)","#7B3FA0"),"map00982":("Drug metab. P450","#D6604D"),
        "map03050":("Proteasome","#11838E"),"map00020":("TCA cycle","#8C510A")}
tx=[]
for k,(lab,col) in HILITE.items():
    if k in ndrug and k in Gtemp:
        xx,yy=ndrug[k],Gtemp[k][0]
        axA.scatter([xx],[yy],s=46,color=col,ec="white",lw=0.7,zorder=5)
        t=axA.text(xx,yy,lab,fontsize=7.6,color=col,fontweight="bold",zorder=6); t.set_path_effects(STROKE); tx.append(t)
if HAVE_AT and tx: adjust_text(tx,ax=axA,arrowprops=dict(arrowstyle="-",color="#999",lw=0.5),expand=(1.6,2.0))
_leg=[("Translation","#1F78B4"),("Respiration","#E6820D"),("TCA cycle","#8C564B"),
      ("Proteasome","#17BECF"),("P450 detox","#E7298A"),("Efflux","#6A3D9A")]
axA.legend(handles=[Line2D([0],[0],marker="o",ls="",mfc=c,mec="white",ms=6,label=l) for l,c in _leg],
           loc="upper left",fontsize=6.8,frameon=False,title="Pathway category",title_fontsize=7.5,
           handletextpad=0.3,labelspacing=0.3,borderpad=0.2)
axA.set_xlabel("Prothioconazole pathway enrichment\n(NES, pooled 15 + 21 °C)"); axA.set_ylabel("Warming pathway enrichment (NES)")
axA.set_xlim(-lim,lim); axA.set_ylim(-lim,lim); axA.set_aspect("equal"); axA.spines[["top","right"]].set_visible(False)

# ---- B drug at 15/21/27 vs warming ----
axB=fig.add_subplot(gs[1])
axB.axhline(0,color="#e6e6e6",lw=0.8); axB.axvline(0,color="#e6e6e6",lw=0.8)
axB.plot([-9,9],[-9,9],ls="--",color="#c4c4c4",lw=1,zorder=1)
xline=np.linspace(-8,8,100); hand=[]
for t in ["15","21","27"]:
    yd=np.array([drug[t][g]["lfc"] for g in common]); c=TEMP[t]
    axB.scatter(xv,yd,s=5,color=c,alpha=0.09,ec="none",zorder=2,rasterized=True)
    b,a=np.polyfit(xv,yd,1); rr=np.corrcoef(xv,yd)[0,1]
    axB.plot(xline,b*xline+a,color=c,lw=3,zorder=4,solid_capstyle="round")
    hand.append(Line2D([0],[0],color=c,lw=3,label=f"Drug at {t} °C   slope={b:+.2f},  r={rr:+.2f}"))
axB.set_xlabel("Warming effect, 27 vs 15 °C (log$_2$FC, no drug)"); axB.set_ylabel("Prothioconazole effect (log$_2$FC)")
axB.set_xlim(-8,8); axB.set_ylim(-7,7); axB.spines[["top","right"]].set_visible(False)
axB.legend(handles=hand,fontsize=9.5,loc="upper left",frameon=False)
# ---- C/D physiology-transcriptome integration (module scores vs traits) ----
import pandas as _pd
_ROOT=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),".."))
_ms=_pd.read_csv(f"{_ROOT}/tables/revision/integration/module_scores_by_condition_extended.csv").set_index("cond")
_ph=_pd.read_csv(f"{_ROOT}/tables/revision/integration/physiology_condition_values.csv").set_index("cond")
_mt=_pd.read_csv(f"{_ROOT}/tables/revision/integration/module_trait_correlations.csv")
def _integ_panel(ax,module,trait,xlab,ylab,rx=0.04):
    x=_ms[module]; y=_ph[trait]
    r=float(_mt[(_mt.module==module)&(_mt.trait==trait)].pearson_r.iloc[0])
    b=np.polyfit(x,y,1); xr=np.linspace(x.min(),x.max(),50)
    ax.plot(xr,np.polyval(b,xr),color="#9a9a9a",lw=1.2,ls="--",zorder=1)
    for cond in _ms.index:
        t,d=cond.split("_"); mk="o" if d=="0" else "^"
        ax.scatter(x[cond],y[cond],s=90,color=TEMP[t],marker=mk,ec="white",lw=0.8,zorder=3)
        ax.annotate(f"{t} °C, {'0' if d=='0' else '2'}",(x[cond],y[cond]),textcoords="offset points",
                    xytext=(7,4),fontsize=7.2,color=TEMP[t])
    ax.text(rx,0.94,f"r = {r:+.2f}",transform=ax.transAxes,fontsize=10,color="#333",va="top")
    ax.set_xlabel(xlab,fontsize=9.5); ax.set_ylabel(ylab,fontsize=9.5)
    ax.tick_params(labelsize=8); ax.spines[["top","right"]].set_visible(False)
axC=fig.add_subplot(gs[1,0]); _integ_panel(axC,"Ribosome","growth","Ribosome module score (VST z)","Growth rate (C C$^{-1}$ h$^{-1}$)")
axD=fig.add_subplot(gs[1,1]); _integ_panel(axD,"AOX","CUE","AOX module score (VST z)","CUE",rx=0.42)
axD.legend(handles=[Line2D([0],[0],marker="o",ls="",mfc="#777",mec="white",ms=8,label="0 mg L$^{-1}$"),
                    Line2D([0],[0],marker="^",ls="",mfc="#777",mec="white",ms=8,label="2 mg L$^{-1}$")],
           loc="lower left",fontsize=8,frameon=False,title="Prothioconazole",title_fontsize=8.5)
panel_letters(fig,[(axA,"A",0.055),(axB,"B",0.05),(axC,"C",0.055),(axD,"D",0.05)])
fig.savefig(f"{FIG_DIR}/Figure6_interaction.png",dpi=300,bbox_inches="tight")
fig.savefig(f"{FIG_DIR}/Figure6_interaction.pdf",bbox_inches="tight")

# ===================== SUPPLEMENTARY: gene-level convergence =================
def peak(g): return max((drug[t][g]["lfc"] for t in ["15","21","27"]),key=abs)
sigany=lambda g: any(drug[t][g]["padj"]<0.05 and abs(drug[t][g]["lfc"])>=1 for t in ["15","21","27"])
yp=np.array([peak(g) for g in common]); sg=np.array([sigany(g) for g in common])
figS,axS=plt.subplots(figsize=(5.4,5.1)); figS.subplots_adjust(left=0.13,right=0.96,top=0.90,bottom=0.12)
axS.axhline(0,color="#e6e6e6",lw=0.8); axS.axvline(0,color="#e6e6e6",lw=0.8); axS.plot([-8,8],[-8,8],ls="--",color="#bbb",lw=1)
axS.scatter(xv[~sg],yp[~sg],s=5,color="#e8e8e8",alpha=0.45,ec="none",rasterized=True)
xs,ys=xv[sg],yp[sg]; shared=np.sign(xs)==np.sign(ys)
axS.scatter(xs[shared&(xs>0)],ys[shared&(xs>0)],s=11,color=UP,alpha=0.55,ec="none")
axS.scatter(xs[shared&(xs<0)],ys[shared&(xs<0)],s=11,color=DN,alpha=0.55,ec="none")
axS.scatter(xs[~shared],ys[~shared],s=11,color=DIV,alpha=0.6,ec="none")
rS=np.corrcoef(xs,ys)[0,1]; conc=100*np.mean(shared)
axS.text(0.03,0.97,f"r = {rS:.2f}\n{conc:.0f}% same direction",transform=axS.transAxes,va="top",fontsize=9.5,color="#333")
txS=[]
for g,lab,col in [(aoxg[0] if aoxg else None,"AOX2 (shared)","#B2182B"),("Mycgr3G50464","CDR1 (drug-specific)","#6A3D9A")]:
    if g and g in common:
        axS.scatter([warm[g]["lfc"]],[peak(g)],s=48,color=col,ec="white",lw=0.9,zorder=5)
        t=axS.text(warm[g]["lfc"],peak(g),lab,fontsize=8,color=col,fontweight="bold"); t.set_path_effects(STROKE); txS.append(t)
if HAVE_AT and txS: adjust_text(txS,ax=axS,arrowprops=dict(arrowstyle="-",color="#999",lw=0.5),expand=(1.5,1.9))
axS.set_xlabel("Warming effect, 27 vs 15 °C (log$_2$FC)"); axS.set_ylabel("Peak prothioconazole effect across 15–27 °C (log$_2$FC)")
axS.set_xlim(-8,8); axS.set_ylim(-8,8); axS.set_aspect("equal"); axS.spines[["top","right"]].set_visible(False)
axS.set_title("Gene-level convergence of drug and warming",loc="left",fontsize=11,fontweight="bold")
axS.legend(handles=[Line2D([0],[0],marker="o",ls="",mfc=UP,mec="none",ms=7,label="shared ↑"),
                    Line2D([0],[0],marker="o",ls="",mfc=DN,mec="none",ms=7,label="shared ↓"),
                    Line2D([0],[0],marker="o",ls="",mfc=DIV,mec="none",ms=7,label="divergent")],
           fontsize=8,loc="lower right",frameon=False)
figS.savefig(f"{FIG_DIR}/FigureS3_gene_convergence.png",dpi=300,bbox_inches="tight")
figS.savefig(f"{FIG_DIR}/FigureS3_gene_convergence.pdf",bbox_inches="tight")
print(f"Figure 6 done — pathway r={rB:.2f} ({len(paths)} pathways);  S3 gene r={rS:.2f} ({conc:.0f}% concordant)")
