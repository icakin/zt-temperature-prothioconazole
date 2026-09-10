import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import gridspec
from scipy import stats
plt.rcParams.update({"font.family":"DejaVu Sans","svg.fonttype":"none","pdf.fonttype":42,
 "ps.fonttype":42,"font.size":9,"axes.labelsize":10,"xtick.labelsize":8.5,
 "ytick.labelsize":8.5,"axes.linewidth":0.8,"xtick.major.width":0.8,
 "ytick.major.width":0.8,"xtick.major.size":3.5,"ytick.major.size":3.5})
C15,C27="#2166AC","#B2182B"; INK="#1B2420"; MUT="#7A8288"; LINE="#A6AEB3"

# ---- geometry: complex box widened ~15%; AOX kept narrower to match its label
LX,BWL = 22,36
RX,BWR = 73,52
BY,BH  = 60,15
UQY    = 86                      # lowered from 90: breathing room at panel top

def draw(axS):
    axS.set_xlim(0,100); axS.set_ylim(0,100); axS.axis("off")
    def box(cx,cy,w,h,title,fc,ec,fs,tc=INK):
        axS.add_patch(FancyBboxPatch((cx-w/2,cy-h/2),w,h,boxstyle="round,pad=0,rounding_size=3",
            fc=fc,ec=ec,lw=1.0,zorder=3))
        axS.text(cx,cy,title,ha="center",va="center",fontsize=fs,fontweight="bold",color=tc,zorder=4)
    def arr(x1,y1,x2,y2,lw=0.9,ms=8):
        axS.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=ms,
            color=LINE,lw=lw,shrinkA=0,shrinkB=0,zorder=2))
    box(50,UQY,36,12,"UQ pool","#EAEFEC","#93A399",10.5)
    axS.plot([50,50],[UQY-6,78],color=LINE,lw=0.9,zorder=2)
    axS.plot([LX,RX],[78,78],color=LINE,lw=0.9,solid_joinstyle="round",zorder=2)
    for cx in (LX,RX): arr(cx,78,cx,BY+BH/2+1.2)
    box(LX,BY,BWL,BH,"AOX","#FBEAE7",C27,11.5,tc=C27)
    box(RX,BY,BWR,BH,"complex III → IV","#E9F0F8","#7FA3C8",9.6)
    for cx in (LX,RX):
        arr(cx,BY-BH/2-1.2,cx,46,ms=7.5,lw=0.9)
        axS.text(cx,43,"O$_2$ → H$_2$O",fontsize=8.6,color=MUT,ha="center",va="top")
    axS.text(LX,36,"no proton\npumping",fontsize=8.6,color=MUT,ha="center",va="top",
             style="italic",linespacing=1.30)
    axS.text(RX,36,"proton pumping\n→ ATP",fontsize=8.6,color=MUT,ha="center",va="top",
             style="italic",linespacing=1.30)
    axS.text(4,86,"SHAM",ha="left",va="center",fontsize=11.5,fontweight="bold",color=C27)
    axS.plot([12,12],[81,BY+BH/2+2.4],color=C27,lw=1.6,zorder=5)
    axS.plot([5,19],[BY+BH/2+2.4]*2,color=C27,lw=2.6,zorder=5,solid_capstyle="butt")
    axS.text(LX,25,"hypothesized greater\ncontribution at 27 °C",ha="center",va="top",
             fontsize=9.5,color=C27,fontweight="bold",linespacing=1.20)

S15={"Control":[0.04686363,0.04698807,0.04592120],"0.0125":[0.03787549,0.03791131,0.03758877],
"0.025":[0.03858113,0.04206666,0.04042717],"0.05":[0.04545486,0.04682765,0.04291770],
"0.1":[0.04347594,0.05082732,0.04153762],"0.2":[0.04094986,0.04368220,0.04079319]}
S27={"Control":[0.056911552,0.02707208,0.0375065883],"0.0125":[0.035679923,0.03455479,0.0293858560],
"0.025":[0.048367115,0.02489322,0.0378213980],"0.05":[0.030323167,0.01517697,0.0252240356],
"0.1":[0.012129864,0.02539396,0.0045534091],"0.2":[0.015641014,0.01677232,0.0137551349]}
doses=["0.0125","0.025","0.05","0.1","0.2"]
xd=np.array([float(d) for d in doses]); ldall=np.log10(xd)

fig=plt.figure(figsize=(7.1,5.0),dpi=300)
outer=gridspec.GridSpec(1,2,figure=fig,width_ratios=[0.505,0.495],wspace=0.26)
axS=fig.add_subplot(outer[0])
right=gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=outer[1],height_ratios=[1,0.90],hspace=0.55)
topb=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=right[0],width_ratios=[0.14,1],wspace=0.13)
axC=fig.add_subplot(topb[0]); ax=fig.add_subplot(topb[1],sharey=axC); axE=fig.add_subplot(right[1])
draw(axS)

xs=np.logspace(np.log10(0.0115),np.log10(0.225),120); slopes={}
for S,col,mk,off,key in [(S15,C15,"o",0.965,"15"),(S27,C27,"^",1.035,"27")]:
    cv=np.array(S["Control"])
    axC.scatter([1.0]*3,cv,s=13,color=col,alpha=0.40,lw=0,zorder=2)
    axC.errorbar([1.0],[cv.mean()],yerr=[cv.std(ddof=1)],fmt=mk,ms=5,color=col,mfc=col,
                 mec="white",mew=0.7,ecolor=col,elinewidth=1.15,capsize=2.1,capthick=1.15)
    ax.axhline(cv.mean(),color=col,lw=0.7,ls=(0,(4,4)),alpha=0.28,zorder=1)
    M=np.array([S[d] for d in doses]); sl=[stats.linregress(ldall,M[:,c]).slope for c in range(3)]
    slopes[key]=sl
    for c in range(3):
        ax.scatter(xd*off,M[:,c],s=11,color=col,alpha=0.32,lw=0,zorder=2)
        lr=stats.linregress(ldall,M[:,c])
        ax.plot(xs,lr.intercept+lr.slope*np.log10(xs),color=col,lw=0.7,alpha=0.30,zorder=3)
    ax.errorbar(xd,M.mean(1),yerr=M.std(1,ddof=1),fmt=mk,ms=5,color=col,mfc=col,mec="white",
                mew=0.7,ecolor=col,elinewidth=1.15,capsize=2.1,capthick=1.15,zorder=5)
    mi=stats.linregress(np.repeat(ldall,3),M.flatten(order="C"))
    ax.plot(xs,mi.intercept+np.mean(sl)*np.log10(xs),color=col,lw=1.7,zorder=4)
axC.set_xlim(0.72,1.28); axC.set_xticks([1.0]); axC.set_xticklabels(["ctrl"],fontsize=8.0)
axC.set_ylim(0,0.062); axC.set_yticks([0,0.02,0.04,0.06])
axC.set_ylabel("$r$ (h$^{-1}$)",fontsize=10)
axC.spines[["top","right"]].set_visible(False); axC.tick_params(labelsize=8.5)
ax.set_xscale("log"); ax.set_xticks(xd)
ax.set_xticklabels(["0.0125","0.025","0.05","0.1","0.2"],fontsize=8.0)
ax.minorticks_off(); ax.set_xlim(0.0110,0.235)
ax.spines[["top","right","left"]].set_visible(False); ax.tick_params(labelleft=False,left=False)
ax.set_xlabel("SHAM (mM)",fontsize=10)
ax.text(0.228,0.0535,"15 °C",color=C15,fontsize=9.0,fontweight="bold",ha="right")
ax.text(0.228,0.0205,"27 °C",color=C27,fontsize=9.0,fontweight="bold",ha="right")
for i,(key,col,mk) in enumerate([("15",C15,"o"),("27",C27,"^")]):
    sl=np.array(slopes[key]); m=sl.mean(); ci=stats.t.interval(0.95,2,loc=m,scale=stats.sem(sl))
    axE.scatter([i]*3,sl,s=16,color=col,alpha=0.42,lw=0,zorder=2)
    axE.errorbar([i],[m],yerr=[[m-ci[0]],[ci[1]-m]],fmt=mk,ms=5.6,color=col,mfc=col,mec="white",
                 mew=0.7,ecolor=col,elinewidth=1.2,capsize=2.7,capthick=1.2,zorder=4)
axE.axhline(0,color="#B0B4B8",lw=0.75,ls=(0,(3,3)),zorder=1)
axE.set_xlim(-0.7,1.7); axE.set_xticks([0,1]); axE.set_xticklabels(["15 °C","27 °C"],fontsize=9)
axE.set_ylim(-0.048,0.012)
axE.set_ylabel("Δ slope\n(h$^{-1}$ per tenfold)",fontsize=10)
axE.spines[["top","right"]].set_visible(False); axE.tick_params(labelsize=8.5)
axE.text(-0.45,-0.0405,"Δ slope = $-$0.0238\nWelch $P$ = 0.020",ha="left",fontsize=8.0,
         color="#333",va="bottom",linespacing=1.5)

# ---- panel letters: one rule for all three (3 mm left, 1.5 mm above the axes box)
fig.canvas.draw()
DX = (5.0/25.4)/fig.get_figwidth(); DY = (1.5/25.4)/fig.get_figheight()
for lab,a in [("a",axS),("b",axC),("c",axE)]:
    ap=a.get_position()
    fig.text(ap.x0-DX, ap.y1+DY, lab, fontsize=12, fontweight="bold", ha="left", va="bottom")

plt.savefig("tables/revision/aox/Fig_SHAM_MAIN.png",dpi=600,bbox_inches="tight")
plt.savefig("tables/revision/aox/Fig_SHAM_MAIN.pdf",bbox_inches="tight")
print("saved")
