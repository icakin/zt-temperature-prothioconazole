#!/usr/bin/env python3
# FigureS13 — identifiability tests: (A) leave-one-temperature-out for the
# CYP51 safety-margin model; (B) AOX blocking under the condition constraints.
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "scripts")
import fig_style as st
st.apply()

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.3))
fig.subplots_adjust(left=0.10, right=0.97, top=0.90, bottom=0.18, wspace=0.34)

# ---- A: LOTO ----
L = pd.read_csv("tables/gem/loto_margin.csv")
TEMPS = L["T_C"].values
cols = dict(zip([15,18,21,24,26,27,28], st.TEMP_RAMP))
for _, r in L.iterrows():
    c = cols[int(r["T_C"])]
    axA.plot([r["EC50_obs"], r["EC50_obs"]], [r["EC50_lower"], r["EC50_upper"]],
             color="#DDDDDD", lw=0, zorder=1)
    axA.errorbar(r["EC50_obs"], r["EC50_pred_LOTO"],
                 xerr=[[r["EC50_obs"]-r["EC50_lower"]], [r["EC50_upper"]-r["EC50_obs"]]],
                 fmt="o", ms=6, color=c, ecolor=c, elinewidth=0.9, capsize=2,
                 mec="white", mew=0.7, zorder=3)
    axA.annotate(f"{int(r['T_C'])}", (r["EC50_obs"], r["EC50_pred_LOTO"]),
                 textcoords="offset points", xytext=(7,4), fontsize=6, color=c)
lim = [0.6, 12]
axA.plot(lim, lim, ls="--", color="#BBBBBB", lw=0.9, zorder=2)
axA.set_xscale("log"); axA.set_yscale("log"); axA.set_xlim(*lim); axA.set_ylim(*lim)
axA.set_xticks([1,2,4,8]); axA.set_xticklabels(["1","2","4","8"])
axA.set_yticks([1,2,4,8]); axA.set_yticklabels(["1","2","4","8"])
axA.set_xlabel("Observed EC$_{50}$ (mg L$^{-1}$, 95% CI)")
axA.set_ylabel("Held-out prediction (mg L$^{-1}$)")
axA.text(0.03, 0.97, "leave-one-temperature-out\nr = −0.14 (log$_{10}$); 4/7 within CI",
         transform=axA.transAxes, va="top", fontsize=6.3, color=st.INK2)

# ---- B: AOX blocking ----
C = pd.read_csv("tables/gem/aox_pfba_cost.csv")
order = ["15_0","15_2","21_0","21_2","27_0","27_2"]
C = C.set_index("cond").loc[order]
x = np.arange(len(order))
colmap = {15: st.TEMP3[15], 21: st.TEMP3[21], 27: st.TEMP3[27]}
bars = [colmap[int(c.split("_")[0])] for c in order]
axB.bar(x, C["pFBA_cost_increase_%"], 0.6, color=bars, ec="white", lw=0.6)
for i, c in enumerate(order):
    axB.text(i, C["pFBA_cost_increase_%"].iloc[i] + 0.12,
             "feasible", ha="center", fontsize=5.6, color=st.INK2)
axB.set_xticks(x)
axB.set_xticklabels([c.replace("_", " °C\n") + " mg L$^{-1}$" for c in order], fontsize=5.4)
axB.set_ylabel("Extra total flux when AOX is blocked (%)")
axB.set_ylim(0, max(C["pFBA_cost_increase_%"])*1.35)
axB.text(0.03, 0.97, "all six conditions remain feasible\nwithout AOX (O$_2$ via cytochrome chain,\nless substrate-level ATP)",
         transform=axB.transAxes, va="top", fontsize=6.3, color=st.INK2)

for ax, L_ in [(axA,"A"), (axB,"B")]:
    pos = ax.get_position()
    fig.text(pos.x0-0.058, pos.y1+0.03, L_, fontsize=12, fontweight="bold")
fig.savefig("figures/FigureS13.png", dpi=600, bbox_inches="tight")
fig.savefig("figures/FigureS13.pdf", bbox_inches="tight")
print("S13 done")
