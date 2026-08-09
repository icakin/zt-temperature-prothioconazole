#!/usr/bin/env python3
# FigureS12 — model robustness: dose x temperature surface, anchoring, medium
import sys, os
import numpy as np, pandas as pd, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "scripts")
import fig_style as st
st.apply()

fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.4))
fig.subplots_adjust(left=0.09, right=0.96, top=0.93, bottom=0.09, wspace=0.34, hspace=0.42)
(axA, axB), (axC, axD) = axes
TEMPS = [15,18,21,24,26,27,28]; DOSES=[0.06,0.12,0.25,0.5,1.0,2.0,4.0]

# ---- A: measured (Hill) vs model relative-growth surfaces (curves per T) ----
hill = pd.read_csv("tables/gem/surface_hill_relative.csv", index_col=0); hill.columns=[float(c) for c in hill.columns]
modr = pd.read_csv("tables/gem/surface_model_relative.csv", index_col=0); modr.columns=[float(c) for c in modr.columns]
for T, col in zip(TEMPS, st.TEMP_RAMP):
    axA.plot(DOSES, [hill.loc[T,c] for c in DOSES], color=col, lw=1.3)
    axA.plot(DOSES, [modr.loc[T,c] for c in DOSES], color=col, lw=1.1, ls=(0,(3,2)))
axA.set_xscale("log"); axA.set_xlabel("Prothioconazole (mg L$^{-1}$)")
axA.set_ylabel("Relative growth  $g(T,c)/g(T,0)$")
axA.plot([],[], color=st.INK2, lw=1.3, label="measured (Hill)")
axA.plot([],[], color=st.INK2, lw=1.1, ls=(0,(3,2)), label="model")
axA.legend(fontsize=6, frameon=False, loc="lower left")

# ---- B: scatter model vs measured ----
obs = np.array([[hill.loc[T,c] for c in DOSES] for T in TEMPS])
prd = np.array([[modr.loc[T,c] for c in DOSES] for T in TEMPS])
for i, (T, col) in enumerate(zip(TEMPS, st.TEMP_RAMP)):
    axB.scatter(obs[i], prd[i], s=22, color=col, ec="white", lw=0.4, label=f"{T} °C")
axB.plot([0,1.05],[0,1.05], ls="--", color="#bbb", lw=0.8)
r = np.corrcoef(obs.ravel(), prd.ravel())[0,1]
axB.text(0.04, 0.95, f"r = {r:.2f}\n(49 conditions)", transform=axB.transAxes, va="top", fontsize=6.5, color=st.INK2)
axB.set_xlabel("Measured relative growth"); axB.set_ylabel("Model relative growth")
axB.legend(fontsize=5.4, frameon=False, loc="lower right", ncol=2)

# ---- C: anchoring sensitivity ----
shifts = {"-4": ("20 °C anchor", "#64A5DE"), "-2": ("22 °C anchor", "#457BBE"),
          "0": ("24 °C anchor (used)", "#1A1A1A"), "+2": ("26 °C anchor", "#D64B21")}
for sh, (lab, col) in shifts.items():
    f = f"tables/gem/anchor_tpc_shift{sh}.csv"
    d = pd.read_csv(f)
    tc, gc = d.columns[0], d.columns[1]
    axC.plot(d[tc], d[gc], color=col, lw=1.5 if sh=="0" else 1.1, label=lab,
             zorder=4 if sh=="0" else 2)
axC.axvline(36.8, color="#C8C8C8", lw=0.8, ls=(0,(2,2)))
axC.text(36.5, 0.19, "CTmax invariant", rotation=90, fontsize=5.8, color=st.INK2, ha="right")
axC.set_xlim(5, 40); axC.set_xlabel("Temperature (°C)"); axC.set_ylabel("Model growth rate (h$^{-1}$)")
axC.legend(fontsize=5.8, frameon=False, loc="upper left")

# ---- D: medium sensitivity ----
ms = pd.read_csv("tables/gem/medium_sensitivity.csv").set_index("cond")
conds = ["15_0","15_2","21_0","21_2","27_0","27_2"]
x = np.arange(len(conds)); w = 0.36
axD.bar(x-w/2, [ms.loc[c,"AOX_share_minimal"]*100 for c in conds], w,
        color="#457BBE", label="sucrose minimal (used)")
axD.bar(x+w/2, [ms.loc[c,"AOX_share_rich"]*100 for c in conds], w,
        color="#C9A227", label="enriched (+38 uptakes)")
axD.set_xticks(x); axD.set_xticklabels([c.replace("_"," °C\n")+" mg L$^{-1}$" for c in conds], fontsize=5.2)
axD.set_ylabel("AOX share of O$_2$ (%)")
axD.legend(fontsize=6, frameon=False, loc="upper left")

for ax, L in [(axA,"A"),(axB,"B"),(axC,"C"),(axD,"D")]:
    pos = ax.get_position()
    fig.text(pos.x0-0.055, pos.y1+0.012, L, fontsize=12, fontweight="bold")
fig.savefig("figures/FigureS12.png", dpi=600, bbox_inches="tight")
fig.savefig("figures/FigureS12.pdf", bbox_inches="tight")
print("S12 done")
