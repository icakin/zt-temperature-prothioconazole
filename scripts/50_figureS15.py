#!/usr/bin/env python3
# =============================================================================
# 50_figureS15.py — Figure S15: climate relevance of the measured thermal bands.
#   A  Mean spray-window hours per season in each assayed band, 1961-1990 vs
#      1996-2025, pooled over ten NW-European winter-wheat sites.
#   B  Annual supra-optimal (>24 C) spray-window hours with a linear trend.
#   C  Per-site change in supra-optimal hours between the two normals.
# Input : tables/revision/climate/ (script 49)
# RUN FROM the master "X0123 copy" folder:  python3 scripts/50_figureS15.py
# =============================================================================
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0, "scripts")
import fig_style as st
st.apply()

TAB, OUT = "tables/revision/climate", "figures"
os.makedirs(OUT, exist_ok=True)
POOL = pd.read_csv(os.path.join(TAB, "climate_pooled_summary.csv"))
ANN = pd.read_csv(os.path.join(TAB, "climate_annual_pooled.csv"))
PS = pd.read_csv(os.path.join(TAB, "climate_per_site_change.csv"))

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(9.6, 3.3))
fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.26, wspace=0.42)

# ---- A: hours per band, early vs recent ------------------------------------
bands = ["synergy_15_20", "neutral_20_24", "antagonism_24_28", "above_assayed"]
lab = ["15–20 °C\nsynergistic", "20–24 °C\nnear T$_{opt}$",
       "24–28 °C\nantagonistic", ">28 °C\nnot assayed"]
cols = ["#27519E", "#64A5DE", "#D64B21", "#A31621"]
P = POOL.set_index("band")
x = np.arange(len(bands))
axA.bar(x - 0.19, [P.loc[b, "early_mean_hours"] for b in bands], 0.36,
        color="#C9C9C9", ec="white", lw=0.6)
axA.bar(x + 0.19, [P.loc[b, "recent_mean_hours"] for b in bands], 0.36,
        color=cols, ec="white", lw=0.6)
for i, b in enumerate(bands):
    d = P.loc[b, "change_hours"]
    y = max(P.loc[b, "early_mean_hours"], P.loc[b, "recent_mean_hours"])
    axA.text(i, y * 1.03, f"{d:+.0f} h", ha="center", va="bottom",
             fontsize=6.0, color=st.INK2)
axA.set_xticks(x); axA.set_xticklabels(lab, fontsize=5.6)
axA.set_ylabel("Spray-window hours per season")
from matplotlib.patches import Patch
axA.legend(handles=[Patch(fc="#C9C9C9", ec="white", label="1961–1990"),
                    Patch(fc="#7A7A7A", ec="white", label="1996–2025")],
           fontsize=6.0, frameon=False, loc="upper right",
           title="bars coloured by band", title_fontsize=5.6)

# ---- B: annual supra-optimal hours with trend -------------------------------
axB.plot(ANN["year"], ANN["supra_optimal"], color="#BBBBBB", lw=0.8, zorder=1)
axB.scatter(ANN["year"], ANN["supra_optimal"], s=13, color="#D64B21", zorder=3)
sl = stats.linregress(ANN["year"], ANN["supra_optimal"])
xx = np.array([ANN["year"].min(), ANN["year"].max()])
axB.plot(xx, sl.intercept + sl.slope * xx, color="#1A1A1A", lw=1.4, zorder=4)
axB.text(0.03, 0.95, f"{sl.slope*10:+.1f} h per decade\np = {sl.pvalue:.1e}",
         transform=axB.transAxes, va="top", fontsize=6.4, color=st.INK2)
axB.set_xlabel("Year")
axB.set_ylabel("Supra-optimal (>24 °C)\nspray-window hours")

# ---- C: per-site change -----------------------------------------------------
PS = PS.copy()
PS["supra_change"] = PS["antagonism_24_28_change"] + PS["above_assayed_change"]
PS = PS.sort_values("supra_change")
yy = np.arange(len(PS))
axC.barh(yy, PS["supra_change"], color="#D64B21", ec="white", lw=0.6, height=0.68)
axC.axvline(0, color="#888888", lw=0.8)
axC.set_yticks(yy)
axC.set_yticklabels(PS["site"], fontsize=5.8)
axC.set_xlabel("Change in supra-optimal hours\n(1996–2025 minus 1961–1990)")

for ax, L in [(axA, "A"), (axB, "B"), (axC, "C")]:
    pos = ax.get_position()
    fig.text(pos.x0 - 0.055, pos.y1 + 0.035, L, fontsize=12, fontweight="bold")

fig.savefig(os.path.join(OUT, "FigureS15.png"), dpi=600, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "FigureS15.pdf"), bbox_inches="tight")
print("Figure S15 done")
