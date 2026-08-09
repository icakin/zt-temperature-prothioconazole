#!/usr/bin/env python3
# =============================================================================
# 48_figureS14.py — Figure S14: model-free test of the temperature x dose
# interaction (no Bliss framework, no reference temperature).
#   A  Interaction residual F(T,c) - [f(T) + h(c)] versus temperature, one
#      line per prothioconazole concentration.
#   B  The same residual as a temperature x concentration surface, directly
#      comparable with the Bliss deviation heatmap of Fig. 2F.
#   C  Out-of-sample evidence: exact leave-one-out error of the separable
#      model versus the interaction model.
# Input : tables/revision/interaction/ (script 47)
# RUN FROM the master "X0123 copy" folder:  python3 scripts/48_figureS14.py
# =============================================================================
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
sys.path.insert(0, "scripts")
import fig_style as st
st.apply()

TAB, OUT = "tables/revision/interaction", "figures"
os.makedirs(OUT, exist_ok=True)
G = pd.read_csv(os.path.join(TAB, "interaction_residual_surface.csv"))
S = pd.read_csv(os.path.join(TAB, "interaction_modelfree_stats.csv")).set_index("stat")["value"]
doses = sorted(G["conc"].unique())

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(9.6, 3.3))
fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.20, wspace=0.44)

# ---- A: residual vs temperature, one line per dose --------------------------
ramp = st.DOSE_RAMP[1:]                       # drop the control gray
for c, col in zip(doses, ramp):
    sub = G[G["conc"] == c].sort_values("Temp")
    axA.plot(sub["Temp"], sub["delta_med"], color=col, lw=1.3, label=f"{c:g}")
axA.axhline(0, color="#BBBBBB", lw=0.9, ls=(0, (3, 2)))
axA.axvline(24, color="#DDDDDD", lw=0.8)
axA.text(24.2, axA.get_ylim()[1] * 0.93, "T$_{opt}$", fontsize=6, color=st.INK2)
axA.set_xlabel("Temperature (°C)")
axA.set_ylabel("Interaction residual\n$F(T,c) - [f(T)+h(c)]$  (log scale)")
leg = axA.legend(title="mg L$^{-1}$", fontsize=5.6, title_fontsize=5.8,
                 loc="lower right", ncol=2, frameon=True, framealpha=0.82,
                 borderpad=0.35, handlelength=1.4, columnspacing=1.0)
leg.get_frame().set_facecolor("white"); leg.get_frame().set_linewidth(0)

# ---- B: the residual as a surface, comparable with Fig. 2F ------------------
piv = G.pivot_table(index="conc", columns="Temp", values="delta_med"
                    ).sort_index(ascending=True)   # row 0 = lowest dose
v = np.nanmax(np.abs(piv.values))
im = axB.imshow(piv.values, aspect="auto", cmap=st.DIV_CMAP, origin="lower",
                norm=TwoSlopeNorm(vcenter=0, vmin=-v, vmax=v),
                extent=[G["Temp"].min(), G["Temp"].max(), -0.5, len(piv) - 0.5])
axB.set_yticks(range(len(piv)))
axB.set_yticklabels([f"{c:g}" for c in piv.index], fontsize=6.0)
axB.set_ylabel("Prothioconazole (mg L$^{-1}$)")
axB.set_xlabel("Temperature (°C)")
axB.axvline(24, color="#333333", lw=0.8, ls=(0, (3, 2)))
cb = fig.colorbar(im, ax=axB, pad=0.02, fraction=0.046)
cb.set_label("Interaction residual (log scale)", fontsize=6.0)
cb.ax.tick_params(labelsize=5.6)
axB.set_title("blue, growth below separable (synergy)\n"
              "red, growth above separable (antagonism)",
              fontsize=5.8, fontweight="normal", color=st.INK2, pad=3)

# ---- C: exact leave-one-out evidence ---------------------------------------
lab = ["separable\n$f(T)+h(c)$", "interaction\n$+\\,ti(T,c)$"]
axC.bar([0, 1], [S["loo_rmse_sep"], S["loo_rmse_int"]], 0.5,
        color=["#C9C9C9", "#27519E"], ec="white", lw=0.6)
axC.set_xticks([0, 1]); axC.set_xticklabels(lab, fontsize=6.2)
axC.set_ylabel("Exact leave-one-out RMSE (log growth)")
axC.set_ylim(0, S["loo_rmse_sep"] * 1.38)
for i, k in enumerate(["loo_rmse_sep", "loo_rmse_int"]):
    axC.text(i, S[k] * 1.02, f"{S[k]:.4f}", ha="center", va="bottom", fontsize=6.2)
axC.text(0.5, 0.97,
         f"paired $\\Delta$ squared LOO error\n"
         f"{S['loo_paired_mean']:.4f} "
         f"[{S['loo_paired_lo']:.4f}, {S['loo_paired_hi']:.4f}], "
         f"$P$ = {S['loo_paired_p']:.3f}\n"
         f"$\\Delta$AIC = {S['dAIC']:.1f};  "
         f"$ti(T,c)$: $P$ = {S['ti_p']:.1e}",
         transform=axC.transAxes, ha="center", va="top", fontsize=6.0,
         color=st.INK2)

for ax, L in [(axA, "A"), (axB, "B"), (axC, "C")]:
    pos = ax.get_position()
    fig.text(pos.x0 - 0.058, pos.y1 + 0.035, L, fontsize=12, fontweight="bold")

fig.savefig(os.path.join(OUT, "FigureS14.png"), dpi=600, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "FigureS14.pdf"), bbox_inches="tight")
print("Figure S14 done")
