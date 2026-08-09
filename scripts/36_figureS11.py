#!/usr/bin/env python3
# =============================================================================
# figureS11_gem_diagnostics.py — Supplementary Figure S11: GEM grounding
# diagnostics (run from project root).
#   A  TemStaPro pseudo-Tm distribution over the 742 model enzymes
#   B  DLTKcat temperature response: per-reaction median log10 kcat(T),
#      with CYP51 (r_0317) highlighted
#   C  Required AOX flux share vs AOX transcript abundance (the rho = 0.99
#      correspondence behind main Figure 6B)
#   D  Bayesian calibration as an identifiability analysis: emergent prior,
#      posterior median and 90% posterior-predictive band vs the measured TPC.
#      Only dTm is curve-constrained (-14.1 K); the posterior overcorrects the
#      falling limb, so downstream analyses use the uncalibrated model.
# =============================================================================
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as st
st.apply()

OUT = "figures"
T3 = st.TEMP3

fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(7.4, 6.2))
fig.subplots_adjust(left=0.09, right=0.97, top=0.94, bottom=0.09,
                    wspace=0.32, hspace=0.36)

# ---- A: Tm distribution ------------------------------------------------------
bp = pd.read_csv("etcGEMs/strains/zt_ipo323/thermal/BestParamsTopt.csv", index_col=0)
tm = bp["Tm"] - 273.15
axA.hist(tm, bins=36, color="#457BBE", ec="white", lw=0.4)
for T, c in T3.items():
    axA.axvline(T, color=c, lw=1.0, ls=(0, (3, 2)))
    axA.text(T, axA.get_ylim()[1]*0.97, f"{T}", color=c, ha="center", va="top", fontsize=6.5)
axA.set_xlabel("TemStaPro pseudo-T$_m$ (°C)")
axA.set_ylabel(f"Enzymes (n = {len(tm)})")

# ---- B: DLTKcat kcat(T) ------------------------------------------------------
dl = pd.read_csv("etcGEMs/strains/zt_ipo323/dltkcat/output.csv")
med = dl.groupby("Temp_C")["pred_log10kcat"].quantile([0.25, 0.5, 0.75]).unstack()
axB.fill_between(med.index, med[0.25], med[0.75], color="#457BBE", alpha=0.20, lw=0)
axB.plot(med.index, med[0.5], color="#27519E", lw=1.4, label="All reactions (median, IQR)")
cy = dl[dl["rxn_id"].str.startswith("r_0317")].groupby("Temp_C")["pred_log10kcat"].median()
axB.plot(cy.index, cy.values, color="#D64B21", lw=1.4, label="CYP51 (r_0317)")
axB.set_xlabel("Temperature (°C)"); axB.set_ylabel("Predicted log$_{10}$ k$_{cat}$ (s$^{-1}$)")
axB.legend(fontsize=6.2, frameon=False, loc="upper left")

# ---- C: AOX flux share vs transcript ----------------------------------------
fx = pd.read_csv("tables/gem/ztGEM_condition_fluxes.csv").set_index("cond")
q = pd.read_csv("tables/gem/mc_uncertainty_quantiles.csv")
q.columns = ["cond", "qq", "AOX_share", "PO", "glut_flux", "CUE"]
aox_med = q[q["qq"] == 0.5].set_index("cond")["AOX_share"] * 100   # same source as main Fig 6B
order = ["15_0", "15_2", "21_0", "21_2", "27_0", "27_2"]
mkr = {0: "o", 2: "^"}
for cond in order:
    t, d = (int(x) for x in cond.split("_"))
    axC.scatter(fx.loc[cond, "AOX_vst"], aox_med[cond],
                s=55, color=T3[t], marker=mkr[d], ec="white", lw=0.7, zorder=3)
    axC.annotate(f"{t} °C, {d}", (fx.loc[cond, "AOX_vst"], aox_med[cond]),
                 textcoords="offset points", xytext=(6, 3), fontsize=6.0, color=T3[t])
rho = aox_med[order].corr(fx.loc[order, "AOX_vst"], method="spearman")
axC.text(0.04, 0.95, f"Spearman ρ = {rho:.2f}", transform=axC.transAxes, va="top",
         fontsize=7, color=st.INK2)
axC.set_xlabel("AOX transcript (VST)")
axC.set_ylabel("AOX share of O$_2$ (%)")

# ---- D: calibration identifiability analysis --------------------------------
pp = np.load("etcGEMs/strains/zt_ipo323/outputs/calibration_phase1/posterior_predictive.npz")
Td, plo, pmed, phi = pp["temps_C"], pp["lo"], pp["med"], pp["hi"]
prior, obsT, obs = pp["prior"], pp["obs_T"], pp["obs"]
msk = Td <= 40
axD.plot(Td[msk], prior[msk], ls=(0, (4, 2)), color="#27519E", lw=1.4,
         label="Emergent (uncalibrated)", zorder=3)
axD.fill_between(Td[msk], plo[msk], phi[msk], color="#D64B21", alpha=0.15, lw=0,
                 zorder=2, label="Calibrated, 90% posterior band")
axD.plot(Td[msk], pmed[msk], color="#D64B21", lw=1.4,
         label="Calibrated (posterior median)", zorder=4)
axD.scatter(obsT, obs, s=30, color="#1A1A1A", zorder=5, label="Measured growth")
axD.set_xlim(5, 40); axD.set_ylim(0, 0.28)
axD.set_xlabel("Temperature (°C)"); axD.set_ylabel("Growth rate (h$^{-1}$)")
axD.legend(fontsize=5.8, frameon=False, loc="upper right")

for ax, L in [(axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")]:
    p = ax.get_position()
    fig.text(p.x0 - 0.055, p.y1 + 0.015, L, fontsize=12, fontweight="bold")

fig.savefig(os.path.join(OUT, "FigureS11.png"), dpi=600, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "FigureS11.pdf"), bbox_inches="tight")
print("S11 done")
