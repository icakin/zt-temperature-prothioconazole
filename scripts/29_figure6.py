#!/usr/bin/env python3
# =============================================================================
# figure6_gem.py — Figure 6: the enzyme- and temperature-constrained GEM.
#   A  Uncalibrated model TPC vs measured growth. Enzyme catalytic optima
#      are anchored at the measured optimum (see Methods); shape, magnitude
#      and high-T failure are the emergent features. Calibration -> S11D.
#   B  Required AOX share of O2 per condition (Monte-Carlo 95% bands) —
#      the model's demand matches AOX transcription (rho = 0.99)
#   C  EC50(T): observed Hill vs the CYP51 safety-margin model and the
#      GEM implementation (Vmax(T) + competitive inhibition)
#   D  Model vs measured CUE per condition (condition-constrained pFBA)
# House style from fig_style.py; canvas sized to print at 183 mm.
# =============================================================================
import os, sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as st
st.apply()

GEM = "tables/gem"
CAL = os.environ.get("CAL_DIR", "etcGEMs/strains/zt_ipo323/outputs/calibration_phase1")
OUT = "figures"

T3 = st.TEMP3
mkr = {0: "o", 2: "^"}

fig = plt.figure(figsize=(9.6, 7.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1, 1],
                      wspace=0.30, hspace=0.38, left=0.09, right=0.965, top=0.94, bottom=0.08)

# ===== A: emergent vs calibrated TPC vs measurement ==========================
axA = fig.add_subplot(gs[0, 0])
npz = os.path.join(CAL, "posterior_predictive.npz")
pp = np.load(npz)
Td, lo, med, hi = pp["temps_C"], pp["lo"], pp["med"], pp["hi"]
prior, obsT, obs = pp["prior"], pp["obs_T"], pp["obs"]
m = Td <= 40
axA.plot(Td[m], prior[m], color="#27519E", lw=1.6,
         label="Uncalibrated model", zorder=3)
axA.scatter(obsT, obs, s=34, color="#1A1A1A", zorder=5, label="Measured growth")
axA.axvline(24.0, color="#C8C8C8", lw=0.8, ls=(0, (2, 2)), zorder=1)
axA.text(23.3, 0.40, "T$_{opt}$ 24.0 °C", ha="right", va="top",
         fontsize=6.2, color=st.INK2)
axA.set_yscale("log"); axA.set_ylim(3e-3, 0.45)
axA.set_xlabel("Temperature (°C)"); axA.set_ylabel("Growth rate (h$^{-1}$)")
axA.legend(loc="lower left", fontsize=6.4, frameon=False)
axA.set_xlim(5, 40)
# Calibration posterior deliberately not shown here: only dTm is
# curve-constrained and the posterior overcorrects the falling limb --
# shown honestly in Figure S11D as an identifiability analysis.

# ===== B: required AOX share per condition (MC bands) ========================
axB = fig.add_subplot(gs[0, 1])
q = pd.read_csv(os.path.join(GEM, "mc_uncertainty_quantiles.csv"))
q.columns = ["cond", "qq", "AOX_share", "PO", "glut_flux", "CUE"]
piv = q.pivot(index="cond", columns="qq", values="AOX_share") * 100
order = ["15_0", "15_2", "21_0", "21_2", "27_0", "27_2"]
fx = pd.read_csv(os.path.join(GEM, "ztGEM_condition_fluxes.csv")).set_index("cond")
for i, cond in enumerate(order):
    t, d = (int(x) for x in cond.split("_"))
    c = T3[t]
    lo_, md_, hi_ = piv.loc[cond, 0.025], piv.loc[cond, 0.5], piv.loc[cond, 0.975]
    axB.vlines(i, lo_, hi_, color=c, lw=2.2, alpha=0.9)
    axB.scatter([i], [md_], s=52, color=c, marker=mkr[d], ec="white", lw=0.7, zorder=4)
rho = pd.Series({c: fx.loc[c, "AOX_share_of_O2"] for c in order}).corr(
      pd.Series({c: fx.loc[c, "AOX_vst"] for c in order}), method="spearman")
axB.text(0.04, 0.95, f"vs AOX transcript:\nSpearman ρ = {rho:.2f}",
         transform=axB.transAxes, va="top", fontsize=6.6, color=st.INK2)
axB.set_xticks(range(6))
axB.set_xticklabels(["15 °C\n0", "15 °C\n2", "21 °C\n0", "21 °C\n2", "27 °C\n0", "27 °C\n2"], fontsize=6.2)
axB.set_ylabel("Required AOX share of O$_2$ (%)")
axB.set_xlabel("Condition (temperature, prothioconazole mg L$^{-1}$)")

# ===== C: EC50(T) — observed vs safety margin vs GEM =========================
axC = fig.add_subplot(gs[1, 0])
hp = pd.read_csv("tables/physiology/growth/03_hill_parameters_by_temperature.csv")
mf = pd.read_csv(os.path.join(GEM, "cyp51_margin_fit.csv")).sort_values("T_C")
axC.plot(mf["T_C"], mf["EC50_pred"], color="#D64B21", lw=1.5,
         label=f"CYP51 safety-margin model (E$_a$ = 0.47 eV)", zorder=3)
axC.errorbar(hp["temperature"], hp["EC50"],
             yerr=[hp["EC50"] - hp["EC50_lower"], hp["EC50_upper"] - hp["EC50"]],
             fmt="o", ms=4.5, color="#1A1A1A", elinewidth=0.9, capsize=2,
             label="Observed (Hill EC$_{50}$)", zorder=4)
gem_pts = {15: 1.93, 21: 2.05, 27: 4.91}   # ecFBA implementation (ztGEM_v04_report §GEM-EC50)
axC.scatter(list(gem_pts), list(gem_pts.values()), s=58, marker="D",
            facecolor="none", edgecolor="#27519E", lw=1.4,
            label="GEM implementation (ecFBA)", zorder=5)
axC.set_yscale("log"); axC.set_yticks([1, 2, 4, 8]); axC.set_yticklabels(["1", "2", "4", "8"])
axC.set_xlabel("Temperature (°C)"); axC.set_ylabel("EC$_{50}$ (mg L$^{-1}$)")
axC.legend(loc="upper left", fontsize=6.4, frameon=False)

# ===== D: effective P/O ratio collapse (predicted, not imposed) ==============
axD = fig.add_subplot(gs[1, 1])
pivP = q.pivot(index="cond", columns="qq", values="PO")
for i, cond in enumerate(order):
    t, d = (int(x) for x in cond.split("_"))
    c = T3[t]
    axD.vlines(i, pivP.loc[cond, 0.025], pivP.loc[cond, 0.975], color=c, lw=2.2, alpha=0.9)
    axD.scatter([i], [pivP.loc[cond, 0.5]], s=52, color=c, marker=mkr[d], ec="white", lw=0.7, zorder=4)
axD.axhline(1.0, ls=(0, (3, 3)), color="#C8C8C8", lw=0.8)
axD.set_xticks(range(6))
axD.set_xticklabels(["15 °C\n0", "15 °C\n2", "21 °C\n0", "21 °C\n2", "27 °C\n0", "27 °C\n2"], fontsize=6.2)
axD.set_ylabel("Effective P/O ratio (ATP per O)")
axD.set_xlabel("Condition (temperature, prothioconazole mg L$^{-1}$)")
axD.text(0.97, 0.95, "energetic efficiency collapses\nunder drug + heat",
         transform=axD.transAxes, ha="right", va="top", fontsize=6.4,
         fontstyle="italic", color=st.INK2)
hand = ([Line2D([0], [0], marker="s", ls="", mfc=T3[t], mec="white", ms=6, label=f"{t} °C") for t in (15, 21, 27)]
        + [Line2D([0], [0], marker=mkr[d], ls="", mfc="#777", mec="white", ms=6,
                  label=f"{d} mg L$^{{-1}}$") for d in (0, 2)])
axD.legend(handles=hand, fontsize=6.2, frameon=False, loc="center left", ncol=2)

# panel letters
for ax, L in [(axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")]:
    p = ax.get_position()
    fig.text(p.x0 - 0.055, p.y1 + 0.012, L, fontsize=13, fontweight="bold")

fig.savefig(os.path.join(OUT, "Figure6.png"), dpi=600, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure6.pdf"), bbox_inches="tight")
print("Figure 6 done | calibration source:", CAL)
