#!/usr/bin/env python3
# =============================================================================
# 56_shift_scale_diagnostics.py — model-free invariants and the diagnostic
# figure for the shift-versus-scale test. Pre-registered in ShiftScale_spec.md.
#
#   scaling invariant   : log s(T,c) - log s_ref(c) has ~zero slope in log c
#                         (non-positive s omitted from THIS diagnostic only)
#   translation invariant: between-temperature horizontal gap at fixed
#                         suppression level is constant across levels
# Figure: (A) scale-aligned collapse, (B) shift-aligned collapse,
#         (C) bootstrap held-out RMSE, scale vs shift, with thresholds.
#
# RUN FROM the master "X0123 copy" folder:  python3 scripts/56_shift_scale_diagnostics.py
# Outputs -> tables/revision/geometry/  and figures/ (diagnostic_shift_scale)
# =============================================================================
import os, json, sys
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "scripts")
import fig_style as st
st.apply()

OUT = "tables/revision/geometry"
os.makedirs(OUT, exist_ok=True)
TEMPS = [15, 18, 21, 24, 26, 27, 28]
POS = np.array([0.06, 0.12, 0.25, 0.5, 1.0, 2.0, 4.0])
DOSES = np.concatenate([[0.0], POS]); nT, nD = len(TEMPS), len(POS)

raw = pd.read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv")
raw["conc"] = raw["Dose"].replace("Control", "0").astype(float)
raw["Temp"] = raw["T"].astype(float)
raw = raw[(raw["growth_C_per_C_h"] > 0) & np.isfinite(raw["growth_C_per_C_h"])]
cm = raw.groupby(["Temp", "conc"])["growth_C_per_C_h"].mean()
G = np.full((nT, nD), np.nan)
for ti, T in enumerate(TEMPS):
    mu0 = cm[(T, 0.0)]
    for h, c in enumerate(POS):
        if (T, c) in cm.index:
            G[ti, h] = cm[(T, c)] / mu0
S = 1.0 - G                                            # fractional suppression


def Hgrid(h, k): return POS ** h / (POS ** h + k ** h)


def fit_full(model):
    obs = S.copy(); m = np.isfinite(obs)
    def pred(th):
        if model == "scale":
            h, k0 = np.exp(th[0]), np.exp(th[1]); a = 1/(1+np.exp(-th[2:9]))
            P = a[:, None] * Hgrid(h, k0)[None, :]
        else:
            h, M = np.exp(th[0]), 1/(1+np.exp(-th[1])); k = np.exp(th[2:9])
            P = M * np.stack([Hgrid(h, k[i]) for i in range(nT)])
        return P
    def resid(th): P = pred(th); return P[m] - obs[m]
    f = least_squares(resid, np.r_[0.0, 0.0, np.zeros(7)], method="lm", max_nfev=8000)
    return f.x, pred(f.x)


th_sc, P_sc = fit_full("scale")
th_sh, P_sh = fit_full("shift")
alpha = 1/(1+np.exp(-th_sc[2:9]))                       # per-T amplitude
h_sh, M_sh, k_sh = np.exp(th_sh[0]), 1/(1+np.exp(-th_sh[1])), np.exp(th_sh[2:9])

# ---- non-positive suppression, listed and omitted from log diagnostic -------
neg = [(TEMPS[ti], POS[h], float(S[ti, h])) for ti in range(nT) for h in range(nD)
       if np.isfinite(S[ti, h]) and S[ti, h] <= 0]
pd.DataFrame(neg, columns=["T", "conc", "s_obs"]).to_csv(
    os.path.join(OUT, "nonpositive_suppression_omitted.csv"), index=False)

# ---- scaling invariant ------------------------------------------------------
logS = np.log(np.where(S > 0, S, np.nan))
ref = np.nanmean(logS, axis=0)                          # pooled common shape (log)
slopes = []
for ti, T in enumerate(TEMPS):
    y = logS[ti] - ref; x = np.log10(POS)
    ok = np.isfinite(y)
    if ok.sum() >= 3:
        b = np.polyfit(x[ok], y[ok], 1)[0]
        slopes.append({"T": T, "scaling_slope": b, "n_points": int(ok.sum())})
SC = pd.DataFrame(slopes)
SC.to_csv(os.path.join(OUT, "scaling_invariant_slopes.csv"), index=False)

# ---- translation invariant --------------------------------------------------
levels = [0.2, 0.3, 0.4, 0.5]
contour = {}
for ti, T in enumerate(TEMPS):
    ok = np.isfinite(S[ti])
    s_t, x_t = S[ti, ok], np.log10(POS[ok])
    order = np.argsort(s_t)
    if s_t[order][0] <= min(levels) and s_t[order][-1] >= max(levels) and ok.sum() >= 3:
        f = interp1d(s_t[order], x_t[order], bounds_error=False)
        contour[T] = {p: float(f(p)) for p in levels}
rows = []
for p in levels:
    xs = {T: contour[T][p] for T in contour if np.isfinite(contour[T][p])}
    if len(xs) >= 2:
        vals = np.array(list(xs.values())); vals -= vals.mean()
        rows.append({"level": p, "n_temps": len(xs), "logdose_spread_sd": float(vals.std(ddof=1))})
TR = pd.DataFrame(rows)
# translation holds if the per-temperature offsets are constant across levels:
common = {}
for T in contour:
    offs = [contour[T][p] for p in levels if np.isfinite(contour[T].get(p, np.nan))]
    if len(offs) >= 2:
        common[T] = np.std(np.array(offs) - np.mean([contour[t][p] for p in levels for t in contour
                                                     if np.isfinite(contour[t].get(p, np.nan))]))
TR.to_csv(os.path.join(OUT, "translation_invariant_contours.csv"), index=False)

print("=== scaling invariant (slope of log-suppression ratio on log dose; ~0 supports scaling) ===")
print(SC.to_string(index=False))
print(f"  median |slope| = {SC['scaling_slope'].abs().median():.3f}")
print(f"\n  non-positive suppression cells omitted from log diagnostic: {len(neg)}")
for t, c, s in neg: print(f"    {t} C, {c} mg/L, s = {s:.3f}")
print("\n=== translation invariant (between-T log-dose spread at each suppression level) ===")
print(TR.to_string(index=False))
print("  translation holds if spread is ~constant across levels (offsets level-independent)")

# =============================== figure ======================================
dec = json.load(open(os.path.join(OUT, "shift_scale_decision.json")))
bs = np.load(os.path.join(OUT, "shift_scale_bootstrap.npz"))
cols = {T: c for T, c in zip(TEMPS, st.TEMP_RAMP)}
x = np.log10(POS)

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(9.6, 3.3))
fig.subplots_adjust(left=0.075, right=0.985, top=0.88, bottom=0.20, wspace=0.40)

# A: scale-aligned collapse — s / alpha(T) should superimpose
for ti, T in enumerate(TEMPS):
    ok = np.isfinite(S[ti])
    axA.plot(x[ok], S[ti, ok] / alpha[ti], "o-", color=cols[T], ms=3, lw=1, label=f"{T}")
axA.plot(x, Hgrid(np.exp(th_sc[0]), np.exp(th_sc[1])), "k--", lw=1.4, label="common shape")
axA.set_xlabel("log$_{10}$ dose (mg L$^{-1}$)")
axA.set_ylabel("suppression / amplitude  $s/\\alpha(T)$")
axA.set_title("A  scale alignment", loc="left", fontsize=8, fontweight="bold")
axA.legend(fontsize=5.2, frameon=False, ncol=2, loc="upper left")

# B: shift-aligned collapse — s vs (log dose - log k(T)) should superimpose
for ti, T in enumerate(TEMPS):
    ok = np.isfinite(S[ti])
    axB.plot(x[ok] - np.log10(k_sh[ti]), S[ti, ok], "o-", color=cols[T], ms=3, lw=1)
xx = np.linspace(-2.2, 1.2, 100)
axB.plot(xx, M_sh / (1 + 10 ** (-h_sh * xx)), "k--", lw=1.4)
axB.set_xlabel("aligned log$_{10}$ dose  (log c $-$ log $k(T)$)")
axB.set_ylabel("suppression $s$")
axB.set_title("B  shift alignment", loc="left", fontsize=8, fontweight="bold")

# C: bootstrap held-out RMSE, scale vs shift
axC.hist(bs["rmse_scale"], bins=40, color="#D64B21", alpha=0.55, label="scale")
axC.hist(bs["rmse_shift"], bins=40, color="#27519E", alpha=0.55, label="shift")
axC.axvline(dec["RMSE_scale"], color="#D64B21", lw=1.4)
axC.axvline(dec["RMSE_shift"], color="#27519E", lw=1.4)
axC.set_xlabel("held-out RMSE (relative growth)")
axC.set_ylabel("bootstrap count")
axC.set_title("C  out-of-sample error", loc="left", fontsize=8, fontweight="bold")
axC.legend(fontsize=6, frameon=False, loc="upper right")
axC.text(0.02, 0.97, f"P(shift<scale)={dec['p_shift']:.2f}\n"
                     f"|ΔRMSE|={dec['median_abs_rmse_diff']:.3f} vs 0.25ν={dec['materiality_floor']:.3f}\n"
                     f"outcome: {dec['outcome'].split('(')[0].strip()}",
         transform=axC.transAxes, va="top", fontsize=5.8, color=st.INK2)

fig.savefig("figures/diagnostic_shift_scale.png", dpi=300, bbox_inches="tight")
fig.savefig("figures/diagnostic_shift_scale.pdf", bbox_inches="tight")
print("\nfigure -> figures/diagnostic_shift_scale.png")
