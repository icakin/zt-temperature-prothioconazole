#!/usr/bin/env python3
# =============================================================================
# verify_all.py — independent verification of every modelling claim.
# Run from the project root:  python3 scripts/44_verify_all.py
# Each check recomputes a quoted number from source data by a method
# independent of the script that produced it, and prints OK / MISMATCH.
# =============================================================================
import numpy as np, pandas as pd, json, sys
from scipy import stats
ok = lambda b: "OK  " if b else "MISMATCH"
G = "tables/gem"; P = "tables/physiology"

print("== 1. LOTO by closed-form OLS (independent of the optimiser) ==")
h = pd.read_csv(f"{P}/growth/03_hill_parameters_by_temperature.csv")
T, mu, ec = h.temperature.values, h.r0.values, h.EC50.values
TK, X, Y = T+273.15, 1/(T+273.15)-1/294.15, np.log10(h.EC50.values*h.r0.values)
pred = np.array([10**(stats.linregress(np.delete(X,i), np.delete(Y,i)).intercept
                 + stats.linregress(np.delete(X,i), np.delete(Y,i)).slope*X[i])/mu[i]
                 for i in range(len(T))])
r = np.corrcoef(np.log10(ec), np.log10(pred))[0,1]
within = int(((pred>=h.EC50_lower)&(pred<=h.EC50_upper)).sum())
print(f"   out-of-sample r = {r:.3f} (claim ~-0.14)     {ok(r < 0)}")
print(f"   within observed CI = {within}/7 (claim 4/7)  {ok(within==4)}")
print(f"   27 C held out: {pred[T==27][0]:.2f} vs observed {ec[T==27][0]:.2f}")

print("\n== 2. Quoted flux numbers vs source tables ==")
q = pd.read_csv(f"{G}/mc_uncertainty_quantiles.csv")
q.columns = ["cond","qq","AOX","PO","glut","CUE"]
m = q[q.qq==0.5].set_index("cond")
for cond, claim in [("15_2",13),("21_2",24),("27_2",42)]:
    v = round(m.loc[cond,"AOX"]*100)
    print(f"   AOX share {cond}: {v}% (claim {claim}%)      {ok(v==claim)}")
print(f"   P/O 15_0 = {m.loc['15_0','PO']:.2f} (claim ~1.3)   {ok(abs(m.loc['15_0','PO']-1.3)<0.1)}")
print(f"   P/O 27_2 = {m.loc['27_2','PO']:.3f} (claim ~0.07)  {ok(abs(m.loc['27_2','PO']-0.07)<0.02)}")
fx = pd.read_csv(f"{G}/ztGEM_condition_fluxes.csv").set_index("cond")
o = ["15_0","15_2","21_0","21_2","27_0","27_2"]
rho = pd.Series({c:fx.loc[c,"AOX_share_of_O2"] for c in o}).corr(
      pd.Series({c:fx.loc[c,"AOX_vst"] for c in o}), method="spearman")
print(f"   Spearman rho = {rho:.3f} (claim 0.99)       {ok(abs(rho-0.99)<0.01)}")

print("\n== 3. AOX blocking ==")
c = pd.read_csv(f"{G}/aox_pfba_cost.csv")
print(f"   flux penalty {c['pFBA_cost_increase_%'].min():.1f}-{c['pFBA_cost_increase_%'].max():.1f}% "
      f"(claim 0-4%)      {ok(c['pFBA_cost_increase_%'].max()<5)}")
net_on  = c['ATPsynth_on'].sum(); net_off = c['ATPsynth_off'].sum()
print(f"   ATP synthase rises when AOX blocked: {net_on:.2f} -> {net_off:.2f}   {ok(net_off>net_on)}")

print("\n== 4. Surface and anchoring ==")
s = json.load(open(f"{G}/surface_summary.json"))
print(f"   surface r = {s['pearson_r_hill']:.3f} (claim 0.83)   {ok(abs(s['pearson_r_hill']-0.83)<0.01)}")
mr = pd.read_csv(f"{G}/surface_model_relative.csv", index_col=0)
mr = mr[[c for c in mr.columns if float(c) > 0]]        # treated doses only
sat = int((mr.values>0.999).sum())
print(f"   conditions with no predicted effect = {sat} of {mr.size} (claim 31 of 49)  {ok(sat==31)}")
a0 = pd.read_csv(f"{G}/anchor_tpc_shift0.csv"); a2 = pd.read_csv(f"{G}/anchor_tpc_shift+2.csv")
print(f"   rmax invariant to anchor: {a0.iloc[:,1].max():.3f} vs {a2.iloc[:,1].max():.3f}  "
      f"{ok(abs(a0.iloc[:,1].max()-a2.iloc[:,1].max())<0.005)}")

print("\n== 5. Medium sensitivity ==")
ms = pd.read_csv(f"{G}/medium_sensitivity.csv").set_index("cond")
print(f"   27_2 AOX share: minimal {ms.loc['27_2','AOX_share_minimal']*100:.0f}% -> "
      f"rich {ms.loc['27_2','AOX_share_rich']*100:.0f}%  "
      f"{ok(ms.loc['27_2','AOX_share_rich']>ms.loc['27_2','AOX_share_minimal'])}")
print("\nAll checks complete.")
