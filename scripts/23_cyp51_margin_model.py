#!/usr/bin/env python3
# cyp51_margin_model.py — an enzyme-kinetic "safety margin" layer for CYP51
# that makes fungicide potency temperature-dependent WITHIN the GEM framework.
#
# Model: CYP51 capacity follows Arrhenius kinetics (optionally with heat
# denaturation, per AEM aem.03965-14: MgCYP51 active at 22C, dead at 37C).
# Ergosterol demand scales with measured control growth mu(T).
# Competitive inhibition: effective capacity = Vmax(T) / (1 + D/Ki).
# Growth is unaffected until effective capacity < demand; growth halves when
# effective capacity = half demand  =>  EC50(T) = Ki * (2*sigma(T) - 1),
# where sigma(T) = Vmax(T)/demand(T) is the safety margin.
#
# Fit to the seven measured EC50s; compare against a binding-only (van't Hoff)
# model; then implement Vmax(T) in ztGEM and recover EC50(T) from the network.
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"sans-serif","font.size":9,"axes.linewidth":0.8})

kB = 8.617333e-5   # eV/K
hill = pd.read_csv("tables/physiology/growth/03_hill_parameters_by_temperature.csv")
pec  = pd.read_csv("tables/physiology/growth/04_pec50_by_temperature.csv")
H = hill.merge(pec[["temperature","pEC50_median","pEC50_sd"]], on="temperature")
T_C = H["temperature"].values.astype(float)
T_K = T_C + 273.15
mu  = H["r0"].values
ec50_obs = H["EC50"].values
# work on log10(EC50); sd from pEC50 sd (same scale)
y = np.log10(ec50_obs); sy = H["pEC50_sd"].values
Tref = 294.15  # 21 C

def sigma(T_K, A, Ea, Tm=None, dk=2.0):
    v = A*np.exp(-Ea/kB*(1/T_K - 1/Tref))
    if Tm is not None:
        v = v / (1 + np.exp((T_K - Tm)/dk))   # 2-state heat denaturation
    return v/mu

## ---- model 1: margin model (Ki, A, Ea) --------------------------------------
def resid_margin(p):
    Ki, A, Ea = np.exp(p[0]), np.exp(p[1]), p[2]
    s = sigma(T_K, A, Ea)
    ec = Ki*np.maximum(2*s - 1, 1e-6)
    return (np.log10(ec) - y)/sy
fit1 = least_squares(resid_margin, x0=[np.log(1.0), np.log(0.08), 0.5],
                     bounds=([-8,-8,0.0],[8,8,2.5]))
Ki1, A1, Ea1 = np.exp(fit1.x[0]), np.exp(fit1.x[1]), fit1.x[2]
rss1 = np.sum(fit1.fun**2)

## ---- model 1b: margin + denaturation (Ki, A, Ea, Tm) ------------------------
def resid_margin_dn(p):
    Ki, A, Ea, Tm = np.exp(p[0]), np.exp(p[1]), p[2], p[3]
    s = sigma(T_K, A, Ea, Tm=Tm)
    ec = Ki*np.maximum(2*s - 1, 1e-6)
    return (np.log10(ec) - y)/sy
fit1b = least_squares(resid_margin_dn, x0=[np.log(1.0), np.log(0.08), 0.5, 304.0],
                      bounds=([-8,-8,0.0,297.0],[8,8,2.5,320.0]))
Ki1b, A1b, Ea1b, Tm1b = np.exp(fit1b.x[0]), np.exp(fit1b.x[1]), fit1b.x[2], fit1b.x[3]
rss1b = np.sum(fit1b.fun**2)

## ---- model 2: binding-only van't Hoff (log EC50 linear in 1/T) --------------
X = np.vstack([np.ones_like(T_K), 1/T_K]).T
W = np.diag(1/sy**2)
beta = np.linalg.solve(X.T@W@X, X.T@W@y)
y2 = X@beta
rss2 = np.sum(((y2 - y)/sy)**2)

## ---- model 3: constant EC50 (null) ------------------------------------------
ybar = np.sum(y/sy**2)/np.sum(1/sy**2)
rss3 = np.sum(((ybar - y)/sy)**2)

def aicc(rss, k, n=7):
    aic = n*np.log(rss/n) + 2*k
    return aic + 2*k*(k+1)/(n-k-1)
print("model comparison (weighted RSS; AICc):")
print(f"  margin (Ki,A,Ea):        RSS={rss1:7.2f}  AICc={aicc(rss1,3):7.2f}  "
      f"Ki={Ki1:.2f} mg/L  Ea={Ea1:.2f} eV")
print(f"  margin+denat (4 par):    RSS={rss1b:7.2f}  AICc={aicc(rss1b,4):7.2f}  "
      f"Ki={Ki1b:.2f}  Ea={Ea1b:.2f} eV  Tm={Tm1b-273.15:.1f} C")
print(f"  van't Hoff binding only: RSS={rss2:7.2f}  AICc={aicc(rss2,2):7.2f}")
print(f"  constant EC50 (null):    RSS={rss3:7.2f}  AICc={aicc(rss3,1):7.2f}")

s1 = sigma(T_K, A1, Ea1)
print("\nfitted safety margin sigma(T) [margin model]:")
for t, s_, e_o, e_p in zip(T_C, s1, ec50_obs, Ki1*(2*s1-1)):
    print(f"  {t:4.0f} C  sigma={s_:5.2f}  EC50 obs {e_o:5.2f}  pred {e_p:5.2f}")

## ---- figure -----------------------------------------------------------------
Tg = np.linspace(14, 29, 200); TgK = Tg + 273.15
from scipy.interpolate import interp1d
mu_i = interp1d(T_C, mu, kind="cubic", fill_value="extrapolate")
def sigma_g(A, Ea, Tm=None, dk=2.0):
    v = A*np.exp(-Ea/kB*(1/TgK - 1/Tref))
    if Tm is not None: v = v/(1+np.exp((TgK-Tm)/dk))
    return v/mu_i(Tg)
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9))
plt.subplots_adjust(left=0.06, right=0.98, top=0.85, bottom=0.16, wspace=0.32)
ax = axes[0]
ax.errorbar(T_C, ec50_obs, yerr=[ec50_obs-10**(y-sy), 10**(y+sy)-ec50_obs],
            fmt="o", color="black", capsize=3, label="measured")
ax.plot(Tg, Ki1*np.maximum(2*sigma_g(A1,Ea1)-1,1e-6), color="#B2182B", label="margin model")
ax.plot(Tg, 10**(beta[0]+beta[1]/TgK), color="#2166AC", ls="--", label="binding-only")
ax.set_xlabel("temperature (°C)"); ax.set_ylabel("EC$_{50}$ (mg L$^{-1}$)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("A  Safety-margin model fits EC$_{50}$(T)", loc="left", fontsize=9.5)
ax = axes[1]
ax.plot(Tg, A1*np.exp(-Ea1/kB*(1/TgK-1/Tref)), color="#1B7837", label="CYP51 capacity (Arrhenius)")
ax.plot(Tg, mu_i(Tg), color="#6A3D9A", label="ergosterol demand (∝ growth)")
ax.set_xlabel("temperature (°C)"); ax.set_ylabel("rate (h$^{-1}$ equivalents)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("B  Capacity vs demand", loc="left", fontsize=9.5)
ax = axes[2]
ax.plot(Tg, sigma_g(A1,Ea1), color="black")
ax.scatter(T_C, s1, color="black", s=30, zorder=3)
ax.axhline(1, color="grey", ls="--", lw=0.8)
ax.set_xlabel("temperature (°C)"); ax.set_ylabel("safety margin $\\sigma$ = capacity/demand")
ax.set_title("C  Margin is smallest when cool → synergy;\nlargest when hot → antagonism",
             loc="left", fontsize=9.5)
fig.savefig("figures/diagnostics/FigureG3_cyp51_margin.png", dpi=200)
out = pd.DataFrame({"T_C": T_C, "mu": mu, "EC50_obs": ec50_obs,
                    "sigma_fit": s1, "EC50_pred": Ki1*(2*s1-1)})
out.to_csv("tables/gem/cyp51_margin_fit.csv", index=False)
import json
json.dump({"Ki_mgL": float(Ki1), "A": float(A1), "Ea_eV": float(Ea1),
           "RSS_margin": float(rss1), "RSS_vantHoff": float(rss2),
           "RSS_null": float(rss3)}, open("tables/gem/cyp51_margin_params.json","w"))
print("\nsaved figure + tables")
