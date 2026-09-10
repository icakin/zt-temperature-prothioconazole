#!/usr/bin/env python3
# =============================================================================
# 40_power_aox_experiment.py — sample-size / power for the AOX inhibitor growth
# arm, using the residual log-growth variance of the EXISTING SensorDish data.
# Vial-level residual is a LOWER BOUND on biological (between-culture) variance.
# =============================================================================
import numpy as np, pandas as pd
from scipy import stats
d=pd.read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv")
d["conc"]=d["Dose"].replace("Control","0").astype(float); d["Tc"]=d["T"].astype(int)
d=d[d.growth_C_per_C_h>0].copy(); d["lg"]=np.log(d.growth_C_per_C_h)
resid=[]
for _,s in d.groupby(["Tc","conc"]):
    if len(s)>=2: resid.extend((s.lg-s.lg.mean()).values)
sd=np.array(resid).std(ddof=1)
print(f"vial-level residual SD(log growth) = {sd:.3f}  (~{100*np.sqrt(np.exp(sd**2)-1):.1f}% CV; lower bound on biological SD)")

def n_for(eff,vfac,power=0.80,alpha=0.05):
    delta=np.log(eff); crit=stats.norm.ppf(1-alpha/2)
    for n in range(2,1000):
        se=np.sqrt(vfac*sd**2/n); ncp=delta/se
        if 1-stats.norm.cdf(crit-ncp)+stats.norm.cdf(-crit-ncp)>=power: return n,se
    return None,None

print("\nMATERIAL-EFFECT detection (difference-in-differences on log growth):")
for eff in [1.10,1.15,1.20]:
    nb,seb=n_for(eff,2.0); nu,seu=n_for(eff,4.0)
    print(f"  {int((eff-1)*100)}% change: paired/blocked n={nb}/arm | unblocked n={nu}/arm")
print("\nEQUIVALENCE <5% (TOST, true effect 0): CI half-width < log(1.05)")
for n in range(2,1000):
    if stats.norm.ppf(0.95)*np.sqrt(2*sd**2/n)<np.log(1.05):
        print(f"  need n={n}/arm (blocked) -> infeasible; P1 reframed as falsification, not equivalence"); break
