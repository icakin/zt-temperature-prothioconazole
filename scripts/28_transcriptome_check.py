#!/usr/bin/env python3
# 28_transcriptome_check.py — decisive first check for the transcriptome-constrained
# etcGEM. Each enzyme's abundance is capped by its MEASURED expression in that
# condition (draw_prot_i upper bound = kappa * counts_i,condition); the shared pool
# is opened so the per-enzyme transcript caps are what bind. One global scale kappa
# is calibrated on the controls. Then DRUG growth is whatever the network achieves
# on the drug-condition proteome - not fitted.
#   Q1: does a drug-induced growth drop emerge at all?
#   Q2: is the drop bigger at 15 C than 27 C? (right sign for temperature-dependent
#       potency; observed EC50 rises with temperature.)
# Thermal kcat layer is OFF here (DLTKcat-21) so the drug/control RATIO at each
# temperature is driven purely by the measured expression differences.
import os, sys, json, logging
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from etc_build import ETCGEM
logging.getLogger("cobra").setLevel(logging.ERROR)
OUT = "tables/revision/v3build"

cm = pd.read_csv("tables/rnaseq_stage/cond_means.csv", index_col=0)  # genes x conditions
g = ETCGEM()
NG0 = json.load(open(os.path.join(OUT,"63_calibration.json")))["NGAM_0"]

# map each enzyme to its draw reaction and to expression
draws = {}
for rxn in g.M.reactions:
    if rxn.id.startswith("draw_prot_"):
        draws[rxn.id[len("draw_prot_"):]] = rxn
g.M.reactions.prot_pool_exchange.bounds = (0.0, 1e6)   # open pool (non-binding)
# reference expression = control 21C, to keep kappa interpretable
ref = cm["21_0"].replace(0, np.nan)

def set_condition(cond, kappa):
    col = cm[cond]
    for enz, rxn in draws.items():
        c = col.get(enz, np.nan)
        if c != c:  # enzyme not measured -> leave generous (won't bind)
            rxn.upper_bound = 1e3
        else:
            rxn.upper_bound = float(kappa * c)
    g.set_maintenance(NG0)

def growth(cond, kappa):
    set_condition(cond, kappa); return g.max_growth()

# observed control + drug growth (dose 0 and 2) at 15/21/27
d = pd.read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv")
d["conc"]=d["Dose"].replace("Control","0").astype(float); d["Tc"]=d["T"].astype(int)
d=d[d.growth_C_per_C_h>0]
obs={}
for T in [15,21,27]:
    for D in [0,2]:
        obs[(T,D)]=d[(d["Tc"]==T)&(d["conc"]==D)]["growth_C_per_C_h"].mean()

# calibrate kappa so control-21 growth == observed control-21
target = obs[(21,0)]
lo,hi=1e-9,1e-2
for _ in range(60):
    mid=np.sqrt(lo*hi); gm=growth("21_0",mid)
    if gm>=target: hi=mid
    else: lo=mid
    if (hi-lo)/hi<1e-3: break
KAP=hi
print(f"calibrated kappa={KAP:.3e}; control-21 growth={growth('21_0',KAP):.4f} (target {target:.4f})")

rows=[]
print("\ncond   model_ctrl model_drug  drop%   | obs_ctrl obs_drug obs_drop%")
for T in [15,21,27]:
    gc=growth(f"{T}_0",KAP); gd=growth(f"{T}_2",KAP)
    drop=1-gd/gc if gc>0 else np.nan
    oc,od=obs[(T,0)],obs[(T,2)]; odrop=1-od/oc
    rows.append(dict(T=T,model_ctrl=gc,model_drug=gd,model_drop=drop,obs_drop=odrop))
    print(f"{T}C    {gc:.4f}    {gd:.4f}   {100*drop:5.1f}%  |  {oc:.4f}  {od:.4f}  {100*odrop:5.1f}%")

R=pd.DataFrame(rows)
print("\n--- verdicts ---")
print("Q1 drug drop emerges (all temps drop>0):", bool((R.model_drop>0).all()))
print("Q2 right sign (drop 15 > drop 27):", bool(R.model_drop.iloc[0] > R.model_drop.iloc[2]),
      f"[model 15={100*R.model_drop.iloc[0]:.1f}% vs 27={100*R.model_drop.iloc[2]:.1f}%]")
print("   observed drop 15={:.1f}% vs 27={:.1f}%".format(100*R.obs_drop.iloc[0],100*R.obs_drop.iloc[2]))
R.to_csv(os.path.join(OUT,"transcriptome_check.csv"),index=False)
json.dump(dict(kappa=KAP,rows=rows),open(os.path.join(OUT,"transcriptome_check.json"),"w"),indent=2,default=float)
print("\nwrote transcriptome_check.csv")
