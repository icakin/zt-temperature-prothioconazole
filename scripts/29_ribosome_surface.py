#!/usr/bin/env python3
# 29_ribosome_surface.py — one etcGEM + one calibration for the FULL surface.
#   control mu_met(T): calibrated ETC-GEM (network-mechanistic, fits TPC).
#   drug: growth-law coupling mu = mu_met(T) * (1 - f(T)*g(D))^beta, where f(T) is
#   the MEASURED ribosomal-regulon repression fraction at 2 mg/L (drug-induced,
#   interpolated across temperature) and g(D) is the dose-response shape (Hill,
#   g(2)=1). The TEMPERATURE dependence of potency is NOT fitted - it comes from the
#   measured ribosome data. Only the dose shape (K,n) and growth-law exponent beta
#   are fitted. Compares to all 56 measured conditions.
import os, sys, json, logging
import numpy as np, pandas as pd
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from etc_build import ETCGEM
logging.getLogger("cobra").setLevel(logging.ERROR)
OUT = "tables/revision/v3build"

# --- measured ribosome repression f(T) at 2 mg/L ---
de=pd.read_csv("tables/rnaseq_stage/DE_proth_at_21C_2vs0.csv")
desc=de["description"].fillna("").str.lower(); nm=de["gene_name"].fillna("").str.lower()
kegg=de["KEGG"].fillna("").astype(str); pfam=de["PFAM"].fillna("").str.lower()
ribo=de[(desc.str.contains("ribosomal protein"))|(kegg.str.contains("map03010"))|
        (pfam.str.contains("ribosomal_"))|(nm.str.match(r"rp[ls]"))]["gene_id"].tolist()
cm=pd.read_csv("tables/rnaseq_stage/cond_means.csv",index_col=0)
ribo=[gname for gname in ribo if gname in cm.index]
Rsum={c:cm.loc[ribo,c].sum() for c in cm.columns}
f_meas={T:1-Rsum[f"{T}_2"]/Rsum[f"{T}_0"] for T in [15,21,27]}
# interpolate repression to all assay temps (hold flat beyond 27)
allT=[15,18,21,24,26,27,28]
fT={T:float(np.interp(T,[15,21,27],[f_meas[15],f_meas[21],f_meas[27]])) for T in allT}
print("measured ribosome repression f(T):",{T:round(fT[T],3) for T in allT})

# --- calibrated etcGEM control mu_met(T) ---
cal=json.load(open(os.path.join(OUT,"63_calibration.json")));s0=json.load(open(os.path.join(OUT,"64_thermal_gate0.json")))
K=cal["kappa_star"];NG0=cal["NGAM_0"];AL=cal["alpha"];Ea=s0["Ea_selected"];DTM=s0["dTm"]
g=ETCGEM(); mu_met={}
for T in allT:
    g.set_state(float(T),Ea,DTM,K); g.set_maintenance(NG0+AL*max(0,T-24)**2); mu_met[T]=g.max_growth()

# --- measured growth surface ---
d=pd.read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv")
d["conc"]=d["Dose"].replace("Control","0").astype(float); d["Tc"]=d["T"].astype(int)
d=d[d.growth_C_per_C_h>0]
obs=d.groupby(["Tc","conc"])["growth_C_per_C_h"].mean().reset_index()
mu_obs0={T:obs[(obs["Tc"]==T)&(obs["conc"]==0)]["growth_C_per_C_h"].mean() for T in allT}

# --- fit dose-shape g(D) (Hill, normalised g(2)=1) and beta ---
def gD(D,Kd,n): return (D**n/(Kd**n+D**n))/(2.0**n/(Kd**n+2.0**n))
def predict(T,D,Kd,n,beta):
    if D==0: return mu_met[T]
    return mu_met[T]*max(1-fT[T]*gD(D,Kd,n),1e-6)**beta
def resid(p):
    Kd,n,beta=p; out=[]
    for x in obs.itertuples():
        if x.conc>0: out.append(predict(x.Tc,x.conc,Kd,n,beta)-x.growth_C_per_C_h)
    return out
sol=least_squares(resid,[1.0,1.5,1.0],bounds=([0.05,0.5,0.3],[10,4,3]))
Kd,n,beta=sol.x
print(f"fitted dose-shape: K={Kd:.3f} mg/L, n={n:.2f}, growth-law beta={beta:.2f}")

pred=[];meas=[]
for x in obs.itertuples():
    pred.append(predict(x.Tc,x.conc,Kd,n,beta)); meas.append(x.growth_C_per_C_h)
pred=np.array(pred);meas=np.array(meas)
r=np.corrcoef(pred,meas)[0,1]; rmse=np.sqrt(np.mean((pred-meas)**2))
# interaction (relative) surface
gp=[predict(x.Tc,x.conc,Kd,n,beta)/mu_met[x.Tc] for x in obs.itertuples() if x.conc>0]
go=[x.growth_C_per_C_h/mu_obs0[x.Tc] for x in obs.itertuples() if x.conc>0]
r_int=np.corrcoef(gp,go)[0,1]
print(f"full surface: abs r={r:.3f} RMSE={rmse:.4f} | interaction r={r_int:.3f} (n={len(pred)})")

# figure
fig,axes=plt.subplots(2,4,figsize=(13,6.4),sharex=True,sharey=True)
dg=np.linspace(0,4,60)
for ax,T in zip(axes.ravel(),allT):
    sub=obs[obs["Tc"]==T].sort_values("conc")
    ax.plot(sub["conc"],sub["growth_C_per_C_h"],"o",color="#1b2631",ms=5,label="measured")
    ax.plot(dg,[predict(T,D,Kd,n,beta) for D in dg],"-",color="#1f7a3f",lw=2,label="etcGEM+ribosome")
    ax.set_title(f"{T} °C  (ribo repr {100*fT[T]:.0f}%)",fontsize=9.5); ax.set_ylim(0,None)
axes.ravel()[-1].axis("off"); axes[0,0].legend(frameon=False,fontsize=9)
axes[1,0].set_xlabel("prothioconazole (mg L⁻¹)"); axes[0,0].set_ylabel("growth (h⁻¹)"); axes[1,0].set_ylabel("growth (h⁻¹)")
fig.suptitle(f"etcGEM + ribosome sector: temperature-dependence from MEASURED ribosome repression, not fitted "
             f"[abs r={r:.2f}, interaction r={r_int:.2f}]",fontsize=10.5)
fig.tight_layout()
for ext in ["png","pdf"]: fig.savefig(os.path.join(OUT,f"fig_ribosome_surface.{ext}"),dpi=170)
json.dump(dict(f_meas=f_meas,fT=fT,Kd=Kd,n=n,beta=beta,r_abs=r,rmse=rmse,r_int=r_int,n_ribo=len(ribo)),
          open(os.path.join(OUT,"ribosome_surface_fit.json"),"w"),indent=2,default=float)
print("wrote fig_ribosome_surface.png/.pdf")
