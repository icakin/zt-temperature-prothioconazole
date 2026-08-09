#!/usr/bin/env python3
# mc_uncertainty.py — Monte-Carlo propagation of measurement + assumption
# uncertainty through the condition-constrained ztGEM.
# Sampled per draw: mu (r0 CI x drug-effect CI), respiration (loglinear-a CI x
# effect CI), biomass carbon fraction FC ~ N(0.48, 0.03), maintenance
# GAM_extra ~ N(79.2, 12). Outputs quantiles of required AOX share, effective
# P/O, glutathione-pathway flux, model CUE per condition.
import cobra, pandas as pd, numpy as np, logging, json, sys
from cobra.flux_analysis import pfba
logging.getLogger("cobra").setLevel(logging.ERROR)
rng = np.random.default_rng(42)
N = 120

M0 = cobra.io.read_sbml_model("models/gem/ztGEM_v03.xml")
hill = pd.read_csv("tables/physiology/growth/03_hill_parameters_by_temperature.csv").set_index("temperature")
resp = pd.read_csv("tables/physiology/respiration/03_loglinear_params_by_temperature.csv").set_index("temperature")
eff  = pd.read_csv("tables/physiology/fungicide_effect_sizes_by_temperature_summary.csv")
eff2 = eff[eff["Dose"].astype(float)==2.0]
def eff_draw(trait, T):
    r = eff2[(eff2["Trait"]==trait)&(eff2["T"].astype(int)==T)].iloc[0]
    sd = (r["q97.5"]-r["q2.5"])/(2*1.96)
    return rng.normal(r["median_effect"], sd)
def ln_draw(med, lo, hi):
    sd = (np.log(hi)-np.log(lo))/(2*1.96)
    return float(np.exp(rng.normal(np.log(med), sd)))

conds = ["15_0","15_2","21_0","21_2","27_0","27_2"]
GROWTH,O2EX,SUC,NGAM = "r_2111","r_1992","r_2058","r_4046"
# glutathione subsystem members (from groups)
glut_rxns = []
for g in M0.groups:
    if "lutathione" in (g.name or ""):
        glut_rxns += [m.id for m in g.members]
glut_rxns = sorted(set(glut_rxns))
print("glutathione rxns:", len(glut_rxns))

# one persistent model per condition (bounds updated per draw)
res = []
Mw = M0.copy()
med = Mw.medium
for rid in ["r_1714","r_4502","r_4504"]:
    if rid in med: med[rid]=0
med[SUC]=1000; med[O2EX]=1000
Mw.medium = med
for cond in conds:
    T, d = (int(x) for x in cond.split("_"))
    n_ok = 0
    for i in range(N):
        FC = float(np.clip(rng.normal(0.48, 0.03), 0.40, 0.56))
        GAMX = float(np.clip(rng.normal(79.2, 12.0), 40, 120))
        mu = ln_draw(hill.loc[T,"r0"], hill.loc[T,"r0_lower"], hill.loc[T,"r0_upper"])
        R  = ln_draw(resp.loc[T,"a"], resp.loc[T,"a_lower"], resp.loc[T,"a_upper"])
        if d:
            mu *= float(np.exp(eff_draw("Growth", T)))
            R  *= float(np.exp(eff_draw("Respiration", T)))
        O2 = R*FC/12.011*1000
        Mw.reactions.get_by_id(NGAM).bounds = (GAMX*mu, GAMX*mu)
        Mw.reactions.get_by_id(GROWTH).bounds = (mu, mu)
        Mw.reactions.get_by_id(O2EX).bounds = (-O2, -O2)
        Mw.reactions.get_by_id(SUC).bounds = (-1000, 0)
        Mw.objective = Mw.reactions.get_by_id(SUC); Mw.objective_direction = "max"
        v = Mw.slim_optimize()
        if v is None or v != v: continue
        Mw.reactions.get_by_id(SUC).bounds = (v, v)
        try:
            s = pfba(Mw)
        except Exception:
            Mw.reactions.get_by_id(SUC).bounds = (-1000, 0); continue
        Mw.reactions.get_by_id(SUC).bounds = (-1000, 0)
        aox = s.fluxes["r_AOX"]; atp = s.fluxes["r_0226"]
        co2 = s.fluxes["r_1672"]; cbio = mu*FC/12.011*1000
        res.append({"cond": cond, "draw": i,
                    "AOX_share": aox*0.5/O2, "PO": atp/(2*O2),
                    "glut_flux": float(np.abs(s.fluxes[glut_rxns]).sum()),
                    "CUE": cbio/(cbio+co2) if (cbio+co2)>0 else np.nan})
        n_ok += 1
    print(cond, "ok draws:", n_ok, flush=True)
D = pd.DataFrame(res)
D.to_csv("tables/gem/mc_uncertainty_draws.csv", index=False)
Q = D.groupby("cond")[["AOX_share","PO","glut_flux","CUE"]].quantile([0.025,0.5,0.975]).round(4)
Q.to_csv("tables/gem/mc_uncertainty_quantiles.csv")
print(Q.to_string())
