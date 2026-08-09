#!/usr/bin/env python3
# Definitive AOX identifiability test.
# Step 1: enumerate TRUE O2 sinks (exclude transport reactions).
# Step 2: with AOX closed, which sink takes the oxygen, and is it physiological?
# Step 3: close AOX + artefactual small-molecule oxidases -> feasible?
import cobra, pandas as pd, numpy as np, logging
from cobra.flux_analysis import pfba, flux_variability_analysis
logging.getLogger("cobra").setLevel(logging.ERROR)
FC, GAM = 0.48, 79.2
M0 = cobra.io.read_sbml_model("models/gem/ztGEM_v03.xml")
phys = pd.read_csv("tables/revision/integration/physiology_condition_values.csv")
phys["O2"] = phys["respiration"]*FC/12.011*1000
GROWTH,O2EX,SUC,NGAM = "r_2111","r_1992","r_2058","r_4046"
O2IDS = {"s_1275","s_1277","s_1278","s_1279","s_2817"}   # O2 in all compartments

# --- true sinks: consume O2 and are NOT pure transport (O2 on both sides) ---
def is_transport(r):
    ins = {m.id for m,v in r.metabolites.items() if v<0}
    outs = {m.id for m,v in r.metabolites.items() if v>0}
    return bool(ins & O2IDS) and bool(outs & O2IDS)
sinks = [r for r in M0.reactions
         if any(m.id in O2IDS and v<0 for m,v in r.metabolites.items())
         and not is_transport(r) and not r.id.startswith("r_1992")]
print(f"true O2-consuming reactions (excluding transport): {len(sinks)}")
for r in sinks[:40]:
    print(f"   {r.id:10s} {r.name[:70]}")

def setup(M, c, closed=()):
    med = M.medium
    for rid in ["r_1714","r_4502","r_4504"]:
        if rid in med: med[rid]=0
    med[SUC]=1000; med[O2EX]=c["O2"]; M.medium=med
    M.reactions.get_by_id(GROWTH).bounds=(c["growth"],c["growth"])
    M.reactions.get_by_id(O2EX).bounds=(-c["O2"],-c["O2"])
    m_=GAM*c["growth"]; M.reactions.get_by_id(NGAM).bounds=(m_,m_)
    for rid in closed:
        if rid in [x.id for x in M.reactions]: M.reactions.get_by_id(rid).bounds=(0,0)
    return M

def run(c, closed=()):
    M = setup(M0.copy(), c, closed)
    M.objective = M.reactions.get_by_id(SUC); M.objective_direction="max"
    v = M.slim_optimize()
    if v != v: return None
    M.reactions.get_by_id(SUC).bounds=(v,v)
    s = pfba(M)
    used = {r.id: s.fluxes[r.id] for r in sinks if abs(s.fluxes[r.id])>1e-6}
    return {"sucrose": -v, "AOX": s.fluxes["r_AOX"], "cytOX": s.fluxes["r_0438"],
            "sinks_used": used}

# artefactual non-respiratory oxidases the optimiser can abuse
ART = ["r_0956"]   # pyridoxine oxidase; extend if others appear
rows=[]
for _, c in phys.iterrows():
    base = run(c)
    noaox = run(c, ("r_AOX",))
    strict = run(c, ("r_AOX",)+tuple(ART))
    rec = {"cond": c["cond"], "AOX_share_base": base["AOX"]*0.5/c["O2"],
           "noAOX_feasible": noaox is not None,
           "noAOX_top_sink": max(noaox["sinks_used"].items(), key=lambda kv: abs(kv[1]))[0] if noaox else "-",
           "noAOX_top_flux": max(noaox["sinks_used"].values(), key=abs) if noaox else np.nan,
           "noAOX_noArtefact_feasible": strict is not None}
    rows.append(rec)
R = pd.DataFrame(rows)
print("\n" + R.round(3).to_string(index=False))
R.to_csv("tables/gem/aox_identifiability.csv", index=False)
