#!/usr/bin/env python3
# =============================================================================
# 45_enzyme_cost_aox.py — Is the AOX solution cheaper in PROTEIN, or only in FLUX?
# Run from the project root:  python3 scripts/45_enzyme_cost_aox.py
#
# Uses the enzyme-constrained model (ecZtGEM), imposes the same measured
# physiology as the paper, and asks two questions:
#   (1) per oxygen consumed, how much protein does each terminal oxidase draw?
#   (2) under a minimum-PROTEIN objective (rather than minimum-flux/pFBA),
#       does the network still route oxygen through AOX?
# =============================================================================
import cobra, pandas as pd, numpy as np, logging
logging.getLogger("cobra").setLevel(logging.ERROR)
FC, GAM = 0.48, 79.2
M0   = cobra.io.read_sbml_model("models/gem/ecZtGEM_full.xml")
phys = pd.read_csv("tables/revision/integration/physiology_condition_values.csv")
phys["O2"] = phys["respiration"]*FC/12.011*1000
pool  = M0.metabolites.get_by_id("prot_pool")
DRAWS = [r for r in M0.reactions if r.id.startswith("draw_prot_")]
GROWTH, O2EX, SUC, NGAM = "r_2111", "r_1992", "r_2058", "r_4046"

print("== what the model charges each terminal oxidase (uniform kcat) ==")
def per_o2(rid):
    r = M0.reactions.get_by_id(rid)
    prot = sum(abs(v)*abs([d for d in m.reactions if d.id.startswith("draw_")][0].metabolites[pool])
               for m, v in r.metabolites.items() if m.id.startswith("prot_") and v < 0)
    nsub = sum(1 for m, v in r.metabolites.items() if m.id.startswith("prot_") and v < 0)
    o2   = sum(abs(v) for m, v in r.metabolites.items() if m.id in ("s_1275","s_1278") and v < 0)
    return nsub, prot/o2 if o2 else np.nan
for rid in ["r_AOX_No1", "r_0438_No1", "r_0438_No2"]:
    n, po = per_o2(rid)
    print(f"   {rid:12s} subunits charged = {n}   protein per O2 = {po:.3e}")
print("   ATP synthase (r_0226) enzyme-charged:",
      "YES" if any(m.id.startswith("prot_") and v < 0
                   for m, v in M0.reactions.get_by_id("r_0226").metabolites.items()) else "NO  <-- limitation")

print("\n== AOX share under two objectives, identical physiological constraints ==")
mc = pd.read_csv("tables/gem/mc_uncertainty_quantiles.csv")
mc.columns = ["cond","qq","AOX","PO","g","CUE"]
minflux = mc[mc.qq == 0.5].set_index("cond")["AOX"]*100
rows = []
for _, c in phys.iterrows():
    M = M0.copy(); med = M.medium
    for rid in ["r_1714","r_4502","r_4504"]:
        if rid in med: med[rid] = 0
    med[SUC] = 1000; med[O2EX] = c["O2"]; M.medium = med
    M.reactions.get_by_id(GROWTH).bounds = (c["growth"], c["growth"])
    M.reactions.get_by_id(O2EX).bounds  = (-c["O2"], -c["O2"])
    m_ = GAM*c["growth"]; M.reactions.get_by_id(NGAM).bounds = (m_, m_)
    M.objective = M.problem.Objective(
        sum(M.reactions.get_by_id(d.id).flux_expression*abs(d.metabolites[pool]) for d in DRAWS),
        direction="min")
    s = M.optimize()
    rows.append({"cond": c["cond"],
                 "AOX_share_minFlux_pct":    round(minflux[c["cond"]], 1),
                 "AOX_share_minProtein_pct": round(s.fluxes["r_AOX_No1"]*0.5/c["O2"]*100, 1),
                 "total_protein_minProtein": round(s.objective_value, 5)})
R = pd.DataFrame(rows)
print(R.to_string(index=False))
R.to_csv("tables/gem/aox_objective_dependence.csv", index=False)
print("\nCONCLUSION: the AOX solution is least-FLUX, not least-PROTEIN. Under a")
print("proteome-minimising objective the network prefers the cytochrome route,")
print("so the model cannot identify the terminal route; the transcriptomic data can.")
