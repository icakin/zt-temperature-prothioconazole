#!/usr/bin/env python3
# medium_sensitivity.py — Task 3: does the AOX requirement survive when the
# model medium is enriched toward complex YMS (amino acids + nucleobases open)?
import cobra, pandas as pd, numpy as np, logging, re, json
from cobra.flux_analysis import pfba, flux_variability_analysis
logging.getLogger("cobra").setLevel(logging.ERROR)

FC = 0.48
M0 = cobra.io.read_sbml_model("models/gem/ztGEM_v03.xml")
phys = pd.read_csv("tables/revision/integration/physiology_condition_values.csv")
phys["O2_mmol_gDW_h"] = phys["respiration"] * FC / 12.011 * 1000
GROWTH, O2EX, SUC, NGAM = "r_2111", "r_1992", "r_2058", "r_4046"
GAM_EXTRA = 79.2

# find amino-acid + nucleobase exchange reactions by metabolite name
AA = ["alanine","arginine","asparagine","aspartate","cysteine","glutamate","glutamine",
      "glycine","histidine","isoleucine","leucine","lysine","methionine","phenylalanine",
      "proline","serine","threonine","tryptophan","tyrosine","valine",
      "adenine","guanine","uracil","cytosine"]
rich_ex = []
for r in M0.exchanges:
    met = list(r.metabolites)[0]
    nm = (met.name or "").lower()
    if any(a in nm for a in AA) and "trna" not in nm and "peptide" not in nm:
        rich_ex.append(r.id)
print(f"rich-medium exchanges opened: {len(rich_ex)}")

def run(medium_mode):
    rows = []
    for _, c in phys.iterrows():
        M = M0.copy()
        med = M.medium
        for rid in ["r_1714","r_4502","r_4504"]:
            if rid in med: med[rid] = 0
        med[SUC] = 1000
        med[O2EX] = c["O2_mmol_gDW_h"]
        if medium_mode == "rich":
            for rid in rich_ex: med[rid] = 0.5     # mmol/gDW/h each, modest supply
        M.medium = med
        M.reactions.get_by_id(GROWTH).bounds = (c["growth"], c["growth"])
        maint = GAM_EXTRA * c["growth"]
        M.reactions.get_by_id(NGAM).bounds = (maint, maint)
        M.reactions.get_by_id(O2EX).bounds = (-c["O2_mmol_gDW_h"], -c["O2_mmol_gDW_h"])
        M.objective = M.reactions.get_by_id(SUC); M.objective_direction = "max"
        v = M.slim_optimize()
        if v != v: rows.append({"cond": c["cond"], "AOX_share": np.nan}); continue
        M.reactions.get_by_id(SUC).bounds = (v, v)
        s = pfba(M)
        aox, o2 = s.fluxes["r_AOX"], -s.fluxes[O2EX]
        # ATP per O for effective P/O: use ATP synthase flux / (2*O2)
        atp_syn = s.fluxes["r_0226"] if "r_0226" in s.fluxes else np.nan
        # min feasible AOX at optimum (the formal requirement)
        fva = flux_variability_analysis(M, ["r_AOX"], fraction_of_optimum=1.0)
        rows.append({"cond": c["cond"], "AOX_share": aox*0.5/o2,
                     "AOX_min_share": fva.loc["r_AOX","minimum"]*0.5/o2,
                     "PO": atp_syn/(2*o2) if atp_syn==atp_syn else np.nan})
    return pd.DataFrame(rows).set_index("cond")

base = run("minimal"); rich = run("rich")
out = base.join(rich, lsuffix="_minimal", rsuffix="_rich").round(3)
print(out.to_string())
out.to_csv("tables/gem/medium_sensitivity.csv")
json.dump({"n_rich_exchanges": len(rich_ex)}, open("tables/gem/medium_sensitivity_meta.json","w"))
print("saved medium_sensitivity.csv")
