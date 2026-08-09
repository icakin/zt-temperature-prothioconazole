#!/usr/bin/env python3
# constrain_ztGEM.py — physiology-constrained analysis of ztGEM v0.3.
# For each of the six RNA-seq conditions: fix measured growth rate and O2
# consumption, run pFBA, and record the AOX flux the network REQUIRES, the
# carbon budget (model CUE), and ATP yields. Then the boss's question: does
# equal fractional CYP51 inhibition cost different growth at 15 vs 27 C?
import cobra, pandas as pd, numpy as np, logging
from cobra.flux_analysis import pfba
logging.getLogger("cobra").setLevel(logging.ERROR)

FC = 0.48         # g carbon per g DW (yeast-like biomass; assumption, flagged)
M0 = cobra.io.read_sbml_model("models/gem/ztGEM_v03.xml")
phys = pd.read_csv("tables/revision/integration/physiology_condition_values.csv")

# unit conversion: respiration (h-1, fgC/fgC biomass/h) -> mmol O2 /gDW/h (RQ=1)
phys["O2_mmol_gDW_h"] = phys["respiration"] * FC / 12.011 * 1000
print(phys[["cond","growth","respiration","O2_mmol_gDW_h","CUE"]].round(3).to_string(index=False))

GROWTH = "r_2111"   # growth pseudo-reaction
O2EX   = "r_1992"   # oxygen exchange
SUC    = "r_2058"   # sucrose exchange
NGAM   = "r_4046"   # non-growth ATP maintenance

rows = []
sols = {}
for _, c in phys.iterrows():
    M = M0.copy()
    med = M.medium
    for rid in ["r_1714","r_4502","r_4504"]:
        if rid in med: med[rid] = 0
    med[SUC] = 1000
    med[O2EX] = c["O2_mmol_gDW_h"]
    M.medium = med
    # fix growth and O2 consumption to measured values
    M.reactions.get_by_id(GROWTH).bounds = (c["growth"], c["growth"])
    M.reactions.get_by_id(O2EX).bounds = (-c["O2_mmol_gDW_h"], -c["O2_mmol_gDW_h"])
    # minimize sucrose uptake first (economic assumption), then pFBA at that uptake
    M.objective = M.reactions.get_by_id(SUC)
    M.objective_direction = "max"        # uptake is negative; max = least uptake
    v_suc = M.slim_optimize()
    if v_suc is None or v_suc != v_suc:
        print(c["cond"], "INFEASIBLE"); continue
    M.reactions.get_by_id(SUC).bounds = (v_suc, v_suc)
    s = pfba(M)
    sols[c["cond"]] = s
    aox = s.fluxes["r_AOX"]
    cyt = s.fluxes["r_0438"]              # cytochrome c oxidase (O2-consuming)
    o2  = -s.fluxes[O2EX]
    suc = -s.fluxes[SUC]
    co2 = s.fluxes["r_1672"] if "r_1672" in s.fluxes else np.nan   # CO2 exchange
    # carbon budget: sucrose C in vs biomass C vs CO2 C
    c_in = suc * 12
    c_bio = c["growth"] * FC / 12.011 * 1000    # mmol C into biomass /gDW/h
    cue_model = c_bio / (c_bio + co2) if co2 == co2 else np.nan
    ngam = s.fluxes[NGAM]
    rows.append({"cond": c["cond"], "mu_fixed": c["growth"], "O2_fixed": o2,
                 "sucrose_uptake": suc, "CO2_out": co2,
                 "AOX_flux": aox, "cytOX_flux": cyt,
                 "AOX_share_of_O2": aox * 0.5 / o2,
                 "model_CUE": cue_model, "measured_CUE": c["CUE"], "NGAM": ngam})
R = pd.DataFrame(rows)
R.to_csv("tables/gem/ztGEM_condition_fluxes.csv", index=False)
print(R.round(3).to_string(index=False))

# correlate required AOX flux with measured AOX expression (VST)
vsd = pd.read_csv("tables/revision/expression/vsd_matrix.csv", index_col=0)
cd = pd.read_csv("tables/revision/expression/coldata.csv")
cd["cond"] = cd["temperature"].astype(int).astype(str)+"_"+cd["prothioconazole"].astype(int).astype(str)
aox_expr = {cc: vsd.loc["Mycgr3G72918", cd.loc[cd["cond"]==cc, "sample"]].mean()
            for cc in R["cond"]}
R["AOX_vst"] = R["cond"].map(aox_expr)
from scipy import stats
ok = R.dropna(subset=["AOX_flux","AOX_vst"])
r_ = stats.pearsonr(ok["AOX_flux"], ok["AOX_vst"])
rs = stats.spearmanr(ok["AOX_flux"], ok["AOX_vst"])
print(f"required AOX flux vs measured AOX expression: Pearson r={r_.statistic:.3f} "
      f"(p={r_.pvalue:.3f}), Spearman rho={rs.statistic:.3f} (p={rs.pvalue:.3f})")
R.to_csv("tables/gem/ztGEM_condition_fluxes.csv", index=False)

# ---- CYP51 inhibition scan under condition constraints (15_0 vs 27_0) -------
print("\nCYP51 inhibition scan (O2 fixed at condition value, growth maximized):")
scan_rows = []
for cond in ["15_0", "27_0"]:
    c = phys[phys["cond"] == cond].iloc[0]
    M = M0.copy()
    med = M.medium
    for rid in ["r_1714","r_4502","r_4504"]:
        if rid in med: med[rid] = 0
    med[SUC] = 1000; med[O2EX] = c["O2_mmol_gDW_h"]
    M.medium = med
    M.reactions.get_by_id(O2EX).bounds = (-c["O2_mmol_gDW_h"], 0)
    mu0 = M.slim_optimize()
    v0 = abs(pfba(M).fluxes["r_0317"])
    for f in [1, 0.5, 0.2, 0.1, 0.05, 0.02]:
        with M:
            M.reactions.r_0317.bounds = (0, v0 * f)
            mu = M.slim_optimize()
            scan_rows.append({"cond": cond, "cyp51_fraction": f,
                              "mu": mu, "rel_growth": mu / mu0})
S = pd.DataFrame(scan_rows)
S.to_csv("tables/gem/ztGEM_cyp51_scan.csv", index=False)
print(S.pivot(index="cyp51_fraction", columns="cond", values="rel_growth").round(3).to_string())
