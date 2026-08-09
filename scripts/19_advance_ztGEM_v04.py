#!/usr/bin/env python3
# advance_ztGEM_v04.py — ztGEM v0.4:
#  (1) maintenance fitted to the four sub-optimal control temperatures
#      (ATP_maint = 79.2 * mu; NGAM ~ 0; NNLS on 15/18/21/24C controls)
#  (2) required AOX shares under fitted maintenance + FVA ranges
#  (3) E-Flux transcriptome-constrained models per condition, calibrated at
#      21_0 only, predicting growth / O2 / AOX share / CUE for the other five
import cobra, pandas as pd, numpy as np, logging, json
from cobra.flux_analysis import pfba, flux_variability_analysis
logging.getLogger("cobra").setLevel(logging.ERROR)

FC = 0.48
GAM_EXTRA = 79.2       # mmol ATP /gDW per unit mu (fitted, maintenance_fit_points.csv)
M0 = cobra.io.read_sbml_model("models/gem/ztGEM_v03.xml")
phys = pd.read_csv("tables/revision/integration/physiology_condition_values.csv")
phys["O2"] = phys["respiration"] * FC / 12.011 * 1000
conds = ["15_0","15_2","21_0","21_2","27_0","27_2"]
GROWTH, O2EX, SUC, NGAM = "r_2111", "r_1992", "r_2058", "r_4046"

def prep(M, o2cap=1000, block_ferm=False):
    med = M.medium
    for rid in ["r_1714","r_4502","r_4504"]:
        if rid in med: med[rid] = 0
    med[SUC] = 1000; med[O2EX] = o2cap
    M.medium = med
    if block_ferm:   # respiratory-culture assumption (RQ = 1 in physiology pipeline)
        ferm = ["ethanol","acetate","glycerol","acetaldehyde","(r)-lactate",
                "(s)-lactate","succinate","pyruvate","formate","ethyl acetate"]
        for r in M.exchanges:
            nm = (list(r.metabolites)[0].name or "").lower()
            if any(nm == f or nm.startswith(f + " ") for f in ferm):
                r.upper_bound = 0.0

## ---- (2) required AOX with fitted maintenance ------------------------------
rows = []
for cond in conds:
    c = phys[phys["cond"]==cond].iloc[0]
    M = M0.copy(); prep(M, c["O2"])
    maint = GAM_EXTRA * c["growth"]
    M.reactions.get_by_id(NGAM).bounds = (maint, maint)
    M.reactions.get_by_id(GROWTH).bounds = (c["growth"], c["growth"])
    M.reactions.get_by_id(O2EX).bounds = (-c["O2"], -c["O2"])
    M.objective = M.reactions.get_by_id(SUC); M.objective_direction = "max"
    v = M.slim_optimize()
    if v is None or v != v:
        rows.append({"cond": cond}); continue
    M.reactions.get_by_id(SUC).bounds = (v, v)
    s = pfba(M)
    fva = flux_variability_analysis(M, [M.reactions.r_AOX], fraction_of_optimum=0.0)
    rows.append({"cond": cond, "AOX_flux": s.fluxes["r_AOX"],
                 "AOX_share": s.fluxes["r_AOX"]*0.5/c["O2"],
                 "AOX_share_min": fva.loc["r_AOX","minimum"]*0.5/c["O2"],
                 "AOX_share_max": fva.loc["r_AOX","maximum"]*0.5/c["O2"],
                 "sucrose": -v})
RF = pd.DataFrame(rows)
print("required AOX (fitted maintenance):")
print(RF.round(3).to_string(index=False))
RF.to_csv("tables/gem/ztGEM_v04_required_AOX_fitted_maintenance.csv", index=False)

## ---- (3) E-Flux predictions -------------------------------------------------
vsd = pd.read_csv("tables/revision/expression/vsd_matrix.csv", index_col=0)
cd = pd.read_csv("tables/revision/expression/coldata.csv")
cd["cond"] = cd["temperature"].astype(int).astype(str)+"_"+cd["prothioconazole"].astype(int).astype(str)
expr = 2 ** vsd
cond_expr = {cc: expr[cd.loc[cd["cond"]==cc, "sample"]].mean(axis=1) for cc in conds}

def split_top(s, kw):
    depth = 0; out = []; cur = []
    for t in s.replace("(", " ( ").replace(")", " ) ").split():
        if t == "(": depth += 1; cur.append(t)
        elif t == ")": depth -= 1; cur.append(t)
        elif t == kw and depth == 0: out.append(" ".join(cur)); cur = []
        else: cur.append(t)
    out.append(" ".join(cur)); return out

def gpr_activity(rule, ev):
    if not rule: return None
    s = rule.strip()
    ors = split_top(s, "or")
    if len(ors) > 1:
        vals = [gpr_activity(x, ev) for x in ors]
        vals = [v for v in vals if v is not None]
        return sum(vals) if vals else None
    ands = split_top(s, "and")
    if len(ands) > 1:
        vals = [gpr_activity(x, ev) for x in ands]
        vals = [v for v in vals if v is not None]
        return min(vals) if vals else None
    s = s.strip()
    if s.startswith("(") and s.endswith(")"): return gpr_activity(s[1:-1], ev)
    return ev.get(s)

def eflux_model(cc, cap):
    ev = cond_expr[cc].to_dict()
    M = M0.copy(); prep(M, block_ferm=True)
    acts = {}
    for rxn in M.reactions:
        if rxn.id in (GROWTH, NGAM) or rxn.boundary: continue
        a = gpr_activity(rxn.gene_reaction_rule, ev)
        if a is not None: acts[rxn.id] = a
    scale = np.percentile(list(acts.values()), 95)
    for rid, a in acts.items():
        r = M.reactions.get_by_id(rid)
        ub = cap * min(a/scale, 2.0)
        r.upper_bound = min(r.upper_bound, ub)
        if r.lower_bound < 0: r.lower_bound = max(r.lower_bound, -ub)
    return M

def eflux_predict(cc, cap):
    M = eflux_model(cc, cap)
    # maintenance stoichiometrically coupled to growth: v_NGAM = GAM_EXTRA * mu
    ng = M.reactions.get_by_id(NGAM); ng.bounds = (0, 1000)
    gr = M.reactions.get_by_id(GROWTH)
    cpl = M.problem.Constraint(
        ng.flux_expression - GAM_EXTRA * gr.flux_expression, lb=0, ub=0)
    M.add_cons_vars(cpl)
    M.objective = gr; M.objective_direction = "max"
    mu = M.slim_optimize()
    if mu is None or mu != mu: return None
    s = pfba(M)
    mu = s.fluxes[GROWTH]; o2 = -s.fluxes[O2EX]; co2 = s.fluxes["r_1672"]
    cbio = mu * FC/12.011*1000
    return {"mu": mu, "O2": o2, "AOX_share": s.fluxes["r_AOX"]*0.5/o2 if o2>1e-9 else np.nan,
            "CUE": cbio/(cbio+co2) if co2==co2 and (cbio+co2)>0 else np.nan}

# calibrate CAP so predicted growth at 21_0 matches measured
target = phys.set_index("cond").loc["21_0","growth"]
lo, hi = 0.5, 400.0
for _ in range(25):
    mid = (lo+hi)/2
    p = eflux_predict("21_0", mid)
    if p is None or p["mu"] < target: lo = mid
    else: hi = mid
CAP = (lo+hi)/2
print(f"\nE-Flux calibrated at 21_0: CAP={CAP:.4f} (pred mu={eflux_predict('21_0',CAP)['mu']:.4f}, obs {target:.4f})")

pred = []
for cc in conds:
    p = eflux_predict(cc, CAP)
    c = phys.set_index("cond").loc[cc]
    pred.append({"cond": cc, "mu_pred": p["mu"] if p else np.nan, "mu_obs": c["growth"],
                 "O2_pred": p["O2"] if p else np.nan, "O2_obs": c["O2"],
                 "AOX_share_pred": p["AOX_share"] if p else np.nan,
                 "CUE_pred": p["CUE"] if p else np.nan, "CUE_obs": c["CUE"]})
PR = pd.DataFrame(pred)
print(PR.round(3).to_string(index=False))
from scipy import stats
oth = PR[PR["cond"] != "21_0"].dropna()
r_mu = stats.pearsonr(oth["mu_pred"], oth["mu_obs"])
r_o2 = stats.pearsonr(oth["O2_pred"], oth["O2_obs"])
r_cue = stats.pearsonr(oth["CUE_pred"], oth["CUE_obs"])
rs_mu = stats.spearmanr(oth["mu_pred"], oth["mu_obs"])
print(f"held-out predictions (5 conds): growth r={r_mu.statistic:.3f} (p={r_mu.pvalue:.3f}, rho={rs_mu.statistic:.2f}), "
      f"O2 r={r_o2.statistic:.3f}, CUE r={r_cue.statistic:.3f}")
PR.to_csv("tables/gem/ztGEM_v04_eflux_predictions.csv", index=False)
json.dump({"GAM_extra": GAM_EXTRA, "NGAM": 0.0, "EFlux_CAP": CAP,
           "calibration_condition": "21_0"}, open("tables/gem/ztGEM_v04_maintenance_fit.json","w"))
print("saved v0.4 outputs")
