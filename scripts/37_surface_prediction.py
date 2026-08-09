#!/usr/bin/env python3
# =============================================================================
# surface_prediction.py — Task 1: can the network + CYP51 competitive
# inhibition reproduce the full measured dose x temperature growth surface?
#
# Implementation (matches manuscript Methods):
#  - baseline thermal capacity per T from the measured CONTROL growth only
#    (relative scaling; the dose dimension is never fitted)
#  - CYP51 capacity: pure Arrhenius, Ea and Ki from the margin fit
#    (fitted to the seven EC50s), scale anchored at 21 C reference
#  - growth(T, c) = FBA max growth with CYP51 ub = cap(T)/(1 + c/Ki)
#  - model relative surface g(T,c)/g(T,0) compared with measured
#    r(T,c)/r(T,0) across all 56 conditions
# =============================================================================
import numpy as np, pandas as pd, cobra, json

kB = 8.617e-5
DOSES = [0.0, 0.06, 0.12, 0.25, 0.5, 1.0, 2.0, 4.0]
TEMPS = [15, 18, 21, 24, 26, 27, 28]

# ---- measured surface (per-condition mean growth) ----------------------------
d = pd.read_csv('tables/physiology/derived_N0_R_results_with_carbon.csv')
d['Dose'] = d['Dose'].replace('Control', '0').astype(float)
meas = d.groupby(['T', 'Dose'])['growth_C_per_C_h'].mean().unstack()
# dose 0 column name check
d0col = 0.0 if 0.0 in meas.columns else sorted(meas.columns)[0]
meas_rel = meas.div(meas[d0col], axis=0)

# ---- margin parameters -------------------------------------------------------
P = json.load(open('tables/gem/cyp51_margin_params.json'))
Ki, A, Ea = P['Ki_mgL'], P['A'], P['Ea_eV']
mf = pd.read_csv('tables/gem/cyp51_margin_fit.csv').set_index('T_C')
sigma21 = mf.loc[21.0, 'sigma_fit']

# ---- model -------------------------------------------------------------------
m = cobra.io.read_sbml_model('models/gem/ztGEM_v03.xml')
sol0 = m.optimize()
g_max = sol0.objective_value
cyp = m.reactions.get_by_id('r_0317')
# demand slope: CYP51 flux per unit growth (network constant)
from cobra.flux_analysis import pfba
fx = pfba(m)
v_cyp_max = abs(fx.fluxes['r_0317'])
dslope = v_cyp_max / g_max
print(f"g_max={g_max:.3f}/h  v_cyp={v_cyp_max:.3e}  demand slope={dslope:.3e}")

mu = mf['mu']                      # measured control growth per T
biom = m.reactions.get_by_id(str(m.objective.expression.args[0].args[1]).replace('1.0*','')) if False else None
biomass_rxn = [r for r in m.reactions if m.objective.expression.has(r.forward_variable)][0]

rows = []
for T in TEMPS:
    TK, T21K = T + 273.15, 294.15
    g_ref = g_max * mu[T] / mu[21.0]                # thermal capacity from control only
    cap21 = sigma21 * dslope * (g_max * 1.0)        # capacity at 21 C anchored to fitted margin
    cap = cap21 * np.exp(-Ea / kB * (1/TK - 1/T21K))
    for c in DOSES:
        with m as mm:
            mm.reactions.get_by_id(biomass_rxn.id).upper_bound = g_ref
            mm.reactions.get_by_id('r_0317').upper_bound = cap / (1 + c / Ki)
            g = mm.slim_optimize()
        rows.append({'T': T, 'Dose': c, 'g_model': g})
S = pd.DataFrame(rows)
piv = S.pivot(index='T', columns='Dose', values='g_model')
mod_rel = piv.div(piv[0.0], axis=0)

# ---- compare -----------------------------------------------------------------
common_doses = [c for c in DOSES if c in meas_rel.columns and c > 0]
obs = np.array([[meas_rel.loc[T, c] for c in common_doses] for T in TEMPS]).ravel()
prd = np.array([[mod_rel.loc[T, c] for c in common_doses] for T in TEMPS]).ravel()
r = np.corrcoef(obs, prd)[0, 1]
rmse = np.sqrt(np.mean((obs - prd) ** 2))
print(f"\nSURFACE: n={len(obs)} treated conditions  Pearson r={r:.3f}  RMSE={rmse:.3f}")

# Bliss sign structure (referenced to 24 C, like the data analysis)
def bliss(rel_grid, absolute_grid):
    g24_0 = absolute_grid.loc[24, 0.0]
    out = {}
    for T in TEMPS:
        for c in common_doses:
            fT = absolute_grid.loc[T, 0.0] / g24_0
            fF = absolute_grid.loc[24, c] / g24_0
            fo = absolute_grid.loc[T, c] / g24_0
            out[(T, c)] = fo - fT * fF
    return out
bl_model = bliss(mod_rel, piv)
bl_meas = bliss(meas_rel, meas)
agree = sum(1 for k in bl_model if np.sign(bl_model[k]) == np.sign(bl_meas[k]))
print(f"BLISS sign agreement: {agree}/{len(bl_model)} conditions")
cool = [k for k in bl_model if k[0] <= 18]; hot = [k for k in bl_model if k[0] >= 26]
print("model mean Δ cool:", np.mean([bl_model[k] for k in cool]).round(3),
      " hot:", np.mean([bl_model[k] for k in hot]).round(3))
print("meas  mean Δ cool:", np.mean([bl_meas[k] for k in cool]).round(3),
      " hot:", np.mean([bl_meas[k] for k in hot]).round(3))

S.to_csv('tables/gem/surface_model_predictions.csv', index=False)
meas_rel.to_csv('tables/gem/surface_measured_relative.csv')
mod_rel.to_csv('tables/gem/surface_model_relative.csv')
json.dump({'pearson_r': float(r), 'rmse': float(rmse), 'n': int(len(obs)),
           'bliss_sign_agreement': f"{agree}/{len(bl_model)}",
           'model_bliss_cool': float(np.mean([bl_model[k] for k in cool])),
           'model_bliss_hot': float(np.mean([bl_model[k] for k in hot]))},
          open('tables/gem/surface_summary.json', 'w'))
print("saved surface tables")

# ================= refined comparison: Hill-characterised surface =============
hp = pd.read_csv('tables/physiology/growth/03_hill_parameters_by_temperature.csv').set_index('temperature')
hill_rel = pd.DataFrame({c: 1.0 / (1.0 + (c / hp['EC50']) ** hp['n']) for c in common_doses})
obs_h = np.array([[hill_rel.loc[T, c] for c in common_doses] for T in TEMPS]).ravel()
r_h = np.corrcoef(obs_h, prd)[0, 1]
rmse_h = np.sqrt(np.mean((obs_h - prd) ** 2))
print(f"vs HILL surface: r={r_h:.3f}  RMSE={rmse_h:.3f}")

# Bliss signs vs the paper's posterior-median deviations
bp = pd.read_csv('tables/physiology/growth/08_bliss_deviation_pointwise.csv')
bp = bp.set_index(['temperature', 'conc'])['delta_median']
keys = [k for k in bl_model if k in bp.index]
agree_p = sum(1 for k in keys if np.sign(bl_model[k]) == np.sign(bp.loc[k]))
# restrict to conditions the paper calls confidently (|delta| CI excluding 0 -> use prob columns)
bp2 = pd.read_csv('tables/physiology/growth/08_bliss_deviation_pointwise.csv')
bp2['confident'] = (bp2['prob_synergy'] > 0.9) | (bp2['prob_antagonism'] > 0.9)
conf = bp2[bp2['confident']].set_index(['temperature', 'conc'])['delta_median']
keys_c = [k for k in bl_model if k in conf.index]
agree_c = sum(1 for k in keys_c if np.sign(bl_model[k]) == np.sign(conf.loc[k]))
print(f"BLISS sign vs paper posterior medians: {agree_p}/{len(keys)} all; "
      f"{agree_c}/{len(keys_c)} where the data are confident (P>0.9)")
json.dump({'pearson_r_raw': float(r), 'pearson_r_hill': float(r_h),
           'rmse_hill': float(rmse_h), 'n': int(len(obs)),
           'bliss_sign_all': f"{agree_p}/{len(keys)}",
           'bliss_sign_confident': f"{agree_c}/{len(keys_c)}"},
          open('tables/gem/surface_summary.json', 'w'))
hill_rel.to_csv('tables/gem/surface_hill_relative.csv')
print("refined summary saved")
