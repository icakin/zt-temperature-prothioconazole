#!/usr/bin/env python3
# Leave-one-temperature-out prediction for the CYP51 safety-margin model.
import numpy as np, pandas as pd, json
from scipy.optimize import least_squares
kB = 8.617e-5; Tref = 294.15

hill = pd.read_csv('tables/physiology/growth/03_hill_parameters_by_temperature.csv')
pec = pd.read_csv('tables/physiology/growth/04_pec50_by_temperature.csv')
H = hill.merge(pec[['temperature','pEC50_sd']], on='temperature')
T_C = H['temperature'].values; T_K = T_C + 273.15
mu = H['r0'].values
y = np.log10(H['EC50'].values); sy = H['pEC50_sd'].values

def fit_predict(train_idx, test_idx):
    Ttr, ytr, sytr, mutr = T_K[train_idx], y[train_idx], sy[train_idx], mu[train_idx]
    def resid(p):
        Ki, A, Ea = np.exp(p[0]), np.exp(p[1]), p[2]
        s = A*np.exp(-Ea/kB*(1/Ttr - 1/Tref))/mutr
        pred = Ki*np.maximum(2*s - 1, 1e-9)
        return (np.log10(pred) - ytr)/sytr
    fit = least_squares(resid, x0=[np.log(1.0), np.log(0.08), 0.5],
                        bounds=([-30,-30,0],[10,10,2]))
    Ki, A, Ea = np.exp(fit.x[0]), np.exp(fit.x[1]), fit.x[2]
    sT = A*np.exp(-Ea/kB*(1/T_K[test_idx] - 1/Tref))/mu[test_idx]
    return Ki*np.maximum(2*sT - 1, 1e-9), Ea

rows = []
for i in range(len(T_C)):
    tr = [j for j in range(len(T_C)) if j != i]
    pred, Ea = fit_predict(tr, [i])
    rows.append({'T_C': T_C[i], 'EC50_obs': 10**y[i], 'EC50_pred_LOTO': pred[0], 'Ea_train': Ea})
L = pd.DataFrame(rows)
L['log10_err'] = np.log10(L['EC50_pred_LOTO']) - np.log10(L['EC50_obs'])
print(L.round(3).to_string(index=False))
r = np.corrcoef(np.log10(L['EC50_obs']), np.log10(L['EC50_pred_LOTO']))[0,1]
print(f"\nLOTO log-EC50 correlation r = {r:.3f}")
h = L[L.T_C >= 26]
print(f"held-out 27C: obs {L[L.T_C==27]['EC50_obs'].iloc[0]:.2f}  pred {L[L.T_C==27]['EC50_pred_LOTO'].iloc[0]:.2f}")
print(f"fold-rise 15->27 obs {L[L.T_C==27]['EC50_obs'].iloc[0]/L[L.T_C==15]['EC50_obs'].iloc[0]:.2f}  "
      f"pred(27 held out) {L[L.T_C==27]['EC50_pred_LOTO'].iloc[0]/L[L.T_C==15]['EC50_obs'].iloc[0]:.2f}")
L.to_csv('tables/gem/loto_margin.csv', index=False)
json.dump({'loto_r_log10': float(r)}, open('tables/gem/loto_summary.json','w'))
