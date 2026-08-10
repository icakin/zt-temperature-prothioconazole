#!/usr/bin/env python3
# =============================================================================
# 55_shift_scale_cv.py — shift-versus-scale geometry test of the growth
# dose-response surface. Pre-registered in ShiftScale_spec.md revision 2.
#
# Two matched 9-parameter families for fractional suppression s = 1 - g:
#   scale (amplitude): s = alpha(T) * H(c; h, k0)      shared h,k0; per-T alpha
#   shift (position) : s = M      * H(c; h, k(T))      shared h,M ; per-T k
#   mixed (ceiling)  : s = alpha(T)* H(c; h, k(T))     diagnostic only, 15 par
# H(c;h,k)=c^h/(c^h+k^h), H(0)=0 so controls anchor g(T,0)=1 with no pseudodose.
#
# Test: leave-one-positive-dose-out CV (7 folds), condition-level RMSE on g,
# equal condition weighting. Uncertainty: block bootstrap over (T, replicate-
# series), B=2000, each series carried whole. Thresholds fixed: win >= 0.90,
# materiality >= 0.25*nu. (numpy hot path; statistics identical to spec.)
#
# RUN FROM the master "X0123 copy" folder:  python3 scripts/55_shift_scale_cv.py
# Outputs -> tables/revision/geometry/
# =============================================================================
import os, json
import numpy as np, pandas as pd
from scipy.optimize import least_squares

OUT = "tables/revision/geometry"
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(20260810)
B = 2000
WIN, MAT = 0.90, 0.25

TEMPS = [15, 18, 21, 24, 26, 27, 28]
POS = np.array([0.06, 0.12, 0.25, 0.5, 1.0, 2.0, 4.0])
DOSES = np.concatenate([[0.0], POS])                 # col 0 = control
nT, nD, nAll = len(TEMPS), len(POS), len(POS) + 1

raw = pd.read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv")
raw["conc"] = raw["Dose"].replace("Control", "0").astype(float)
raw["Temp"] = raw["T"].astype(float)
raw = raw[(raw["growth_C_per_C_h"] > 0) & np.isfinite(raw["growth_C_per_C_h"])]
raw = raw.rename(columns={"growth_C_per_C_h": "mu"})
assert len(raw) == 156, len(raw)

# T x series x dose tensor (NaN where a vial is absent). 3 series per T.
REPS = ["R1", "R2", "R3"]
Mten = np.full((nT, 3, nAll), np.nan)
for _, r in raw.iterrows():
    ti = TEMPS.index(int(r["Temp"])); ri = REPS.index(r["Replicate"])
    di = int(np.where(np.isclose(DOSES, r["conc"]))[0][0])
    Mten[ti, ri, di] = r["mu"]
assert np.isfinite(Mten).sum() == 156


def g_from_tensor(M):
    """condition-mean relative growth g[T,pos-dose]; NaN where absent."""
    cm = np.nanmean(M, axis=1)                       # nT x nD
    ctrl = cm[:, 0:1]
    g = cm[:, 1:] / ctrl                             # nT x nD-1
    return g                                          # NaN propagates


G_OBS = g_from_tensor(Mten)
mask_obs = np.isfinite(G_OBS)
print(f"positive-dose conditions with data: {int(mask_obs.sum())} of 49")

# ---- Hill and model predictions on the fixed POS grid -----------------------
def Hgrid(h, k):
    ch = POS ** h
    return ch / (ch + k ** h)                        # length nD


def predict_surface(model, theta):
    """return nT x nD-1 predicted s for the whole positive grid."""
    if model == "scale":
        h, k0 = np.exp(theta[0]), np.exp(theta[1])
        alpha = 1 / (1 + np.exp(-theta[2:9]))
        return alpha[:, None] * Hgrid(h, k0)[None, :]
    if model == "shift":
        h, M = np.exp(theta[0]), 1 / (1 + np.exp(-theta[1]))
        k = np.exp(theta[2:9])
        return M * np.stack([Hgrid(h, k[i]) for i in range(nT)])
    if model == "mixed":
        h = np.exp(theta[0])
        alpha = 1 / (1 + np.exp(-theta[1:8]))
        k = np.exp(theta[8:15])
        return alpha[:, None] * np.stack([Hgrid(h, k[i]) for i in range(nT)])


NPAR = {"scale": 9, "shift": 9, "mixed": 15}
X0 = {"scale": np.r_[0.0, 0.0, np.zeros(7)],
      "shift": np.r_[0.0, 0.0, np.zeros(7)],
      "mixed": np.r_[0.0, np.zeros(7), np.zeros(7)]}


def fit(model, g, train_col, x0):
    """fit on training dose columns (indices into POS); g is nT x nD-1."""
    tr = np.array(train_col)
    obs = 1.0 - g[:, tr]
    m = np.isfinite(obs)

    def resid(theta):
        s = predict_surface(model, theta)[:, tr]
        return (s[m] - obs[m])
    f = least_squares(resid, x0, method="lm", max_nfev=4000)
    return f.x


# observed warm-start thetas per (model, held dose) ---------------------------
WARM = {m: {} for m in NPAR}
for m in NPAR:
    for h in range(nD):
        train = [j for j in range(nD) if j != h]
        WARM[m][h] = fit(m, G_OBS, train, X0[m])


def loto_err(model, g, warm):
    """leave-one-dose-out held-out errors on relative growth; nT x nD-1, NaN off-fold cells kept NaN."""
    E = np.full((nT, nD), np.nan)
    for h in range(nD):
        train = [j for j in range(nD) if j != h]
        theta = fit(model, g, train, warm[model][h])
        s = predict_surface(model, theta)
        gp = 1.0 - s[:, h]
        col = g[:, h]
        ok = np.isfinite(col)
        E[ok, h] = gp[ok] - col[ok]
    return E


def pooled_rmse(E):
    e = E[np.isfinite(E)]
    return float(np.sqrt(np.mean(e ** 2)))


# =========================== observed-data result ============================
E_obs = {m: loto_err(m, G_OBS, WARM) for m in NPAR}
rmse = {m: pooled_rmse(E_obs[m]) for m in NPAR}
mae = {m: float(np.mean(np.abs(E_obs[m][np.isfinite(E_obs[m])]))) for m in NPAR}
pdr = {m: {POS[h]: float(np.sqrt(np.nanmean(E_obs[m][:, h] ** 2)))
           for h in range(nD)} for m in NPAR}

print("\n=== pooled held-out RMSE / MAE on relative growth ===")
for m in ["scale", "shift", "mixed"]:
    print(f"  {m:6s}  RMSE={rmse[m]:.4f}  MAE={mae[m]:.4f}  (params {NPAR[m]})")

# =========================== block bootstrap ================================
rmse_b = {m: np.zeros(B) for m in NPAR}
nu_accum = np.zeros((nT, nD)); nu_sq = np.zeros((nT, nD)); nu_n = np.zeros((nT, nD))
for b in range(B):
    idx = RNG.integers(0, 3, size=(nT, 3))            # resample series per T
    Mb = np.stack([Mten[t, idx[t], :] for t in range(nT)])
    gb = g_from_tensor(Mb)                             # nT x nD (positive only)
    fin = np.isfinite(gb)
    nu_accum[fin] += gb[fin]
    nu_sq[fin] += gb[fin] ** 2
    nu_n[fin] += 1
    for m in NPAR:
        rmse_b[m][b] = pooled_rmse(loto_err(m, gb, WARM))
    if (b + 1) % 400 == 0:
        print(f"  bootstrap {b+1}/{B}")

# nu = RMS over conditions of bootstrap SD of g
var_cells = (nu_sq - nu_accum ** 2 / np.maximum(nu_n, 1)) / np.maximum(nu_n - 1, 1)
var_cells = var_cells[nu_n > 1]
nu = float(np.sqrt(np.mean(var_cells)))
floor = MAT * nu

d = rmse_b["shift"] - rmse_b["scale"]                 # + favours scale
p_scale = float(np.mean(rmse_b["scale"] < rmse_b["shift"]))
p_shift = 1.0 - p_scale
med_adv = float(np.median(np.abs(d)))

best_single = "scale" if rmse["scale"] < rmse["shift"] else "shift"
d_mix = rmse_b[best_single] - rmse_b["mixed"]
p_mix = float(np.mean(rmse_b["mixed"] < rmse_b[best_single]))
med_mix = float(np.median(d_mix))

print("\n=== decision quantities ===")
print(f"nu = {nu:.4f}   materiality floor 0.25*nu = {floor:.4f}")
print(f"p_scale = {p_scale:.3f}   p_shift = {p_shift:.3f}")
print(f"observed RMSE: scale {rmse['scale']:.4f}  shift {rmse['shift']:.4f}  "
      f"(shift - scale = {rmse['shift']-rmse['scale']:+.4f})")
print(f"median |RMSE difference| = {med_adv:.4f}   vs floor {floor:.4f}")

scale_clear = (rmse["scale"] < rmse["shift"]) and (p_scale >= WIN) and (med_adv >= floor)
shift_clear = (rmse["shift"] < rmse["scale"]) and (p_shift >= WIN) and (med_adv >= floor)
if scale_clear:
    outcome = "A (scaling clearly wins)"
elif shift_clear:
    outcome = "B (translation clearly wins)"
else:
    outcome = "C1 (mixed geometry needed)" if (med_mix >= floor and p_mix >= WIN) \
        else "C2 (underpowered; potency vs amplitude not separable)"
print(f"\n=== PRE-REGISTERED OUTCOME: {outcome} ===")
print(f"mixed diagnostic: best single = {best_single}; "
      f"median(RMSE_best - RMSE_mixed) = {med_mix:+.4f}; P(mixed better) = {p_mix:.3f}")

print("\n=== per-dose held-out RMSE (edge = extrapolation) ===")
for h in range(nD):
    tag = "extrapolation" if POS[h] in (0.06, 4.0) else "interpolation"
    print(f"  dose {POS[h]:>4g}  scale {pdr['scale'][POS[h]]:.4f}  "
          f"shift {pdr['shift'][POS[h]]:.4f}  mixed {pdr['mixed'][POS[h]]:.4f}  [{tag}]")

# ---- outputs ---------------------------------------------------------------
pd.DataFrame([{"model": m, "pooled_RMSE": rmse[m], "pooled_MAE": mae[m],
               "n_params": NPAR[m], "boot_RMSE_median": float(np.median(rmse_b[m])),
               "boot_RMSE_p05": float(np.percentile(rmse_b[m], 5)),
               "boot_RMSE_p95": float(np.percentile(rmse_b[m], 95))}
              for m in NPAR]).to_csv(os.path.join(OUT, "shift_scale_cv_summary.csv"), index=False)
pd.DataFrame([{"dose": POS[h], "RMSE_scale": pdr["scale"][POS[h]],
               "RMSE_shift": pdr["shift"][POS[h]], "RMSE_mixed": pdr["mixed"][POS[h]],
               "fold_type": "extrapolation_edge" if POS[h] in (0.06, 4.0) else "interpolation_interior"}
              for h in range(nD)]).to_csv(os.path.join(OUT, "shift_scale_per_dose_rmse.csv"), index=False)
rows = []
for ti, T in enumerate(TEMPS):
    for h in range(nD):
        if np.isfinite(G_OBS[ti, h]):
            rows.append({"T": T, "conc": POS[h], "g_obs": G_OBS[ti, h],
                         "err_scale": E_obs["scale"][ti, h], "err_shift": E_obs["shift"][ti, h],
                         "err_mixed": E_obs["mixed"][ti, h]})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "shift_scale_cell_errors.csv"), index=False)
np.savez(os.path.join(OUT, "shift_scale_bootstrap.npz"),
         rmse_scale=rmse_b["scale"], rmse_shift=rmse_b["shift"], rmse_mixed=rmse_b["mixed"])
json.dump({"nu": nu, "materiality_floor": floor, "win_threshold": WIN,
           "p_scale": p_scale, "p_shift": p_shift, "RMSE_scale": rmse["scale"],
           "RMSE_shift": rmse["shift"], "RMSE_mixed": rmse["mixed"],
           "median_abs_rmse_diff": med_adv, "best_single": best_single,
           "p_mixed_better": p_mix, "median_mixed_gain": med_mix,
           "outcome": outcome, "B": B}, open(os.path.join(OUT, "shift_scale_decision.json"), "w"), indent=2)
print(f"\nwrote -> {OUT}")
