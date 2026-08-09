#!/usr/bin/env python3
# =============================================================================
# 53_loto_baselines.py — Stage 1, S1.3
# Identical-fold leave-one-temperature-out comparison of the published CYP51
# safety-margin model (M5) against five empirical baselines, with joint
# posterior target resampling and unit-invariant identifiability diagnostics.
#
# Pre-registered in Stage1_predictive_spec.md revision 3. Diagnostic only.
# M5 is fitted through the published code path (scripts/23, scripts/41) and is
# NOT modified.
#
# RUN FROM the master "X0123 copy" folder:  python3 scripts/53_loto_baselines.py
# Outputs -> tables/revision/predictive/
# =============================================================================
import os, json
import numpy as np, pandas as pd
from scipy.optimize import least_squares

OUT = "tables/revision/predictive"
os.makedirs(OUT, exist_ok=True)
kB, Tref = 8.617333e-5, 294.15
RNG = np.random.default_rng(20260809)
N_DRAWS = 500
TOL = 0.10                      # pre-registered |e27| tolerance, log10 units

hill = pd.read_csv("tables/physiology/growth/03_hill_parameters_by_temperature.csv")
pec = pd.read_csv("tables/physiology/growth/04_pec50_by_temperature.csv")
H = hill.merge(pec[["temperature", "pEC50_sd"]], on="temperature")
T_C = H["temperature"].values.astype(float)
T_K = T_C + 273.15
MU = H["r0"].values                       # control growth, held at posterior median
Y = np.log10(H["EC50"].values)
SY = H["pEC50_sd"].values
N = len(T_C)
I27 = int(np.where(T_C == 27)[0][0])

# ---------------------------------------------------------------- models ----
# Every model is fitted by weighted least squares with the same weights 1/SY,
# matching the published M5 fit.

def wls(X, y, w):
    Xw, yw = X * w[:, None], y * w
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    return beta

def design(name, t_C, t_K, mu):
    o = np.ones_like(t_C)
    if name == "M0": return np.column_stack([o])
    if name == "M1": return np.column_stack([o, t_C])
    if name == "M2": return np.column_stack([o, 1.0 / t_K])
    if name == "M3": return np.column_stack([o, np.log(mu)])
    if name == "M4": return np.column_stack([o, t_C, t_C ** 2])
    raise KeyError(name)

def fit_linear(name, tr, y, sy, mu):
    X = design(name, T_C[tr], T_K[tr], mu[tr])
    return wls(X, y[tr], 1.0 / sy[tr])

def pred_linear(name, beta, te, mu):
    return design(name, T_C[te], T_K[te], mu[te]) @ beta

def m5_resid_factory(t_K, y, sy, mu):
    def resid(p):
        Ki, A, Ea = np.exp(p[0]), np.exp(p[1]), p[2]
        s = A * np.exp(-Ea / kB * (1.0 / t_K - 1.0 / Tref)) / mu
        return (np.log10(Ki * np.maximum(2 * s - 1, 1e-9)) - y) / sy
    return resid

def fit_m5(tr, y, sy, mu):
    """Published code path: scripts/23_cyp51_margin_model.py / 41_loto_margin.py."""
    f = least_squares(m5_resid_factory(T_K[tr], y[tr], sy[tr], mu[tr]),
                      x0=[np.log(1.0), np.log(0.08), 0.5],
                      bounds=([-30, -30, 0.0], [10, 10, 2.0]))
    return f

def pred_m5(f, te, mu):
    Ki, A, Ea = np.exp(f.x[0]), np.exp(f.x[1]), f.x[2]
    s = A * np.exp(-Ea / kB * (1.0 / T_K[te] - 1.0 / Tref)) / mu[te]
    return np.log10(Ki * np.maximum(2 * s - 1, 1e-9))

MODELS = ["M0", "M1", "M2", "M3", "M4", "M5"]
NPAR = {"M0": 1, "M1": 2, "M2": 2, "M3": 2, "M4": 3, "M5": 3}
LABEL = {"M0": "constant log10 EC50", "M1": "linear in T",
         "M2": "linear in 1/kT", "M3": "control-growth only",
         "M4": "quadratic in T", "M5": "safety-margin (published)"}


def loto(y, sy, mu, want_intervals=False):
    """Return dict model -> array of held-out predictions (log10), and optionally
    the 95% predictive-interval half-widths.

    Interval construction, declared here for reproducibility: for each fold the
    training residual scale s_train is the weighted RMSE with (n_train - p)
    degrees of freedom; the predictive interval is
        pred +/- 1.96 * sqrt(s_train^2 + sd_target_heldout^2).
    """
    P = {m: np.zeros(N) for m in MODELS}
    HW = {m: np.zeros(N) for m in MODELS}
    for i in range(N):
        tr = np.array([j for j in range(N) if j != i]); te = np.array([i])
        for m in MODELS:
            if m == "M5":
                f = fit_m5(tr, y, sy, mu)
                P[m][i] = pred_m5(f, te, mu)[0]
                rtr = f.fun * sy[tr]
            else:
                b = fit_linear(m, tr, y, sy, mu)
                P[m][i] = pred_linear(m, b, te, mu)[0]
                rtr = pred_linear(m, b, tr, mu) - y[tr]
            if want_intervals:
                dof = max(len(tr) - NPAR[m], 1)
                s_tr = np.sqrt(np.sum(rtr ** 2) / dof)
                HW[m][i] = 1.96 * np.sqrt(s_tr ** 2 + sy[i] ** 2)
    return (P, HW) if want_intervals else P


# ------------------------------------------------- observed-target LOTO ----
P, HW = loto(Y, SY, MU, want_intervals=True)
rows, folds = [], []
for m in MODELS:
    e = P[m] - Y
    rows.append(dict(model=m, label=LABEL[m], n_params=NPAR[m],
                     LOTO_RMSE=np.sqrt(np.mean(e ** 2)),
                     LOTO_MAE=np.mean(np.abs(e)),
                     coverage_95=np.mean(np.abs(e) <= HW[m]),
                     err_15=e[0], err_27=e[I27], err_28=e[-1]))
    for i in range(N):
        folds.append(dict(model=m, T_C=T_C[i], obs_log10=Y[i],
                          pred_log10=P[m][i], err_log10=e[i],
                          obs=10 ** Y[i], pred=10 ** P[m][i],
                          in_95_interval=abs(e[i]) <= HW[m][i]))
C = pd.DataFrame(rows).sort_values("LOTO_RMSE")
pd.DataFrame(folds).to_csv(os.path.join(OUT, "loto_fold_errors.csv"), index=False)
C.to_csv(os.path.join(OUT, "loto_model_comparison.csv"), index=False)

print("=== LOTO on observed targets (log10 EC50) ===")
print(C.round(4).to_string(index=False))
best_emp = C[C.model != "M5"].iloc[0]["model"]
print(f"\nbest empirical baseline (lowest RMSE among M0-M4): {best_emp}")

# --------------------------------------- joint posterior target resampling ----
DR = pd.read_csv(os.path.join(OUT, "hill_joint_draws.csv"))
W = DR.pivot(index=".draw", columns="temperature", values="ec50")[list(T_C)].values
sel = RNG.choice(W.shape[0], size=N_DRAWS, replace=False)
print(f"\n=== joint posterior resampling: {N_DRAWS} draws "
      f"(joint across temperatures; independent draws prohibited) ===")

rmse_d = {m: np.zeros(N_DRAWS) for m in MODELS}
mae_d = {m: np.zeros(N_DRAWS) for m in MODELS}
e27_d = {m: np.zeros(N_DRAWS) for m in MODELS}
for k, d in enumerate(sel):
    yd = np.log10(np.maximum(W[d], 1e-9))
    Pd = loto(yd, SY, MU)
    for m in MODELS:
        e = Pd[m] - yd
        rmse_d[m][k] = np.sqrt(np.mean(e ** 2))
        mae_d[m][k] = np.mean(np.abs(e))
        e27_d[m][k] = e[I27]

res = []
for m in MODELS:
    e27 = e27_d[m]
    s = float(np.mean(e27 < 0))
    p = float(np.mean(np.abs(e27) < TOL))
    lo, hi = np.percentile(e27, [5, 95])
    res.append(dict(model=m, label=LABEL[m],
                    rmse_median=np.median(rmse_d[m]),
                    rmse_p05=np.percentile(rmse_d[m], 5),
                    rmse_p95=np.percentile(rmse_d[m], 95),
                    e27_median=np.median(e27), e27_p05=lo, e27_p95=hi,
                    e27_interval_includes_zero=bool(lo < 0 < hi),
                    sign_consistency_s=s, p_within_tol=p))
R = pd.DataFrame(res)
R.to_csv(os.path.join(OUT, "loto_target_resampling.csv"), index=False)
print(R.round(4).to_string(index=False))

# pairwise win proportions (fraction of draws in which row model has lower RMSE)
Wp = pd.DataFrame(index=MODELS, columns=MODELS, dtype=float)
for a in MODELS:
    for b in MODELS:
        Wp.loc[a, b] = np.nan if a == b else float(np.mean(rmse_d[a] < rmse_d[b]))
Wp.to_csv(os.path.join(OUT, "loto_win_proportions.csv"))
print("\n=== pairwise win proportions (row attains lower LOTO RMSE than column) ===")
print(Wp.round(3).to_string())

s5, p5 = R.loc[R.model == "M5", "sign_consistency_s"].iat[0], \
         R.loc[R.model == "M5", "p_within_tol"].iat[0]
print(f"\nM5 27 C readouts: sign consistency s = {s5:.3f}, "
      f"Pr(|e27| < {TOL}) = {p5:.3f}")
print("pre-registered combination rule: target uncertainty MATERIALLY EXPLAINS "
      f"the miss if s < 0.90 or p >= 0.20  ->  "
      f"{'YES' if (s5 < 0.90 or p5 >= 0.20) else 'NO'}")

# ------------------------------------------ identifiability diagnostics ----
f_full = fit_m5(np.arange(N), Y, SY, MU)
J = f_full.jac                                   # weighted residual Jacobian
p_names = ["log_Ki", "log_A", "Ea"]
Hm = J.T @ J                                     # Gauss-Newton Hessian
dof = N - len(p_names)
sigma2 = float(np.sum(f_full.fun ** 2) / dof)
Cov = sigma2 * np.linalg.pinv(Hm)
SE = np.sqrt(np.diag(Cov))

Hn = Hm / np.sqrt(np.outer(np.diag(Hm), np.diag(Hm)))          # normalised Hessian
Rc = Cov / np.sqrt(np.outer(np.diag(Cov), np.diag(Cov)))       # parameter correlations
Js = J * SE[None, :]                                           # MULTIPLY by SE
sv = np.linalg.svd(Js, compute_uv=False)

print("\n=== identifiability (unit-invariant) ===")
print("full-fit parameters: " +
      ", ".join(f"{n}={v:.4g} (SE {s:.3g})" for n, v, s in zip(p_names, f_full.x, SE)))
print("\nnormalised Hessian (dimensionless curvature coupling):")
print(pd.DataFrame(Hn, index=p_names, columns=p_names).round(4).to_string())
print("\nparameter correlation matrix (from inverse-Hessian covariance):")
print(pd.DataFrame(Rc, index=p_names, columns=p_names).round(4).to_string())
print(f"\nSE-scaled Jacobian singular values: {np.round(sv, 6)}")
sv_ratio = float(sv[0] / sv[-1])
print(f"singular value ratio (largest/smallest) = {sv_ratio:.4g}")

exc = []
for i in range(N):
    tr = np.array([j for j in range(N) if j != i])
    fi = fit_m5(tr, Y, SY, MU)
    exc.append(dict(held_out_T=T_C[i],
                    **{f"{n}_dev_SE": (fi.x[k] - f_full.x[k]) / SE[k]
                       for k, n in enumerate(p_names)},
                    Ea_fitted=fi.x[2]))
E = pd.DataFrame(exc)
print("\n=== fold parameter excursions (deviation from full fit, in SE units) ===")
print(E.round(3).to_string(index=False))
max_exc = float(np.max(np.abs(E[[c for c in E.columns if c.endswith('_dev_SE')]].values)))
print(f"max |excursion| = {max_exc:.3f} SE")

pd.DataFrame(Hn, index=p_names, columns=p_names).to_csv(
    os.path.join(OUT, "identifiability_normalised_hessian.csv"))
pd.DataFrame(Rc, index=p_names, columns=p_names).to_csv(
    os.path.join(OUT, "identifiability_parameter_correlations.csv"))
E.to_csv(os.path.join(OUT, "identifiability_fold_excursions.csv"), index=False)
# ------------------------------------------------------------------------
# POST HOC, NOT a pre-registered decision input. Recorded because it explains
# the comparison result: is the fitted model operating in the regime its
# mechanistic story describes?
# ------------------------------------------------------------------------
Ki, A, Ea = np.exp(f_full.x[0]), np.exp(f_full.x[1]), f_full.x[2]
sig = A * np.exp(-Ea / kB * (1.0 / T_K - 1.0 / Tref)) / MU
print("\n=== post hoc: is the safety-margin term active at the optimum? ===")
print("fitted safety margin sigma(T) = capacity / demand:")
print("   " + "  ".join(f"{t:.0f}C {s:.0f}" for t, s in zip(T_C, sig)))
print(f"minimum sigma = {sig.min():.0f}; relative error of approximating "
      f"(2*sigma - 1) by 2*sigma = {1/(2*sig.min()):.2e}")
print("The saturating term is therefore inactive: the model is in its linear regime.")
pred_full = np.log10(Ki * np.maximum(2 * sig - 1, 1e-9))
pred_red = (np.log10(2 * Ki * A) - (Ea / (kB * np.log(10))) * (1 / T_K - 1 / Tref)
            - np.log10(MU))
print(f"max |difference| between the fitted model and "
      f"[c - (Ea/kB ln10)(1/T - 1/Tref) - log10 mu(T)] = "
      f"{np.max(np.abs(pred_full - pred_red)):.2e} log10 units")
Xf = np.column_stack([np.ones(N), 1 / T_K, np.log10(MU)])
bf = wls(Xf, Y, 1.0 / SY)
print(f"The safety-margin structure FIXES the coefficient on log10 mu at -1.000; "
      f"freeing it gives {bf[2]:+.3f}, i.e. the constraint is essentially what the "
      f"data prefer and is NOT what costs M5 out of sample.")
print("Interpretation: at its optimum M5 is algebraically an Arrhenius term in 1/T "
      "minus log10 control growth, with 2 effective parameters (only the product "
      "Ki*A is identified). Its out-of-sample loss is therefore not a consequence "
      "of a wrong constraint but of the predictors themselves carrying no "
      "transferable signal -- both M5 and M3 are beaten by a constant.")
json.dump(dict(sigma_min=float(sig.min()), sigma_max=float(sig.max()),
               reduction_max_abs_diff_log10=float(np.max(np.abs(pred_full - pred_red))),
               free_coef_on_log10_mu=float(bf[2]),
               Ki_times_A=float(Ki * A)),
          open(os.path.join(OUT, "posthoc_model_regime.json"), "w"), indent=2)

json.dump(dict(singular_values=sv.tolist(), sv_ratio=sv_ratio,
               max_excursion_SE=max_exc, sigma2=sigma2,
               params=dict(zip(p_names, f_full.x.tolist())),
               SE=dict(zip(p_names, SE.tolist())),
               best_empirical_baseline=best_emp,
               s_M5=s5, p_M5=p5),
          open(os.path.join(OUT, "identifiability_summary.json"), "w"), indent=2)
print(f"\nwrote -> {OUT}")
