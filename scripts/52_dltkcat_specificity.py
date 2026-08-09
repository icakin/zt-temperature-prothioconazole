#!/usr/bin/env python3
# =============================================================================
# 52_dltkcat_specificity.py — Stage 1, S1.2
# Does DLTKcat's prediction for CYP51 (r_0317) carry ENZYME-SPECIFIC temperature
# information, or is it effectively the generic response shared across the
# model's enzymes?
#
# NOT a test of whether DLTKcat agrees with the fitted 0.47 eV. Agreement is not
# required and disagreement is not disqualifying (spec section 3).
#
# Pre-registered window: 15-28 C. The DLTKcat grid is 5-45 C in 4 C steps, so
# the grid points strictly inside the window are 17, 21 and 25 C. That set is
# the PRIMARY analysis, exactly as specified. The smallest bracketing set
# (13-29 C) is reported as a labelled sensitivity check only.
#
# RUN FROM the master "X0123 copy" folder:  python3 scripts/52_dltkcat_specificity.py
# Outputs -> tables/revision/predictive/
# =============================================================================
import os, sys
import numpy as np, pandas as pd

OUT = "tables/revision/predictive"
os.makedirs(OUT, exist_ok=True)
kB = 8.617333e-5          # eV/K
CYP51 = "r_0317"
WIN_PRIMARY = (15.0, 28.0)
WIN_SENS = (13.0, 29.0)

D = pd.read_csv("etcGEMs/strains/zt_ipo323/dltkcat/output.csv")
D["rxn"] = D["rxn_id"].str.split("_No").str[0]
print(f"reactions with DLTKcat predictions: {D['rxn'].nunique()}")
print(f"temperature grid: {sorted(D['Temp_C'].unique())}\n")


def apparent_Ea(sub):
    """Slope of ln kcat vs -1/kT  ->  apparent activation energy in eV."""
    if len(sub) < 2:
        return np.nan
    x = -1.0 / (kB * (sub["Temp_C"].values + 273.15))
    y = np.log(10.0) * sub["pred_log10kcat"].values          # ln kcat
    return np.polyfit(x, y, 1)[0]


def analyse(win, label):
    lo, hi = win
    W = D[(D["Temp_C"] >= lo) & (D["Temp_C"] <= hi)]
    grid = sorted(W["Temp_C"].unique())
    print(f"--- {label} window {lo}-{hi} C, grid points used: {grid} ---")

    # per-reaction apparent Ea (median over isozyme/substrate arms)
    ea = (W.groupby(["rxn", "rxn_id"]).apply(apparent_Ea, include_groups=False)
            .rename("Ea_eV").reset_index()
            .groupby("rxn")["Ea_eV"].median())
    ea = ea[np.isfinite(ea)]

    # normalised kcat(T) curves, each divided by its own value at 21 C
    piv = (W.groupby(["rxn", "Temp_C"])["pred_log10kcat"].median().unstack())
    if 21.0 not in piv.columns:
        raise RuntimeError("21 C reference not on the grid")
    norm = piv.sub(piv[21.0], axis=0)              # log10 scale => ratio
    norm = norm.dropna()
    med_curve = norm.median(axis=0)
    dev = np.sqrt(((norm - med_curve) ** 2).mean(axis=1))     # RMS deviation

    if CYP51 not in ea.index or CYP51 not in dev.index:
        raise RuntimeError("CYP51 (r_0317) missing from DLTKcat outputs")

    ea_c, dev_c = ea[CYP51], dev[CYP51]
    ea_pct = 100.0 * (ea < ea_c).mean()
    dev_pct = 100.0 * (dev < dev_c).mean()

    print(f"  reactions analysed: {len(ea)}")
    print(f"  apparent Ea distribution: median {ea.median():.3f} eV, "
          f"10th {ea.quantile(0.10):.3f}, 90th {ea.quantile(0.90):.3f} eV")
    print(f"  CYP51 apparent Ea = {ea_c:.3f} eV  ->  percentile {ea_pct:.1f}")
    print(f"  normalised-curve RMS deviation distribution: median "
          f"{dev.median():.4f}, 75th {dev.quantile(0.75):.4f}")
    print(f"  CYP51 curve deviation = {dev_c:.4f}  ->  percentile {dev_pct:.1f}")
    print(f"  CYP51 normalised curve: "
          + ", ".join(f"{t:.0f}C {10**norm.loc[CYP51, t]:.3f}" for t in grid))
    print(f"  median normalised curve: "
          + ", ".join(f"{t:.0f}C {10**med_curve[t]:.3f}" for t in grid))

    # pre-registered gate (spec section 6, Outcome C / E)
    gate_ea = (ea_pct <= 10.0) or (ea_pct >= 90.0)
    gate_dev = dev_pct >= 75.0
    veto = dev_pct < 50.0
    print(f"  gate: Ea outside central 80%? {gate_ea} | "
          f"curve deviation above 75th percentile? {gate_dev}")
    print(f"  VETO (curve deviation below median -> information is generic): {veto}")
    print(f"  => specificity gate {'PASSES' if (gate_ea and gate_dev and not veto) else 'FAILS'}\n")

    return dict(window=label, lo=lo, hi=hi, n_reactions=len(ea),
                cyp51_Ea_eV=ea_c, cyp51_Ea_percentile=ea_pct,
                dist_Ea_median=ea.median(), dist_Ea_p10=ea.quantile(0.10),
                dist_Ea_p90=ea.quantile(0.90),
                cyp51_curve_deviation=dev_c, cyp51_curve_dev_percentile=dev_pct,
                dist_dev_median=dev.median(), dist_dev_p75=dev.quantile(0.75),
                gate_Ea_outside_central80=gate_ea,
                gate_dev_above_p75=gate_dev, veto_generic=veto,
                specificity_gate_passes=bool(gate_ea and gate_dev and not veto))


rows = [analyse(WIN_PRIMARY, "PRIMARY (pre-registered)"),
        analyse(WIN_SENS, "SENSITIVITY (bracketing, not pre-registered)")]
R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "dltkcat_cyp51_specificity.csv"), index=False)

print("Context only, carrying no decision weight:")
print(f"  fitted safety-margin activation energy = 0.472 eV")
print(f"  DLTKcat CYP51 apparent Ea (primary)    = {rows[0]['cyp51_Ea_eV']:.3f} eV")
print(f"\nwrote -> {OUT}/dltkcat_cyp51_specificity.csv")
