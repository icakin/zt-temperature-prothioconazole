"""Turn temperature-dependent kcat predictions (DLTKcat) into MMRT parameters.

DLTKcat predicts kcat as a function of temperature for an enzyme-substrate pair.
This module fits an MMRT curve to those kcat-vs-T points to recover the two
interpretable knobs the rest of the pipeline uses -- ``Topt`` and ``dCp`` -- and
then either

  * applies them in place to a GECKO provider (recommended for eciML1515: keeps
    the model's calibrated kcat costs, overrides only the thermal shape), or
  * writes a csv-provider table (rxn_id, mw, kcat, Topt, dCp, group, T0) for the
    plain-GEM route.

The MMRT fit is *linear* in (dH, dS, dCp):

    ln k(T) - ln(kB*T/h)
        = dH * (-1/(R T)) + dS * (1/R)
          + dCp * ( -(T-T0)/(R T) + ln(T/T0)/R )

so it's an ordinary least-squares solve -- fast and robust, no initial guess.

Workflow
--------
1. ``export_targets(model, "targets.csv")`` -> a long table of
   (rxn_id, enzyme_id, temp_C, kcat) with kcat blank. Fill the kcat column with
   DLTKcat predictions (you supply each enzyme's sequence + substrate SMILES +
   the temperatures to DLTKcat itself).
2. ``fit_predictions(df)`` -> per-enzyme Topt/dCp/kcat_T0/r2.
3. ``apply_fits_to_provider(pm, fits)`` or ``write_csv_table(...)``.
"""
from __future__ import annotations

import argparse
from typing import Dict, Optional

import numpy as np
import pandas as pd

from .mmrt import R, _ln_prefactor


# ---------------------------------------------------------------------------
# Core fit
# ---------------------------------------------------------------------------
# DLTKcat global skill: log10(kcat) RMSE ~0.9 -> a floor on the fit residual
# sigma (in ln units) so optimistic per-enzyme residuals don't understate the
# parameter uncertainty. ln = log10 * ln(10).
import math as _math
_SIGMA_FLOOR_LN = 0.9 * _math.log(10.0)


def _topt_of_grid(g, T0, dH, dCp, dS):
    """Argmax-Topt of an MMRT curve on grid g (vectorised over draws if dH/... are
    columns)."""
    lnk = (_ln_prefactor(g) + (-dH - dCp * (g - T0)) / (R * g)
           + (dS + dCp * np.log(g / T0)) / R)
    return g[np.argmax(lnk, axis=-1)]


def fit_mmrt(temps_C, kcats, T0_C: float = 30.0,
             grid=(273.15, 353.15, 4001), with_uncertainty: bool = True,
             n_draws: int = 64, sigma_floor_ln: float = _SIGMA_FLOOR_LN,
             seed: int = 0) -> Dict[str, float]:
    """Fit MMRT to kcat-vs-temperature points by linear least squares.

    Returns dict with Topt_C, dCp, kcat_T0, dH, dS, r2, n, ok, and (if
    with_uncertainty) Topt_sd, dCp_sd. The fit is linear in (dH, dS, dCp), so it
    has covariance Cov = sigma^2 (X^T X)^-1; sd_Topt/sd_dCp come from sampling
    (dH,dS,dCp) ~ N(coef, Cov) and computing Topt/dCp per draw (Topt is nonlinear
    in the coefficients, so we sample rather than linearise). sigma is floored by
    DLTKcat's global skill.
    ``ok`` is False if there are too few points, non-positive kcats, a
    non-negative dCp (no thermal peak), or a poor fit.
    """
    T = np.asarray(temps_C, float) + 273.15
    k = np.asarray(kcats, float)
    m = np.isfinite(T) & np.isfinite(k) & (k > 0)
    T, k = T[m], k[m]
    out = dict(Topt_C=np.nan, dCp=np.nan, kcat_T0=np.nan, dH=np.nan, dS=np.nan,
               r2=np.nan, n=int(T.size), ok=False, Topt_sd=np.nan, dCp_sd=np.nan)
    if T.size < 4 or np.unique(np.round(T, 3)).size < 4:
        return out
    T0 = T0_C + 273.15
    y = np.log(k) - _ln_prefactor(T)
    X = np.column_stack([
        -1.0 / (R * T),                                     # dH
        np.full_like(T, 1.0 / R),                           # dS
        -(T - T0) / (R * T) + np.log(T / T0) / R,           # dCp
    ])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    dH, dS, dCp = map(float, coef)

    # fitted curve on a fine grid -> numeric Topt
    g = np.linspace(*grid)
    lnk_g = _ln_prefactor(g) + (-dH - dCp * (g - T0)) / (R * g) + \
        (dS + dCp * np.log(g / T0)) / R
    Topt = float(g[int(np.argmax(lnk_g))])
    kcat_T0 = float(np.exp(_ln_prefactor(T0) - dH / (R * T0) + dS / R))

    yhat = X @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    peaked = dCp < 0 and grid[0] + 1 < Topt < grid[1] - 1
    out.update(Topt_C=Topt - 273.15, dCp=dCp, kcat_T0=kcat_T0, dH=dH, dS=dS,
               r2=r2, ok=bool(peaked and (np.isnan(r2) or r2 > 0.8)))

    if with_uncertainty:
        try:
            dof = max(T.size - 3, 1)
            sigma2 = max(ss_res / dof, sigma_floor_ln ** 2)
            cov = sigma2 * np.linalg.inv(X.T @ X)
            draws = np.random.default_rng(seed).multivariate_normal(coef, cov, size=n_draws)
            dH_d, dS_d, dCp_d = draws[:, 0], draws[:, 1], draws[:, 2]
            lnk = (_ln_prefactor(g)[None, :]
                   + (-dH_d[:, None] - dCp_d[:, None] * (g[None, :] - T0)) / (R * g[None, :])
                   + (dS_d[:, None] + dCp_d[:, None] * np.log(g / T0)[None, :]) / R)
            Topt_draws = g[np.argmax(lnk, axis=1)]
            out["Topt_sd"] = float(np.std(Topt_draws))
            out["dCp_sd"] = float(np.std(dCp_d))
        except np.linalg.LinAlgError:
            pass
    return out


def fit_predictions(df: pd.DataFrame, key: str = "rxn_id", T0_C: float = 30.0,
                    temp_col: str = "temp_C", kcat_col: str = "kcat") -> pd.DataFrame:
    """Fit MMRT per group in a long-format prediction table.

    df columns: <key>, temp_col, kcat_col (e.g. rxn_id/enzyme_id, temp_C, kcat).
    """
    rows = []
    for kval, sub in df.groupby(key):
        fit = fit_mmrt(sub[temp_col].values, sub[kcat_col].values, T0_C)
        fit[key] = kval
        rows.append(fit)
    cols = [key, "Topt_C", "dCp", "kcat_T0", "r2", "n", "ok", "dH", "dS",
            "Topt_sd", "dCp_sd"]
    return pd.DataFrame(rows)[cols]


# ---------------------------------------------------------------------------
# GECKO enzyme map / target export
# ---------------------------------------------------------------------------
def enzyme_map(model_path: str, T0_C: float = 30.0, **gecko_kw) -> pd.DataFrame:
    """Per-reaction (rxn_id, enzyme_id, mw, kcat_ref_s, group) from a GECKO model."""
    from .providers import from_gecko
    pm = from_gecko(model_path, T0=T0_C + 273.15, **gecko_kw)
    return pd.DataFrame([
        dict(rxn_id=e.rxn_id, enzyme_id=e.enzyme_id, mw=e.mw,
             kcat_ref_s=e.kcat_ref, group=e.group)
        for e in pm.ec.table
    ])


def export_targets(model_path: str, out_csv: str, tmin: float = 5.0,
                   tmax: float = 55.0, n: int = 11, T0_C: float = 30.0,
                   **gecko_kw) -> pd.DataFrame:
    """Write a long template (rxn_id, enzyme_id, mw, group, temp_C, kcat=blank)
    to fill with DLTKcat predictions."""
    emap = enzyme_map(model_path, T0_C=T0_C, **gecko_kw)
    temps = np.linspace(tmin, tmax, n)
    rows = []
    for _, r in emap.iterrows():
        for Tc in temps:
            rows.append(dict(rxn_id=r.rxn_id, enzyme_id=r.enzyme_id,
                             mw=r.mw, group=r.group, temp_C=round(float(Tc), 2),
                             kcat=""))
    out = pd.DataFrame(rows)
    out.to_csv(out_csv, index=False)
    print(f"[dltkcat] wrote {len(out)} target rows "
          f"({emap.shape[0]} reactions x {n} temps) -> {out_csv}")
    return out


# ---------------------------------------------------------------------------
# DLTKcat input preparation
# ---------------------------------------------------------------------------
# Currency / cofactor metabolites that are never "the substrate" for kcat.
CURRENCY = {
    "h2o", "h", "co2", "o2", "nh4", "hco3", "pi", "ppi", "pppi",
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "ctp", "cdp", "cmp",
    "utp", "udp", "ump", "itp", "idp", "imp", "ttp", "dttp",
    "datp", "dadp", "dgtp", "dctp", "dutp", "dttp",
    "nad", "nadh", "nadp", "nadph", "fad", "fadh2", "fmn", "fmnh2",
    "coa", "so3", "so4", "h2s", "h2o2", "no2", "no3", "n2",
    "na1", "k", "cl", "mg2", "ca2", "fe2", "fe3", "mn2", "zn2",
    "cu2", "cu", "cobalt2", "mobd", "ni2", "cd2", "q8", "q8h2",
    "mql8", "mqn8", "2dmmql8", "2dmmq8", "actp",
}


def _carbon_count(formula: Optional[str]) -> int:
    import re
    if not formula:
        return 0
    m = re.search(r"C(\d*)(?![a-z])", formula)  # 'C' not followed by lowercase (avoid Cl, Ca)
    if not m:
        return 0
    return int(m.group(1)) if m.group(1) else 1


def _base_bigg(met) -> str:
    """Compartment-free BiGG base id, e.g. 'atp_c[c]' / 'atp_c' -> 'atp'."""
    b = met.annotation.get("bigg.metabolite") if met.annotation else None
    if isinstance(b, list):
        b = b[0]
    if b:
        return str(b)
    s = met.id
    for suf in ("[c]", "[e]", "[p]", "[m]"):
        s = s.replace(suf, "")
    return s.rsplit("_", 1)[0] if "_" in s else s


def _clean_name(name: str) -> str:
    """'Cytidine C9H13N3O5 [cytosol]' -> 'Cytidine' (best effort)."""
    import re
    s = re.sub(r"\s*\[[^\]]*\]\s*$", "", name).strip()      # drop trailing [compartment]
    s = re.sub(r"\s+C\d+H\d+[A-Za-z0-9]*\s*$", "", s).strip()  # drop trailing formula
    return s


def _is_pmet(met) -> bool:
    """GECKO per-reaction arm pseudo-metabolite (pmet_<RXN>), no real chemistry."""
    return met.id.startswith("pmet") or (met.formula is None and _base_bigg(met) == "pmet")


def select_substrate(rxn, model=None, currency=CURRENCY, _depth=0):
    """Pick the primary substrate metabolite of a reaction (reactant side).

    Excludes enzymes and currency metabolites and prefers the largest carbon
    count. GECKO splits many reactions into an arm reaction (producing a
    ``pmet_<RXN>`` pseudo-metabolite) plus isozyme reactions that consume it; if
    the reactant is such a pmet we trace back through the arm reaction (needs
    ``model``) to recover the real substrate.
    """
    cands = []
    for mm, coef in rxn.metabolites.items():
        if coef >= 0 or mm.id.startswith("prot_"):
            continue
        if _is_pmet(mm):
            if model is not None and _depth < 3:
                for arm in mm.reactions:
                    if arm.id != rxn.id and arm.metabolites.get(mm, 0) > 0:
                        s = select_substrate(arm, model, currency, _depth + 1)
                        if s is not None:
                            cands.append(s)
            continue
        if _base_bigg(mm) in currency:
            continue
        cands.append(mm)
    if not cands:
        return None
    cands.sort(key=lambda m: _carbon_count(m.formula), reverse=True)
    return cands[0]


def build_dltkcat_input(model_path: str, out_csv: str, tmin: float = 5.0,
                        tmax: float = 55.0, n: int = 11, T0_C: float = 30.0,
                        **gecko_kw) -> pd.DataFrame:
    """Write the CSV DLTKcat's `convert_input` consumes.

    One row per (reaction, temperature) with the enzyme UniProt id and the
    chosen substrate *name* (DLTKcat resolves SMILES + sequence from these).
    Also carries rxn_id/bigg/mnx columns so predictions can be mapped back and
    any unresolved substrates fixed by hand. Reactions with no non-currency
    substrate are reported and skipped.
    """
    from .providers import from_gecko
    pm = from_gecko(model_path, T0=T0_C + 273.15, **gecko_kw)
    model = pm.ec.model
    temps = np.linspace(tmin, tmax, n)
    rows, skipped = [], 0
    for e in pm.ec.table:
        rxn = model.reactions.get_by_id(e.rxn_id)
        sub = select_substrate(rxn, model=model)
        if sub is None or not e.enzyme_id:
            skipped += 1
            continue
        mnx = sub.annotation.get("metanetx.chemical", "") if sub.annotation else ""
        if isinstance(mnx, list):
            mnx = mnx[0]
        for Tc in temps:
            rows.append(dict(rxn_id=e.rxn_id, enz=e.enzyme_id,
                             sub=_clean_name(sub.name), bigg=_base_bigg(sub),
                             mnx=mnx, Temp_C=round(float(Tc), 2),
                             Temp_K=round(float(Tc) + 273.15, 2)))
    out = pd.DataFrame(rows)
    out.to_csv(out_csv, index=False)
    nrxn = out["rxn_id"].nunique() if len(out) else 0
    print(f"[dltkcat] wrote {len(out)} rows for {nrxn} reactions x {n} temps "
          f"({skipped} reactions skipped: no enzyme/substrate) -> {out_csv}")
    print("[dltkcat] next: DLTKcat convert_input(enz_col='enz', sub_col='sub') "
          "-> normalize temps -> predict.py; then parse_dltkcat_output(...)")
    return out


def parse_dltkcat_output(pred, kcat_col: str = "kcat", key: str = "rxn_id",
                         temp_col: str = "Temp_C", T0_C: float = 30.0,
                         log10: bool = True) -> pd.DataFrame:
    """Fit MMRT from a DLTKcat prediction table back to per-key Topt/dCp.

    DLTKcat outputs log10(kcat); set log10=True to exponentiate. Requires the
    <key>, temp_col and kcat_col columns (keep rxn_id through the DLTKcat run).
    """
    df = pred if isinstance(pred, pd.DataFrame) else pd.read_csv(pred)
    df = df.copy()
    df["kcat"] = 10.0 ** df[kcat_col] if log10 else df[kcat_col]
    return fit_predictions(df, key=key, T0_C=T0_C, temp_col=temp_col, kcat_col="kcat")


# ---------------------------------------------------------------------------
# Apply fits
# ---------------------------------------------------------------------------
def apply_fits_to_provider(pm, fits: pd.DataFrame, key: str = "rxn_id",
                           only_ok: bool = True) -> int:
    """Override Topt/dCp on a provider's enzyme table from fitted parameters.

    key = 'rxn_id' (per-reaction) or 'enzyme_id' (per-protein, broadcast to all
    reactions that enzyme catalyses). kcat/base_cost are left untouched so the
    model's calibrated enzyme costs are preserved. Returns number of entries
    updated; call happens in place and refreshes the constraint.
    """
    if only_ok and "ok" in fits.columns:
        fits = fits[fits["ok"]]
    lut = {row[key]: row for _, row in fits.iterrows()}
    n = 0
    for e in pm.ec.table:
        k = e.rxn_id if key == "rxn_id" else e.enzyme_id
        if k in lut:
            row = lut[k]
            e.Topt = float(row["Topt_C"]) + 273.15
            e.dCp = float(row["dCp"])
            # also feed the two-state/unfolding transition-state heat capacity
            # (dCpt, J/mol/K) so DLTKcat grounds the envelope curvature in
            # unfolding mode (harmless in mmrt mode, which ignores dCpt).
            e.dCpt = float(row["dCp"]) * 1000.0
            n += 1
    pm.ec.refresh_params()
    print(f"[dltkcat] applied fits to {n}/{len(pm.ec.table)} enzymes (key={key})")
    return n


def write_csv_table(fits: pd.DataFrame, model_path: str, out_csv: str,
                    key: str = "rxn_id", T0_C: float = 30.0,
                    default_Topt_offset: float = 8.0, default_dCp: float = -8.0,
                    **gecko_kw) -> pd.DataFrame:
    """Write a from_kcat_csv table (rxn_id, mw, kcat, Topt, dCp, group, T0)
    merging fitted Topt/dCp with the model's mw/kcat/group. Reactions without a
    good fit fall back to the default thermal knobs."""
    emap = enzyme_map(model_path, T0_C=T0_C, **gecko_kw)
    ok = fits[fits["ok"]] if "ok" in fits.columns else fits
    lut = {row[key]: row for _, row in ok.iterrows()}
    rows = []
    for _, r in emap.iterrows():
        k = r.rxn_id if key == "rxn_id" else r.enzyme_id
        f = lut.get(k)
        Topt = float(f["Topt_C"]) + 273.15 if f is not None else T0_C + 273.15 + default_Topt_offset
        dCp = float(f["dCp"]) if f is not None else default_dCp
        rows.append(dict(rxn_id=r.rxn_id, mw=r.mw, kcat=r.kcat_ref_s,
                         Topt=Topt, dCp=dCp, group=r.group, T0=T0_C + 273.15))
    out = pd.DataFrame(rows)
    out.to_csv(out_csv, index=False)
    print(f"[dltkcat] wrote provider table for {len(out)} reactions "
          f"({len(lut)} with fitted thermal params) -> {out_csv}")
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="DLTKcat -> MMRT (Topt, dCp) tools")
    sub = ap.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("targets", help="export a prediction template")
    t.add_argument("--model", required=True)
    t.add_argument("--out", required=True)
    t.add_argument("--tmin", type=float, default=5.0)
    t.add_argument("--tmax", type=float, default=55.0)
    t.add_argument("--n", type=int, default=11)
    t.add_argument("--t0", type=float, default=30.0)

    f = sub.add_parser("fit", help="fit MMRT to a filled prediction table")
    f.add_argument("--pred", required=True)
    f.add_argument("--key", default="rxn_id")
    f.add_argument("--t0", type=float, default=30.0)
    f.add_argument("--out", required=True)

    c = sub.add_parser("csv", help="fit + write a from_kcat_csv provider table")
    c.add_argument("--pred", required=True)
    c.add_argument("--model", required=True)
    c.add_argument("--key", default="rxn_id")
    c.add_argument("--t0", type=float, default=30.0)
    c.add_argument("--out", required=True)

    p = sub.add_parser("prep", help="write DLTKcat convert_input CSV (enz, sub, temp)")
    p.add_argument("--model", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--tmin", type=float, default=5.0)
    p.add_argument("--tmax", type=float, default=55.0)
    p.add_argument("--n", type=int, default=11)
    p.add_argument("--t0", type=float, default=30.0)

    pr = sub.add_parser("parse", help="fit MMRT from a DLTKcat prediction table")
    pr.add_argument("--pred", required=True)
    pr.add_argument("--kcat-col", default="kcat")
    pr.add_argument("--key", default="rxn_id")
    pr.add_argument("--temp-col", default="Temp_C")
    pr.add_argument("--t0", type=float, default=30.0)
    pr.add_argument("--no-log10", action="store_true", help="predictions are raw kcat, not log10")
    pr.add_argument("--out", required=True)

    a = ap.parse_args(argv)
    if a.cmd == "targets":
        export_targets(a.model, a.out, a.tmin, a.tmax, a.n, a.t0)
    elif a.cmd == "fit":
        df = pd.read_csv(a.pred)
        fits = fit_predictions(df, key=a.key, T0_C=a.t0)
        fits.to_csv(a.out, index=False)
        print(f"[dltkcat] fitted {int(fits['ok'].sum())}/{len(fits)} usable -> {a.out}")
    elif a.cmd == "csv":
        df = pd.read_csv(a.pred)
        fits = fit_predictions(df, key=a.key, T0_C=a.t0)
        write_csv_table(fits, a.model, a.out, key=a.key, T0_C=a.t0)
    elif a.cmd == "prep":
        build_dltkcat_input(a.model, a.out, a.tmin, a.tmax, a.n, a.t0)
    elif a.cmd == "parse":
        fits = parse_dltkcat_output(a.pred, kcat_col=a.kcat_col, key=a.key,
                                    temp_col=a.temp_col, T0_C=a.t0,
                                    log10=not a.no_log10)
        fits.to_csv(a.out, index=False)
        print(f"[dltkcat] fitted {int(fits['ok'].sum())}/{len(fits)} usable -> {a.out}")


if __name__ == "__main__":
    main()
