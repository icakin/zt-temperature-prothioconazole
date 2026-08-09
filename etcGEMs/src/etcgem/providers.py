"""Model providers: turn a genome-scale model into (cobra model, cost table).

Three routes, all returning a ``ProvidedModel``:

* ``toy_ecoli_core``  -- self-contained, no downloads. e_coli_core (ships with
  cobrapy) plus synthetic MW / kcat / Topt / dCp per reaction. Use it to smoke-
  test the whole pipeline and to develop new analyses offline.

* ``from_gecko``      -- read a real GECKO enzyme-constrained model (SBML/.mat/
  .json) and extract per-reaction kcats from its protein-usage stoichiometry.
  Assigns MMRT knobs (Topt, dCp) from defaults you can override per enzyme.

* ``from_kcat_csv``   -- attach a table of kcats (e.g. DLKcat / DLTKcat output)
  to any base GEM. Columns: rxn_id, mw, kcat[, Topt, dCp, group, T0].

The GECKO and CSV routes are the ones you point at ecYeastGEM on your machine;
the toy route is what the sandbox tests run on.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .enzyme_cost import EnzymeConstrainedModel, EnzymeCostTable, EnzymeEntry


@dataclass
class ProvidedModel:
    ec: EnzymeConstrainedModel
    T0: float                     # reference temperature (K) of the kcats
    biomass_rxn: str
    name: str
    closed_free_o2_sinks: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Model-quality correction: uncosted free-energy (O2) sink reactions
# ---------------------------------------------------------------------------
# Base reaction IDs (GECKO adds No<n> isozyme and _REV suffixes) of O2-consuming
# side reactions that carry ~zero enzyme cost and short-circuit the electron
# transport chain, creating a futile O2 cycle (Parsa's gas-flux audit: at LB/40 C,
# ~75% of gross O2 turnover was this cycle). Closing them routes O2 through the
# genuine, enzyme-costed cytochrome oxidase (CYTBO3) at ~-0.3% growth. NB: catalase
# (CATNo1) and superoxide dismutase (SPODMNo1) are REAL enzymes (kcat 1e5-1e9/s) and
# are deliberately NOT in this list. Extend per-organism as other models are audited.
FREE_O2_SINK_REACTIONS = [
    "QMO2",    # 2 O2 + ubiquinol -> 2 superoxide + ubiquinone      (ygiN / P0ADU2; ~free)
    "QMO3",    # 2 O2 + menaquinol -> 2 superoxide + menaquinone    (same enzyme; ~free)
    "MOX",     # malate + O2 -> H2O2 + oxaloacetate                 (no enzyme-cost entry)
    "CU1Opp",  # 4 Cu+ + O2 -> 4 Cu2+ + 2 H2O                       (no enzyme-cost entry)
]


def close_free_energy_sinks(model, bases: Optional[List[str]] = None) -> List[str]:
    """Close uncosted O2-sink reactions (Parsa's audit) so flux routes through the
    genuine enzyme-costed respiratory chain. Matches each base ID plus any GECKO
    ``No<n>`` isozyme / ``_REV`` suffix, and sets the forward (and reverse, if the
    reaction is reversible) bound to 0. Returns the list of closed reaction IDs.
    Reusable/auditable: pass ``bases`` to close a different set per organism. Does
    NOT touch catalase / SOD (real kcat)."""
    bases = list(bases) if bases is not None else FREE_O2_SINK_REACTIONS
    pat = re.compile(r"^(" + "|".join(re.escape(b) for b in bases) + r")(No\d+)?(_REV(No\d+)?)?$")
    closed = []
    for r in model.reactions:
        if pat.match(r.id):
            if r.upper_bound > 0.0:
                r.upper_bound = 0.0
            if r.lower_bound < 0.0:
                r.lower_bound = 0.0
            closed.append(r.id)
    return closed


# ---------------------------------------------------------------------------
# Budget calibration
# ---------------------------------------------------------------------------
def calibrate_budget(model, table: EnzymeCostTable, T0: float, biomass_rxn: str,
                     target_fraction: float = 0.6) -> float:
    """Pick a pool budget so constrained growth at T0 is ~target_fraction of the
    unconstrained optimum (guarantees the enzyme constraint actually binds).

    Bisection on the budget; cheap because each step is one LP.
    """
    with model:
        model.objective = biomass_rxn
        uncon = model.slim_optimize()
    if not uncon or uncon <= 0 or not np.isfinite(uncon):
        raise RuntimeError("Unconstrained model does not grow; check biomass reaction.")
    target = target_fraction * uncon

    # An upper bound on useful budget: cost if every enzyme ran at unit flux.
    hi = sum(e.cost(T0) for e in table) or 1.0
    lo = 0.0
    # Calibrate on a disposable copy so we don't leave a pool constraint behind.
    ecm = EnzymeConstrainedModel(model.copy(), table, default_budget=hi)
    ecm.model.objective = biomass_rxn
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        ecm.set_budget(mid)
        ecm.set_temperature(T0)
        g = ecm.model.slim_optimize()
        g = 0.0 if (g is None or not np.isfinite(g)) else g
        if g < target:
            lo = mid
        else:
            hi = mid
    return hi


# ---------------------------------------------------------------------------
# Toy provider (offline)
# ---------------------------------------------------------------------------
def toy_ecoli_core(T0: float = 303.15, seed: int = 0,
                   topt_mean_offset: float = 8.0, topt_sd: float = 6.0,
                   dCp_mean: float = -8.0, dCp_sd: float = 2.0,
                   n_groups: int = 4, target_fraction: float = 0.6) -> ProvidedModel:
    """E. coli core with synthetic enzyme kinetics. Deterministic given seed.

    Enzyme optima are drawn a little above T0 with spread ``topt_sd`` so the
    organismal TPC (which folds in the proteome budget) peaks below the mean
    enzyme Topt -- the usual empirical pattern.
    """
    from cobra.io import load_model
    rng = np.random.default_rng(seed)
    model = load_model("textbook")
    biomass = str(model.objective.expression).split()[0]  # fallback below
    biomass = _find_biomass(model)

    entries = []
    for rxn in model.reactions:
        if rxn.id.startswith("EX_") or rxn.id == biomass:
            continue
        if rxn.id in ("ATPM",):
            continue
        if not rxn.metabolites:
            continue
        mw = float(np.clip(rng.lognormal(mean=np.log(40), sigma=0.5), 10, 150))     # kDa
        kcat = float(np.clip(rng.lognormal(mean=np.log(50), sigma=1.0), 1, 800))    # 1/s
        Topt = T0 + topt_mean_offset + rng.normal(0, topt_sd)
        dCp = float(np.clip(rng.normal(dCp_mean, dCp_sd), -20, -2))
        group = f"grp{rng.integers(0, n_groups)}"
        entries.append(EnzymeEntry(rxn.id, mw, kcat, Topt, dCp, T0, group))
    table = EnzymeCostTable(entries)

    budget = calibrate_budget(model, table, T0, biomass, target_fraction)
    group_budgets = _group_budgets_from_reference(table, T0, budget)
    ec = EnzymeConstrainedModel(model, table, default_budget=budget,
                                group_budgets=group_budgets)
    ec.model.objective = biomass
    return ProvidedModel(ec=ec, T0=T0, biomass_rxn=biomass, name="toy_ecoli_core")


def _group_budgets_from_reference(table, T0, total_budget, slack=1.5):
    """Per-group sub-budgets, set generous (slack x each group's proportional
    share) so group caps only bind when a sweep tightens allocation."""
    per_group = {}
    for g, ents in table.by_group().items():
        share = sum(e.cost(T0) for e in ents)
        per_group[g] = share
    tot = sum(per_group.values()) or 1.0
    return {g: slack * total_budget * v / tot for g, v in per_group.items()}


def _find_biomass(model) -> str:
    # objective reaction name is most reliable
    for r in model.reactions:
        if r.objective_coefficient != 0:
            return r.id
    for r in model.reactions:
        if "biomass" in r.id.lower():
            return r.id
    raise RuntimeError("Could not identify biomass reaction.")


# ---------------------------------------------------------------------------
# GECKO ecModel extractor
# ---------------------------------------------------------------------------
def load_enzyme_thermal_params(path: str):
    """Read a per-enzyme thermal parameter table (MRes / Li 2021 format) indexed
    by UniProt id, with columns Topt, Tm (K), Length, T90, dCpt (J/mol/K).

    Returns a DataFrame indexed by UniProt id. The first (unnamed) column is the
    UniProt id in the MRes tables."""
    import pandas as pd
    df = pd.read_csv(path)
    idcol = df.columns[0]
    df = df.set_index(df[idcol].astype(str))
    return df


def apply_thermal_params(table, params_df, key: str = "enzyme_id",
                         use_dCpt: bool = True):
    """Set each entry's grounded Topt/Tm/T90/length/dCpt from a params table,
    joined by UniProt id. Topt/Tm are converted to K assumption already met (the
    MRes tables store K). ``use_dCpt=False`` leaves dCpt unset (the emergent model
    uses an independent literature MMRT dCp prior instead of the table's tuned
    per-enzyme values). Returns (n_matched, n_total) coverage."""
    import numpy as np
    n_matched = 0
    for e in table.entries:
        uid = getattr(e, key, None)
        if uid is None or uid not in params_df.index:
            continue
        row = params_df.loc[uid]
        if hasattr(row, "iloc") and getattr(row, "ndim", 1) > 1:
            row = row.iloc[0]   # duplicate ids -> take first
        def _get(col):
            v = row.get(col) if hasattr(row, "get") else None
            return float(v) if v is not None and np.isfinite(v) else None
        topt = _get("Topt")
        if topt is not None:
            e.Topt = topt
        e.Tm = _get("Tm")
        e.T90 = _get("T90")
        e.length = _get("Length")
        if use_dCpt:
            e.dCpt = _get("dCpt")
        n_matched += 1
    return n_matched, len(table.entries)


def from_gecko(model_path: str, T0: float = 303.15,
               default_Topt_offset: float = 7.0, default_dCp: float = -12.0,
               prot_prefix: str = "prot_", pool_id: str = "prot_pool",
               biomass_rxn: Optional[str] = None,
               target_fraction: Optional[float] = None,
               pool_scale: float = 1.0,
               thermal_model: str = "mmrt", ngam_temperature: bool = False,
               ngam_rxn: Optional[str] = None,
               enzyme_params: Optional[str] = None,
               budget_override: Optional[float] = None,
               enzyme_params_use_dCpt: bool = True,
               dcp_prior_kJ: float = -4.0,
               close_free_o2_sinks: bool = True) -> ProvidedModel:
    """Extract an enzyme cost table from a GECKO-style ecModel.

    Handles two encodings:
      (A) 'full' GECKO: draw reactions convert `prot_pool` -> `prot_<id>` with
          coefficient = MW; metabolic reactions consume `prot_<id>` with
          coefficient = 1/kcat. We recover MW and kcat per reaction.
      (B) 'short'/sMOMENT: reactions consume `prot_pool` directly with
          coefficient MW/kcat. We keep that as base_cost (temperature scaling is
          invariant to the MW/kcat split).

    Topt/dCp are not stored in GECKO models, so we assign defaults; override them
    afterwards via the returned table (e.g. from DLTKcat predictions).
    """
    import cobra

    model = _load_any(model_path)
    # Model-quality correction (default ON): close the uncosted O2-sink reactions
    # (Parsa's audit) so O2 routes through the real respiratory chain, not a futile cycle.
    closed_sinks = []
    if close_free_o2_sinks:
        closed_sinks = close_free_energy_sinks(model)
        if closed_sinks:
            print(f"[from_gecko] closed {len(closed_sinks)} uncosted O2-sink reaction(s) "
                  f"(Parsa's audit): {closed_sinks}")
    if biomass_rxn is None:
        biomass_rxn = _find_biomass(model)

    # Auto-detect the pool metabolite id (GECKO adds a compartment suffix, e.g.
    # 'prot_pool[c]', so the bare 'prot_pool' won't match get_by_id).
    all_met_ids = {m.id for m in model.metabolites}
    if pool_id not in all_met_ids:
        cand = [i for i in all_met_ids if "prot_pool" in i or i.startswith(pool_id)]
        if cand:
            pool_id = sorted(cand, key=len)[0]
    # Enzyme pseudo-metabolites = prot_ prefix, excluding the pool itself.
    prot_mets = {m.id for m in model.metabolites
                 if m.id.startswith(prot_prefix) and m.id != pool_id}

    # MW per enzyme from draw reactions (route A: pool -> single enzyme).
    mw_of: Dict[str, float] = {}
    draw_rxns = set()
    for rxn in model.reactions:
        prods = [m for m in rxn.metabolites if m.id in prot_mets and rxn.metabolites[m] > 0]
        cons_pool = any(m.id == pool_id and rxn.metabolites[m] < 0 for m in rxn.metabolites)
        if cons_pool and len(prods) == 1:
            draw_rxns.add(rxn.id)
            enz = prods[0].id
            for m, c in rxn.metabolites.items():
                if m.id == pool_id:
                    mw_of[enz] = abs(c)   # |pool coeff| == MW (kDa)

    pool_met = model.metabolites.get_by_id(pool_id) if pool_id in all_met_ids else None
    entries = []
    for rxn in model.reactions:
        if rxn.id == biomass_rxn or rxn.id.startswith("EX_") or rxn.id in draw_rxns:
            continue
        enz_consumed = {m.id: -rxn.metabolites[m] for m in rxn.metabolites
                        if m.id in prot_mets and rxn.metabolites[m] < 0}
        pool_consumed = rxn.metabolites.get(pool_met, 0.0) if pool_met is not None else 0.0

        if enz_consumed:  # route A: coefficient == 1/kcat[1/h]
            # if a complex, use the limiting (largest coefficient == smallest kcat)
            enz_id, inv_kcat_h = max(enz_consumed.items(), key=lambda kv: kv[1])
            kcat_s = 1.0 / (inv_kcat_h * 3600.0) if inv_kcat_h > 0 else 50.0
            mw = mw_of.get(enz_id, 40.0)
            entries.append(EnzymeEntry(rxn.id, mw, kcat_s,
                                       T0 + default_Topt_offset, default_dCp, T0,
                                       group=_subsystem(rxn),
                                       enzyme_id=_clean_uniprot(enz_id, prot_prefix)))
        elif pool_consumed < 0:  # route B: base_cost = MW/kcat directly
            entries.append(EnzymeEntry(rxn.id, mw=40.0, kcat_ref=50.0,
                                       Topt=T0 + default_Topt_offset, dCp=default_dCp,
                                       T0=T0, group=_subsystem(rxn),
                                       base_cost=abs(pool_consumed)))
    if not entries:
        raise RuntimeError(
            "No enzyme-linked reactions found. Check prot_prefix/pool_id match "
            "this model's naming (inspect metabolite ids).")

    table = EnzymeCostTable(entries)
    # Budget: reuse the model's existing pool bound if present, else calibrate.
    if budget_override is not None:
        # EMERGENT magnitude: the pool budget is derived from independent data
        # (P_total x f_metab x sigma), NOT the growth-calibrated GECKO bound.
        budget = float(budget_override)
    else:
        budget = _existing_pool_bound(model, pool_id)
        if budget is None:
            if target_fraction is None:
                target_fraction = 0.6
            budget = calibrate_budget(model, table, T0, biomass_rxn, target_fraction)
        # DEPRECATED knob: pool_scale<1 tuned the pool to make peak growth respond.
        # For the emergent model pool_scale=1.0 (nothing tuned to the growth curve).
        budget *= pool_scale
    group_budgets = _group_budgets_from_reference(table, T0, budget)

    # Unfolding mode: overlay grounded per-enzyme Topt/Tm/dCpt/length before the
    # model precomputes its unfolding thermodynamics; report coverage.
    if thermal_model == "unfolding" and enzyme_params:
        params_df = load_enzyme_thermal_params(enzyme_params)
        n_match, n_tot = apply_thermal_params(table, params_df, key="enzyme_id",
                                              use_dCpt=enzyme_params_use_dCpt)
        print(f"[from_gecko] unfolding mode: matched grounded Topt/Tm for "
              f"{n_match}/{n_tot} enzymes ({100*n_match/max(1,n_tot):.1f}%); "
              f"the rest use dataset-mean fallbacks")

    ec = EnzymeConstrainedModel(model, table, default_budget=budget,
                                group_budgets=group_budgets,
                                thermal_model=thermal_model,
                                ngam_temperature=ngam_temperature, ngam_rxn=ngam_rxn,
                                unfold_means={"dCpt": dcp_prior_kJ * 1000.0})
    ec.model.objective = biomass_rxn
    # P1c: reconcile the two redundant enzyme-mass pools into ONE proteome budget.
    # The etcgem sMOMENT/sector pool (P_total x f_metab x sigma, medium- & growth-law-
    # aware, scaled by kcat_scale) is the sole proteome accounting; relax the leftover
    # GECKO base pool supply (prot_pool_exchange) so it can no longer cap growth
    # independently (the magnitude levers kcat_scale/sigma/P_total do not touch it).
    _relax_base_protein_pool(ec.model, pool_id)
    return ProvidedModel(ec=ec, T0=T0, biomass_rxn=biomass_rxn,
                         name=f"gecko:{model.id}", closed_free_o2_sinks=closed_sinks)


# A GEM-compatible LB (rich) medium: amino acids, nucleosides/bases, vitamins,
# ions and glucose. Component EX_ ids from the MResProject data/media.csv definition
# (reused; consistent with the Machado et al. 2018 LB composition).
_CARBON_BASES = ["glc__D", "ac", "glyc", "succ", "lac__D", "lac__L", "xyl__D",
                 "gal", "fru", "man", "malt", "sucr", "pyr", "etoh", "cit"]
LB_COMPONENTS = [
    "adn", "ala__L", "amp", "arg__L", "aso3", "asp__L", "ca2", "cbl1", "cd2", "cl",
    "cmp", "cobalt2", "cro4", "cu2", "cys__L", "dad_2", "dcyt", "fe2", "fe3", "fol",
    "glc__D", "glu__L", "gly", "gmp", "gsn", "h2o", "h2s", "h", "hg2", "his__L",
    "hxan", "ile__L", "ins", "k", "leu__L", "lipoate", "lys__L", "met__L", "mg2",
    "mn2", "mobd", "na1", "nac", "nh4", "ni2", "o2", "phe__L", "pheme", "pi",
    "pnto__R", "pro__L", "pydx", "ribflv", "ser__L", "so4", "thm", "thr__L",
    "thymd", "trp__L", "tyr__L", "ump", "ura", "uri", "val__L", "zn2"]


def set_medium(pm, medium="glucose_minimal", carbon="glc__D", aerobic=True,
               uptake_ub=1000.0, lb_media_csv=None, bhi_media_csv=None):
    """Set the growth medium as AVAILABILITY, not pinned uptake rates: open the
    `EX_<met>_e_REV` uptakes for the medium's components and close other carbon
    sources; the enzyme-constrained model then determines actual uptake.

    medium="glucose_minimal" (default): a single carbon source (+ O2 if aerobic)
    on top of the model's minimal-salt defaults. medium="LB"/"BHI": a rich medium
    opening the amino-acid / nucleoside / vitamin / ion component uptakes
    (:data:`LB_COMPONENTS`, or the given ``lb_media_csv`` / ``bhi_media_csv``).
    BHI is not encoded in silico in the literature (even the Rothia GEM used
    LB/TSB as defined rich stand-ins), so it is approximated by a curated rich
    component list and mapped to the LB medium-matched sector allocation.
    Returns (n_opened, n_missing)."""
    model = pm.ec.model if hasattr(pm, "ec") else pm
    # switch the medium-matched sector allocation, if wired. BHI has no measured
    # proteome (DeyuWang is LB/Glucose/Glycerol), so it uses the LB (rich) curve.
    sector_medium = "LB" if medium == "BHI" else medium
    alloc = getattr(getattr(pm, "ec", None), "_alloc_from_data", None)
    if alloc is not None and hasattr(alloc, "set_active_medium"):
        alloc.set_active_medium(sector_medium)
    if medium in ("LB", "BHI"):
        comps = LB_COMPONENTS
        csv = bhi_media_csv if medium == "BHI" else lb_media_csv
        if csv:
            import pandas as pd
            df = pd.read_csv(csv)
            comps = [str(n)[3:-2] for n in df["Name"] if str(n).startswith("EX_")]
        # close non-medium carbon sources, then open every component uptake present
        rich_set = set(comps)
        for base in _CARBON_BASES:
            rev = f"EX_{base}_e_REV"
            if rev in model.reactions and base not in rich_set:
                model.reactions.get_by_id(rev).upper_bound = 0.0
        opened, missing = 0, 0
        for base in comps:
            rev = f"EX_{base}_e_REV"
            if rev in model.reactions:
                model.reactions.get_by_id(rev).upper_bound = uptake_ub
                opened += 1
            else:
                missing += 1
        model.solver.update()
        return opened, missing
    # glucose_minimal (single carbon source)
    for base in _CARBON_BASES:
        rev = f"EX_{base}_e_REV"
        if rev in model.reactions:
            model.reactions.get_by_id(rev).upper_bound = uptake_ub if base == carbon else 0.0
    o2 = "EX_o2_e_REV"
    if o2 in model.reactions:
        model.reactions.get_by_id(o2).upper_bound = uptake_ub if aerobic else 0.0
    model.solver.update()
    return 1 + int(aerobic), 0


def _existing_pool_bound(model, pool_id):
    for rxn in model.reactions:
        # the pool exchange/supply reaction produces prot_pool
        prod = [m for m in rxn.metabolites if m.id == pool_id and rxn.metabolites[m] > 0]
        if prod and rxn.upper_bound < 1e6:
            return float(rxn.upper_bound)
    return None


def _relax_base_protein_pool(model, pool_id, big=1e6):
    """Relax the GECKO base protein-pool supply reaction (which produces the pool
    metabolite with a finite upper bound) so it no longer independently caps growth.
    The etcgem sMOMENT/sector pool becomes the sole enzyme-mass constraint. Returns
    the original bound that was relaxed (or None)."""
    for rxn in model.reactions:
        prod = [m for m in rxn.metabolites if m.id == pool_id and rxn.metabolites[m] > 0]
        if prod and rxn.upper_bound < 1e6:
            old = float(rxn.upper_bound)
            rxn.upper_bound = float(big)
            model.solver.update()
            print(f"[from_gecko] reconciled proteome pools: relaxed base pool supply "
                  f"'{rxn.id}' ({old:.4g} -> inf); the etcgem sMOMENT/sector pool is now "
                  f"the sole enzyme-mass budget.")
            return old
    return None


def _subsystem(rxn) -> str:
    s = getattr(rxn, "subsystem", None)
    return s if s else "default"


def _clean_uniprot(prot_met_id: str, prot_prefix: str) -> str:
    """'prot_P0A8F4[c]' -> 'P0A8F4'."""
    s = prot_met_id
    if s.startswith(prot_prefix):
        s = s[len(prot_prefix):]
    return s.split("[")[0].strip()


# ---------------------------------------------------------------------------
# CSV kcat loader (DLKcat / DLTKcat)
# ---------------------------------------------------------------------------
def from_kcat_csv(model_path: str, csv_path: str, T0: float = 303.15,
                  default_Topt_offset: float = 8.0, default_dCp: float = -8.0,
                  biomass_rxn: Optional[str] = None,
                  target_fraction: float = 0.6) -> ProvidedModel:
    """Attach a kcat table to a plain GEM.

    CSV columns (header required): rxn_id, mw, kcat  and optionally
    Topt, dCp, group, T0. Missing optional columns fall back to defaults.
    """
    model = _load_any(model_path)
    if biomass_rxn is None:
        biomass_rxn = _find_biomass(model)
    entries = []
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            rid = row["rxn_id"].strip()
            if rid not in model.reactions:
                continue
            t0 = float(row.get("T0") or T0)
            entries.append(EnzymeEntry(
                rxn_id=rid,
                mw=float(row.get("mw") or 40.0),
                kcat_ref=float(row["kcat"]),
                Topt=float(row.get("Topt") or (t0 + default_Topt_offset)),
                dCp=float(row.get("dCp") or default_dCp),
                T0=t0,
                group=(row.get("group") or "default").strip() or "default",
            ))
    if not entries:
        raise RuntimeError("No CSV rows matched reactions in the model.")
    table = EnzymeCostTable(entries)
    budget = calibrate_budget(model, table, T0, biomass_rxn, target_fraction)
    group_budgets = _group_budgets_from_reference(table, T0, budget)
    ec = EnzymeConstrainedModel(model, table, default_budget=budget,
                                group_budgets=group_budgets)
    ec.model.objective = biomass_rxn
    return ProvidedModel(ec=ec, T0=T0, biomass_rxn=biomass_rxn,
                         name=f"csv:{model.id}")


def _read_sbml_safe(path: str):
    """Read an SBML model, sanitising COBRA id-encodings that break the LP
    backend. GECKO SBML encodes spaces as ``__32__`` which cobra decodes to a
    literal space in the id (e.g. 'protein pseudoreaction'); optlang/GLPK reject
    whitespace in variable names. We decode normally then replace spaces."""
    import cobra
    from cobra.io.sbml import F_REPLACE, F_REACTION, F_SPECIE, F_GENE
    base = dict(F_REPLACE)
    _wrap = lambda f: (lambda s: f(s).replace(" ", "_"))
    fr = dict(base)
    for key in (F_REACTION, F_SPECIE, F_GENE):
        if key in fr:
            fr[key] = _wrap(base[key])
    return cobra.io.read_sbml_model(path, f_replace=fr)


def _load_any(path: str):
    import cobra
    p = path.lower()
    if p.endswith((".xml", ".sbml", ".xml.gz")):
        return _read_sbml_safe(path)
    if p.endswith(".json"):
        return cobra.io.load_json_model(path)
    if p.endswith(".mat"):
        return cobra.io.load_matlab_model(path)
    if p.endswith(".yml") or p.endswith(".yaml"):
        return cobra.io.load_yaml_model(path)
    raise ValueError(f"Unsupported model format: {path}")
