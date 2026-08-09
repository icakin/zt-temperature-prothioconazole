"""Enzyme-constraint layer: a temperature-dependent proteome pool budget.

We represent enzyme demand abstractly as a *cost table* -- one entry per
catalysed reaction giving molecular weight, a reference kcat and its MMRT
temperature-response knobs (Topt, dCp) plus an allocation group. This keeps the
engine organism-agnostic: the same code drives a toy E. coli core model, a
GECKO ecYeastGEM, or a table of DLKcat/DLTKcat predictions.

The constraint is the sMOMENT total-protein pool:

    sum_i  c_i(T) * (v_i_fwd + v_i_rev)   <=  P_budget

with per-flux protein cost

    c_i(T) = MW_i / ( kcat_i(T) * 3600 )        [ g protein / (mmol gDW^-1 h^-1) ]

Units: MW in kDa == g/mmol, kcat in 1/s (x3600 -> 1/h), flux in mmol/gDW/h, so
c_i*v_i is g enzyme / gDW and P_budget is the g protein / gDW allocated to the
modelled enzymes (the classic phi * sigma * f mass fraction). Optional per-group
sub-budgets let you probe proteome *allocation* between pathways, not just the
total.

Temperature enters only through kcat_i(T); raising T past an enzyme's Topt
lowers its kcat, inflating c_i, tightening the budget and eventually starving
growth -- the mechanistic origin of the organismal thermal performance curve.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional

import numpy as np

from .mmrt import MMRTParams

KCAT_S_TO_H = 3600.0


@dataclass
class EnzymeEntry:
    rxn_id: str
    mw: float          # kDa (== g/mmol)
    kcat_ref: float    # 1/s at T0
    Topt: float        # K
    dCp: float         # kJ/mol/K (negative)
    T0: float          # K, reference temperature of kcat_ref
    group: str = "default"
    enzyme_id: Optional[str] = None     # UniProt id (for keying DLTKcat fits)
    base_cost: Optional[float] = None   # override g/(mmol gDW^-1 h^-1); for GECKO short models
    # Two-state unfolding parameters (only used when thermal_model="unfolding";
    # None -> fall back to dataset means). Tm/T90 in K, length in residues, dCpt
    # (transition-state heat capacity) in J/mol/K -- see unfolding.py.
    Tm: Optional[float] = None
    T90: Optional[float] = None
    length: Optional[float] = None
    dCpt: Optional[float] = None

    def __post_init__(self):
        # Per-flux protein cost at T0. relative_kcat(T) is invariant to kcat_ref,
        # so base_cost fully determines the temperature-scaled cost.
        if self.base_cost is None:
            self.base_cost = self.mw / (self.kcat_ref * KCAT_S_TO_H)

    def mmrt(self, Topt: Optional[float] = None, dCp: Optional[float] = None) -> MMRTParams:
        return MMRTParams(
            kcat_ref=self.kcat_ref,
            T0=self.T0,
            Topt=self.Topt if Topt is None else Topt,
            dCp=self.dCp if dCp is None else dCp,
        )

    def cost(self, T: float, Topt: Optional[float] = None, dCp: Optional[float] = None) -> float:
        """Per-flux protein cost c_i(T) = base_cost / relative_kcat(T)."""
        rel = float(self.mmrt(Topt, dCp).relative_kcat(T))
        rel = max(rel, 1e-6)  # guard against numerical zero at extreme T
        return self.base_cost / rel


@dataclass
class EnzymeCostTable:
    """Container of EnzymeEntry with convenience accessors."""
    entries: List[EnzymeEntry] = field(default_factory=list)

    def __len__(self):
        return len(self.entries)

    def __iter__(self):
        return iter(self.entries)

    @property
    def groups(self) -> List[str]:
        return sorted({e.group for e in self.entries})

    def by_group(self) -> Dict[str, List[EnzymeEntry]]:
        out: Dict[str, List[EnzymeEntry]] = {}
        for e in self.entries:
            out.setdefault(e.group, []).append(e)
        return out


# ---------------------------------------------------------------------------
# Perturbation spec applied to the whole table during a sweep.
# ---------------------------------------------------------------------------
@dataclass
class Perturbation:
    """A single point in sensitivity space.

    Attributes
    ----------
    dTopt        : uniform shift (K) added to every enzyme's Topt
    topt_scale   : multiplies each enzyme's (Topt - T0) spread about T0 -> widens
                   or compresses the *heterogeneity* of optima across enzymes
    dCp_scale    : multiplies every enzyme's dCp (curvature / thermal breadth)
    dTm          : uniform shift (K) added to every enzyme's melting temperature Tm
                   (unfolding mode only) -> moves the denaturation collapse / CTmax
    kappa_scale  : multiplier on the effective in-vivo translation efficiency in the
                   biosynthesis cap: (translation_coeff / kappa_scale) * v_bio <=
                   f_bio * P_total, so kappa_scale > 1 relaxes the translation cap and
                   raises the achievable r_max (sectors mode only). The honest
                   magnitude lever -- a genuinely uncertain/borrowed quantity.
    kcat_scale   : global multiplier on metabolic enzyme turnover kcat -- scales the
                   *level* of kcat(T) (equivalently divides every per-flux cost), so a
                   value > 1 lowers enzyme demand and raises the pool-limited flux
                   ceiling / r_max. It scales the level, NOT the shape, so it leaves
                   E_a / T_opt / CT_max unchanged. The in-vitro->in-vivo kcat gap is the
                   physical prior; a single global scalar (per-enzyme kcat corrections
                   are unidentifiable from one growth curve).
    tm_scale     : multiplies each enzyme's (Tm - mean_Tm) about the distribution mean
                   (mirrors topt_scale for optima) -> controls how synchronised the
                   unfolding collapse is (narrow: sharp cliff; broad: gradual shoulder).
    ngam_scale     : multiplier on the maintenance (NGAM) amplitude.
    ngam_steepness : multiplier on how fast maintenance rises with temperature (the
                     peak-rounding lever). Both no-ops at 1.
    sigma_sat    : effective average in-vivo enzyme saturation sigma (the GECKO/
                   sMOMENT saturation coefficient in the total-protein budget
                   P_metab = sigma * f_metab * P_total). None -> use the model's
                   nominal sigma (identity). Because sigma multiplies the whole
                   enzyme-mass budget it raises BOTH the metabolic pool AND the
                   biosynthesis/ribosome cap (unlike kcat_scale, which only relaxes
                   the metabolic per-flux cost), so it is the effective magnitude
                   lever once the metabolic pool is cheap and the ribosome budget
                   binds. Applied as the multiplier sigma_sat / sigma_nom on both
                   sector caps (sectors mode only).
    budget       : total proteome pool P (g/gDW); None keeps model default
    group_alloc  : per-group multiplier on that group's sub-budget (allocation)
    """
    dTopt: float = 0.0
    topt_scale: float = 1.0
    dCp_scale: float = 1.0
    dTm: float = 0.0
    tm_scale: float = 1.0
    kappa_scale: float = 1.0
    kcat_scale: float = 1.0
    ngam_scale: float = 1.0
    ngam_steepness: float = 1.0
    budget: Optional[float] = None
    group_alloc: Dict[str, float] = field(default_factory=dict)
    # Proteome-sector allocation (opt-in; None -> use the scalar pool / set_budget
    # path, i.e. exactly the pre-sector behaviour). Only meaningful when the model
    # has sectors wired (see sectors.add_proteome_sectors).
    f_metab: Optional[float] = None
    f_maint: Optional[float] = None
    maint_to_bio: Optional[float] = None   # shift maintenance->biosynthesis at fixed f_metab
    sigma_sat: Optional[float] = None      # effective in-vivo saturation (scales both caps)

    def uses_allocation(self) -> bool:
        return (self.f_metab is not None or self.f_maint is not None
                or self.maint_to_bio is not None or self.sigma_sat is not None)

    def topt_of(self, e: EnzymeEntry) -> float:
        return e.T0 + self.topt_scale * (e.Topt - e.T0) + self.dTopt

    def dCp_of(self, e: EnzymeEntry) -> float:
        return e.dCp * self.dCp_scale


class EnzymeConstrainedModel:
    """Attach and drive a temperature-dependent pool constraint on a cobra model."""

    POOL = "enzyme_pool"

    def __init__(self, model, table: EnzymeCostTable, default_budget: float,
                 group_budgets: Optional[Dict[str, float]] = None,
                 thermal_model: str = "mmrt", ngam_temperature: bool = False,
                 ngam_rxn: Optional[str] = None,
                 unfold_means: Optional[Dict[str, float]] = None):
        self.model = model
        self.table = EnzymeCostTable(
            [e for e in table.entries if e.rxn_id in model.reactions]
        )
        dropped = len(table) - len(self.table)
        if dropped:
            print(f"[enzyme_cost] {dropped} table entries had no matching reaction; skipped.")
        self.default_budget = default_budget
        self.group_budgets = dict(group_budgets or {})
        self._pool = None
        self._group_cons: Dict[str, object] = {}
        self._sectors = None   # populated by sectors.add_proteome_sectors (opt-in)
        self._alloc_from_data = None   # proteome_alloc.TemperatureAllocation (opt-in)
        # Thermal model: "mmrt" (default, peak-normalised MMRT -- unchanged) or
        # "unfolding" (two-state denaturation keyed on per-enzyme Tm, after
        # Li 2021 / the MRes; see unfolding.py). ngam_temperature adds the
        # temperature-dependent maintenance term (unfolding only).
        self.thermal_model = thermal_model
        self.ngam_temperature = bool(ngam_temperature)
        um = unfold_means or {}
        self._unfold_mean_Tm = float(um.get("Tm", 273.15 + 55.6))   # thesis mean 55.6 C
        self._unfold_mean_len = float(um.get("length", 300.0))
        self._unfold_mean_dCpt = float(um.get("dCpt", -4000.0))     # J/mol/K
        # Precompute parameter arrays and variable handles for a vectorised
        # temperature update (the hot path in a sweep).
        ents = self.table.entries
        self._base = np.array([e.base_cost for e in ents], float)
        self._T0 = np.array([e.T0 for e in ents], float)
        self._Topt = np.array([e.Topt for e in ents], float)
        self._dCp = np.array([e.dCp for e in ents], float)
        rxns = [model.reactions.get_by_id(e.rxn_id) for e in ents]
        self._fwd = [r.forward_variable for r in rxns]
        self._rev = [r.reverse_variable for r in rxns]
        self._group_pos: Dict[str, np.ndarray] = {}
        for g in {e.group for e in ents}:
            self._group_pos[g] = np.array([i for i, e in enumerate(ents) if e.group == g])
        self._ngam_rxn = self._detect_ngam(ngam_rxn) if self.ngam_temperature else None
        if self.thermal_model == "unfolding":
            self._build_unfolding()
        self._build()

    def _detect_ngam(self, ngam_rxn):
        cands = [ngam_rxn] if ngam_rxn else ["ATPM", "NGAM"]
        for c in cands:
            if c and c in self.model.reactions:
                return self.model.reactions.get_by_id(c)
        return None

    def _build_unfolding(self):
        """Precompute per-enzyme two-state unfolding thermodynamics (dHTH, dSTS,
        dCpu in J) from Tm/T90/length, plus the transition-state dCpt, applying
        dataset-mean fallbacks for enzymes without measured values."""
        from . import unfolding as U
        ents = self.table.entries
        hth, sts, cpu, cpt, tms = [], [], [], [], []
        for e in ents:
            Tm = e.Tm if e.Tm is not None and np.isfinite(e.Tm) else self._unfold_mean_Tm
            length = e.length if e.length is not None and np.isfinite(e.length) else self._unfold_mean_len
            dHTH, dSTS, dCpu = U.unfolding_params(Tm, e.T90, length)
            hth.append(dHTH); sts.append(dSTS); cpu.append(dCpu); tms.append(Tm)
            cpt.append(e.dCpt if e.dCpt is not None and np.isfinite(e.dCpt) else self._unfold_mean_dCpt)
        self._uHTH = np.array(hth, float)
        self._uSTS = np.array(sts, float)
        self._uCpu = np.array(cpu, float)
        self._uCpt = np.array(cpt, float)
        self._Tm = np.array(tms, float)   # per-enzyme Tm (K), for the dTm reference scale

    # -- construction -------------------------------------------------------
    def _build(self):
        m = self.model
        self._pool = m.problem.Constraint(0, lb=0, ub=self.default_budget, name=self.POOL)
        cons = [self._pool]
        for g, b in self.group_budgets.items():
            c = m.problem.Constraint(0, lb=0, ub=b, name=f"enzyme_pool_{g}")
            self._group_cons[g] = c
            cons.append(c)
        m.add_cons_vars(cons)
        m.solver.update()
        self.set_temperature(self._ref_T())  # initialise coefficients

    def _ref_T(self) -> float:
        return float(np.median([e.T0 for e in self.table])) if len(self.table) else 303.15

    def refresh_params(self):
        """Rebuild the vectorised parameter arrays from the table entries.

        Call after mutating entry Topt/dCp/base_cost in place (e.g. applying
        DLTKcat fits), then set_temperature to push new coefficients."""
        ents = self.table.entries
        self._base = np.array([e.base_cost for e in ents], float)
        self._T0 = np.array([e.T0 for e in ents], float)
        self._Topt = np.array([e.Topt for e in ents], float)
        self._dCp = np.array([e.dCp for e in ents], float)
        if self.thermal_model == "unfolding":
            self._build_unfolding()
        self.set_temperature(self._ref_T())

    # -- drivers ------------------------------------------------------------
    def _costs(self, T: float, pert: Perturbation) -> np.ndarray:
        """Vectorised per-flux cost array c_i(T) over all enzymes.

        The MMRT shape is normalised to its own peak, so the reference kcat is
        treated as each enzyme's *maximum* turnover (at Topt) and any temperature
        deviation only raises cost. This keeps the proteome pool binding across
        the whole curve (no super-efficiency plateau) while preserving T0 as the
        pivot for the topt_scale perturbation.
        """
        if self.thermal_model == "unfolding":
            return self._costs_unfolding(T, pert)
        from .mmrt import relative_kcat_vec
        Topt_eff = self._T0 + pert.topt_scale * (self._Topt - self._T0) + pert.dTopt
        dCp_eff = self._dCp * pert.dCp_scale
        rel = relative_kcat_vec(T, self._T0, Topt_eff, dCp_eff)
        peak = relative_kcat_vec(Topt_eff, self._T0, Topt_eff, dCp_eff)  # value at each Topt
        s = np.nan_to_num(rel / peak, nan=1e-6, posinf=1e6, neginf=1e-6)
        s = np.clip(s, 1e-6, 1e6)             # normalised shape, <=1 at Topt; finite guard
        return self._base / (s * pert.kcat_scale)

    def _costs_unfolding(self, T: float, pert: Perturbation) -> np.ndarray:
        """Two-state-unfolding per-flux cost (Li 2021 / MRes):

            cost_i(T) = base_cost_i / ( rel_kcat_i(T) * f_N_i(T) )

        where rel_kcat_i is the transition-state turnover factor anchored at Topt
        (NOT peak-normalised) and f_N_i is the native fraction (-> 0 above Tm), so
        denaturation sets the falling limb. The envelope knobs still apply:
        dTopt/topt_scale shift the enzyme optima, dCp_scale scales the
        transition-state heat capacity, and dTm shifts the denaturation temperature
        Tm (a +dTm shift of every Tm is applied by evaluating the native fraction at
        T - dTm, following the reference etc.py ``Tadj``; so +dTm pushes CTmax up).
        base_cost is the cost at Topt, as in the reference model."""
        from . import unfolding as U
        Topt_eff = self._T0 + pert.topt_scale * (self._Topt - self._T0) + pert.dTopt
        dCpt_eff = self._uCpt * pert.dCp_scale
        rk = U.rel_kcat(T, self._uHTH, self._uSTS, self._uCpu, dCpt_eff, Topt_eff)
        # dTm shifts every Tm uniformly; tm_scale stretches/compresses the Tm spread
        # about its mean (narrow -> synchronised sharp collapse; broad -> gradual
        # shoulder). Both applied as a per-enzyme shift of the evaluation temperature.
        Tm_shift = pert.dTm + (pert.tm_scale - 1.0) * (self._Tm - np.mean(self._Tm))
        fN = U.native_fraction(T - Tm_shift, self._uHTH, self._uSTS, self._uCpu)
        # numerical guard: clamp the turnover*fold product to a finite, strictly
        # positive band so no enzyme becomes free (denom=inf -> cost=0) or NaN under
        # extreme perturbations, which would leave the LP degenerate and hang the
        # solver. A no-op in the physical regime (rk*fN in [1e-6, 1]).
        prod = np.nan_to_num(rk * fN, nan=1e-6, posinf=1e6, neginf=1e-6)
        denom = np.clip(prod, 1e-6, 1e6)
        return self._base / (denom * pert.kcat_scale)

    def set_temperature(self, T: float, pert: Optional[Perturbation] = None):
        """Recompute all pool coefficients for temperature T (K)."""
        pert = pert or Perturbation()
        c = self._costs(T, pert)
        pool_coef = {}
        for v_f, v_r, ci in zip(self._fwd, self._rev, c):
            ci = float(ci)
            pool_coef[v_f] = ci
            pool_coef[v_r] = ci
        self._pool.set_linear_coefficients(pool_coef)
        for g, con in self._group_cons.items():
            pos = self._group_pos[g]
            gc = {}
            for i in pos:
                ci = float(c[i])
                gc[self._fwd[i]] = ci
                gc[self._rev[i]] = ci
            con.set_linear_coefficients(gc)
        # Temperature-dependent maintenance (unfolding mode). Li/MRes call set_NGAMT
        # at EVERY temperature, so maintenance rises with T and rounds the peak --
        # apply it both without AND with the sector model wired.
        if (self.thermal_model == "unfolding" and self.ngam_temperature
                and self._ngam_rxn is not None):
            from . import unfolding as U
            if self._sectors is None:
                val = U.ngam_T(T, scale=pert.ngam_scale, steepness=pert.ngam_steepness)
                self._ngam_rxn.lower_bound = val
                self._ngam_rxn.upper_bound = max(val, self._ngam_rxn.upper_bound)
            else:
                # relative T-factor about the 25 C anchor (scale cancels in the ratio,
                # so it is applied separately as an amplitude), on top of the
                # sector-owned ATPM bound that set_allocation scales by f_maint.
                T0 = 273.15 + 25.0
                r = (U.ngam_T(T, steepness=pert.ngam_steepness)
                     / U.ngam_T(T0, steepness=pert.ngam_steepness))
                self._ngam_T_factor = float(r) * float(pert.ngam_scale)
                s = self._sectors
                atpm = s["atpm_rxn"]
                if atpm is not None and s["f_maint_nom"] > 0:
                    fmaint = getattr(self, "_last_fmaint", None) or s["f_maint_nom"]
                    val = s["atpm_nom_lb"] * (fmaint / s["f_maint_nom"]) * self._ngam_T_factor
                    atpm.lower_bound = val
                    atpm.upper_bound = max(val, atpm.upper_bound)
        if pert.budget is not None or pert.group_alloc:
            self.set_budget(pert.budget, pert.group_alloc)

    def set_budget(self, budget: Optional[float] = None,
                   group_alloc: Optional[Dict[str, float]] = None):
        if budget is not None:
            self._pool.ub = budget
        if group_alloc:
            for g, mult in group_alloc.items():
                if g in self._group_cons:
                    base = self.group_budgets[g]
                    self._group_cons[g].ub = base * mult

    def set_allocation(self, f_metab: Optional[float] = None,
                       f_maint: Optional[float] = None,
                       kappa_scale: float = 1.0,
                       sigma_sat: Optional[float] = None):
        """Set the three-sector proteome partition (Basan/Scott) in place.

        f_bio = 1 - f_metab - f_maint. Updates the metabolic pool bound
        (f_metab*P_total), the biosynthesis cap (f_bio*P_total) and, if present,
        the maintenance-ATP lower bound (scaled by f_maint / f_maint_nominal).
        ``kappa_scale`` (>0) multiplies the biosynthesis cap bound, equivalent to
        dividing the translation coefficient by kappa_scale -- the in-vivo
        translation-efficiency lever on r_max. ``sigma_sat`` (>0) overrides the
        effective in-vivo saturation: it multiplies BOTH sector caps by
        sigma_sat / sigma_nom (the whole enzyme-mass budget scales with the average
        saturation), so unlike kcat_scale it also relaxes the ribosome/biosynthesis
        cap. None -> nominal sigma (identity). Requires sectors wired via
        sectors.add_proteome_sectors."""
        s = self._sectors
        if s is None:
            raise RuntimeError("proteome sectors not enabled; call "
                               "sectors.add_proteome_sectors(pm, cfg) first")
        fm = s["f_metab_nom"] if f_metab is None else float(f_metab)
        fmaint = s["f_maint_nom"] if f_maint is None else float(f_maint)
        fbio = 1.0 - fm - fmaint
        eps = 1e-9
        if min(fm, fmaint, fbio) < -eps or fm + fmaint > 1.0 + eps:
            raise ValueError(f"invalid sector simplex: f_metab={fm}, "
                             f"f_maint={fmaint}, f_bio={fbio}")
        P = s["P_total"]
        # effective-saturation multiplier on the whole enzyme-mass budget (both caps)
        sig = 1.0 if sigma_sat is None else float(sigma_sat) / float(s.get("sigma_nom", 0.45))
        if s.get("growth_law"):
            # coupled growth law: f_maint drives the split; the biosynthesis intercept
            # f_bio_0 is fixed and f_metab_0 = 1 - f_maint - f_bio_0 (conserving). The
            # measured f_metab arg is ignored (the mu-dependent law sets it via the
            # v_bio coefficients wired at build). Both caps are LINEAR in mu.
            f_bio_0 = s["f_bio_0"]
            f_metab_0 = 1.0 - fmaint - f_bio_0
            self._pool.ub = f_metab_0 * P * sig
            s["bio_constraint"].ub = f_bio_0 * P * float(kappa_scale) * sig
        else:
            self._pool.ub = fm * P * sig
            s["bio_constraint"].ub = fbio * P * float(kappa_scale) * sig
        self._last_fmaint = fmaint
        atpm = s["atpm_rxn"]
        if atpm is not None and s["f_maint_nom"] > 0:
            # include the T-dependent maintenance factor set by set_temperature
            nf = getattr(self, "_ngam_T_factor", 1.0)
            val = s["atpm_nom_lb"] * (fmaint / s["f_maint_nom"]) * nf
            atpm.lower_bound = val
            atpm.upper_bound = max(val, atpm.upper_bound)

    # -- diagnostics --------------------------------------------------------
    def enzyme_mass(self, solution, T: float, pert: Optional[Perturbation] = None) -> float:
        """Total enzyme mass used (g/gDW) in a solution at temperature T."""
        pert = pert or Perturbation()
        total = 0.0
        for e in self.table:
            v = abs(solution.fluxes.get(e.rxn_id, 0.0))
            total += e.cost(T, Topt=pert.topt_of(e), dCp=pert.dCp_of(e)) * v
        return total
