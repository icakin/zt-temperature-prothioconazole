# Claude Code prompt — P1c: reconcile the two redundant enzyme-mass pools into a single proteome budget (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). Model-structure fix + emergent diagnostic. It
does NOT re-run the Bayesian (that is the next P2 re-run) and does NOT rewrite the report. Growth only.

NOTE TO USER: launch in an auto-approving mode.

CONTEXT (from the P2c ceiling diagnostic): the ~1.55 "structural" ceiling is NOT structural — it is a
PROTEOME double-accounting. The model carries TWO redundant enzyme-mass pools that both bind:
- our sMOMENT SECTOR pool (P_total*f_metab*sigma etc.) — binds ~1.44;
- the leftover GECKO base pool `prot_pool_exchange` (upper bound ~0.091 g/gDW) — binds ~1.55.
`kcat_scale` scales only our sMOMENT per-flux cost, NOT the GECKO base pool, so the magnitude levers
saturate against a pool they cannot touch. Relaxing ALL proteome jumps growth to 7.6 (uptake/ATPM/
GAM/precursors all had slack). Fix: ONE proteome budget, with the magnitude levers acting on the true
binding pool. Outputs of the diagnostic: strains/eciML1515/outputs/ceiling_diagnostic/ (summary.json,
NOTE.md) — read them.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
strains/eciML1515/outputs/ceiling_diagnostic/{summary.json,NOTE.md}, src/etcgem/{enzyme_cost.py
(the sMOMENT sector pool + biosynthesis cap; how base_cost/kcat_scale enter),sectors.py (the sector
partition + pool budget; how it relates to the GECKO prot_pool_exchange),providers.py (the GECKO
ecModel prot_pool_exchange bound; set_medium BHI),config.py,tpc.py,calibration.py}, and
strains/eciML1515/strain.yaml. Do NOT re-run the Bayesian or rewrite the report.

PART A - reconcile to a SINGLE proteome budget
- The GECKO base pool (`prot_pool_exchange`, ~0.091 g/gDW) and the sMOMENT sector (metabolic) pool
  both account for the metabolic-enzyme mass — redundantly. Make the sector / growth-law pool the
  SOLE proteome accounting:
  * relax/remove the GECKO base `prot_pool_exchange` upper bound as an independent binding constraint
    (set it non-binding: to +inf, or exactly to the sector total budget) so it can no longer cap
    growth on its own; and
  * confirm the sector partition (metabolic + biosynthesis + maintenance + chaperone summing to
    P_total, with the growth-law coupling) is then the only proteome limit.
- Ensure the magnitude levers (kcat_scale, sigma, P_total) act on the resulting SINGLE binding pool
  (kcat_scale must effectively relax the true binding constraint, not a redundant one).
- SANITY: check the reconciled single budget is sensible in magnitude (the sector metabolic pool
  P_total*f_metab*sigma ~ 0.1 g/gDW is the same order as the GECKO 0.091 — they are two computations
  of the same quantity; keep the sector one, which is medium- and growth-law-aware). Note the value.

PART B - re-check emergent behaviour + the (now movable) ceiling
- At the rich BHI optimum, report the reconciled emergent rmax, and show that kcat_scale / sigma now
  MOVE the ceiling (bracket a couple of values) — i.e. the levers act on the true pool.
- Report what enzyme-budget scaling (equivalently sigma, or kcat_scale) is needed to reach the
  observed ~2.4, and state honestly whether it is WITHIN the literature sigma range (0.4-0.5) or
  ABOVE it (~0.77) — i.e. how much of the rich-medium gap is the fixable redundancy vs a genuinely
  higher in-vivo capacity (also note the observed 2.4 is a digitized high-end value).

PART C - outputs (no re-tune, no report rewrite)
- Regenerate the emergent reference TPC + emergent-vs-Van-Derlinden validation on the reconciled
  model; report the descriptors. Save strains/eciML1515/outputs/pool_reconciliation/NOTE.md (what
  changed; the single budget; the sigma/kcat needed for 2.4 vs literature) for P3.
- Do NOT re-run the Bayesian (the P2 re-run follows) and do NOT rewrite the report.

VERIFY (report all)
1. Single proteome budget: the GECKO base `prot_pool_exchange` no longer independently binds; the
   sector/growth-law pool is the sole proteome limit; the reconciled budget value.
2. kcat_scale / sigma now move the rich-optimum ceiling (bracketed); the sigma (or budget scaling)
   needed to reach ~2.4, and whether it is within (0.4-0.5) or above (~0.77) the literature range.
3. Reconciled emergent descriptors (rmax, Topt, CTmax, Ea) at rich BHI; validation vs Van Derlinden.
4. NOTE.md written; Bayesian NOT re-run; report NOT rewritten.

CONSTRAINTS
- Model fix (pool reconciliation) + emergent diagnostic only. No re-tune, no report rewrite.
- Keep ONE proteome budget = the sector/growth-law pool; the magnitude levers must act on it.
- Report ACTUAL numbers; be honest about the residual (fixable redundancy vs above-literature capacity).
- Autonomous; commit in parts: "reconcile redundant enzyme pools -> single sector/growth-law budget",
  "magnitude levers act on the true pool + emergent diagnostics (ceiling now movable)".
```
