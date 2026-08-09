# Claude Code prompt — allocation vs envelope Shapley variance decomposition (H1.3)

Paste the block below into Claude Code, run from the pipeline project folder.

---

```
Add a new analysis to this enzyme- and temperature-constrained GEM (etc-GEM) TPC
pipeline: decompose the variation a single genome can generate in its thermal
performance curve (TPC) into an ALLOCATION-set part, an ENVELOPE-set part, and
their INTERACTION, using a two-group functional-ANOVA / Shapley variance
decomposition. This operationalises hypothesis H1.3 of the project (a TPC has a
genome-set envelope and an allocation-set magnitude, separable). Do NOT change
any existing scientific/numerical code — reuse it.

ORIENT FIRST
- The project may be in one of two layouts: a package `tpc_pipeline` (flat) or,
  if a reorg has been applied, `src/etcgem`. Detect which, and add the new module
  and CLI entry in the same location and style as the existing sensitivity.py /
  dltkcat.py modules. Read README.md, sensitivity.py, enzyme_cost.py (the
  Perturbation dataclass), tpc.py (compute_tpc, TPC.descriptors) and __main__.py /
  cli.py before writing anything. Reuse: Perturbation, compute_tpc,
  TPC.descriptors, and the LHS helper in sensitivity.py (_lhs). Match code style.

THE TWO PARAMETER GROUPS (both already fields of Perturbation)
- ENVELOPE (genome-set thermal params): dTopt, topt_scale, dCp_scale.
- ALLOCATION (proteome): budget_scale (-> Perturbation.budget = default_budget *
  budget_scale) and any per-group allocation multipliers (group_alloc).
- Nominal / identity point: dTopt=0, topt_scale=1, dCp_scale=1, budget_scale=1,
  group_alloc={}. At nominal, compute_tpc reproduces the base TPC.

DESIGN (crossed, so the ANOVA is exact and balanced)
- Draw M allocation samples via LHS over the allocation ranges (envelope at
  nominal), and N envelope samples via LHS over the envelope ranges (allocation
  at nominal). Use sensitivity._lhs; separate seeds per group.
- Evaluate the model on the FULL crossed grid: for every (allocation_i,
  envelope_j) build one Perturbation combining allocation_i's allocation-field
  values with envelope_j's envelope-field values, run compute_tpc, and store each
  TPC descriptor. This yields an M x N matrix per descriptor
  (Topt_C, rmax, CTmin_C, CTmax_C, niche_width_C, B80_C, Ea_eV, skewness).
- Also store the MARGINAL curve ensembles for plotting the achievable ranges:
  allocation-only = the M curves at nominal envelope; envelope-only = the N curves
  at nominal allocation. (Store full growth curves for these two only, not the
  whole grid.)

DECOMPOSITION MATH (per descriptor, on the M x N matrix f_ij; drop NaN cells and
report how many valid)
- grand mean:      mu   = mean(f_ij)
- alloc main:      a_i  = mean_j f_ij - mu           (allocation-only effect)
- envelope main:   e_j  = mean_i f_ij - mu           (envelope-only effect)
- interaction:     g_ij = f_ij - mu - a_i - e_j
- variances:       V_A = mean_i a_i^2 ; V_E = mean_j e_j^2 ; V_AE = mean_ij g_ij^2
- total:           V = V_A + V_E + V_AE
- fractions (Sobol grouped, sum to 1): S_A=V_A/V, S_E=V_E/V, S_AE=V_AE/V
- total-effect indices: T_A = S_A + S_AE, T_E = S_E + S_AE
- Shapley effects (exact for two groups — split interaction evenly):
    phi_A = S_A + S_AE/2 ,  phi_E = S_E + S_AE/2   (phi_A + phi_E = 1)
Report all of V_A,V_E,V_AE,V, S_A,S_E,S_AE, T_A,T_E, phi_A,phi_E, n_valid.
(For balanced grids with imbalance from dropped NaNs, use the available cells and
note it; if a whole row/column is NaN, drop it before the ANOVA.)

CONFIG (add a `decomposition:` block; keep everything else as-is)
  decomposition:
    n_allocation: 24
    n_envelope: 24
    seed: 1
    allocation_params:
      budget_scale: [0.6, 1.1]
    envelope_params:
      dTopt:      [-6.0, 6.0]
      topt_scale: [0.7, 1.4]
      dCp_scale:  [0.5, 2.0]
  Reuse the same provider/temperature_grid/solver_timeout/crit_frac as the sweep.

MODULE API (new file, e.g. decomposition.py)
- run_decomposition(pm, temps_C, allocation_ranges, envelope_ranges,
    n_alloc, n_env, seed, crit_frac) -> a result object holding: the allocation
    and envelope sample tables, the per-descriptor M x N grids, the two marginal
    curve ensembles (+ temps), and a per-descriptor decomposition table (the math
    above). Include a .save(out_dir).
- Provide a runnable entry consistent with the repo: if there is a unified cli.py,
  add a `decompose` subcommand (e.g. `... decompose --strain eciML1515` or
  `--config`); otherwise add `python -m tpc_pipeline.decompose --config PATH`
  mirroring the existing resume.py / __main__.py pattern (set solver timeout from
  config; support --no-plots). Output dir: <output_dir>/decomposition/.

OUTPUTS (in <output_dir>/decomposition/)
- decomposition_table.csv : one row per descriptor with all the quantities above.
- grids.npz               : the M x N descriptor grids + allocation/envelope samples.
- allocation_only_curves.npy, envelope_only_curves.npy, temps_C.npy.
- figures (matplotlib, match plotting.py style):
  * achievable_ranges.png : two panels (allocation-only fan | envelope-only fan),
    each with median +/- IQR band and the nominal curve, so the width of each fan
    at each temperature is visible (this is the H1.3 achievable range).
  * variance_partition.png : for the key descriptors (Topt_C, rmax, CTmax_C,
    niche_width_C, Ea_eV), a stacked bar of S_A / S_E / S_AE (allocation /
    envelope / interaction fractions).
  * shapley_effects.png : grouped bar of phi_A vs phi_E per descriptor.
- A short summary.json: for each descriptor, which axis dominates and the fractions.

VERIFY (do these and report)
1. Toy model first (fast): run on configs/example_toy.yaml (or the toy provider)
   with small M=N (e.g. 12). Assert per descriptor that S_A+S_E+S_AE == 1 (within
   1e-6), all fractions in [0,1], phi_A+phi_E == 1. Print the table.
2. Check the expected H1.3 pattern and report whether it holds: allocation should
   dominate the variance of the MAGNITUDE descriptor (rmax / B0), while the
   envelope should dominate the TEMPERATURE descriptors (Topt_C, CTmax_C,
   niche_width_C). This is the scientific sanity check, not a hard assertion.
3. eciML1515 smoke run with a reduced grid (e.g. M=N=12, coarse temperature grid)
   to confirm it runs end-to-end and writes all outputs + figures. Use a per-LP
   solver timeout so no single sample hangs.
4. Update README with a short "Allocation vs envelope decomposition (H1.3)"
   section: the two groups, the crossed design, the Shapley-for-two-groups split,
   the CLI command, and the caveat below.

CAVEAT TO DOCUMENT
- The variance fractions are defined relative to the chosen input distributions
  (uniform over the configured ranges), so state the ranges with the result; they
  are a structural/in-silico decomposition of what the model can generate, not a
  claim about real cells. Keep allocation ranges in the pool-binding regime (e.g.
  budget_scale <= ~1.1) or the allocation axis will look artificially inert.

CONSTRAINTS
- Reuse Perturbation, compute_tpc, TPC.descriptors, _lhs; do not duplicate or
  modify the MMRT / enzyme-cost / provider logic.
- Write the descriptor extraction so it can later target respiration/CUE curves,
  not only growth (e.g. take the descriptors table as the input to the ANOVA), but
  do not implement respiration/CUE here.
- Keep it a self-contained, reviewable addition; commit as "add allocation/envelope
  Shapley variance decomposition (H1.3)".
```
