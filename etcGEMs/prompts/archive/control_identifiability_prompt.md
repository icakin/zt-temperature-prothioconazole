# Claude Code prompt — per-enzyme thermal control coefficients + identifiability (H1.1/H1.2)

Self-contained. Run from the project root (`.../MICROADAPT/etcGEMs`). Adds a new
per-enzyme sensitivity/identifiability analysis to the etc-GEM pipeline. Do NOT
change any existing scientific/numerical code — reuse it.

---

```
Add a per-enzyme control and identifiability analysis to this enzyme- and
temperature-constrained GEM (etc-GEM) TPC pipeline. It answers: which individual
enzymes' thermal parameters set the organismal thermal-performance-curve (TPC)
envelope (Topt, CT_max), which enzymes' capacity sets the baseline rate B0, and —
the identifiability flip side — which parameters the growth TPC is insensitive to
and therefore cannot be inferred from growth data alone (they need proteome / flux
data). This operationalises H1.1 (sequence-predictable envelope set by a few
thermostability-limiting enzymes), informs H1.2 (allocation/rate control), and
produces an experimental-design ranking of which enzymes most need accurate
predicted/measured parameters.

ORIENT FIRST
- Detect the layout/package (etcgem under src/, or tpc_pipeline) and add the module
  + CLI entry in the same place/style as sensitivity.py / dltkcat.py / (if present)
  decomposition.py. Read: providers.py (from_gecko), enzyme_cost.py (EnzymeEntry
  fields Topt/dCp/base_cost, EnzymeConstrainedModel.refresh_params, _costs), tpc.py
  (compute_tpc, TPC.descriptors), mmrt.py (relative_kcat_vec), and how
  dltkcat.apply_fits_to_provider mutates single entries then calls refresh_params —
  reuse that per-entry mutate-then-refresh pattern.
- This is a STRAIN-level diagnostic: it needs only a strain (no sweep experiment).

PER-ENZYME PARAMETERS (already on each EnzymeEntry)
- Topt_i  : the enzyme's optimum temperature (envelope).
- dCp_i   : MMRT curvature (breadth).
- kcat_i / base_cost_i : the enzyme's per-flux proteome cost = MW_i/(kcat_i*3600);
            lowering base_cost_i == raising kcat_i (capacity / allocation weight).

CONTROL COEFFICIENTS (define exactly)
Let D be an organismal TPC descriptor (Topt_C, CT_max_C, niche_width_C, rmax=B0,
Ea_eV). Compute by CENTRAL finite difference, perturbing ONE enzyme entry, calling
refresh_params(), recomputing, then restoring the entry.
1. Thermal control coefficient (envelope):
     CC[D, Topt_i] = ( D(Topt_i+dT) - D(Topt_i-dT) ) / (2*dT)      [dT default 1 K]
     CC[D, dCp_i]  via a fractional step on dCp_i (default 10%).
   For temperature descriptors these are ~dimensionless (K per K). Report the
   SUMMATION check: sum_i CC[Topt_org, Topt_i] should be order 1 and positive-
   dominated (organismal Topt is ~a control-weighted mix of enzyme optima).
2. Rate / flux control coefficient (B0), at a fixed temperature T:
     FCC_i(T) = dln(mu) / dln(kcat_i)  ~ ( ln mu(kcat_i*(1+f)) - ln mu(kcat_i*(1-f)) )
                 / (2f)                                              [f default 2%]
   where mu is growth at T. By the summation theorem sum_i FCC_i(T) ~ 1 at a
   pool-limited optimum — report the sum as a sanity check. Only enzymes with
   nonzero usage at T can have nonzero FCC, so restrict the finite differences to
   used enzymes.

TWO-STAGE DESIGN (so it is tractable on ~2500 enzymes)
Stage A - cheap screen over ALL enzymes, no per-enzyme re-solve:
  Solve the nominal TPC once. At each analysis temperature T, from the solution
  get enzyme usage share u_i(T) = cost_i(T)*|v_i(T)| / pool_budget. Compute the
  analytic thermal sensitivity of each enzyme's cost, s_i(T) = d ln cost_i/dT =
  -d ln(relative_kcat_i)/dT (from mmrt). Screen scores:
    rate_screen_i   = u_i(T)                       (who dominates the proteome at T)
    thermal_screen_i= u_i(T_high) * max(s_i(T_high), 0)   (heavily used AND losing
                       kcat fast on the hot flank -> candidate envelope setter)
  Rank; take the top `screen_top_k` (default 100) enzymes for Stage B.
Stage B - targeted finite differences:
  Compute the full thermal CCs (Topt_i, dCp_i vs Topt_org/CT_max/niche_width) for
  the top thermal candidates, and FCC_i(T) at the analysis temperatures for the
  top rate candidates (or all used enzymes if that is still cheap). Use a solver
  timeout so no single solve hangs; reuse warm starts where possible.

TEMPERATURES
- Default the three analysis temperatures to sub-optimal, optimal and supra-optimal
  taken from the nominal TPC (e.g. Topt-10, Topt, min(CT_max-2, Topt+8)), or accept
  an explicit list. These mirror the proposal's sub/opt/supra sampling.

IDENTIFIABILITY MAP
- A parameter p_i is identifiable from the growth TPC to the extent the descriptors
  respond to it: ident_i = max_D | CC_norm[D, p_i] | where CC_norm normalises each
  descriptor's CC by that descriptor's spread across the enzymes (so descriptors are
  comparable). Flag p_i "identifiable from growth TPC" if ident_i > threshold
  (config), else "requires omics (proteome/flux)". Enzymes never limiting at any T
  have ~0 control on every descriptor -> non-identifiable from growth alone; this is
  the etc-GEM under-determination that motivates the proteome/exometabolome/13C-MFA
  tiers. Report per-parameter and aggregate fractions.
- Caveat to document: this is a first-order, control-magnitude proxy for
  identifiability, not a full Fisher-information / profile-likelihood analysis.

CONFIG / EXPERIMENT
- If the repo uses experiment overlays, add experiments/control.yaml with kind:
  control and a `control:` block; run via `etcgem control --strain NAME
  --experiment control`. If not, read a `control:` block from the strain config or
  accept CLI flags. Schema:
    control:
      temperatures_C: null           # null -> derive sub/opt/supra from nominal TPC
      perturb: {Topt_K: 1.0, dCp_frac: 0.1, kcat_frac: 0.02}
      screen_top_k: 100
      descriptors: [Topt_C, CT_max_C, niche_width_C, rmax, Ea_eV]
      identifiable_threshold: 0.05

MODULE API (new file, e.g. control.py)
- run_control(pm, temps_C, params, ...) -> result holding: the usage-share table by
  temperature; Stage-A screen scores; Stage-B thermal CC table and FCC table; the
  identifiability table; and summary rankings. Include .save(out_dir).
- CLI subcommand `control` (strain-level). Output dir: strains/NAME/outputs/control/.
  Set solver timeout from config; dump resolved_config.yaml.

OUTPUTS (in .../outputs/control/)
- thermal_control.csv    : per enzyme (id, uniprot/enzyme_id, group), CC[Topt_org,
  Topt_i], CC[CT_max, Topt_i], CC[.., dCp_i], + screen score + rank.
- rate_control.csv       : per enzyme, FCC_i at each analysis temperature.
- usage_by_temperature.csv : u_i(T) for all used enzymes at each analysis T (shows
  which enzyme limits where, and how the limiting set switches with temperature).
- identifiability.csv    : per enzyme x parameter (Topt_i, dCp_i, kcat_i) the
  normalised control, the identifiable flag, and which descriptor it acts on most.
- summary.json           : top-N thermal-determinant enzymes (set Topt/CT_max),
  top-N rate-limiting enzymes at sub/opt/supra, the two summation checks, and the
  identifiable vs requires-omics fractions.
- figures (match plotting.py style):
  * thermal_control_bar.png : top-K |CC[Topt_org, Topt_i]| and |CC[CT_max, Topt_i]|.
  * bottleneck_vs_temperature.png : u_i(T) vs T for the top few limiting enzymes
    (the limiting-enzyme switch across temperature).
  * identifiability_hist.png : distribution of ident_i with the threshold marked.

VERIFY (do these; report)
1. Toy model first (configs toy / strains/_toy): run with small screen_top_k.
   Assert the analysis completes and both summation checks are order 1
   (sum_i FCC_i(T) ~ 1 at the optimum within a tolerance; sum of thermal CCs
   positive and order 1). Print the top enzymes.
2. Report the expected H1 pattern: thermal control should be CONCENTRATED in a
   small number of enzymes (envelope set by few -> consistent with sequence-
   predictability, H1.1), and the identifiable-from-growth fraction should be a
   minority of parameters (motivating omics, H1.2). State whether it holds.
3. eciML1515 smoke run (reduced screen_top_k, coarse TPC) end-to-end writing all
   outputs + figures with a per-LP solver timeout.
4. Update README with a short "Per-enzyme control & identifiability" section: the
   two control-coefficient families, the two-stage screen, the identifiability
   proxy + caveat, and the CLI command.

CONSTRAINTS
- Reuse EnzymeConstrainedModel per-entry mutation + refresh_params, compute_tpc,
  TPC.descriptors, mmrt; do NOT modify the MMRT / enzyme-cost / provider internals.
- Only enzymes with nonzero usage can have nonzero rate control — prune to those.
- Write descriptor extraction so it can later target respiration/CUE, not only
  growth, but do not implement those here.
- Keep it a self-contained, reviewable addition; commit as "add per-enzyme thermal
  control coefficients + identifiability (H1.1/H1.2)".
```
