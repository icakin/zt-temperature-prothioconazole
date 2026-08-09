# Claude Code prompt — three fixes: proteome-wide identifiability, dCp→Ea calibration, sector-based H1.3 decomposition (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`).

NOTE TO USER: launch in an auto-approving mode (accept edits + allow commands, e.g.
`--dangerously-skip-permissions`) so it runs unattended. This does code changes,
a recalibration, several pipeline re-runs, and re-renders the report; it will take a
while. Everything is committed in git so it is revertible.

---

```
Work AUTONOMOUSLY end to end. Do not pause for approval between steps; make
reasonable decisions and commit your own work. Only stop if a step hard-fails in a
way you cannot resolve, and then report exactly what blocked you. Print a concise
summary at the end. Read the relevant modules before editing: src/etcgem/control.py,
providers.py, tpc.py (Ea/descriptor computation), sensitivity.py (_make_perturbation),
decomposition.py, config.py, cli.py; strains/eciML1515/strain.yaml; experiments/*.yaml;
reports/etcgem/{report.qmd,assemble.py}. Do NOT change the core MMRT / cost math.

These fixes address a review of the current report:
- (A) The identifiability table shows 44/45 "identifiable from growth", but that is
  computed only over the top-K screened (high-control) enzymes, so it contradicts
  the text and would undercut H1.2. Make identifiability proteome-wide.
- (B) Nominal Ea = 1.26 eV is ~2x the conserved ~0.65 eV; it is a pure consequence of
  the default dCp = -12 (chosen for curve width, not Ea). Calibrate dCp to a target Ea.
- (C) The H1.3 decomposition uses a single scalar allocation axis (budget_scale), so
  the near-perfect separability is partly by construction. Re-run it with the
  mechanistic sector allocation axis (f_metab/f_maint).

=====================================================================
FIX A - proteome-wide identifiability (src/etcgem/control.py)
=====================================================================
Currently identifiability.csv covers only the top-K finite-difference enzymes. Change
it to cover ALL enzymes using the cheap first-order screen as the identifiability
proxy, refined by finite difference where available:
- For every enzyme (not just top-K), compute the first-order control/identifiability
  proxy from the existing screen quantities (usage share u_i(T) x analytic parameter
  sensitivity), normalised comparably to the current `ident` metric. Enzymes with ~0
  usage get ~0 control -> not identifiable from growth. Where a top-K finite-
  difference control coefficient exists, use it (more accurate) and mark the row
  `refined = True`; otherwise `refined = False` (screen-only proxy).
- Write identifiability.csv over ALL enzymes x parameters (Topt_i, dCp_i, kcat_i),
  keeping columns rxn_id, enzyme_id, parameter, ident, top_descriptor,
  identifiable_from_growth, refined.
- In control's summary.json report BOTH: the PROTEOME-WIDE identifiable fraction
  (expected small) and, separately, the fraction among the top-K control enzymes
  (expected high). Keep the top-K thermal_control.csv ranking unchanged.
- Document in the code/README that identifiability here is a first-order,
  control-magnitude proxy (not a full Fisher-information analysis).
Update reports/etcgem/report.qmd section 5 so tbl-ident reports the proteome-wide
fraction (and can mention the top-K fraction), consistent with the prose that few
enzymes dominate control while most parameters are weakly identifiable (H1.1/H1.2).
Test on the toy strain (`etcgem control --strain _toy --experiment control_quick` or
equivalent) and confirm identifiability.csv now has ~one row per enzyme-parameter and
the proteome-wide identifiable fraction is small. Commit: "identifiability: compute
proteome-wide, not just top-K".

=====================================================================
FIX B - calibrate default dCp to a target Ea
=====================================================================
Add a calibration that chooses default_dCp so the nominal rising-limb Ea hits a
target (Ea increases monotonically with |dCp|, so use bisection):
- Add calibrate_dCp_to_Ea(build_provider_fn/cfg, target_Ea_eV, lo=-20, hi=-3,
  tol=0.02) that bisects dCp, rebuilding the provider and computing the nominal TPC
  Ea each step (reuse tpc.compute_tpc + TPC.descriptors), returning the dCp giving
  Ea≈target. Put it where providers/config can reach it.
- Add a CLI `etcgem calibrate-dcp --strain NAME --target-ea 0.65` that runs it,
  PRINTS the calibrated dCp and achieved Ea, and updates strains/NAME/strain.yaml
  provider.default_dCp in place. IMPORTANT: the target is configurable and defaults
  to 0.65 eV (metabolic-theory value) — the PI should set it to the measured
  bacterial growth-TPC Ea; note this in the CLI help and README.
- A shallower dCp broadens the TPC and may push CT_max beyond the temperature grid.
  After calibrating, check the nominal CT_max; if CT_max > stop_C - 2, raise the
  strain's temperature_grid stop_C (e.g. to 70) and n proportionally so CT_max stays
  resolved, and report that you did so.
Run it for eciML1515 (target 0.65 eV unless told otherwise), updating strain.yaml.
Commit: "calibrate default dCp to target Ea (eciML1515 -> ~0.65 eV)".

=====================================================================
FIX C - sector-based H1.3 decomposition
=====================================================================
- Verify sensitivity._make_perturbation maps f_metab/f_maint allocation samples onto
  Perturbation (calling set_allocation); extend it if not, so run_decomposition can
  take sector allocation params. (The sectors sweep already works, so the plumbing
  likely exists.)
- Create experiments/decomposition_sectors.yaml: kind: decompose, with
    proteome_sectors: {enabled: true, atpm_reaction: null}   # sector override
    decomposition:
      n_allocation: 16
      n_envelope: 16
      seed: 1
      allocation_params: {f_metab: [0.4, 0.6], f_maint: [0.05, 0.25]}
      envelope_params:  {dTopt: [-6,6], topt_scale: [0.7,1.4], dCp_scale: [0.5,2.0]}
  If experiment-level proteome_sectors override does not activate sectors, temporarily
  enable in strain.yaml for this run then revert (as in the sectors sweep).
- The allocation axis in the decomposition is now the sector split, so the H1.3
  separability claim rests on mechanistic allocation degrees of freedom, not a scalar
  budget. Keep the ANOVA/Shapley math unchanged.

=====================================================================
RE-RUN + RE-RENDER (after A, B, C are in)
=====================================================================
With the recalibrated dCp now in strain.yaml, regenerate the report's source runs
WITH figures (do not pass --no-plots), because dCp changed the model:
  etcgem sweep     --strain eciML1515 --experiment default
  etcgem decompose --strain eciML1515 --experiment decomposition_sectors
  etcgem control   --strain eciML1515 --experiment control      (full, not _quick)
  etcgem sweep     --strain eciML1515 --experiment calibrated
  etcgem sweep     --strain eciML1515 --experiment dltkcat_ext   (if that is the dltkcat run)
  etcgem sweep     --strain eciML1515 --experiment sectors       (sectors enabled as above)
Then update reports/etcgem/assemble.py so the `decompose` source points at the new
decompose_decomposition_sectors run, and the control source at the full control run;
run `python reports/etcgem/assemble.py`; then `cd reports/etcgem && quarto render`.
Update report.qmd section 4 to note the allocation axis is now the proteome sectors
(f_metab/f_maint), and section 1/9 to note Ea is calibrated to the target. Leave the
author placeholders.

VERIFY + SUMMARISE (report all)
1. Proteome-wide identifiable fraction (should be small) and top-K fraction (high).
2. Calibrated dCp and achieved nominal Ea (~target); confirm CT_max still resolved.
3. Sector-based decomposition: report phi_A/phi_E for rmax and Topt_C and whether the
   allocation/envelope separation persists with the mechanistic allocation axis.
4. Report re-renders to _output/report.pdf with the updated tables/figures.

CONSTRAINTS
- Autonomous; commit each fix separately; do not change core MMRT/cost math.
- Never pass --no-plots for the report's source runs.
- Leave strains/eciML1515/strain.yaml with proteome_sectors disabled at the end
  (sectors active only within the sector experiments/runs).
- Identifiability is a first-order proxy; Ea target is configurable; state both.
```
