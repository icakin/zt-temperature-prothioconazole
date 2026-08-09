# Claude Code prompt — P2: Bayesian tuning on Van Derlinden (MG1655, rich/BHI), unified param set incl. a new kcat_scale knob (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). AFTER P1. This is the pivotal stage: tune
the model to the exact-strain rich curve and learn what corrections the data demand. Additive
inverse analysis; emergent model unchanged; **no report edits and no sensitivity/decomposition
re-run** (those are P3). Growth only.

NOTE TO USER: launch in an auto-approving mode. This run is much faster with Gurobi (see PART A0):
install `gurobipy` AND a free academic licence FIRST (the bundled trial licence is size-limited and
will reject the ~7,700-reaction ecModel). With Gurobi the run drops from ~3 h to ~30-45 min and the
solves are exact (no GLPK stalls); without it the prompt still runs on GLPK, just slowly.

CONTEXT (from P1): at the rich (BHI) operating point the emergent model predicts Van Derlinden's
shape reasonably (Topt 36.7 vs 40) but under-predicts magnitude ~2.3x (rmax 1.03 vs 2.40),
over-predicts the upper limit (CTmax 50 vs ~46), and gives Ea 0.82 vs a digitized-observed 0.64.
CAUTION: the observed Ea (0.64) is a figure-digitized cloud-centre and may be artefactually low
(the model's 0.82 is nearer the ~0.85 bacterial benchmark) — do NOT let the fit chase it; the
broad dCp prior + discrepancy term should absorb it.

TWO IMPORTANT SET-UP POINTS:
1) The shared default operating point is STILL glucose-minimal (default_medium_proteome: Glucose).
   P2 MUST explicitly fit at the RICH (BHI) operating point — set_medium(pm,"BHI") with the LB
   sector allocation — not rely on the default.
2) Add a kcat_scale knob (does not exist yet) — the biologically-plausible metabolic-magnitude
   lever, replacing budget_scale (they are degenerate for the pool constraint).

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first:
src/etcgem/{enzyme_cost.py (Perturbation; base_cost / rel_kcat in _costs_unfolding; kappa_scale),
calibration.py (Simulator, load curve, priors, emcee, log_likelihood_sd, trusted-curve loader),
validation.py (vanderlinden_bhi CurveSpec + _set_medium -> BHI),providers.py (set_medium BHI),
sectors.py,cli.py}, the Van Derlinden CSV
(strains/eciML1515/thermal/sources/vanderlinden2012_intfoodmicro/vanderlinden2012_mg1655_bhi_tpc.csv).
Do NOT modify the emergent-model defaults, the report, or the sensitivity/decomposition code.

PART 0 - housekeeping: Van Derlinden ONLY (remove the Erdos cross-check)
- Remove the erdos_LB secondary cross-check from the validation (drop it from validation.py's
  SECONDARY list and regenerate the validation outputs/figure so only vanderlinden_bhi appears).
  The whole pipeline — validation and tuning — is Van Derlinden only. (Noll is already gone.)

RE-RUN CONTEXT (this supersedes the first P2 run): P1b fixed the model — the peak now rounds
(T-dependent maintenance restored) and the coupled proteome-conserving growth law is available as
a toggle. kcat_scale already exists (from the first P2); tm_scale, ngam_scale, ngam_steepness were
added in P1b. KEY: turn the growth law ON. In the first P2 the biosynthesis cap bound at the peak,
so kappa relaxed it (->1.56) but kcat_scale was wasted (it acts on the metabolic pool, which had
slack). With the growth law ON the binding constraint shifts to the coupled metabolic sector, so
kcat_scale becomes the EFFECTIVE magnitude lever. So this re-run tests: with kcat unlocked, what
kcat_scale is demanded to reach ~2.4, and is it physically plausible (a few-fold in-vivo/in-vitro
gap) or implausibly large (=> genuinely structural).

PART A0 - SPEEDUPS (same posterior, faster): fast solver + efficient sampling
These are legitimate accelerations only — they do NOT change the model, likelihood, priors, or the
target posterior; they change how fast/cleanly it is computed. Do all three.
(a) FAST LP SOLVER with graceful fallback + PRE-FLIGHT:
    - At start-up, try to switch every provider/worker model to Gurobi:
      `pm.ec.model.solver = "gurobi"` (and the same in _worker_init). Wrap in try/except; if Gurobi
      is unavailable or unlicensed, fall back to the current GLPK. PRINT loudly which solver is in
      use ("[solver] gurobi" vs "[solver] GLPK - slow; install gurobipy + academic licence").
    - STOP-AND-FLAG (do NOT silently proceed on GLPK): Gurobi is expected to be installed and
      licensed in this venv. If the switch to Gurobi FAILS and the run would fall back to GLPK,
      ABORT immediately with a clear message: "[solver] Gurobi NOT active - the venv running this
      job is not seeing Gurobi. Expected [solver] gurobi. Stopping so you can fix the environment
      (activate the venv with gurobipy + ~/gurobi.lic) rather than run ~3 h on GLPK." Do NOT start
      the multi-hour chain in this case. (Only proceed on GLPK if the user has explicitly overridden
      with an ALLOW_GLPK flag.)
    - PRE-FLIGHT before launching the chain: build the model at the rich BHI operating point and the
      optimal T and solve ONCE. If on Gurobi, confirm the FULL-SIZE model solves (catch the
      "size-limited licence" / model-too-large error and ABORT with a clear message telling the user
      to install the free academic licence) so a licence problem fails NOW, not 2 h in. Report the
      single-solve wall-time for gurobi vs the known GLPK time.
    - SOLVER TIMEOUT: the existing 2 s cap is a GLPK band-aid for degenerate stalls and can silently
      return suboptimal (wrong) likelihoods. On Gurobi, raise it (e.g. 30 s; Gurobi returns the true
      optimum in ms) so no eval is truncated. Keep the finite/optimality guards either way; if any
      solve returns non-optimal status, return -inf for that theta (do not accept a stalled value).
(b) CONVERGENCE-BASED EARLY STOP (instead of a fixed n_steps): every ~200 steps compute the
    integrated autocorrelation time (sampler.get_autocorr_time(tol=0)); STOP once the chain is
    longer than ~50*max(tau) AND the min n_eff across parameters is >= 400. Keep a hard n_steps_max
    cap. This saves the tail when converged and loses nothing when not. Report the stopping reason.
(c) WARM-START to shrink burn-in: before sampling, run a quick optimiser to the posterior mode
    (scipy.optimize; e.g. a short differential_evolution or L-BFGS on the negative log-posterior),
    then initialise the walkers in a TIGHT Gaussian ball around that mode (still covering the prior
    widths enough to explore). This lets n_burn drop (~500 -> ~150) because the ensemble no longer
    has to migrate from the emergent point. If the optimiser fails, fall back to the current
    emergent-point initialisation. Walker count: keep it a multiple of the process-pool size for
    full core use.

PART A - confirm knobs + turn the growth law ON
- Do NOT re-add knobs: kcat_scale (global metabolic-turnover level lever; leaves Ea/Topt/CTmax
  unchanged) exists from the first P2; tm_scale, ngam_scale, ngam_steepness exist from P1b. Just
  confirm they are wired and no-op at their identity values.
- Turn the coupled growth-law partition ON (biosynthesis_growth_law: true) at the rich operating
  point, so kcat_scale acts on the binding (metabolic) constraint.

PART B - unified free-parameter set + rich operating point + provenance priors
- Fit at the RICH (BHI) operating point with the GROWTH LAW ON: explicitly set_medium(pm,"BHI")
  (LB sector allocation) + biosynthesis_growth_law=true; do NOT rely on the glucose-minimal default.
- Free-parameter set (the SAME canonical set the sensitivity/decomposition will use in P3):
  envelope {dTopt, topt_scale, dCp_scale, dTm, tm_scale} + magnitude {kcat_scale, kappa_scale} +
  allocation {f_metab, f_maint} + maintenance {ngam_scale, ngam_steepness} + a noise/discrepancy
  term sigma_disc. DROP budget_scale (degenerate with kcat_scale). Hold P_total FIXED.
  FREE sigma (the in-vivo saturation / total enzyme budget) — this is the EFFECTIVE magnitude
  lever now that the pools are reconciled (P1c): kcat_scale SATURATES at ~1.72 (once metabolic
  enzymes are cheap the biosynthesis/ribosome budget binds, which kcat cannot relax), whereas
  scaling the whole budget via sigma raises both caps and reaches ~2.4. sigma and kcat_scale are
  partly degenerate for the metabolic pool — report the combined magnitude correction and the ridge.
  The growth-law slope stays a CONFIG value (degenerate with kappa_scale via kappa_eff, so free
  only kappa_scale, not the slope).
- PROVENANCE-BASED priors:
  * measured allocation f_metab, f_maint: TIGHT Normal priors centred on the measured (rich/LB)
    values with a small SD (measurement-level wiggle only — we are not free-fitting the proteome).
  * kcat_scale: broad LogNormal about 1 (the in-vitro->in-vivo gap is real and several-fold);
  * sigma (in-vivo saturation / enzyme budget): prior centred on the literature ~0.45, BROAD but
    BOUNDED to (0, 1) (physical); the posterior sigma vs the 0.4-0.5 literature range IS the
    magnitude headline (how much of the rich-medium capacity is above literature);
  * kappa_scale: broad LogNormal about 1;
  * dCp_scale: LogNormal about 1 (sd ~0.30 on log);
  * dTopt ~ N(0, ~4 K); topt_scale ~ LogNormal about 1 (modest); dTm ~ N(0, ~6 K);
    tm_scale ~ LogNormal about 1 (modest — the falling-shoulder spread);
  * ngam_scale ~ LogNormal about 1 (broad — the Li maintenance amplitude is borrowed from yeast);
    ngam_steepness ~ LogNormal about 1 (broad — the maintenance T-dependence, borrowed);
  * sigma_disc: HalfNormal (there is NO measured SD for Van Derlinden -> free discrepancy term).

PART C - fit Van Derlinden on raw ABSOLUTE rates
- Load the Van Derlinden curve (7-46 C, 17 rate points, no SD) via the shared trusted-curve
  loader; fit on raw absolute rate (1/h) with the Gaussian discrepancy likelihood (var = sigma_disc^2).
  The wide 7-46 C span constrains the envelope (dTopt, dCp_scale, dTm) that Noll could not.
- emcee, gradient-free, parallelised, warm-started at the posterior mode (PART A0c); run until the
  PART A0b convergence criterion (n_eff >~ 400 per parameter) or n_steps_max; keep the finite/
  optimality guards. Report solver used, acceptance, autocorrelation, n_eff, stopping reason, and
  wall-time (and the gurobi-vs-GLPK single-solve time from the pre-flight).

PART D - outputs (NO report edits)
- Save under strains/eciML1515/outputs/calibration_vanderlinden_v3/ (do NOT overwrite the earlier
  runs; v1 = flat-top/kcat-wasted, v2 = growth-law-on but redundant pool. v3 = reconciled single
  pool with sigma freed):
  * chain + summary.json (per-param posterior median/90% CI, convergence, curve, medium=BHI,
    wall-time).
  * PRIOR-vs-POSTERIOR TPC on raw absolute rate: emergent/prior curve, posterior-predictive band
    (median + 90% CI), and the Van Derlinden points, 7-46 C.
  * demanded_corrections.csv: per free parameter, prior vs posterior + interpretable correction.
    EXPECTED (report whether borne out): magnitude up (kcat_scale and/or kappa_scale > 1 to reach
    ~2.4), dTm DOWN ~4 K (pull CTmax 50 -> ~46), dTopt small, dCp weakly constrained (do not
    over-tune to the possibly-flattened observed Ea), f_metab/f_maint near their measured values.
  * corner plot (expect a magnitude-degeneracy ridge between kcat_scale and kappa_scale).
  * console SUMMARY answering the headline questions: on the reconciled single pool with sigma
    freed, does the posterior reach ~2.4? THE MAGNITUDE HEADLINE = the posterior sigma vs the
    0.4-0.5 literature range: how much of the rich-medium capacity is above literature (the
    "genuinely higher in-vivo capacity" statement), and how sigma/kcat_scale share the magnitude
    (the degeneracy ridge). Does the growth-law coupling keep the good shape (Topt ~40, Ea ~0.64),
    and is CTmax corrected by dTm? Note the observed 2.4 is a figure-digitized high-end value.
- Keep emergent vs tuned distinct; nothing here is relabelled as prediction. Do NOT edit the report.

VERIFY (report all)
0. Solver in use printed: MUST be "[solver] gurobi". If Gurobi is not active, the run ABORTED early
   with the stop-and-flag message (did NOT grind on GLPK). Pre-flight full-size solve passed (or
   aborted with a clear licence message); single-solve time reported. Early-stop criterion +
   warm-start active; these are confirmed NOT to change the posterior (same model/likelihood/
   priors), only the speed.
1. kcat_scale works (raises rmax; leaves Ea/Topt/CTmax unchanged); it is global; budget_scale dropped.
2. Fit is at the rich BHI operating point (not the glucose default); free set + provenance priors as above.
3. Posterior: per-parameter median/90% CI + demanded corrections; does it reach ~2.4; demanded
   kcat_scale and whether it's physically plausible; dTm correction to CTmax; which magnitude
   lever binds.
4. Convergence (n_eff >~400); prior-vs-posterior TPC + corner exist; emergent model + report unchanged.

CONSTRAINTS
- Additive inverse analysis only; no changes to emergent defaults, core math (beyond the kcat_scale
  hook), or the report. Exact likelihood + gradient-free emcee (not ABC, not Stan/brms).
- Rich (BHI) operating point, explicitly set. Measured allocation tight-prior'd (not free-fit).
  Magnitude expressed via kcat_scale/kappa_scale (pool held at measured; the physical-range
  interpretation discriminates pool vs kcat).
- Fit RAW ABSOLUTE rates; single Van Derlinden curve (joint/per-curve fits are later).
- Autonomous; commit in parts: "validation: Van Derlinden only (drop Erdos cross-check)",
  "calibration: Gurobi solver with GLPK fallback + pre-flight; autocorr early-stop; warm-start init",
  "add kcat_scale (global metabolic turnover) knob",
  "calibration: unified param set + rich BHI operating point + provenance priors",
  "emcee tuning on Van Derlinden (long chains)",
  "calibration_vanderlinden outputs: prior-vs-posterior TPC, corrections, corner".
```
