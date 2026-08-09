# Claude Code prompt — P2c: diagnose the ~1.55 structural growth ceiling (what binds when the proteome is relaxed) (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`). A small, fast diagnostic (a handful of FBA
solves + shadow prices/FVA) to IDENTIFY the non-proteome constraint that caps rich-medium growth at
~1.55 h-1. Not a re-tune, not a big report change. Growth only.

NOTE TO USER: launch in an auto-approving mode.

CONTEXT (from the P2 re-run, calibration_vanderlinden_v2): with the growth law ON, kcat_scale
unlocked (x2.56, plausible) and rmax reached only 1.44 against a ~1.55 ceiling that is a NON-PROTEOME
network constraint — the residual to the observed ~2.4 is structural. This diagnostic NAMES that
ceiling so P3 can report it precisely (and tells us whether it is a fixable artefact or genuinely
structural). Operating point: rich BHI, at the optimal temperature (~40 C).

---

```
Work AUTONOMOUSLY; commit in parts; print a summary. Read first: src/etcgem/{enzyme_cost.py (the
metabolic-pool + biosynthesis-cap constraints, ATPM/NGAM),sectors.py,providers.py (set_medium BHI),
tpc.py,calibration.py (load the v2 posterior medians)}, strains/eciML1515/outputs/
calibration_vanderlinden_v2/ (tuned params), strains/eciML1515/strain.yaml, and the base iML1515
model (GAM/NGAM/ATPM, biomass reaction, exchange bounds). Do NOT re-tune or change model defaults.

PART A - reproduce and bracket the ceiling (rich BHI, optimal T ~40 C)
- Build the model at the rich (BHI) operating point (set_medium(pm,"BHI"), growth law ON) with the
  v2 posterior-median tuned parameters. Confirm the tuned rmax (~1.44).
- Now RELAX the enzyme/proteome limits and re-solve, in two steps, reporting rmax at each:
  (i) remove ONLY the biosynthesis (translation) cap (keep the metabolic pool);
  (ii) remove BOTH the biosynthesis cap AND the metabolic pool (i.e. the underlying iML1515 FBA at
       this medium/maintenance, proteome-unconstrained).
  If (ii) is ~1.55, the ceiling is genuinely NON-PROTEOME (structural). If (ii) is much higher, the
  ceiling was proteome after all — flag and reconcile with the P2 finding.

PART B - identify WHAT binds in the proteome-unconstrained solution
- With proteome removed (case ii), maximise growth and inspect the binding constraints:
  * SHADOW PRICES / duals: list the metabolites/constraints with non-zero shadow price (the active
    limits). Highlight the ATP metabolite's shadow price (is ATP the limiting currency?).
  * ACTIVE EXCHANGE BOUNDS: which uptake reactions sit at their upper bound — carbon (glucose and
    the BHI amino acids/components), O2, ions? Report each uptake flux vs its bound. Confirm the
    BHI medium actually opened enough carbon/O2 (not an artificially low uptake_ub).
  * MAINTENANCE: the ATPM/NGAM(T) flux and its shadow price at ~40 C.
  * GAM: the growth-associated maintenance (ATP coefficient in the biomass reaction) — is ATP
    regeneration the bottleneck?
  * Optionally a targeted FVA on the top flux-carrying / suspected-bottleneck reactions.
- Explicitly test and report each usual suspect: (1) carbon uptake bound, (2) O2 uptake bound,
  (3) ATPM/NGAM, (4) GAM / ATP balance, (5) a specific biomass-precursor or cofactor pathway.

PART C - classify + confirm by relaxing the culprit
- Classify the ceiling: FIXABLE ARTEFACT (e.g. an uptake_ub set too low, an O2 cap, a GAM value) vs
  GENUINELY STRUCTURAL (biomass composition / ATP yield / a network gap).
- CONFIRM the identification: relax the identified binding constraint (e.g. raise the offending
  uptake bound, or lower GAM to a literature value) and report whether the ceiling lifts toward
  ~2.4. State by how much, and whether the magnitude gap is closable.

PART D - outputs (small; no report rewrite)
- Save under strains/eciML1515/outputs/ceiling_diagnostic/: a summary.json (rmax at each relaxation
  step; the binding constraints + shadow prices; the classification; the relax-the-culprit result)
  and a short console SUMMARY that NAMES the ceiling in one or two sentences.
- Write a 3-5 line note (ceiling_diagnostic/NOTE.md) stating what the ~1.55 ceiling is and whether
  it is fixable or structural, for P3 to fold into the interpretation. Do NOT rewrite the report here.

VERIFY (report all)
1. rmax at: tuned (~1.44); biosynthesis-cap removed; both proteome constraints removed (~1.55?). Is
   the ceiling non-proteome (structural) or proteome (reconcile)?
2. The binding constraint(s) NAMED with shadow prices; the checked suspects (carbon/O2 uptake, ATPM,
   GAM, biomass precursor) each reported.
3. Classification (fixable vs structural) and the relax-the-culprit test (does the ceiling lift to
   ~2.4?).
4. summary.json + NOTE.md written; report NOT rewritten.

CONSTRAINTS
- Diagnostic only: no re-tune, no model-default changes (relaxations are temporary/in-solve), no
  report rewrite (P3 folds in the NOTE).
- Rich BHI operating point, optimal T, tuned (v2) parameters for the baseline.
- Report ACTUAL numbers (rmax, shadow prices, fluxes); name the ceiling concretely.
- Autonomous; commit in parts: "ceiling diagnostic: bracket rmax under proteome relaxation",
  "ceiling diagnostic: identify binding constraint (shadow prices / active bounds) + classify".
```
