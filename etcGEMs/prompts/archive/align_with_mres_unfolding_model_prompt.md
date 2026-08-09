# Claude Code prompt — align the etcGEM with the MRes/Li thermal model (two-state unfolding + grounded per-enzyme Topt/Tm)

Run from the project root (`.../MICROADAPT/etcGEMs`). Autonomous. This adds the
data-grounded, denaturation-based thermal model from the MRes (Madkaikar 2023) /
Li et al. 2021 as a new mode, so the falling limb / thermal limits come from
per-enzyme melting temperature (Tm) rather than from a fictional dCp distribution.
No BRENDA, no TOMER needed for E. coli — reuse the MRes's assembled parameter tables.

NOTE TO USER: launch in an auto-approving mode so it runs unattended.

---

```
Work AUTONOMOUSLY end to end; commit each part; print a summary. This is a MODEL
addition; keep it BACKWARD-COMPATIBLE (a new thermal mode; the existing MMRT-only
peak-normalised path must stay byte-for-byte identical when the new mode is off).

BACKGROUND / SOURCES
- We currently model kcat(T) with MMRT + peak-normalisation and a homogeneous dCp,
  which couples breadth to Ea and makes the falling limb depend on a hard-to-ground
  dCp. The MRes (github.com/adimadkaikar/MResProject) and Li et al. 2021 instead use
  a two-state native<->unfolded equilibrium keyed on per-enzyme Tm (so denaturation
  sets the high-temperature collapse), MMRT for kcat(T), and a temperature-dependent
  maintenance (NGAM). Their sensitivity analysis showed Tm/denaturation drives the
  high-T fall, Topt sets the peak, and dCp barely matters.

STEP 0 - get the MRes assets
- git clone https://github.com/adimadkaikar/MResProject (outside the etcGEMs tree).
- READ code/etcpy/etc.py (the temperature-dependence code: native fraction, kcat(T),
  NGAM) and code/etcpy/tempDepCode.py (parameter sampling); these implement the
  thesis equations (native fraction sigmoid in dGu(T); kcat(T) via transition-state/
  MMRT; NGAM Arrhenius-like). Use them as the reference implementation.
- Use data/model_enzyme_params*.csv (initial per-enzyme Topt, Tm, dCp) and
  BestParamsTopt*.csv (tuned) as the per-enzyme parameter source. These are keyed by
  enzyme/gene/UniProt for E. coli iJO1366.

STEP 1 - add per-enzyme Tm + a two-state unfolding thermal mode
- Extend EnzymeEntry with Tm (melting temperature, K) and the unfolding thermodynamic
  parameters needed by the two-state model (dHu/dCp_u or the protein-length
  approximations from the thesis Eq 2.14 when Tm/T90 are missing).
- Add a strain/provider option `thermal_model: mmrt | unfolding` (default mmrt =
  current behaviour). In `unfolding` mode, the effective per-flux enzyme cost is
      cost_i(T) = base_cost_i / ( rel_kcat_i(T) * f_N_i(T) )
  where rel_kcat_i(T) is the MMRT temperature factor and f_N_i(T) in [0,1] is the
  NATIVE (folded) fraction from the two-state model (sigmoid in dGu(T); -> 0 above Tm),
  so denaturation inflates cost sharply near Tm and produces the falling limb. Follow
  etcpy/etc.py for the exact functional forms and constants.
- Add the temperature-dependent NGAM (maintenance) term (rising with T, with a basal
  floor) as in the thesis Eq 2.3 / etcpy; expose it as an option (detect the ATPM/NGAM
  reaction). Keep it off in mmrt mode.
- Do NOT peak-normalise in unfolding mode (kcat(Topt) is anchored as in the MRes/Li
  model). Keep peak-normalisation for the mmrt mode.

STEP 2 - ingest grounded per-enzyme Topt/Tm and map to eciML1515
- Build a loader that reads the MRes parameter table(s) and joins them to the
  eciML1515 enzymes by UniProt id (our EnzymeEntry.enzyme_id). Prefer BestParamsTopt
  (tuned) if present, else model_enzyme_params (initial). Report COVERAGE (how many of
  our enzymes matched). For unmatched enzymes, fall back to the dataset means from the
  thesis: Tm ~ 55.6 C (sd 7.6), Topt ~ 37 C, dCp ~ -4 kJ/mol/K. (Note in the README
  that sequence predictors — DeepET/TOMER for Topt/Tm — are the route for OTHER strains
  in the library where no measured values exist; not needed for E. coli here.)
- Wire this into from_gecko: in unfolding mode with a params table configured, set each
  enzyme's Topt/Tm/dCp from the joined table (fallbacks for gaps).

STEP 3 - verify against the MRes and check decoupling
- Build eciML1515 in unfolding mode with the grounded parameters; compute the nominal
  TPC. Confirm: a realistic E. coli curve (Topt ~37-40 C, steep decline above ~46 C,
  CT_max ~48-50 C set by Tm, not by the grid), and that Ea (rising limb) and
  breadth/CT_max are now DECOUPLED (Ea from kcat(T)/Topt; the falling limb from Tm).
- If feasible, score the nominal curve against the MRes empirical data
  (ExpGrowth.csv from the repo) and report R^2 / how it compares to the MRes's 0.94.
- BACKWARD-COMPAT GATE: with thermal_model=mmrt, a toy and an eciML1515 nominal TPC are
  byte-identical (1e-9) to current code.

STEP 4 - re-run report + document
- Re-run the report's source runs in unfolding mode (default sweep, decompose_sectors,
  control, calibrated) WITH figures; `python reports/etcgem/assemble.py`;
  `cd reports/etcgem && quarto render`.
- Update the report: note the thermal limits are now data-grounded via per-enzyme Tm
  (Leuenberger/Meltome) and Topt (Li-Engqvist), the falling limb is set by denaturation
  not dCp (so breadth and Ea decouple), and this follows the MRes [@Madkaikar2023] and
  Li et al. [@Li2021] model. Add Madkaikar2023 to references.bib (MRes thesis, Imperial
  College London, 2023) and any missing data-source refs (Leuenberger2017 melting
  proteome; LiEngqvist2019 already ~ ref for enzyme Topt).

VERIFY (report all)
1. Coverage: how many eciML1515 enzymes got measured Topt/Tm from the MRes table vs
   fell back to means.
2. Decoupling + realism: nominal Topt, CT_max, Ea in unfolding mode; confirm CT_max is
   set by Tm and Ea is independent of it.
3. Backward-compat gate (mmrt mode identical).
4. Report re-renders.

CONSTRAINTS
- New mode is additive; mmrt mode unchanged (gate is a hard requirement).
- Reuse EnzymeConstrainedModel/refresh_params, compute_tpc, TPC.descriptors; follow the
  MRes etcpy equations for the unfolding/NGAM forms; do not invent thermodynamics.
- Autonomous; commit parts separately: "add per-enzyme Tm + two-state unfolding thermal
  mode", "ingest MRes grounded Topt/Tm for eciML1515", "unfolding-mode re-run + report".
```
