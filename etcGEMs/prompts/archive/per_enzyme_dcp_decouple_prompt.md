# Claude Code prompt — per-enzyme thermal heterogeneity (dCp/Topt) to decouple breadth and Ea (autonomous)

Run from the project root (`.../MICROADAPT/etcGEMs`).

NOTE TO USER: launch in an auto-approving mode so it runs unattended. This is a model
change plus a recalibration and report re-runs; committed in git, so revertible.

---

```
Work AUTONOMOUSLY end to end; commit each part; print a summary. Read first:
src/etcgem/providers.py (from_gecko — where every EnzymeEntry currently gets Topt =
T0+default_Topt_offset and dCp = default_dCp), enzyme_cost.py (EnzymeEntry, refresh_
params), tpc.py (Ea + CT_max descriptors), thermal_sampling.py (nominal_thermal /
apply_thermal_sample / restore_thermal), dltkcat.py (apply_fits_to_provider sets
per-enzyme Topt/dCp), cli.py (cmd_calibrate_dcp, calibrate_dCp_to_Ea), and
strains/eciML1515/strain.yaml. Do NOT change the core MMRT / cost math.

PROBLEM (verified): the NOMINAL model is homogeneous — from_gecko assigns the SAME
Topt (T0+7) and SAME dCp to all ~2500 enzymes — so the organismal TPC is one MMRT
curve and a single dCp couples the rising-limb Ea to the falling-limb breadth/CT_max
(calibrating Ea=0.65 eV forces CT_max ~75.6 C). Fix: make the nominal thermal
parameters PER-ENZYME and heterogeneous (mean curvature sets Ea; the spread/tail sets
breadth/CT_max), so the two decouple. Source per-enzyme values from DLTKcat where
available and from a distribution otherwise.

PART 1 - per-enzyme thermal heterogeneity in the nominal build
- Extend the strain provider config with an optional thermal-distribution block
  (keep default_Topt_offset / default_dCp as the MEANS for backward compatibility):
    provider:
      ...
      Topt_sd: 0.0        # K, spread of enzyme optima about (T0 + default_Topt_offset)
      dCp_sd:  0.0        # spread of per-enzyme MMRT curvature about default_dCp
      thermal_seed: 0
      use_dltkcat_fits: false   # if true, override per-enzyme Topt/dCp from
                                #   strains/<name>/dltkcat/fits.csv where a good fit exists
- In from_gecko, when Topt_sd>0 or dCp_sd>0, draw each enzyme's Topt and dCp
  independently from Normal(mean, sd) (seeded by thermal_seed), clipped so dCp stays
  negative (< -1e-3) and Topt stays in a sane range; when both sd==0, behaviour is
  IDENTICAL to now (homogeneous) — this must be a byte-for-byte backward-compat path.
- If use_dltkcat_fits, after the draw, OVERRIDE per-enzyme Topt/dCp from the strain's
  dltkcat/fits.csv for enzymes with a usable fit (reuse dltkcat.apply_fits_to_provider
  or its per-entry logic; match on rxn_id/enzyme_id). Report how many enzymes were set
  from DLTKcat vs from the distribution. (DLTKcat currently covers few enzymes, so the
  distribution supplies heterogeneity for the rest; broadening DLTKcat coverage is a
  complementary data step, not required here.)

PART 2 - joint calibration that decouples Ea and CT_max
- Add calibrate_envelope(cfg, target_Ea_eV, target_CTmax_C, ...) that tunes TWO knobs
  to hit TWO targets using the (largely monotone, largely separable) relationships:
    * dCp_mean (= provider.default_dCp) controls the rising-limb Ea (steeper -> higher Ea);
    * Topt_sd controls breadth / CT_max (more spread -> lower CT_max, because the
      low-Topt tail denatures first).
  Implement as coordinate descent: bisect dCp_mean to hit Ea (at current Topt_sd), then
  bisect Topt_sd to hit CT_max (at current dCp_mean); iterate ~4-6 times to convergence.
  Each evaluation builds the provider and reads the nominal TPC Ea and CT_max
  (tpc.compute_tpc + TPC.descriptors). Keep dCp_sd fixed (e.g. proportional to |dCp_mean|)
  or expose it; guard dCp<0, Topt_sd>=0; extend the temperature grid if CT_max would run
  off it.
- Add a CLI `etcgem calibrate-envelope --strain NAME --target-ea 0.65 --target-ctmax 48`
  that runs it, PRINTS the calibrated (default_dCp, Topt_sd) and the achieved (Ea,
  CT_max), and updates strains/NAME/strain.yaml (default_dCp, Topt_sd, dCp_sd). Keep the
  existing scalar `calibrate-dcp` for backward compatibility.
- Run it for eciML1515 with target Ea 0.65 eV and a realistic target CT_max (e.g. 48-50 C;
  make it configurable — the PI can set it from the measured panel). Update strain.yaml.

PART 3 - re-run + report
- With the calibrated heterogeneous thermal parameters in strain.yaml, regenerate the
  report's source runs WITH figures (no --no-plots): default sweep, decompose
  (decomposition_sectors), control, calibrated, sectors; then
  `python reports/etcgem/assemble.py` and `cd reports/etcgem && quarto render`.
- Update the report's "# Interpretation and caveats" bullet that currently says one
  MMRT dCp couples breadth and Ea: state that per-enzyme heterogeneity now DECOUPLES
  them (Ea from the mean curvature, breadth/CT_max from the Topt spread / low-Topt
  tail), report the achieved Ea and CT_max, and note per-enzyme dCp is sourced from
  DLTKcat where available and from a distribution otherwise. Do not overclaim; keep it
  a structural/in-silico result.

VERIFY (report all)
1. BACKWARD-COMPAT GATE: with Topt_sd=0 and dCp_sd=0, the nominal TPC descriptors are
   byte-identical (within 1e-9) to the current code. Do not proceed if this fails.
2. DECOUPLING: after calibrate-envelope, the nominal TPC achieves Ea ~= target AND
   CT_max ~= target SIMULTANEOUSLY (which the single-scalar model could not) — report
   both numbers and the calibrated (default_dCp, Topt_sd).
3. Report how many enzymes got per-enzyme Topt/dCp from DLTKcat vs the distribution.
4. Report re-renders to _output/report.pdf with the updated caveat.

CONSTRAINTS
- Backward compatible: sd==0 -> identical numbers (item 1 is a gate).
- Reuse EnzymeEntry/refresh_params, compute_tpc, TPC.descriptors, apply_fits_to_provider,
  the calibration bisection pattern; do not change core MMRT/cost math.
- Autonomous; commit parts separately: "per-enzyme thermal heterogeneity in nominal
  build", "joint Ea/CT_max envelope calibration", "re-run + report update".
```
